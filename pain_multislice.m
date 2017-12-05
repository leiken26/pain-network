function [CI_C, CI_P, Q_C, Q_P] = pain_multislice(myomega, mygamma, thresh, shuffled, num_reps, num_slice_reps)
%
% Run multislice modularity computation on pain and control subjects and
% store the results in a .mat file.
%
% Input parameters:
% - myomega: connection strength for multislice modularity (e.g. 0.1)
% - mygamma: modularity resolution parameter (e.g. 1.5)
% - thresh: (proportional) threshold for generating adjacency matrices (e.g. 0.1)
% - shuffled: boolean flag to indicate shuffling of subjects among classes
%   (e.g. false)
% - num_reps: number of total repetitions, i.e., random selection of 
%   subjects (e.g. 1000)
% - num_slice_reps: number of repetitions for same subject set (e.g. 10)
%
% Required software (besides Matlab):
% - Brain Connectivity Toolbox: https://sites.google.com/site/bctnet/
% - GenLouvain: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
%
% KL, last update: 2017/11.
% leibnitz@nict.go.jp

num_slices = 70;       % number of subjects (slices) for each class

% Load ROI information from file
load('data/BSA_AAL_composite_140ROIs_Extended.mat', 'BrodmannArea', ...
    'Names', 'Networks', 'XYZ_MNIcoordinates');
rois = Names.BSA_SulcusBasedLong;   % ROI names
num_rois = size(rois,1);            % number of ROIs

% generate adjacency matrices from correlation files
[AdjsC, AdjsP] = generate_adjacencies(thresh);

CI_C = zeros(num_rois, num_reps*num_slice_reps*num_slices);
CI_P = zeros(num_rois, num_reps*num_slice_reps*num_slices);
Q_C  = zeros(num_reps, num_slice_reps);
Q_P  = zeros(num_reps, num_slice_reps);

% loop for each repetition
for rep = 1:num_reps % outer loop   
    if (shuffled == true)
        [AdjsCrep, AdjsPrep] = shuffle_subs(AdjsC, AdjsP);
        shuffle_str = 'mix';
    else
        AdjsCrep = AdjsC;
        AdjsPrep = AdjsP;
        shuffle_str = 'C_P';
    end
    CI_Crep = zeros(num_rois, num_slice_reps*num_slices);
    CI_Prep = zeros(num_rois, num_slice_reps*num_slices);
    fprintf('%s %3d/%3d: ',shuffle_str, rep,num_reps);
    for slice_rep = 1:num_slice_reps % inner loop (just repeats modularity)
        fprintf('o');
        
        % randomly permute C subjects, truncate to num_slices
        rc = randperm(num_slices);
        % run genlouvain and store resulting modules in CI_Crep
        [Sc,Qc] = wrapper_genlouvain(AdjsCrep(rc),myomega,mygamma);
        CI_Crep(:,(slice_rep-1)*num_slices+1:slice_rep*num_slices) = Sc;
        Q_C(rep,slice_rep) = Qc;
        
        % now do the same with P
        rp = randperm(num_slices);
        % run genlouvain and store resulting modules in CI_Crep
        [Sp,Qp] = wrapper_genlouvain(AdjsPrep(rp),myomega,mygamma);
        CI_Prep(:,(slice_rep-1)*num_slices+1:slice_rep*num_slices) = Sp;
        Q_P(rep,slice_rep) = Qp;
    end
    % concatenate the computed modules from each repetition CI_Crep and
    % CI_Prep to CI_C and CI_P, respectively
    CI_C(:,(rep-1)*num_slice_reps*num_slices+1:rep*num_slice_reps*num_slices) = CI_Crep;
    CI_P(:,(rep-1)*num_slice_reps*num_slices+1:rep*num_slice_reps*num_slices) = CI_Prep;
    fprintf(' done\n');
end

%% Helper functions
    function [S,Q] = wrapper_genlouvain(A, myomega, mygamma)
        % wrapper function for GenLouvain code obtained from:
        %  http://netwiki.amath.unc.edu/GenLouvain/GenLouvain
        N=length(A{1});
        T=length(A);
        B=spalloc(N*T,N*T,(N+T)*N*T);
        twomu=0;
        for s=1:T
            k=sum(A{s});
            twom=sum(k);
            twomu=twomu+twom;
            indx=(1:N)+(s-1)*N;
            B(indx,indx)=A{s}-mygamma*(k'*k)/twom; %#ok<SPRIX>
        end
        twomu=twomu+T*myomega*N*(T-1);
        all2all = N*[(-T+1):-1,1:(T-1)];
        B = B + myomega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
        [S,Q] = genlouvain(B, 10000, 0);
        Q = Q/twomu;
        S = reshape(S,N,T);
    end

    function [AdjsC, AdjsP] = generate_adjacencies(thr)
        % load correlation matrices from files and generate adjacency
        % matrices from them by using threshold thr.
        % Data is from 3 sites (JP, UK, US) and has 2 classes (P, C)
        folders   = { 'UK_Cam_V5', 'UK_Cam_V5', 'US_OP_V5', 'US_OP_V5', ...
            'JP_CiN_V5', 'JP_CiN_V5' };
        fileflg   = { 'P', 'C' };    % indicates pain/control subject groups
        num_sites = length(folders); % number of subject classes
        
        % load correlation matrices, convert to adjacencies for each subject
        allAdj = cell(1,num_sites); % stores all adjacency matrices
        for i = 1:num_sites
            myflag = fileflg{1+mod(i,2)}; % odd folder is C, even is P
            mypref = sprintf('data/%s',folders{i}); % folder prefix
            mydir  = dir(sprintf('%s/ROICorrelation_%s*.mat',mypref,myflag));
            if (isempty(mydir) == true)
                error('ROICorrelation matrices incomplete');
            end
            num_subs = length(mydir); % number of subjects of this class
            
            siteAdj  = zeros(num_rois,num_rois,num_subs);
            % load each subject's file from site i and convert to
            % adjacency, store in siteAdj array
            for j = 1:num_subs
                load(sprintf('%s/%s',mypref,mydir(j).name), 'ROICorrelation');
                Adj = weight_conversion(threshold_proportional(ROICorrelation,...
                    thr),'binarize'); % generate binary adjacency matrix
                siteAdj(:,:,j) = Adj;
            end
            allAdj{i} = siteAdj; %
        end
        % concatenate all adjacencies from same class (C or P) into one big
        % cell array of all adjacency matrices
        AdjsC = squeeze(num2cell(cat(3,allAdj{1}, allAdj{3}, allAdj{5}), [1,2]));
        AdjsP = squeeze(num2cell(cat(3,allAdj{2}, allAdj{4}, allAdj{6}), [1,2]));
    end

    function [AdjsC, AdjsP] = shuffle_subs(AdjsC, AdjsP)
        % shuffle the order of adjacency matrices within each group,
        % then recombine them to contain half from C and half from P
        shuffle_C_idx = randperm(num_slices);
        shuffle_P_idx = randperm(num_slices);
        AdjsC1 = AdjsC(shuffle_C_idx(1:num_slices/2));
        AdjsC2 = AdjsC(shuffle_C_idx(num_slices/2+1:end));
        AdjsP1 = AdjsP(shuffle_P_idx(1:num_slices/2));
        AdjsP2 = AdjsP(shuffle_P_idx(num_slices/2+1:end));
        
        AdjsC = [AdjsC1; AdjsP1];
        AdjsP = [AdjsC2; AdjsP2];
    end
end