function run_pain_multislice(myomega, mygamma)
%
% Load data files obtained by 'pain_multislice' with modularity resolution
% gamma and multislice coupling strength omega, then compute and compare
% agreement difference matrices.
%
% Input parameters:
% - myomega: connection strength for multislice modularity (e.g. 0.1)
% - mygamma: modularity resolution parameter (e.g. 1.5)
%
% KL, last update: 2017/11.
% leibnitz@nict.go.jp

% these parameters don't need to be modified
thr = 0.1;      % correlation threshold remains fixed at 0.1
num_rois = 140; % number of ROIs
num_subs = 70;  % number of subjects in subject set

% next 2 lines should match the value that were set in pain_multislice.m
%num_slice_reps = 10; % number of repetitions per subject set
%num_reps = 1000;     % number of shuffles with new subject set
num_slice_reps = 10
num_reps = 1000

num_slices = num_subs * num_slice_reps; % number of slices per subject set

% Load the data from saved files. This requires that the results already
% had been computed by running pain_multislice with the shuffled parameter
% set to both true and false.
files1 = sprintf('results/agreement_cp_omega%.2f_gamma%.2f_thresh%.2f_rep%d_srep%d.mat', ...
    myomega, mygamma, thr, num_reps, num_slice_reps);
if (exist(files1, 'file') ~= 2)
    % files1 does not exist yet, so we must first generate it
    [CI_C, CI_P, Q_C, Q_P] = pain_multislice(myomega, mygamma, ...
        thr, false, num_reps, num_slice_reps);
    save(files1, 'CI_C', 'CI_P', 'Q_C', 'Q_P');
else
    load(files1, 'CI_C', 'CI_P', 'Q_C', 'Q_P');
end
CI_C_cp = CI_C; % modularity partitions group 1 (controls)
CI_P_cp = CI_P; % modularity partitions group 2 (pain)
Q_C_cp  = Q_C;  % modularity value of group 1 (controls)
Q_P_cp  = Q_P;  % modularity value of group 2 (pain)

% Do the same for the shuffled groups
files2 = sprintf('results/agreement_mix_omega%.2f_gamma%.2f_thresh%.2f_rep%d_srep%d.mat', ...
    myomega, mygamma, thr, num_reps, num_slice_reps);
if (exist(files2, 'file') ~= 2)
    % files2 does not exist yet, so we must first generate it
    [CI_C, CI_P, Q_C, Q_P] = pain_multislice(myomega, mygamma, ...
        thr, true, num_reps, num_slice_reps);
    save(files2, 'CI_C', 'CI_P', 'Q_C', 'Q_P');
else
    load(files2, 'CI_C', 'CI_P', 'Q_C', 'Q_P');
end
CI_mix1 = CI_C; % modularity partitions group 3 (mix1)
CI_mix2 = CI_P; % modularity partitions group 4 (mix2)
Q_mix1  = Q_C;  % modularity value of group 3 (mix1)
Q_mix2  = Q_P;  % modularity value of group 4 (mix2)

%% case 1: aggrement across all repetitions
Agr_C = agreement(CI_C_cp)/size(CI_C_cp,2);
Agr_P = agreement(CI_P_cp)/size(CI_P_cp,2);
AgrDiff_PC = Agr_P - Agr_C;
Agr_mix1 = agreement(CI_mix1)/size(CI_mix1,2);
Agr_mix2 = agreement(CI_mix2)/size(CI_mix2,2);
AgrDiff_mix = Agr_mix2 - Agr_mix1;

% figure 1: agreements and their differences for P, C, mix1, mix2
figure(1);
subplot(2,3,1); imagesc(Agr_C, [-1,1]); title('C'); axis('image');
subplot(2,3,2); imagesc(Agr_P, [-1,1]); title('P'); axis('image');
subplot(2,3,3); imagesc(AgrDiff_PC, [-1,1]); title('P-C'); axis('image'); 
subplot(2,3,4); imagesc(Agr_mix1, [-1,1]); title('mix1'); axis('image');
subplot(2,3,5); imagesc(Agr_mix2, [-1,1]); title('mix2'); axis('image');
subplot(2,3,6); imagesc(AgrDiff_mix, [-1,1]); title('mix2-mix1'); axis('image');

%% case 2: aggregation over each run
posval_CP = zeros(num_reps,num_rois);
negval_CP = zeros(num_reps,num_rois);
for i=1:num_reps
    CI_Ci = CI_C_cp(:,(i-1)*num_slices+1:i*num_slices);
    CI_Pi = CI_P_cp(:,(i-1)*num_slices+1:i*num_slices);
    Agr_Ci = agreement(CI_Ci)/size(CI_Ci,2);
    Agr_Pi = agreement(CI_Pi)/size(CI_Pi,2);
    AgrDiffi = Agr_Pi - Agr_Ci;
    posval_CP(i,:) = sum((AgrDiffi > 0) .* AgrDiffi,1);
    negval_CP(i,:) = sum((AgrDiffi < 0) .* AgrDiffi,1);
end
mn_posval_CP = mean(posval_CP,1);
mn_negval_CP = mean(negval_CP,1);
mn_abssum_CP = mn_posval_CP + abs(mn_negval_CP);

posval_mix = zeros(num_reps,num_rois);
negval_mix = zeros(num_reps,num_rois);
for i=1:num_reps
    CI_Ci = CI_mix1(:,(i-1)*num_slices+1:i*num_slices);
    CI_Pi = CI_mix2(:,(i-1)*num_slices+1:i*num_slices);
    Agr_Ci = agreement(CI_Ci)/size(CI_Ci,2);
    Agr_Pi = agreement(CI_Pi)/size(CI_Pi,2);
    AgrDiffi = Agr_Pi - Agr_Ci;
    posval_mix(i,:) = sum((AgrDiffi > 0) .* AgrDiffi,1);
    negval_mix(i,:) = sum((AgrDiffi < 0) .* AgrDiffi,1);
end
absval_mix = posval_mix + abs(negval_mix);

load('data/BSA_AAL_composite_140ROIs_Extended.mat','Names','XYZ_MNIcoordinates');
pvals = zeros(1,num_rois);
names = cell(1,num_rois);
xyz = zeros(num_rois,3);
posneg = [ mn_posval_CP', mn_negval_CP' ];
for i = 1:num_rois
    pcval    = mn_abssum_CP(i);
    mixvals  = absval_mix(:,i);
    pvals(i) = length(find(mixvals > pcval))/length(mixvals);
    names(i) = Names.BSA_SulcusBasedLong(i);
    xyz(i,:) = XYZ_MNIcoordinates(i,:); %#ok<IDISVAR,NODEF>
end
% T = table( sorted_CP(:,1), pvals', sorted_CP(:,2), string(names)' );
% T.Properties.VariableNames = { 'AgrDiff_Sum','p_value','ROI_ID','ROI_name' };
[pvals_sorted, idx] = sort(pvals, 'ascend');
pi = pvals_sorted < 0.05;
pvals_sorted = pvals_sorted(pi);
idx = idx(pi);

%% plot figures
% figure 1
figure(1);
clf
hold on
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.');
for i = 1:length(idx)
    fprintf('%.3f (%.3f, %.3f), %3d, %.3f, %s (%d,%d,%d)\n', ...
        sum(abs(posneg(idx(i),:))), posneg(idx(i),1), posneg(idx(i),2), ...
        idx(i), pvals_sorted(i), string(names(idx(i))), ...
        xyz(idx(i),1), xyz(idx(i),2), xyz(idx(i),3) );
    if (pvals_sorted(i) < 0.01)
        plot3(xyz(idx(i),1), xyz(idx(i),2), xyz(idx(i),3), 'r*');
    end
end
axis('equal');
%disp([ idx', pvals(idx)', string(names(idx))' ]);

%% figure 2: bar plots of C-P and mix1-mix2
figure(2);
clf
subplot(3,1,1); 
hold on
bar(1:num_rois,mean(posval_CP,1),'r'); 
bar(1:num_rois,mean(negval_CP,1),'b'); 
title('pos./neg. P-C'); 
ylim([-10,10]);
box on; 
subplot(3,1,2); 
hold on
bar(1:num_rois,mean(posval_mix,1),'r'); 
bar(1:num_rois,mean(negval_mix,1),'b'); 
title('pos./neg. mix'); 
ylim([-10,10]);
box on

subplot(3,1,3);
hold on
stairs(1:num_rois,mean(mn_abssum_CP,1),'m'); 
stairs(1:num_rois,mean(absval_mix,1),'k'); 
ylim([0,20]);
legend('control/pain', 'null model');
xlabel('ROIs')
box on;

%% plot Q-value histograms
figure(4);
clf
hold on
xx = linspace(min(Q_P_cp(:)),max(Q_C_cp(:)),100);
title(sprintf('Q-values: omega=%.2f, gamma=%.2f', myomega, mygamma));
[y1,x1] = hist(Q_C_cp(:), xx); 
stairs(x1, y1/sum(y1), 'b-');
[y2,x2] = hist(Q_P_cp(:), xx); 
stairs(x2, y2/sum(y2), 'r-');
[y3,x3] = hist(Q_mix1(:), xx); 
stairs(x3, y3/sum(y3), 'g-');
[y4,x4] = hist(Q_mix2(:), xx); 
stairs(x4, y4/sum(y4), 'c-');
box on
legend('C', 'P', 'mix1', 'mix2');
xlabel('Q-value');
ylabel('probability');

fprintf('mean(Q_C_cp) = %.3f, mean(Q_P_cp) = %.3f\n',mean(Q_C_cp(:)), mean(Q_P_cp(:)));
