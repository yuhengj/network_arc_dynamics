%% compute differences in ICD graphs
% extension of ICD calculations 
%  Created by Yuheng Jiang
%  11 Oct 2020

%requires function export_fig
addpath('...export_fig')

outputpath = '';
figspath = [outputpath filesep 'figure' filesep];

matrixpath = [outputpath filesep 'matrix' filesep];
matrixfiles = dir([matrixpath '*_neuronall_matrix.mat']);

%load alpha and beta values obtained from ICD_calc.m
%
icd = @(alpha, x_tau, beta)...
    exp(-alpha.*x_tau.^(beta));


icd_diff_list = cell(length(alpha_values)/3,3);
ifile = 0;

for ii = 1:length(alpha_values)
    if rem(ii,3) == 1
       icd_ini = @(x)icd(alpha_values(ii),x,beta_values(ii));
       ifile  = ifile +1;
       icd_diff_list{ifile,1} = icd_ini;
    elseif rem(ii,3) == 0
        icd_fin = @(x)icd(alpha_values(ii),x,beta_values(ii));
        icd_diff = @(x)icd_fin(x)-icd_ini(x);
        icd_diff_list{ifile,2} = icd_fin;
        icd_diff_list{ifile,3} = icd_diff;
    else
    end
end

%%
for idxn = 1:length(alpha_values)/3
saveWellFrame = matrixfiles(idxn*3).name;
saveWellFrame = saveWellFrame(1:17);
fprintf([saveWellFrame '\n'])

g = figure;
fplot(icd_diff_list{idxn,1}, [0,1],'LineWidth',2)
hold on 
fplot(icd_diff_list{idxn,2}, [0,1],'LineWidth',2)

xlabel('Correlation threshold')
ylabel('Normalised degree')

set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
legend('Baseline','After 4BF')

filename = [figspath saveWellFrame 'icd_diff_plot'];
export_fig(g, [filename '.jpg'], '-jpg', '-native')
close all

end

%%
h = figure;
for idxnn = 1:length(alpha_values)/3
fplot(icd_diff_list{idxnn,3}, [0,1],'LineWidth',1.5)
hold on
end

xlabel('Correlation threshold')
ylabel('Difference in normalised degree')

set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
filename_a = [figspath 'icd_diff'];


export_fig(h, [filename_a '.jpg'], '-jpg', '-native')
%close all

save ([outputpath filesep 'icd_vals.mat'], 'icd_diff_list')