%% Intrinsic connectivity distribution (ICD) survival curve
% icd.m is the survival curve function according to equation 2
% method for deriving ICD is adapted from Scheinost et al., 2012, for more
% information, see references.
% input: connectivity matrix
% output: variables of the survival curve for each culture either at
% baseline, during or after 4BF
%  Created by Yuheng Jiang
%  11 Oct 2020

%% Declare directories and function paths
outputpath = '';
matrixpath = [outputpath filesep 'matrix' filesep];
matrixfiles = dir([matrixpath '*_neuronall_matrix.mat']);
figspath = [outputpath filesep 'figure' filesep];

%requires function export_fig
addpath('...export_fig')

%% Initialise variables
iter = 1000;
alpha_values = zeros(length(matrixfiles),1);
beta_values = zeros(length(matrixfiles),1);

%% Run files
for ifile = 1:length(matrixfiles)

load([matrixpath matrixfiles(ifile).name])
net_str = sum(conn_matrix_o,2)/size(conn_matrix_o,1);

[x_tau, y_norm_deg] = normdeg(net_str, iter);

ft = fittype('exp(-alpha*x_tau^(beta))',...
    'dependent',{'y_norm_deg'},'independent',{'x_tau'},...
    'coefficients',{'alpha','beta'});
f = fit( x_tau, y_norm_deg, ft, 'StartPoint', [1, 0]);
coeff = coeffvalues(f);
alpha_values(ifile,1) = coeff(1);
beta_values(ifile,1) = coeff(2);

g = figure;
plot(f, x_tau, y_norm_deg)
xlabel('Correlation threshold')
ylabel('Normalised degree')

set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])

saveWellFrametime = matrixfiles(ifile).name;
saveWellFrametime = saveWellFrametime(1:20);
filename = [figspath saveWellFrametime 'icd_plot'];
export_fig(g, [filename '.jpg'], '-jpg', '-native')
close all
end

%% Local Functions
function [x_tau, y_norm_deg] = normdeg(net_str, iter)
x_tau = zeros(iter+1,1);
y_norm_deg = x_tau;

for ii = 1:iter+1
    tau = (ii-1)*0.001;
x_tau(ii,1) = tau;
norm_deg = sum(net_str>=tau)/size(net_str,1);
y_norm_deg(ii,1) = norm_deg;
end

end
