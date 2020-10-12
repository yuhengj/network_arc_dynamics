%% Deg/Flow distribution, Hub neuron determination and flow through hubs
% determine hubs and flow through hubs
%  Created by Yuheng Jiang
%  11 Oct 2020

%%
outputpath = '';
idxpath = [outputpath filesep 'idx' filesep];
figspath = [outputpath filesep 'figure' filesep];
matrixpath = [outputpath filesep 'matrix' filesep];
matrixfiles = dir([matrixpath '*_neuronall_matrix.mat']);

hub_thresh = 0.85;
nbins = 10;

flow_tot = zeros(length(matrixfiles),1);
mean_flow_tot = flow_tot;
mean_flow_hub = flow_tot;
nhub_all = flow_tot;
%%
for ifile = 1:length(matrixfiles)
load([matrixpath matrixfiles(ifile).name])
net_deg = sum(conn_matrix_o);
[fi,xi] = ksdensity(net_deg);
h = figure;
histogram(net_deg, nbins, 'Normalization','pdf');
hold on
plot(xi,fi, '-r','LineWidth',1.5)
xlim([0 Inf])
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
hold off

mean_flow_tot(ifile) = mean(net_deg)/max(net_deg);
flow_tot(ifile) = sum(net_deg);


[freq, deg] = ecdf(net_deg); %calculate hubs based on degree distribution 90%
deg_thresh = interp1(freq,deg,hub_thresh);
idx_hub = find(net_deg >= deg_thresh);
n_hub = length(idx_hub); %get number of hubs
nhub_all (ifile) = n_hub;

net_deg_hub = net_deg(idx_hub);
[fhub,xhub] = ksdensity(net_deg_hub);
g = figure;
histogram(net_deg_hub, nbins,'Normalization','pdf');
hold on
plot(xhub,fhub, '-r','LineWidth',1.5)
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
hold off

str_hub = mean(net_deg_hub)/max(net_deg_hub); %normalised strength, divide by max strength in the network
mean_flow_hub (ifile) = str_hub;


saveWellFrame_all = matrixfiles(ifile).name;
saveWellFrame = saveWellFrame_all(1:20);
fprintf([saveWellFrame '\n'])

filename_one = [figspath saveWellFrame '_degree_all'];
export_fig(h, [filename_one '.jpg'], '-jpg', '-native')
filename_two = [figspath saveWellFrame '_degree_hub'];
export_fig(g, [filename_two '.jpg'], '-jpg', '-native')

close all
end