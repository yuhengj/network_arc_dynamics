%% STTC calculations and generation of connectivity matrix
% takes spiketimes of individual neurons to compute connetivity via the
% STTC (spike time tiling coefficeient) 
% run_sttc.m is the funtion that is adapted from Cutts and Eglen (2014),
% for more information regarding STTC, see references
% output: connectivity matrix (per culture)
%  Created by Yuheng Jiang
%  11 Oct 2020

outputpath = '';
datainfopath = [outputpath filesep 'datainfo' filesep];
datainfofiles = dir([datainfopath '*.mat']);
spikespath = [outputpath filesep 'spikes' filesep];
spikesfiles = dir([spikespath '*.mat']);
idxpath = [outputpath filesep 'idx' filesep];
figspath = [outputpath filesep 'figure' filesep];
matrixpath = [outputpath filesep 'matrix' filesep];

%requires function run_sttc()
addpath('')
%requires function weight_conversion() in brain connectivity toolbox
addpath('...2017_01_15_BrainConnectivityToolbox')
%requires function export_fig
addpath('...export_fig')


dt = 0.5; %time window for calculating synchrony 
avgsyn = zeros(length(spikesfiles),1);
avgsyn_topten = avgsyn;

for ifile = 1:length(spikesfiles)
load([spikespath spikesfiles(ifile).name])
load([datainfopath datainfofiles(ifile).name])
Time = [0 dataInfo.recordingLength];

saveWellFrame_all = datainfofiles(ifile).name;
saveWellFrame = saveWellFrame_all(1:17);
fprintf([saveWellFrame '\n'])
load([idxpath saveWellFrame '_idxarcfos.mat'])
group_names = fieldnames(neuron_id);

neu_neg = length(getfield(neuron_id, group_names{4}));
neuronlist = zeros(length(spikeTimes)-neu_neg,1);
endidx = 0; 
 
for group_id = 1:3
    if group_id ==1
        startidx=1;
    else
        startidx = endidx+1;
    end
    list = getfield(neuron_id, group_names{group_id});
     endidx = endidx + length(list);
    neuronlist(startidx:endidx,1) = list;
end

neuronlist = sort(neuronlist);
spikeTimes = spikeTimes(neuronlist);
ncell = length(spikeTimes);
ncell_list = 1:ncell;
conn_matrix = zeros(ncell);
cross_corr_combinations = zeros(ncell,2);

for idxn = 1:ncell
    temp_sttc = zeros(ncell,1);
    if sum(isnan(spikeTimes(idxn).data))==0 ...
            && ~isempty(spikeTimes(idxn).data)
        cross_corr_combinations(:,1) = ncell_list(idxn);
        cross_corr_combinations(:,2) = ncell_list;
        for j = 1 : size(cross_corr_combinations,1)
            x = spikeTimes(cross_corr_combinations(j,1)).data;
            if sum(isnan(spikeTimes(cross_corr_combinations(j,2)).data)) == 0 ...
                    && ~isempty(spikeTimes(cross_corr_combinations(j,2)).data)
                y = spikeTimes(cross_corr_combinations(j,2)).data;
                sttc = run_sttc(length(x), length(y), dt, Time, x, y);
                temp_sttc(j,1) = sttc;
            else
                temp_sttc(j,1) = NaN;
            end
        end
        conn_matrix(idxn,:) = temp_sttc';
        
    else
        conn_matrix(idxn,:) = NaN;
        fprintf([saveWellFrame num2str(idxn) ' is not valid/is empty.\n']);
    end
end

conn_matrix_o = weight_conversion(abs(conn_matrix), 'autofix');

XDisplay = strings(1,ncell);
YDisplay = XDisplay;
h = figure;
heatmap(conn_matrix_o,'GridVisible', 'off', 'CellLabelColor', 'none',...
    'Colormap',parula, 'XDisplayLabels', XDisplay, 'YDisplayLabels', YDisplay);
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])

saveWellFrametime = datainfofiles(ifile).name;
saveWellFrametime = saveWellFrametime(1:20);
filename = [figspath saveWellFrametime 'neuron_all_connmatrix'];
export_fig(h, [filename '.jpg'], '-jpg', '-native')

nsttc = sum(conn_matrix_o,1);
avgsyn(ifile) = mean(nsttc);

topsttc = nsttc(nsttc>=prctile(nsttc(~isnan(nsttc)),90));
avgsyn_topten(ifile) = mean(topsttc);

savename = saveWellFrame_all(1:20);
save([matrixpath savename 'neuronall_matrix.mat'], 'conn_matrix_o', 'conn_matrix')

clearvars cross_corr_combinations idxn j temp_sttc
end

