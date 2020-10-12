%% Initial Processing III
% generate overlay of arc postive nuclei over Arc microscopy images 
% generate overlay of fos (neun) postive nuclei over fos (neun) microscopy images 
% visual check for 1) quality of cultures 2) level of induction 
%  Created by Yuheng Jiang
%  11 Oct 2020

clearvars; close all
% Declare paths and imput parameters
outputpath = '';
%requires function export_fig 
addpath('')

% requires CircularHough() and MidpointCircle()
addpath('')

% requires bfmatlab files from BioFormats function: bfopen()
addpath('')

inputspecs = struct('nframe', [], 'dapi',[], 'txred', []);
inputspecs.nframe = 1;
inputspecs.dapi = 1;
inputspecs.txred = 2;

idxpath = [outputpath filesep 'idx' filesep];
if ~exist(idxpath, 'dir'); mkdir(idxpath); end

figspath = [outputpath filesep 'figure' filesep];
paramspath = [outputpath filesep 'params' filesep];

thresholdpath = [outputpath filesep 'threshold' filesep];
thresholdfiles = dir([thresholdpath '*thresholds.mat']);
threshold_idx = sort(repmat(1:length(thresholdfiles),[1 inputspecs.nframe]));

neuronpath = [outputpath filesep 'neuron' filesep];
neuronfiles = dir([neuronpath '*neuron.mat']);
labelpath = [outputpath filesep 'label' filesep];
labelfiles = dir([labelpath '*label.mat']);

%%
for ifile = 1:length(neuronfiles)

thresholdfile = [thresholdpath thresholdfiles(threshold_idx(ifile)).name];
load(thresholdfile)
neuronfilepath = [neuronpath neuronfiles(ifile).name];
load(neuronfilepath)
labelfilepath = [labelpath labelfiles(ifile).name];
load(labelfilepath)
saveWellFrame = neuronfiles(ifile).name;
saveWellFrame = saveWellFrame(1:17);
fprintf([saveWellFrame '\n'])

neuron_id = computeidx(thresholds, neuron);
save([idxpath filesep saveWellFrame '_idxarcfos.mat'], 'neuron_id')

arc_logi = ismember(label, neuron_id.arcpos_all); 
arc_label = double(label).* arc_logi;
arcfile = [paramspath num2str(threshold_idx(ifile),'%02d') '_ARC.mat'];
load (arcfile)
mat_idx = ifile + inputspecs.nframe - (inputspecs.nframe*threshold_idx(ifile));
arc_mat = MAT_ARC(mat_idx).matrix;
h = overlayseg(arc_mat, arc_label, saveWellFrame);
export_fig(h, [figspath saveWellFrame 'OverlayArc.jpg'], '-jpg', '-m2')

fos_logi = ismember(label, neuron_id.fospos_all); 
fos_label = double(label).* fos_logi;
fosfile = [paramspath num2str(threshold_idx(ifile),'%02d') '_CFOS.mat'];
load (fosfile)
fos_mat = MAT_CFOS(mat_idx).matrix;
g = overlayseg(fos_mat, fos_label, saveWellFrame);
export_fig(g, [figspath saveWellFrame 'OverlayFos.jpg'], '-jpg', '-m2')

end
%% Local Functions

function neuron_id = computeidx(thresholds, neuron)

neuron_cell = struct2cell(neuron);
Arc = cell2mat(neuron_cell(3,:));
Fos = cell2mat(neuron_cell(4,:));
neuron_idx = cell2mat(neuron_cell(1,:));

arc_pos_idx = neuron_idx(Arc>=thresholds.arc);
fos_pos_idx = neuron_idx(Fos>=thresholds.fos);

neuron_id = struct('arcpos_fospos', [], 'arcpos_fosneg', [] , 'arcneg_fospos', [], ...
    'arcneg_fosneg', [], 'arcpos_all', [], 'arcneg_all', [], 'fospos_all', [], 'fosneg_all', []);

neuron_id.arcpos_fospos = intersect(arc_pos_idx, fos_pos_idx);
neuron_id.arcpos_fosneg = arc_pos_idx(~ismember(arc_pos_idx,neuron_id.arcpos_fospos));
neuron_id.arcneg_fospos = fos_pos_idx(~ismember(fos_pos_idx,neuron_id.arcpos_fospos));
neuron_id.arcneg_fosneg = neuron_idx(~ismember(neuron_idx, unique(horzcat(arc_pos_idx,fos_pos_idx))));
neuron_id.arcpos_all = arc_pos_idx;
neuron_id.arcneg_all = neuron_idx(~ismember(neuron_idx,arc_pos_idx));
neuron_id.fospos_all = fos_pos_idx;
neuron_id.fosneg_all = neuron_idx(~ismember(neuron_idx,fos_pos_idx));

end

function h = overlayseg(matfile,label, saveWellFrame)
L_nucleus_perim = bwperim(label);
overlay5 = imadjust(matfile);
overlay5(L_nucleus_perim) = max(max(overlay5));

h = figure('Visible','off');
imagesc(overlay5); colormap('gray')
hold on
himage = imagesc(label2rgb(label, 'jet','w','shuffle'));
set(himage, 'AlphaData', 0.5);
varname = @(x) inputname(1);
titlestr = [saveWellFrame '\' varname(matfile) '\_Overlay'];
title(titlestr)
end
