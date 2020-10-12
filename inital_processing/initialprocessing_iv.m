%% Initial Processing IV
% to update shift in position during recording of movies vs staining post
% calcium imaging 
% running the code will generate neuronal nuclei overlays over calcium
% imaging movie
% manually input values for shifting the overlays (dy, dx)
% apply labels and detect calcium transients, get spiketimes
% update neuron structs with spiketimes{}
% required fast_oopsi() from FluoroSNNAP package codes
% output: calcium movies, spiketimes
%  Created by Yuheng Jiang
%  11 Oct 2020

close all;
clearvars;

dy = []; %down=positive
dx = []; %right=positive

ntime = [];
ifile = [];

outputpath = '';

%requires function export_fig 
addpath('/export_fig')
% requires fast_oopsi from FluoroSNNAP
addpath('.../FluoroSNNAP15.04.18/FluoroSNNAP15.04.18/FluoroSNNAP_code/oopsi-master')
% requires bfmatlab files from BioFormats function: bfopen()
addpath('.../bfmatlab')


inputspecs = struct('nframe', [], 'dapi',[], 'txred', [], 'ntime', [], 'imageType',[]);
inputspecs.nframe = 1;
inputspecs.dapi = 1;
inputspecs.txred = 2;
inputspecs.imageType = 'raw';
inputspecs.ntime = ntime;

moviespecs = struct('gamma',[],'frameRate',[],'amp', [], 'lowerThresh',[]);
moviespecs.gamma = 3.0;



spikespath = [outputpath filesep 'spikes' filesep];
if ~exist(spikespath, 'dir'); mkdir(spikespath); end


neuronpath = [outputpath filesep 'neuron' filesep];
neuronfiles = dir([neuronpath '*neuron.mat']);

figspath = [outputpath filesep 'figure' filesep];
paramspath = [outputpath filesep 'params' filesep];
idxpath = [outputpath filesep 'idx' filesep];

thresholdpath = [outputpath filesep 'threshold' filesep];
thresholdfiles = dir([thresholdpath '*thresholds.mat']);
threshold_idx = sort(repmat(1:length(thresholdfiles),[1 inputspecs.nframe]));

labelpath = [outputpath filesep 'label' filesep];
labelfiles = dir([labelpath '*label.mat']);
rawcapath = [outputpath filesep 'rawca' filesep];

moviepath = [outputpath filesep 'movie' filesep];
if ~exist(moviepath, 'dir'); mkdir(moviepath); end




labelfilepath = [labelpath labelfiles(ifile).name];
load(labelfilepath)
saveWellFrame = labelfiles(ifile).name;
saveWellFrame = saveWellFrame(1:17);
fprintf([saveWellFrame '\n'])

load([idxpath filesep saveWellFrame '_idxarcfos.mat'])

%ca movie
datapath = [rawcapath saveWellFrame num2str(inputspecs.ntime,'%02d') filesep];
load([datapath 'dataInfo.mat'])
frames = 1:dataInfo.nFrames;
M = loadrawdata(datapath,frames);
moviespecs.frameRate = dataInfo.nFrames/dataInfo.recordingLength;
neuronLabels = shiftLabels(label, dy, dx);

arcpos_fospos = double(ismember(neuronLabels, neuron_id.arcpos_fospos))*4;
arcpos_fosneg = double(ismember(neuronLabels, neuron_id.arcpos_fosneg))*3;
arcneg_fospos = double(ismember(neuronLabels, neuron_id.arcneg_fospos))*2;
arcneg_fosneg = double(ismember(neuronLabels, neuron_id.arcneg_fosneg));
segLabels = arcpos_fospos + arcpos_fosneg + arcneg_fospos + arcneg_fosneg;

moviespecs.amp = double(dataInfo.maxAmp);
moviespecs.lowerThresh = moviespecs.amp*1/50;
moviespecs.frameRate = dataInfo.meanFrameRate;


h = figure;
set(h,'Position',[10 10 1000 500])
savename = [saveWellFrame num2str(inputspecs.ntime,'%02d') '_with Neuron Labels.avi'];
mvobj = VideoWriter([moviepath savename]);
set(mvobj,'FrameRate', moviespecs.frameRate); %20
open(mvobj);
makemoviewithlabels(h, mvobj,M,segLabels,moviespecs.gamma,...
    moviespecs.amp, moviespecs.lowerThresh, dataInfo)
close(mvobj)

%
neuronfile = [neuronpath neuronfiles(ifile).name];
load(neuronfile)

%%
V = struct('dt',[],'est_a', [], 'est_b',[], 'est_gam', [], 'est_lam', [], 'est_sig',[]);
V.est_a = 1; %1
V.est_b = 1; %0
V.est_gam = 1;
V.est_lam = 1; %0.1
V.est_sig = 1; %0.1
%V.fast_poiss = 0;
%   est_sig:    1 to estimate sig
%   est_lam:    1 to estimate lam
%   est_gam:    1 to estimate gam
%   est_b:      1 to estimate b
%   est_a:      1 to estimate a

%   fast_poiss:     1 if F_t ~ Poisson, 0 if F_t ~ Gaussian
%   fast_nonlin:    1 if F_t is a nonlinear f(C_t), and 0 if F_t is a linear f(C_t)
%   fast_plot:      1 to plot results after each pseudo-EM iteration, 0 otherwise
%   fast_thr:       1 if thresholding inferred spike train before estiamting {a,b}
%   fast_iter_max:  max # of iterations of pseudo-EM  (1 to use default initial parameters)
%   fast_ignore_post: 1 to keep iterating pseudo-EM even if posterior is not increasing, 0 otherwise

P = struct('gam', []);
%   a:      spatial filter
%   b:      background fluorescence
%   sig:    standard deviation of observation noise
%   gam:    decayish, ie, tau=dt/(1-gam)
%   lam:    firing rate-ish, ie, expected # of spikes per frame

%%
thresholdspecs = struct('lowpct', [], 'threshold', []);
thresholdspecs.lowpct = 3; % percentile to be used as lowca, background variation

%%
thresholdspecs.threshold = 3; %
dt = dataInfo.recordingLength/dataInfo.nFrames;

%%
[catrace, timevec] = computeraw(neuronLabels, M);

%figure; plot(timevec, catrace)

%%
lowpct = thresholdspecs.lowpct; 
[dFcorrected, idx_low, sum_ca, pct_low, catrace_low] = computecorrection(catrace,lowpct,timevec);
%figure; plot(timevec, dFcorrected)

%%
V.dt = dt;
P.gam = 1-V.dt/1.5;
threshold = thresholdspecs.threshold;

%%
spikeTimes = computespikes(dFcorrected, V, P, threshold, sum_ca, pct_low, timevec);

%%
for in = 1:length(spikeTimes)
    neuron(in).spikeTimes{inputspecs.ntime} = spikeTimes(in).data;
end

%%
rasterplot = plotraster(spikeTimes);

%%
export_fig(rasterplot, [figspath saveWellFrame num2str(inputspecs.ntime,'%02d') '_Raster.jpg'], '-jpg', '-m2')


%%
save([neuronpath saveWellFrame 'neuron.mat'], 'neuron')
save([spikespath saveWellFrame num2str(inputspecs.ntime,'%02d') '_spikes.mat'], 'spikeTimes')
save([paramspath saveWellFrame num2str(inputspecs.ntime,'%02d') '_dF_F.mat'],'dFcorrected')
save([paramspath saveWellFrame num2str(inputspecs.ntime,'%02d') '_caparams.mat'], 'datapath',...
    'inputspecs', 'V', 'P', 'thresholdspecs', 'dx', 'dy' )
save([idxpath saveWellFrame num2str(inputspecs.ntime,'%02d') '_idxlow.mat'], 'idx_low')

fprintf([saveWellFrame 'done! \n'])


%% Local Functions
function Mstruct = loadrawdata(savepath,frames)
allfiles = dir([savepath filesep 'raw*.mat']);
maxFrame = max(frames);
Mstruct = [];
idx = 0;
sprintf('Loading data, please wait. ')
while length(Mstruct) < maxFrame
    idx = idx + 1;
    load([savepath filesep allfiles(idx).name])
    Mstruct = [Mstruct; M];
end

Mstruct = Mstruct(frames);
sprintf('Done!')
end

function shiftedLabels = shiftLabels(nucleusLabel, downshift, rightshift)
% function shiftedLabels = shiftLabels(nucleusLabel, downshift, rightshift)
%
% Shift ROI (nucleus labels) by [downshift, rightshift] pixels, replacing
% shifted pixels with 0, producing final matrix of same size as input
% inputs:
% nucleusLabels = from labeled matrix of segmented objects
% downshift = number of pixels to move down (negative: up)
% rightshift = number of pixels to move right (negative: left)
%
% created 2/11/2016

[m, n] = size(nucleusLabel);
L_shift = circshift (nucleusLabel, [downshift, rightshift]);

if downshift >0
    L_shift ((1:downshift),:) = 0;
elseif downshift <0
    L_shift(((m-downshift+1):m),:) = 0;
end

if rightshift >0
    L_shift (:,(1:downshift)) = 0;
elseif downshift <0
    L_shift(:,((n-downshift+1):n)) = 0;
end

shiftedLabels = L_shift;

end


function makemoviewithlabels(h, mvobj,M,neuronLabels,gamma,...
    amp,lowerThresh, dataInfo)
% function makemoviewithdynamiclabels(savepath,catrace,neuronLabels,gamma,...
%     amp,lowerThres,frameRate)
% Make movies with neurons label overlays
% camax = max(max(neuron));
% camin = min(min(neuron));
% caamp = (camax - camin);
% % cafactor = 100/caamp; % image(M) with single elements in M, colors are 6-bit
% % % I'm using 100 to scale in order to make everything brighter and easier to
% % % see.
% cafactor = 64/caamp;
for im = 1:length(M)
    m = double(M(im).matrix);
    localamp = max(max(m));
    localampgray = ceil(localamp*255/amp);
    %     m(m<lowerThres) = 0;
    %     imgamma = gamma_correction(m,[0 amp],[0 255],gamma);
    imgamma = gamma_correction(m,[lowerThresh amp],[0 localampgray],gamma);
    %     imgamma(imgamma<lowerThresgray) = 0;
    image(imgamma); colormap('gray')
    hold on
    himage = imagesc(label2rgb(neuronLabels, 'jet'));
    set(himage, 'AlphaData', 0.2);
    title(['data taken at ' num2str(dataInfo.meanFrameRate) ...
        ' fps; current time = ' num2str(M(im).timestamp,'%5.2f') ' secs']);
    axis equal; axis off;
    
    %     segm = zeros(size(m));
    %     for ineuron = 1:4
    %         segm(neuronLabels == ineuron) = ceil((neuron(ineuron,im)-camin) * cafactor);
    %     end
    %     subplot(1,2,2); image(segm);
    %     axis equal; axis off; colormap('bone')
    frame = getframe(h);
    writeVideo(mvobj,frame);
    cla
    
    % L_nucleus_perim = bwperim(label);
    % overlay5 = imadjust(matfile);
    % overlay5(L_nucleus_perim) = max(max(overlay5));
    % nucleus_no = numel(unique(label))-1;
    %
    % h = figure('Visible','off');
    % imagesc(overlay5); colormap('gray')
    % hold on
    % himage = imagesc(label2rgb(label, 'jet','w','shuffle'));
    % set(himage, 'AlphaData', 0.2);
    % title_str = sprintf('Nucleus Detected, Well %d Frame %d, Total number %d', ...
    %     well_no, frame_no, nucleus_no);
    % title(title_str)
end

end

function [catrace, timevec] = computeraw(neuronLabels, M)
totcells = max(max(neuronLabels));
catrace = zeros(totcells,length(M));
timevec = zeros(1,length(M));
offset = M(1).timestamp;
for iframe = 1:length(M)
    stats = regionprops(neuronLabels,M(iframe).matrix,'MeanIntensity');
    for ineuron = 1:max(max(neuronLabels))
        % if neuronlabel ineuron ~exist, value is NaN
        catrace(ineuron,iframe) = stats(ineuron).MeanIntensity;
    end
    timevec(1,iframe) = M(iframe).timestamp - offset;
end
end

function [dFcorrected, idx_low, sum_ca, pct_low, catrace_low] = computecorrection(catrace,lowpct,timevec)
%%% avgerage of silent neurons

% cannot use low ca in some cases, have to use low variability
sum_ca = max(catrace,[],2) - min(catrace,[],2);
% sum(catrace,2);
% prctile automatically removes NaN in calculations
pct_low = prctile(sum_ca,lowpct);
idx_low = find(sum_ca<pct_low);
catrace_low = catrace(sum_ca<pct_low,:);
catrace_low_avg = mean(catrace_low);

delta_t = 5; %first 5 secs of datachange time, NaN traces remain as NaN throughout
t0idx = find(timevec < delta_t,1,'last');
pct_1 = catrace_low_avg(timevec < delta_t);
F0 = repmat(prctile(pct_1,10),size(catrace(:,1:t0idx)));
dF0 = catrace(:,1:t0idx) - F0;
dF0_F0 = dF0./F0;
dFcorrected = zeros(size(catrace));
dFcorrected(:,1:t0idx) = dF0_F0;

for it = (t0idx+1):length(timevec) %
    t = timevec(it);
    tidx = find(timevec > t-delta_t,1); %%
    pct = catrace_low_avg(tidx:it-1);
    F = repmat(prctile(pct,10),size(catrace(:,it)));
    dF = catrace(:,it) - F;
    dF_F = dF./F;
    dFcorrected(:,it) = dF_F;
end

end

function spikeTimes = computespikes(dFcorrected, V, P, threshold, sum_ca, pct_low, timevec)

inferredspikes = zeros(size(dFcorrected));
for in = 1:size(dFcorrected,1)
    if sum(isnan(dFcorrected(in,:)))>0 %if any value is NaN, trace is not computed
        inferredspikes(in,:) = NaN;
    else
        inferredspikes(in,:) = fast_oopsi(dFcorrected(in,:),V,P);
    end
end

low_spikes = inferredspikes(sum_ca<pct_low,:);%idx of 5th pct
thresh = mean(std(low_spikes,0,2))*threshold;
%exclude traces that are zero

spikeTimes = struct('data',[], 'unit',[]);
spikeTimes = repmat(spikeTimes, size(inferredspikes,1),1);
for in = 1:size(inferredspikes,1)
    data = inferredspikes(in,:);
    if sum(isnan(data))>0
        spikeTimes(in).data = NaN;%NaN trace should be different from a quiet cell
    else
        %idx = find(data > thresh);
        %diffidx = [0; diff(idx)];
        %idx(diffidx == 1) = [];
        spikeTimes(in).data = timevec(data > thresh); %timevec(idx)
    end
    spikeTimes(in).unit = in;
end
end

function i = plotraster(spikeTimes)
% Plot spike times from all neurons; raster plot
i = figure('Visible', 'on');
for in = 1:length(spikeTimes)
    hold on;
    plot(spikeTimes(in).data,in*ones(size(spikeTimes(in).data)),'b.')
end
title('Deconvoluted spikes')
xlabel('[s]')
ylabel('neuron ID')
end





