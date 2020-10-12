%% Determination of calcium peaks, network bursting

outputpath = '';
dfpath = [outputpath filesep 'params' filesep];
dffiles = dir([dfpath '*_dF_F.mat']);
figspath = [outputpath filesep 'figure' filesep];
addpath('...export_fig')

nwaves = zeros(length(dffiles),1);

for ifile = 1:length(dffiles)
load([dfpath dffiles(ifile).name])

saveWellFrame = dffiles(ifile).name;
saveWellFrame = saveWellFrame(1:20);

dF_all = nansum(dFcorrected,1);

[MPH,~] = findpeaks(ksdensity(dF_all),'npeaks',1);
MPH = 0.75 * MPH;
h = figure;
findpeaks(dF_all,'MinPeakHeight', MPH, 'MinPeakProminence', 3);
[PKS,~] = findpeaks(dF_all,'MinPeakHeight', MPH, 'MinPeakProminence', 3);
nwaves(ifile,1) = length(PKS);

set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
filename = [figspath saveWellFrame 'dF_F_peaks'];
export_fig(h, [filename '.jpg'], '-jpg', '-native')

end