%% Calculation pf firing rates
% compute firing rates 
%  Created by Yuheng Jiang
%  11 Oct 2020

outputpath = '';
datainfopath = [outputpath filesep 'datainfo' filesep];
datainfofiles = dir([datainfopath '*.mat']);
spikespath = [outputpath filesep 'spikes' filesep];
spikesfiles = dir([spikespath '*.mat']);
idxpath = [outputpath filesep 'idx' filesep];


avgfr = zeros(length(spikesfiles),1);
avgfr_topten = avgfr;

for ifile = 1:length(spikesfiles)
load([spikespath spikesfiles(ifile).name])
load([datainfopath datainfofiles(ifile).name])

saveWellFrame = datainfofiles(ifile).name;
saveWellFrame = saveWellFrame(1:17);
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

spikeTimes = spikeTimes(neuronlist);
ncell = length(spikeTimes);
nspike = zeros(ncell,1);
for ii = 1:ncell
    if isnan(spikeTimes(ii).data) == 0
        nspike(ii) = length(spikeTimes(ii).data);
    else
        nspike(ii) = NaN;
    end
end

fr = mean(nspike(~isnan(nspike)))/dataInfo.recordingLength;
avgfr(ifile) = fr;

topspikes = nspike(nspike>=prctile(nspike(~isnan(nspike)),90));
avgfr_topten(ifile) = mean(topspikes)/dataInfo.recordingLength;

end