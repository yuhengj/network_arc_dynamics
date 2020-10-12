%% Initial Processing II
% process calcium movies in nd2 to mat
%  Created by Yuheng Jiang
%  11 Oct 2020

close all
clearvars 
ntime = []; %to input time baseline=1, during=2, after=3

inputpath = '';
outputpath = '';

% requires bfmatlab files from BioFormats function: bfopen()
addpath('')

inputspecs = struct('ntime', [], 'imageType',[]);
inputspecs.ntime = ntime; %
inputspecs.imageType = 'raw';

%
rawcapath = [outputpath filesep 'rawca' filesep];
if ~exist(rawcapath, 'dir'); mkdir(rawcapath); end
figspath = [outputpath filesep 'figure' filesep];
if ~exist(figspath, 'dir'); mkdir(figspath); end
spikespath = [outputpath filesep 'spikes' filesep];
if ~exist(spikespath, 'dir'); mkdir(spikespath); end
paramspath = [outputpath filesep 'params' filesep];
if ~exist(paramspath, 'dir'); mkdir(paramspath); end
idxpath = [outputpath filesep 'idx' filesep];
if ~exist(idxpath, 'dir'); mkdir(idxpath); end

neuronpath = [outputpath filesep 'neuron' filesep];
neuronfiles = dir([neuronpath '*neuron.mat']);
labelpath = [outputpath filesep 'label' filesep];
labelfiles = dir([labelpath '*label.mat']);

inputfiles = dir([inputpath '*.nd2']);

neuronfile = [neuronpath neuronfiles(ifile).name];
load(neuronfile)
labelfile = [labelpath labelfiles(ifile).name];
load(labelfile)
saveWellFrame = neuronfiles(ifile).name;
saveWellFrame = saveWellFrame(1:17);

savepath = [rawcapath filesep saveWellFrame num2str(inputspecs.ntime,'%02d')]; 
mkdir(savepath)
datapath = [inputpath inputfiles(ifile).name];
fixedpointimageprocessing(datapath, savepath, inputspecs.imageType)

fprintf([saveWellFrame 'done! \n'])


%% Local Functions 

function fixedpointimageprocessing(datapath, savepath, imageType)
% function fixedpointimageprocessing(datapath)
% This is usually used to analyze movies, like calcium imaging, where the
% location is fixed and several images are taken over time.
%   Inputs: datapath
%           savepath
%           imageType - 'raw' or 'gray' to save the original color depth or
%           gray images (which, of course, saves a lot of space).
%   Outputs: raw or gray imagesc are saved as matrices in .mat format
% Frances Yeh 29/6/2016

% if isempty(strfind(datapath,'.nd2'))
%     datapath = [datapath '.nd2'];
% end
data = bfopen(datapath);
javahashtable = data{2};
maxAmp = 0;
maxAmpG = 0;
nPoints = size(data{1},1);
M = struct('idx',0,'matrix',zeros(size(data{1}{1,1})));
s = whos('M');
structLength = floor(0.5*10^9/s.bytes);
nFiles = ceil(nPoints/structLength);
% limit each file size to 0.5 GB so it doesn't overload the
% memory for other people. Also faster to load.
for is = 1:nFiles-1
    M = struct('idx',0,'matrix',zeros(size(data{1}{1,1})));
    M = repmat(M,structLength,1);
    startidx = (is-1) * structLength;
    for ipoint = 1:structLength
        rawM = data{1}{ipoint+startidx,1};
        maxAmp = max(max(max(rawM)),maxAmp);
        if strcmpi(imageType,'gray')
            normM = double(rawM)/2^16;
            grayM = ceil(normM*255);
            M(ipoint).matrix = uint8(grayM);
            maxAmpG = max(max(max(rawM)),maxAmpG);
        else
            M(ipoint).matrix = rawM;
        end
        M(ipoint).idx = ipoint+startidx;
        str = ['timestamp #' num2str(ipoint+startidx,'%04d')];
        if isempty(javahashtable.get(str))
            str = ['timestamp #' num2str(ipoint+startidx,'%03d')];
        end
        M(ipoint).timestamp = javahashtable.get(str); 
        if ipoint ==1 && startidx ==0
        starttime = M(1).timestamp;
        end
    end
    if strcmpi(imageType,'gray')
        filestr = ['gray' num2str(is-1,'%03d') '.mat'];
    else
        filestr = ['raw' num2str(is-1,'%03d') '.mat'];
    end
   
    save([savepath filesep filestr],'M','-v7.3')
end
is = nFiles;
startidx = (is-1) * structLength;
M = struct('idx',0,'matrix',zeros(size(data{1}{1,1})));
M = repmat(M,length(startidx+1:nPoints),1);
for ipoint = 1:length(startidx+1:nPoints)
    rawM = data{1}{ipoint+startidx,1};
    maxAmp = max(max(max(rawM)),maxAmp);
    if strcmpi(imageType,'gray')
        normM = double(rawM)/2^16;
        grayM = ceil(normM*255);
        M(ipoint).matrix = uint8(grayM);
        maxAmpG = max(max(max(rawM)),maxAmpG);
    else
        M(ipoint).matrix = rawM;
    end
    M(ipoint).idx = ipoint+startidx;
    str = ['timestamp #' num2str(ipoint+startidx,'%04d')];
    if isempty(javahashtable.get(str))
        str = ['timestamp #' num2str(ipoint+startidx,'%03d')];
    end
    M(ipoint).timestamp = javahashtable.get(str);
end
if strcmpi(imageType,'gray')
    filestr = ['gray' num2str(is-1,'%03d') '.mat'];
else
    filestr = ['raw' num2str(is-1,'%03d') '.mat'];
end
dataInfo.maxAmp = maxAmp;
dataInfo.nFrames = nPoints;
dataInfo.recordingLength = M(end).timestamp - starttime;
dataInfo.meanFrameRate = dataInfo.nFrames/dataInfo.recordingLength;
dataInfo.screenResolution = size(data{1}{1,1});
save([savepath filesep filestr],'M','-v7.3')
save([savepath filesep 'dataInfo.mat'],'dataInfo')

end