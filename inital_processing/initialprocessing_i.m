%% Inital Processing I
% generation of neuron structs with nucleus location, arc, fos (or neun), and spike
% times{empty}; generate list of exp for thresholding
%  Created by Yuheng Jiang
%  11 Oct 2020

clearvars; close all
% Declare paths and imput parameters
inputpath = ''; %input microscopy images

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

segspecs = struct('backgroundsize', [], 'meshsize',[], 'gammafactor', [], ...
    'disksize', [], 'houghbounds',[], 'radius', [], 'sizelow', [],...
    'sizehigh', [],'intensity', []);
segspecs.backgroundsize = 200;
segspecs.meshsize = 19;
segspecs.gammafactor = 0.8;
segspecs.disksize = 3; 
segspecs.houghbounds = [5 10];
segspecs.radius = 2; 
segspecs.sizelow = 20; 
segspecs.sizehigh = 2000;
segspecs.intensity = 250; 

thresholdspecs = struct('sigma', [], 'options', []); %times of sigma
thresholdspecs.sigma = 3;
thresholdspecs.options = statset('MaxIter',10000, 'MaxFunEvals',100000);


labelpath = [outputpath filesep 'label' filesep];
if ~exist(labelpath, 'dir'); mkdir(labelpath); end
neuronpath = [outputpath filesep 'neuron' filesep];
if ~exist(neuronpath, 'dir'); mkdir(neuronpath); end
figspath = [outputpath filesep 'figure' filesep];
if ~exist(figspath, 'dir'); mkdir(figspath); end
thresholdpath = [outputpath filesep 'threshold' filesep];
if ~exist(thresholdpath, 'dir'); mkdir(thresholdpath); end
paramspath = [outputpath filesep 'params' filesep];
if ~exist(paramspath, 'dir'); mkdir(paramspath); end


%%
allData = dir([inputpath '*.nd2']);

for ifile = 1:length(allData)
    ND2file = [inputpath allData(ifile).name];
    [MAT_DAPI, MAT_ARC, MAT_CFOS] = convertND2toMAT(ND2file, inputspecs);
    save([paramspath num2str(ifile,'%02d') '_ARC.mat'],'MAT_ARC')
    save([paramspath num2str(ifile,'%02d') '_CFOS.mat'],'MAT_CFOS')
    save([paramspath num2str(ifile,'%02d') '_DAPI.mat'],'MAT_DAPI')
    
    avg_exp = struct('arc',[],'fos',[]);
    avg_exp = repmat(avg_exp, length(MAT_ARC),1);
    arc_exp=cell(1,inputspecs.nframe);
    fos_exp=cell(1,inputspecs.nframe);
    
    for ii = 1:inputspecs.nframe
        
        %nucleus segmentation
        dapi_mat = MAT_DAPI(ii).matrix;
        label = detect_nucleus(dapi_mat, segspecs);
        saveName = ['Well ' num2str(ifile,'%02d') ' Frame ' num2str(ii,'%02d')];
        save([labelpath saveName '_label.mat'], 'label')
        
        %export nucleus overlay
        h = overlayseg(dapi_mat,label, ifile, ii);
        export_fig(h, [figspath saveName '_dapi'], '-jpg', '-m2')
        sprintf(['Done: ' saveName 'nucleus'])
        
        % nucleus expression
        arc_mat = MAT_ARC(ii).matrix;
        fos_mat = MAT_CFOS(ii).matrix;
        [neuron, avg_intensity_arc, avg_intensity_fos, ~] = ...
            nuclearExpression(label, arc_mat, fos_mat);
        save([neuronpath saveName '_neuron.mat'], 'neuron')
        sprintf(['Done: ' saveName 'expression'])
        
        arc_exp{ii} = avg_intensity_arc;
        fos_exp{ii} = avg_intensity_fos;
    end
    arc_exp = vertcat(arc_exp{:});
    fos_exp = vertcat(fos_exp{:});
    thresholds = struct('arc', [], 'fos', []);
    thresholds.arc = mle_threshold(arc_exp, thresholdspecs);
    thresholds.fos = mle_threshold(fos_exp, thresholdspecs);
    save([thresholdpath num2str(ifile,'%02d') '_thresholds.mat'], 'thresholds')
    
    sprintf('Done: Well %d...', ifile)
    save([paramspath num2str(ifile,'%02d') '_exp.mat'], 'avg_exp', 'arc_exp', 'fos_exp')
end

save([paramspath 'userparams.mat'], 'inputspecs', 'inputpath', 'segspecs', 'thresholdspecs')

%% Local Functions

function [MAT_DAPI, MAT_ARC, MAT_CFOS] = convertND2toMAT(ND2file, inputspecs)
% Covert ND2 file to MAT array

    data = bfopen(ND2file);
    MAT_ARC = struct('matrix',zeros(1002,1004,'uint16'));
    MAT_ARC = repmat(MAT_ARC,inputspecs.nframe,1);
    MAT_CFOS = struct('matrix',zeros(1002,1004,'uint16'));
    MAT_CFOS = repmat(MAT_CFOS,inputspecs.nframe,1);
    MAT_DAPI = struct('matrix',zeros(1002,1004,'uint16'));
    MAT_DAPI = repmat(MAT_DAPI,inputspecs.nframe,1);
    
    for nf = 1:inputspecs.nframe
        for cid = 1:3
            if cid == inputspecs.dapi
                rawDAPI =  data{nf,1}{cid,1};
                M_dapi = uint16(rawDAPI);
                MAT_DAPI(nf).matrix = M_dapi;
            elseif cid == inputspecs.txred
                rawARC = data{nf,1}{cid,1};
                M_arc = uint16(rawARC);
                MAT_ARC(nf).matrix = M_arc;
            else
                rawFOS = data{nf,1}{cid,1};
                M_fos = uint16(rawFOS);
                MAT_CFOS(nf).matrix = M_fos;
            end
        end
    end
end

function label = detect_nucleus(dapi_mat, segspecs)
% Detection of nucleus via segmentation and CircularHough method
% Input: file path containing grayDAPI mat file 
% Output: nucleus_detected.mat in the same directory; jpg files for each
% nd2 picture with detected nucleus outlined
% variables to be optimized: gammafactor; CircularHough parameters;
% lower_cutoff; higher_cutoff; intensity_cutoff

% load params from input
background_size = segspecs.backgroundsize;
mesh_size = segspecs.meshsize;
gammafactor = segspecs.gammafactor;
disksize = segspecs.disksize;
houghbounds = segspecs.houghbounds;
radius = segspecs.radius;
lower_cutoff = segspecs.sizelow;
higher_cutoff = segspecs.sizehigh;
intensity_cutoff = segspecs.intensity;

I_a = dapi_mat;
I_b = double(I_a)/2^16;
I = ceil(I_b*255);

background = imopen(I,strel('disk',background_size));
background_I = I - background;

mesh = imopen(background_I,strel('disk',mesh_size));
mesh_I = background_I - mesh;

open_I = gamma_correction(mesh_I,[0 65],[0 255],gammafactor); 
open_I(open_I < 50) = 0;

bw = im2bw(open_I, graythresh(open_I));
h = fspecial('disk', disksize);
bw2 = imfilter(bw,h);

[~, circen, ~] = CircularHough_Grd(open_I, houghbounds);

mask_em = zeros(size(I));
 
for icir = 1:length(circen)
    if (circen(icir,1) - radius > 0 && circen(icir,1) -  radius > 0)
    mask_em = MidpointCircle(mask_em, radius, circen(icir,1), circen(icir,2),255);
    end
end
mask_em = logical(mask_em); 

I_eq_c = imcomplement(open_I);
I_mod = imimposemin(I_eq_c, ~bw2 | mask_em); 

L = watershed(I_mod);
L = imclearborder(L, 4);

cc = bwconncomp(L, 4);
nucleusdata = regionprops(cc,'basic');
nucleus_areas = [nucleusdata.Area];

skewfactor = mean(nucleus_areas)/median(nucleus_areas);
if skewfactor > 1
    P = mean(nucleus_areas) - (3/skewfactor)*std(nucleus_areas);
    Q = mean(nucleus_areas) + (3*skewfactor)*std(nucleus_areas);
elseif skewfactor == 1
    P = mean(nucleus_areas) - 3*std(nucleus_areas);
    Q = mean(nucleus_areas) + 3*std(nucleus_areas);
else
    P = mean(nucleus_areas) - (3*skewfactor)*std(nucleus_areas);
    Q = mean(nucleus_areas) + (3/skewfactor)*std(nucleus_areas);
end

label_matrix = labelmatrix(cc);

if P >= lower_cutoff && Q <= higher_cutoff
    L_fil = ismember(label_matrix, find(nucleus_areas >= P & nucleus_areas <= Q));
elseif P <= lower_cutoff && Q <= higher_cutoff
    L_fil = ismember(label_matrix, find(nucleus_areas >= lower_cutoff & nucleus_areas <= Q));
elseif P >= lower_cutoff && Q >= higher_cutoff
    L_fil = ismember(label_matrix, find(nucleus_areas >= P & nucleus_areas <= higher_cutoff));
else 
    L_fil = ismember(label_matrix, find(nucleus_areas >= lower_cutoff & nucleus_areas <= higher_cutoff));
end


L_fil2=zeros(size(L));
L_fil2(L_fil==1)=L(L_fil==1);

L_fil = L_fil2;

cc_fil = bwconncomp(L_fil, 4);

final_no = numel(unique(L_fil))-1;
nucleus_mask = struct('idx',0,'mask',zeros(1002,1004),'matrix',zeros(1002,1004));
nucleus_mask = repmat(nucleus_mask, final_no,1);
L_nucleus = labelmatrix(cc_fil);

avg_intensity = zeros(final_no,1);

remove_index = zeros(final_no,1);
id = 1;
for inucleus = 1:final_no
    nucleus_mask(inucleus).idx = inucleus;
   nucleus_mask(inucleus).mask = L_nucleus==inucleus;
   logi_mask = logical(nucleus_mask(inucleus).mask);
   nucleus_mask(inucleus).matrix(logi_mask ==1) = open_I(logi_mask==1);
   avg_intensity(inucleus) = mean(nucleus_mask(inucleus).matrix(nucleus_mask(inucleus).matrix~=0));
   if final_no > 50
   if avg_intensity(inucleus) > intensity_cutoff
       remove_index(id) = inucleus;
       id = id+1;
   end
   end
end
logical_mask = ~logical(ismember(L_nucleus, remove_index));
L_nucleus2 = zeros(size(L_nucleus));
L_nucleus2(logical_mask==1) = L_nucleus(logical_mask==1);

logicL = logical(L_nucleus2);
cc = bwconncomp(logicL, 4);
label = labelmatrix(cc);
end

function h = overlayseg(matfile,label, well_no, frame_no)
L_nucleus_perim = bwperim(label);
overlay5 = imadjust(matfile);
overlay5(L_nucleus_perim) = max(max(overlay5));
nucleus_no = numel(unique(label))-1;

h = figure('Visible','off');
imagesc(overlay5); colormap('gray')
hold on
himage = imagesc(label2rgb(label, 'jet','w','shuffle'));
set(himage, 'AlphaData', 0.2);
title_str = sprintf('Nucleus Detected, Well %d Frame %d, Total number %d', ...
    well_no, frame_no, nucleus_no);
title(title_str)
end

function [neuron, avg_intensity_arc, avg_intensity_fos, nucleus_mask] =...
    nuclearExpression(label, arc_mat, fos_mat)
% Quantification of fluorecence intensity within nucleus areas detected

        I_arc = arc_mat;
        I_fos = fos_mat;
        L = label;

        final_no = numel(unique(L))-1;
        nucleus_mask = struct('idx',0,'mask',zeros(1002,1004), 'arc_intensity',...
            zeros(1002,1004), 'fos_intensity',zeros(1002,1004));
        nucleus_mask = repmat(nucleus_mask, final_no,1);
        
        avg_intensity_arc = zeros(final_no,1);
        avg_intensity_fos = zeros(final_no,1);
        
        neuron = struct('Nucleus_Idx', [], 'Centroid',[], 'Arc_Exp',[],'Fos_Exp',[],'spikeTimes',[]);
        neuron = repmat(neuron, final_no, 1);
        center_coord = regionprops(L, 'Centroid');
        
        for inucleus = 1:final_no
            nucleus_mask(inucleus).idx = inucleus;
            nucleus_mask(inucleus).mask = L==inucleus;
            logi_mask = logical(nucleus_mask(inucleus).mask);
            nucleus_mask(inucleus).arc_intensity(logi_mask ==1) = I_arc(logi_mask==1);
            nucleus_mask(inucleus).fos_intensity(logi_mask ==1) = I_fos(logi_mask==1);
            avg_intensity_arc(inucleus) = ...
                mean(nucleus_mask(inucleus).arc_intensity(nucleus_mask(inucleus).arc_intensity~=0));
            avg_intensity_fos(inucleus) = ...
                mean(nucleus_mask(inucleus).fos_intensity(nucleus_mask(inucleus).fos_intensity~=0));
            
            neuron(inucleus).Nucleus_Idx = inucleus;
            neuron(inucleus).Centroid = center_coord(inucleus).Centroid;
            neuron(inucleus).Arc_Exp = avg_intensity_arc(inucleus);
            neuron(inucleus).Fos_Exp = avg_intensity_fos(inucleus);
        end
end

function thresholds = mle_threshold(exp_vals, thresholdspecs)
% Use a mixed normal distribution model to estimate underlying populations
% of positive and negative cells.
% threshold for positive cells is set to be at the 99.7%, ie mu + 3*sigma
% level of the negative populations.
x = exp_vals;
times_sigma = thresholdspecs.sigma;
pdf_normixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
%setup custom function of mixed normal population, to estimate p and q and 
%individual population of mu and sigma
%define initial starting parameters for mle run
pStart = .5;
muStart = quantile(x,[.25 .75]);
sigmaStart = sqrt(var(x) - .25*diff(muStart).^2);
start = [pStart muStart sigmaStart sigmaStart];

lb = [0 -Inf -Inf 0 0];
ub = [1 Inf Inf Inf Inf];

options = thresholdspecs.options;
%statset('mlecustom')
%set max number of iterations to ensure convergence
try
paramEsts = mle(x, 'pdf',pdf_normixture, 'start',start, ...
    'lower',lb, 'upper',ub, 'options',options);
%threshold set to be xxx*sigma of the negative distribution
    thresholds = paramEsts(2)+ times_sigma*paramEsts(4); 
catch ME
    disp(ME.message)
    thresholds = inf; 
end
end
