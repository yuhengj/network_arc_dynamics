%% Arc positive percentage
% calculate Arc positive percentags of each culture
%  Created by Yuheng Jiang
%  11 Oct 2020

outputpath = '';
idxpath = [outputpath filesep 'idx' filesep];
matrixpath = [outputpath filesep 'matrix' filesep];
matrixfiles = dir([matrixpath '*_neuronall_matrix.mat']);

group_no = zeros(length(matrixfiles),4);

%%
for ifile = 1:length(matrixfiles)
saveWellFrame_all = matrixfiles(ifile).name;
saveWellFrame = saveWellFrame_all(1:17);
fprintf([saveWellFrame '\n'])
load([idxpath saveWellFrame '_idxarcfos.mat'])
group_names = fieldnames(neuron_id);

for group_id = 1:4
    list = getfield(neuron_id, group_names{group_id});
    if group_id ==1
        group_no(ifile,1) = length(list);
    elseif group_id ==2
        group_no(ifile,2) = length(list);
    elseif group_id ==3
        group_no(ifile,3) = length(list);
    else
        group_no(ifile,4) = length(list);
    end
    
end
end