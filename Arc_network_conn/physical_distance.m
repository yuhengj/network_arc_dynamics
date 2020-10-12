%% Distance measures
% compute phyical distances between each neuronal pair and obtain the
% distance matrix 
% plot relevant figures

%  Created by Yuheng Jiang
%  11 Oct 2020
%% declare paths and variables 
outputpath = '/Users/yuheng/Desktop/data_fnconn';
% neuronpath = [outputpath filesep 'neuron' filesep];
% neuronfiles = dir([neuronpath '*_neuron.mat']);
idxpath = [outputpath filesep 'neuron_id' filesep];
idxfiles = dir([idxpath '*.mat']);
labelpath = [outputpath filesep 'label' filesep];
labelfiles = dir([labelpath '*.mat']);

distpath = [outputpath filesep 'dist' filesep];
if ~exist(distpath, 'dir'); mkdir(distpath); end

% list of neun files
neun_idx = [];

%% compute dist matrix
for ifile = 1:length(labelfiles)
labelfile = [labelpath labelfiles(ifile).name];
load(labelfile)
% neuronfilepath = [neuronpath neuronfiles(ifile).name];
% load(neuronfilepath)
idxfile = [idxpath idxfiles(ifile).name];
load(idxfile)
group_names = fieldnames(neuron_id);

L = label;
center_coord = regionprops(L, 'Centroid');
 
sort_id = zeros(length(center_coord),1);
id_interest = [];
endidx = 0;


for group_id = 1:4
    if group_id ==1
        startidx=1;
    else
        startidx = endidx+1;
    end
    id_interest = getfield(neuron_id, group_names{group_id});
    endidx = endidx + length(id_interest);
    sort_id(startidx:endidx,1) = id_interest;
end

dist_matrix = zeros(length(sort_id));

for idxn = 1:length(sort_id)
    ineuron = sort_id(idxn);
    cross_corr_combinations = zeros(length(sort_id),2);
    temp_dist = zeros(size(cross_corr_combinations,1),1);
     cross_corr_combinations(:,1) = ineuron;
        cross_corr_combinations(:,2) = sort_id;
        
        for j = 1 : size(cross_corr_combinations,1)
            x =  center_coord(cross_corr_combinations(j,1)).Centroid;
                y =  center_coord(cross_corr_combinations(j,2)).Centroid;
                dist = sqrt((x(2)-y(2))^2/(x(1)-y(1))^2);
                temp_dist(j,1) = dist;
        end
        dist_matrix(idxn,:) = temp_dist';
end

save([distpath num2str(ifile,'%02d') '_dist.mat'], 'dist_matrix')

end

%% dist values according to group
idxpath = [outputpath filesep 'idx' filesep];
idxfiles = dir([idxpath '*.mat']);

distpath = [outputpath filesep 'dist' filesep];
distfiles = dir([distpath '*.mat']);

dist_arc = [];
dist_neun = [];
dist_neun_arc = [];

for idx = 1:length(neun_idx)
ifile = neun_idx(idx);

distfile = [distpath distfiles(ifile).name];
load(distfile)

idxfile = [idxpath idxfiles(ifile).name];
load(idxfile)

cm_fix = weight_conversion(dist_matrix, 'autofix');

arc_arc = cm_fix(1:group_id_list(1),1:group_id_list(1));
arc_neun = cm_fix(1:group_id_list(1),(group_id_list(2)+1):group_id_list(3));
neun_neun = cm_fix((group_id_list(2)+1):group_id_list(3),...
    (group_id_list(2)+1):group_id_list(3));

for ia = 1:size(arc_arc)
dist_arc_i = diag(arc_arc,ia);
dist_arc = vertcat(dist_arc,dist_arc_i);
end

for in = 1:size(neun_neun)
dist_neun_i = diag(neun_neun,in);
dist_neun = vertcat(dist_neun,dist_neun_i);
end

for ian = 1:size(neun_neun)
dist_neun_arc_i = diag(arc_neun,ian);
dist_neun_arc = vertcat(dist_neun_arc,dist_neun_arc_i);
end

end
%% Figure 4.5 A

f1 = plot(ones(length(dist_arc),1)+rand(length(dist_arc),1),dist_arc,'kx');
hold on
f2 = plot(ones(length(dist_neun_arc),1)*4+rand(length(dist_neun_arc),1),dist_neun_arc,'kx');
f3 = plot(ones(length(dist_neun),1)*7+rand(length(dist_neun),1),dist_neun,'kx');
hold off

xlim([0 8.5])
ylim([0 3000])

%% Figure 4.5 B
% load values of cnnectivity conn_arc, conn_neun, conn_neun_arc

g1 = plot(conn_arc, dist_arc,'ro');
hold on
g2 = plot(conn_neun, dist_neun,'kx');
g3 = plot(conn_neun_arc, dist_neun_arc,'b*');
hold off

ylim([0 20])
xlim([0.8 1])
set(axes1,'FontSize',16,'LineWidth',2);

% %%
% all_conn = vertcat(conn_arc,conn_neun,conn_neun_arc);
% all_dist = vertcat(dist_arc,dist_neun,dist_neun_arc);
% plot(all_conn, all_dist,'kx');
% set(axes1,'FontSize',16,'LineWidth',2);

%% Figure 4.5 D
totsteps = 1000;
start_thresh = 0.5;
inc = (1-start_thresh)/totsteps;

thresh = zeros(totsteps,1);
percent_potentiate = zeros(totsteps,1);
percent_depress = zeros(totsteps,1);

dist_pot = zeros(totsteps,1);
dist_dep = zeros(totsteps,1);
for ii = 1:totsteps
   dist_pot(ii) = mean(all_dist(all_conn>thresh(ii)&all_conn_ini<0.5));
   dist_dep(ii) = mean(all_dist(all_conn_ini>thresh(ii)&all_conn<0.5));
end

plot(thresh,dist_pot,'-r');
hold on
plot(thresh,dist_dep,'-b');
hold off
