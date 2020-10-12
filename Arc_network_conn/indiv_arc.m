%% individual Arc stats
% get matrix of ifile before during and after (matrix mat)
% get idx of arc pos neurons vs all neurons (idx mat)

%  Created by Yuheng Jiang
%  11 Oct 2020
%% Declare directories and function paths
outputpath = '/Users/yuheng/Desktop/data_fnconn';
matrixpath = [outputpath filesep 'matrix' filesep 'final' filesep];
matrixfiles = dir([matrixpath '*_matrix.mat']);
idxpath = [outputpath filesep 'idx' filesep];
idxfiles = dir([idxpath '*.mat']);
labelpath = [outputpath filesep 'label' filesep];
labelfiles = dir([labelpath '*.mat']);

figspath = [outputpath filesep 'figure' filesep];
if ~exist(figspath, 'dir'); mkdir(figspath); end

%requires brain connectivity toolbox
addpath('/Users/yuheng/Google Drive/matlab_scripts/2017_01_15_BrainConnectivityToolbox/')
%requires function export_fig
addpath('/Users/yuheng/Dropbox/matlab_scripts/export_fig')

% list of neun files
neun_idx = [];

%% get connection strengths according to groups
conn_arc = [];
conn_neun = [];
conn_neun_arc = [];

for idx = 1:length(neun_idx)
ifile = neun_idx(idx);

matrixfile = [matrixpath matrixfiles(ifile).name];
load(matrixfile)

idxfile = [idxpath idxfiles(ifile).name];
load(idxfile)

cm_fix = weight_conversion(conn_matrix, 'autofix');

arc_arc = cm_fix(1:group_id_list(1),1:group_id_list(1));
arc_neun = cm_fix(1:group_id_list(1),(group_id_list(2)+1):group_id_list(3));
neun_neun = cm_fix((group_id_list(2)+1):group_id_list(3),...
    (group_id_list(2)+1):group_id_list(3));

for ia = 1:size(arc_arc)
conn_arc_i = diag(arc_arc,ia);
conn_arc = vertcat(conn_arc,conn_arc_i);
end

for in = 1:size(neun_neun)
conn_neun_i = diag(neun_neun,in);
conn_neun = vertcat(conn_neun,conn_neun_i);
end

for ian = 1:size(neun_neun)
conn_neun_arc_i = diag(arc_neun,ian);
conn_neun_arc = vertcat(conn_neun_arc,conn_neun_arc_i);
end

end

%% Figure 4.4 A
h1 = histogram(nonzeros(conn_arc));
hold on

h3 = histogram(nonzeros(conn_neun_arc));
h2 = histogram(nonzeros(conn_neun));

h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h1.EdgeColor = 'none';
h1.EdgeAlpha = 0.4;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
h2.EdgeColor = 'none';
h2.EdgeAlpha = 0.4;
h3.Normalization = 'probability';
h3.BinWidth = 0.01;
h3.EdgeColor = 'none';
h3.EdgeAlpha = 0.4;
hold off

ylim([0 0.002])
xlim([0.8 1])

%% Figure 4.4 B
conn_arc_ini = conn_arc;
conn_neun_ini = conn_neun;
conn_neun_arc_ini = conn_neun_arc;

g1 = scatter(conn_arc_ini, conn_arc,'.','MarkerEdgeColor',[0 0.4470 0.7410],...
    'MarkerEdgeAlpha',1);
hold on
g2 = scatter(conn_neun_ini, conn_neun,'.','MarkerEdgeColor',[0.8500 0.3250 0.0980],...
    'MarkerEdgeAlpha',1);
g3 = scatter(conn_neun_arc_ini, conn_neun_arc,'.','MarkerEdgeColor',[0.9290 0.6940 0.1250],...
    'MarkerEdgeAlpha',1);
plot([-0.4 1],[-0.4 1],':k', 'LineWidth',1.5)
hold off

xlim([-0.4 1])
ylim([0.5 1])

%% group stats
results_arc = zeros(length(neun_idx),5);
results_neun = zeros(length(neun_idx),5);
col1 = [];
col2 = [];
col3 = [];
col4 = [];

for idx = 1:length(neun_idx)
ifile = neun_idx(idx);


matrixfile = [matrixpath matrixfiles(ifile).name];
load(matrixfile)

idxfile = [idxpath idxfiles(ifile).name];
load(idxfile)

cm_fix = weight_conversion(conn_matrix, 'autofix');

arc_arc = mean(cm_fix(1:group_id_list(1),1:group_id_list(1)),2);
arc_neun = mean(cm_fix(1:group_id_list(1),(group_id_list(2)+1):group_id_list(3)),2);

neun_arc = mean(cm_fix((group_id_list(2)+1):group_id_list(3),1:group_id_list(1)),2);
neun_neun = mean(cm_fix((group_id_list(2)+1):group_id_list(3),...
    (group_id_list(2)+1):group_id_list(3)),2);

results_arc(idx,1) = group_id_list(1);
results_arc(idx,2) = mean(arc_arc);
results_arc(idx,3) = std(arc_arc);
results_arc(idx,4) = mean(arc_neun);
results_arc(idx,5) = std(arc_neun);

results_neun(idx,1) = group_id_list(3)-(group_id_list(2)+1);
results_neun(idx,2) = mean(neun_arc);
results_neun(idx,3) = std(neun_arc);
results_neun(idx,4) = mean(neun_neun);
results_neun(idx,5) = std(neun_neun);

col1 = vertcat(col1,arc_arc);
col2 = vertcat(col2,arc_neun);
col3 = vertcat(col3,neun_arc);
col4 = vertcat(col4,neun_neun);
end

%% 
plot(ones(length(results_arc),1),results_arc(:,2), 'ko')
hold on
plot (ones(length(results_neun),1)*2,results_neun(:,4), 'ko')
hold off
xlim([0.5 2.5])

%%
% save('during_arc_neun.mat', 'results_arc', 'results_neun', 'col1','col2',...
%    'col3','col4')

%% Figure 4.3 B 
during_arc = results_arc(:,2);
during_neun = results_neun(:,4);

before_arc = results_arc(:,2);
before_neun = results_neun(:,4);

final_arc = results_arc(:,2);
final_neun = results_neun(:,4);

for ixn = 1:length(before_arc)
   plot([1,2,3],[before_arc(ixn),during_arc(ixn),final_arc(ixn)],'-kx','LineWidth',2,...
    'MarkerSize',10) 
   hold on
end
hold off
xlim([0.5 3.5])

%% Figure 4.3 C
for ixn = 1:length(before_neun)
   plot([1,2,3],[before_neun(ixn),during_neun(ixn),final_neun(ixn)],'-ko','LineWidth',2,...
    'MarkerSize',10) 
   hold on
end
hold off
xlim([0.5 3.5])

%% Figure 4.4 C
% all connections threshold potentiated vs depressed
% make vector of all_conn group
all_conn = vertcat(conn_arc,conn_neun,conn_neun_arc);
all_conn_ini = vertcat(conn_arc_ini,conn_neun_ini,conn_neun_arc_ini);

group_id = vertcat (ones(length(conn_arc),1), ones(length(conn_neun),1)*2, ...
    ones(length(conn_neun_arc),1)*3);

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

plot(thresh,percent_potentiate,'-r');
hold on
plot(thresh,percent_depress,'-b');
hold off

%% Figure 4.4 B
X1 = thresh;
YMatrix1 = horzcat(percent_potentiate, percent_depress);
figure1 = figure('Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',3);
set(plot1(1),'DisplayName','Potentiated connections',...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(plot1(2),'DisplayName','Depressed connections',...
    'Color',[0 0.447058826684952 0.74117648601532]);

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 0.95]);
% Set the remaining axes properties
set(axes1,'FontSize',16,'LineWidth',2);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.17 0.81 0.35 0.094]);

