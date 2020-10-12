%% cfos analysis 
% determine the connectivity values of each different group 
% plot ksdensity 
% Figure 4.7

%  Created by Yuheng Jiang
%  11 Oct 2020

%% declare paths and set variables
outputpath = '';
idxpath = [outputpath filesep 'idx' filesep];
idxfiles = dir([idxpath '*.mat']);
matrixpath = [outputpath filesep 'matrix' filesep 'before' filesep];
matrixfiles = dir([matrixpath '*_matrix.mat']);


% list of cfos 
cfos_idx = [];
idx_select = cfos_idx; 

%% get connectivity values from each group

conn_one = [];
conn_two = [];
conn_three = [];
conn_four = [];
conn_five = [];
conn_six = [];
conn_seven = [];
conn_eight = [];
conn_nine = [];
conn_ten = [];

%
venn_vals = zeros(length(idx_select),4);
%
results_perculture = zeros(length(idx_select),10);


for idx = 1:length(idx_select)
ifile = idx_select(idx);

matrixfile = [matrixpath matrixfiles(ifile).name];
load(matrixfile)

idxfile = [idxpath idxfiles(ifile).name];
load(idxfile)

venn_vals(idx,1) = group_id_list(1);
venn_vals(idx,2) = group_id_list(2);
venn_vals(idx,3) = group_id_list(3);
venn_vals(idx,4) = group_id_list(4);

cm_fix = weight_conversion(conn_matrix, 'autofix');

grp_one = cm_fix(1:group_id_list(1),1:group_id_list(1));
grp_two = cm_fix(1:group_id_list(1),group_id_list(1)+1:group_id_list(2));
grp_three = cm_fix(1:group_id_list(1),group_id_list(2)+1:group_id_list(3));
grp_four = cm_fix(1:group_id_list(1),group_id_list(3)+1:group_id_list(4));

grp_five = cm_fix(group_id_list(1)+1:group_id_list(2),...
    group_id_list(1)+1:group_id_list(2));
grp_six = cm_fix(group_id_list(1)+1:group_id_list(2),...
    group_id_list(2)+1:group_id_list(3));
grp_seven = cm_fix(group_id_list(1)+1:group_id_list(2),...
    group_id_list(3)+1:group_id_list(4));

grp_eight = cm_fix(group_id_list(2)+1:group_id_list(3),...
    (group_id_list(2)+1):group_id_list(3));
grp_nine = cm_fix(group_id_list(2)+1:group_id_list(3),...
    (group_id_list(3)+1):group_id_list(4));

grp_ten = cm_fix(group_id_list(3)+1:group_id_list(4),...
    (group_id_list(3)+1):group_id_list(4));

for ione = 1:size(grp_one,1)
conn_one_i = diag(grp_one,ione);
conn_one = vertcat(conn_one,conn_one_i);
end
for itwo = 1:size(grp_two,1)
conn_two_i = diag(grp_two,itwo);
conn_two = vertcat(conn_two,conn_two_i);
end
for ithree = 1:size(grp_three,1)
conn_three_i = diag(grp_three,ithree);
conn_three = vertcat(conn_three,conn_three_i);
end
for ifour = 1:size(grp_four,1)
conn_four_i = diag(grp_four,ifour);
conn_four = vertcat(conn_four,conn_four_i);
end

if size(grp_five,1)==0
    conn_five_i = [];
else
for ifive = 1:size(grp_five,1)
conn_five_i = diag(grp_five,ifive);
conn_five = vertcat(conn_five,conn_five_i);
end
end

if size(grp_six,1)==0
    conn_six_i = [];
else
for isix = 1:size(grp_six,1)
conn_six_i = diag(grp_six,isix);
conn_six = vertcat(conn_six,conn_six_i);
end
end

if size(grp_seven,1)==0
    conn_seven_i = [];
else
for iseven = 1:size(grp_seven,1)
conn_seven_i = diag(grp_seven,iseven);
conn_seven = vertcat(conn_seven,conn_seven_i);
end
end

for ieight = 1:size(grp_eight,1)
conn_eight_i = diag(grp_eight,ieight);
conn_eight = vertcat(conn_eight,conn_eight_i);
end
for inine = 1:size(grp_nine,1)
conn_nine_i = diag(grp_nine,inine);
conn_nine = vertcat(conn_nine,conn_nine_i);
end
for iten = 1:size(grp_ten,1)
conn_ten_i = diag(grp_ten,iten);
conn_ten = vertcat(conn_ten,conn_ten_i);
end

results_perculture(idx,1) = mean(conn_one_i);
results_perculture(idx,2) = mean(conn_two_i);
results_perculture(idx,3) = mean(conn_three_i);
results_perculture(idx,4) = mean(conn_four_i);
results_perculture(idx,5) = mean(conn_five_i);
results_perculture(idx,6) = mean(conn_six_i);
results_perculture(idx,7) = mean(conn_seven_i);
results_perculture(idx,8) = mean(conn_eight_i);
results_perculture(idx,9) = mean(conn_nine_i);
results_perculture(idx,10) = mean(conn_ten_i);

end

%% get the connectivity distribution of each group at baseline, during and after 4BF
[f_one,x_one] = ksdensity(conn_four_ini,'Bandwidth',0.05);
[f_two,x_two] = ksdensity(conn_four_dur,'Bandwidth',0.05);
[f_three,x_three] = ksdensity(conn_four_fin,'Bandwidth',0.05);


figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(x_one,f_one, 'LineWidth', 3, 'DisplayName','Baseline',...
    'Color',[0 0 0]);
plot(x_two,f_two, 'LineWidth', 3, 'DisplayName','During 4BF',...
    'Color',[0.5 0.5 0.5]);
plot(x_three,f_three, 'LineWidth', 3, 'DisplayName','After 4BF',...
    'Color',[0.85 0.85 0.85]);


% Set the remaining axes properties
set(axes1,'FontSize',16,'LineWidth',3);
% Create legend
legend(axes1,'show');
xlim([-0.5 1])

%% compare the connectivity distribution of different groups 
[f_one,x_one] = ksdensity(conn_one_ini,'Bandwidth',0.05);
[f_two,x_two] = ksdensity(conn_five_ini,'Bandwidth',0.05);
[f_three,x_three] = ksdensity(conn_eight_ini,'Bandwidth',0.05);
[f_four,x_four] = ksdensity(conn_ten_ini,'Bandwidth',0.05);


figure1 = figure('Color',[1 1 1]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');

plot(x_one,f_one, 'LineWidth', 3, 'DisplayName','Arc+ cFos+');
plot(x_two,f_two, 'LineWidth', 3, 'DisplayName','Arc+ cFos-');
plot(x_three,f_three, 'LineWidth', 3, 'DisplayName','Arc- cFos+');
plot(x_four,f_four, 'LineWidth', 3, 'DisplayName','Arc- cFos-');

% Set the remaining axes properties
set(axes1,'FontSize',16,'LineWidth',3);
% Create legend
legend(axes1,'show');
xlim([-0.5 1])
