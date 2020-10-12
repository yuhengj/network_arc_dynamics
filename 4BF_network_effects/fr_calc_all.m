%% Plot differences of firing rates 
% fr obtained from fr_calc.m
%  Created by Yuheng Jiang
%  11 Oct 2020

fr_ini = zeros(length(fr)/3,1);
fr_fin = zeros(length(fr)/3,1);
idx = 1;

for ii = 1:length(fr)
   if rem(ii,3) == 1
       fr_ini(idx) = fr(ii);
   elseif rem(ii,3) == 0
       fr_fin(idx) = fr(ii);
       idx = idx+1;
   else
   end
end

fr_diff = fr_fin - fr_ini;
pos_change_idx = find(fr_diff>0);
pos_change_idx = setdiff(pos_change_idx,exclude);
pos_change_arc = arc_pct(pos_change_idx);
mean(pos_change_arc) 

neg_change_idx = find(fr_diff<0);
neg_change_idx = setdiff(neg_change_idx,exclude);
neg_change_arc = arc_pct(neg_change_idx);
mean(neg_change_arc) 

h = figure;

fr_diff = fr_fin - fr_ini;
fr_diff_pos = find(fr_diff>0);

for ii = 1:length(fr_fin)
    if any(ii==fr_diff_pos)==1
    plot([1 2],[fr_ini(ii), fr_fin(ii)],'ro-','LineWidth',1.5);
    hold on
    else
        plot([1 2],[fr_ini(ii), fr_fin(ii)],'bo-','LineWidth',1.5);
    end        
end

xlim([0.75 2.25])
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])
%%
g = figure;
fr_fin = fr_fin(neun_idx);
fr_ini = fr_ini(neun_idx);

fr_diff = fr_fin - fr_ini;
fr_diff_pos = find(fr_diff>0);

for ii = 1:length(fr_fin)
    if any(ii==fr_diff_pos)==1
    plot([1 2],[fr_ini(ii), fr_fin(ii)],'ro-','LineWidth',1.5);
    hold on
    else
        plot([1 2],[fr_ini(ii), fr_fin(ii)],'bo-','LineWidth',1.5);
        hold on
    end        
end

xlim([0.75 2.25])
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])

mean(fr_diff)
std(fr_diff)
std(fr_diff)./sqrt(length(fr_diff))
