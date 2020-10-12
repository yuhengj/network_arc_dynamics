%% plot STTC graph
% sttc obtained from sttc_conn_matrix.m
%  Created by Yuheng Jiang
%  11 Oct 2020

sttc_ini = zeros(length(sttc)/3,1);
sttc_fin = zeros(length(sttc)/3,1);
idx = 1;

for ii = 1:length(sttc)
   if rem(ii,3) == 1
       sttc_ini(idx) = sttc(ii);
   elseif rem(ii,3) == 0
       sttc_fin(idx) = sttc(ii);
       idx = idx+1;
   else
   end
end

sttc_fin = sttc_fin(neun_idx);
sttc_ini = sttc_ini(neun_idx);

sttc_diff = sttc_fin - sttc_ini;
sttc_diff_pos = find(sttc_diff>0);

for idxn = 1:length(sttc_fin)
    if any(idxn==sttc_diff_pos)==1
    plot([1 2],[sttc_ini(idxn), sttc_fin(idxn)],'ro-','LineWidth',1.5);
    hold on
    else
        plot([1 2],[sttc_ini(idxn), sttc_fin(idxn)],'bo-','LineWidth',1.5);
        hold on
    end        
end

xlim([0.75 2.25])
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])

mean(sttc_diff)
std(sttc_diff)
std(sttc_diff)/sqrt(length(sttc_diff));

T =table (fr_diff, sttc_fin, arc_pct_neun);
mdl = stepwiselm(T);
plotAdded(mdl)


