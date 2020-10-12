%% Plot bursts 
% bursts = calcium wave peaks obtained from df_f_peaks.m
%  Created by Yuheng Jiang
%  11 Oct 2020

burst_ini = zeros(length(burst)/3,1);
burst_fin = zeros(length(burst)/3,1);
idx = 1;

for ii = 1:length(burst)
   if rem(ii,3) == 1
       burst_ini(idx) = burst(ii);
   elseif rem(ii,3) == 0
       burst_fin(idx) = burst(ii);
       idx = idx+1;
   else
   end
end

burst_fin = burst_fin(neun_idx);
burst_ini = burst_ini(neun_idx);

burst_diff = burst_fin - burst_ini;
burst_diff_pos = find(burst_diff>0);

for idxn = 1:length(burst_fin)
    if any(idxn==burst_diff_pos)==1
    plot([1 2],[burst_ini(idxn), burst_fin(idxn)],'ro-','LineWidth',1.5);
    hold on
    else
        plot([1 2],[burst_ini(idxn), burst_fin(idxn)],'bo-','LineWidth',1.5);
        hold on
    end        
end

xlim([0.75 2.25])
set(gcf, 'Position',[100 100 600 500])
set(gcf, 'Color', 'w')
set(gca, 'Position',[0.1 0.08 0.7 0.84])

mean(burst_diff)
std(burst_diff)
std(burst_diff)./sqrt(length(burst_diff))
