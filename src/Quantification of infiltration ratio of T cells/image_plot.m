%configuration
clear;
clear all;
RawImage_loadpath = './image_analysis';
Tiff_length = 2048;
Tiff_width = 2048;
z_serial = 31;
well_number = 21;
filter = [-2 -1 0; -1 0 1; 0 1 2 ];



%% %plot time-curves
time = [0, 3, 6];
date = '20210516-2';
wells = [61:81];
replicate_num = 3;
count = zeros(length(wells),1);
area = zeros(length(wells),1);
groups_num = length(wells)/replicate_num;
infiltration = zeros(2*groups_num, replicate_num*length(time));
infiltration_mean = zeros(2*groups_num, length(time)); %Dark and Light
infiltration_std = zeros(2*groups_num, length(time));
infiltration_mean_new = zeros(2*groups_num, length(time));
infiltration_std_new = zeros(2*groups_num, length(time));
save_path = './3D_figure';
        
for T = 1:length(time)
    for i = 1:2
        if i ~= 1
            condition = 'l';
        else
            condition = 'd';
        end
        if strcmp(condition, 'd')
            count_name = sprintf('%s/ZZ_3D_tumor_0516_Jurkat_%dh_dark_count', date, time(T));
            area_name = sprintf('%s/ZZ_3D_tumor_0516_Jurkat_%dh_dark_area', date, time(T));
        else
            count_name = sprintf('%s/ZZ_3D_tumor_0516_Jurkat_%dh_count', date, time(T));
            area_name = sprintf('%s/ZZ_3D_tumor_0516_Jurkat_%dh_area', date, time(T));   
        end

        m = 1;

        for ii = wells
            csv_c = readtable(sprintf('%s/results/%s_%d.csv', RawImage_loadpath, count_name, ii));
            count(m) = table2array(csv_c(1, 4));
            csv_a = readtable(sprintf('%s/results/%s_%d.csv', RawImage_loadpath, area_name, ii));
            area(m) = table2array(csv_a(1, 2));
            m = m+1;
        end

    %infiltration ratio: mean(number of infiltrated T cell/ area of target site)
        for iii = 1: groups_num
            infiltration(groups_num*(i-1)+iii, [1:replicate_num]+replicate_num*(T-1)) = (count([1:replicate_num]+replicate_num*(iii-1))./area([1:replicate_num]+replicate_num*(iii-1)));
            %infiltration(groups_num*(i-1)+iii, [1:replicate_num]+replicate_num*(T-1)) = (count([1:replicate_num]+replicate_num*(iii-1)));
        end
    end

end

initial_val = zeros(2*groups_num,1);
for i = 1:2*groups_num
    initial_val(i) = mean(infiltration(i,1:replicate_num));
end

infiltration_new = infiltration./initial_val; %fold-change

for i = 1:2*groups_num
    for ii = 1:length(time)
    infiltration_mean_new(i,ii) = mean(infiltration_new(i,(ii-1)*replicate_num+[1:replicate_num]));
    infiltration_std_new(i,ii) = std(infiltration_new(i,(ii-1)*replicate_num+[1:replicate_num]));
    end
end


infiltration_mean_new_2 = infiltration_mean_new([7, 5, 3, 11], :);
infiltration_std_new_2 = infiltration_std_new([7, 5, 3, 11], :);

L = zeros(4, length(time));

figure(1);

%bar graph
c = 1:4;
bar(c,infiltration_mean_new_2);
hold on; %error bars for three time points
range_1 = c;
errorbar(range_1, infiltration_mean_new_2(:,2), L(:,2), infiltration_std_new_2(:,2), 'k', 'Linestyle', 'None');  %x坐标[1 2 3]，对应中间一组均值和标准差，此为默认值
hold on;
range_2 = 0.225+range_1;
errorbar(range_2, infiltration_mean_new_2(:,3), L(:,3), infiltration_std_new_2(:,3), 'k', 'Linestyle', 'None'); 
% hold on;
% range_3 = -0.225+range_1;
% errorbar(range_3, infiltration_mean_new_2(:,1), L(:,1), infiltration_std_new_2(:,1), 'k', 'Linestyle', 'None');  

legend(sprintf('%d hour', time(1)),sprintf('%d hour', time(2)),sprintf('%d hour', time(3)),'Location','northeastoutside');
set(gca,'XTick',[1:1:4]);
% set(gca,'XTickLabel',{'P12+STO(T-HA)','P12+STO(T-HA) dark','P12+STO(T-HA/2.5%FBS)','P12+STO(T-HA/2.5%FBS) dark','P12+STO(T-HA)','P12+STO(T-HA) dark',...
%     'P12+STO(T-MC)','P12+STO(T-MC) dark','STO(T-HA)','STO(T-HA) dark','STO(T-MC)','STO(T-MC) dark','P815+STO(T-HA)','P815+STO(T-HA) dark','P12(T-HA)','P12(T-HA) dark'})

% set(gca,'XTickLabel',{'P815+T-HA','P815+T-HA dark','P12+T-MC','P12+T-MC dark','P12+T-HA','P12+T-HA dark'})

% set(gca,'XTickLabel',{'P815-CXCL10+T','P815-CXCL10+T-CXCR3','P815-CXCL10+T-CXCR3 dark'})

set(gca, 'XTickLabel',{'P815+STO(T-CXCR3)', 'STO(T-CXCR3)', 'P815-CXCL10+STO(T-CXCR3)', 'P815-CXCL10+STO(T-CXCR3) dark'});

set(gca, 'XTickLabelRotation', 45,'fontsize',12);
%xlabel('group number');
ylabel('Infiltration ratio (fold-change)');
title('3D tumor spheroid model');
%ylim([0 10*10^(-3)])
saveas(gcf, sprintf('%s/ZZ_3D_tumor_Jurkat_Summary_%s_foldchange_3.png',save_path, date));


