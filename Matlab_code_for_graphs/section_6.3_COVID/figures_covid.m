clear, clc, close all

load ('data_figures.mat')

% specification
fontsz = 8;
paper_unit = 'inches'; % 'inches', 'normalized', 'centimeters', 'points'
paper_size = [3.2 2.3]; % [width height]
paper_position = [0 0 3.2 2.3]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'
color_graph = [0.2 0.5 1.0];
color_mark = [0.2 0.5 1.0];

%impulse response
% Column 1: period	
% Column 2: shutdown firms
% Column 3: output	
% Column 4: hours	
% Column 5: firm debt
% Column 6: firm entry
% Column 7: firm exit
% Column 8: consumption
% Column 9: investment	
% Column 10: shock psi
% Column 11: labor productivity	
% Column 12: zero

% Firm shutdown shock and disutility shock
x1 = aggregates_covid_joint(:,1);
y1 = aggregates_covid_joint(:,2);
y2 = aggregates_covid_joint(:,10);
figure
p1=plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'), hold on
p2 = plot(x1,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'),
plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':')
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2020,2024])
xticks([2020 2021 2022 2023 2024])
ylim([0,20])
hleglines = [p1(1) p2(1)];
legend(hleglines,'firm shutdown','labor disutility','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_covid_firm_shutdown','-dpng','-r300')

% GDP
x1 = aggregates_covid_joint(:,1);
x2 = aggregates_data(:,1);
x3 = aggregates_covid_joint(:,1);
x4 = aggregates_data(:,1);
y1 = aggregates_covid_joint(:,3);
y2 = aggregates_data(:,2);
y3 = aggregates_covid_joint(:,4);
y4 = aggregates_data(:,3);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x1,y3,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle','--');
p4 = plot(x4,y4,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');

%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2020,2024])
ylim([-16,0])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'GDP model','GDP data','Hours model','Hours data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_covid_gdp','-dpng','-r300')

% Entry
x1 = aggregates_covid_joint(:,1);
x2 = aggregates_covid_no_operation_only(:,1);
x3 = aggregates_covid_psi_only(:,1);
y1 = aggregates_covid_joint(:,6);
y2 = aggregates_covid_no_operation_only(:,6);
y3 = aggregates_covid_psi_only(:,6);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2020,2022])
xticks ( [ 2020 2021 2022] )
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1) ];
legend(hleglines,'operation + labor disutility','operation','labor disutility','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_covid_entry','-dpng','-r300')

% Exit
x1 = aggregates_covid_joint(:,1);
x2 = aggregates_covid_no_operation_only(:,1);
x3 = aggregates_covid_psi_only(:,1);
y1 = aggregates_covid_joint(:,7);
y2 = aggregates_covid_no_operation_only(:,7);
y3 = aggregates_covid_psi_only(:,7);


figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2020,2022])
xticks ( [ 2020 2021 2022] )
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'operation + labor disutility','operation','labor disutility','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_covid_exit','-dpng','-r300')





