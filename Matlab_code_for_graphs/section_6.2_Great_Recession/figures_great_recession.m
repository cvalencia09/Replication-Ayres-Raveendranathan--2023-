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
% Column 2: z
% Column 3: output	
% Column 4: hours	
% Column 5: firm debt
% Column 6: firm entry
% Column 7: firm exit
% Column 8: consumption
% Column 9: investment	
% Column 10: theta
% Column 11: labor productivity	
% Column 12: zero

% Productivity shock
x1 = aggregates_joint(:,1);
y1 = aggregates_joint(:,2);
y2 = aggregates_joint(:,10);
figure
p1=plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'), hold on
p2=plot(x1,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--')
plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':')
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-16,4])
hleglines = [p1(1) p2(1)];
legend(hleglines,'productivity','collateral constraint','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_productivity','-dpng','-r300')

% Firm debt
x1 = aggregates_joint(:,1);
x2 = aggregates_data(:,1);
y1 = aggregates_joint(:,5);
y2 = aggregates_data(:,6);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-20,10])
hleglines = [p1(1) p2(1)];
legend(hleglines,'productivity + credit','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_firm_debt','-dpng','-r300')


% Entry
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = firm_entry_exit_great_recession(:,1);
y1 = aggregates_joint(:,6);
y2 = aggregates_theta(:,6);
y3 = aggregates_z(:,6);
y4 = firm_entry_exit_great_recession(:,2);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2007.75,2014])
ylim([-2,2])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_entry','-dpng','-r300')

% Exit
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = firm_entry_exit_great_recession(:,1);
y1 = aggregates_joint(:,7);
y2 = aggregates_theta(:,7);
y3 = aggregates_z(:,7);
y4 = firm_entry_exit_great_recession(:,3);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2007.75,2014])
ylim([-2,2])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_exit','-dpng','-r300')

% Exit by age theta and z
x1 = firm_exit_by_age_joint(:,1);
y1 = firm_exit_by_age_joint(:,2);
y2 = firm_exit_by_age_joint(:,3);
y3 = firm_exit_by_age_joint(:,4);
y4 = firm_exit_by_age_joint(:,5);
y5 = firm_exit_by_age_joint(:,6);
y6 = firm_exit_by_age_joint(:,7);
y7 = firm_exit_by_age_joint(:,8);
y8 = firm_exit_by_age_joint(:,9);
y9 = firm_exit_by_age_joint(:,10);

figure
p1 = plot(x1,y1,'Color',[0 0 0.8],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x1,y2,'Color',[0.2 0.5 0.8],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x1,y3,'Color',[1.0 0.5 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,y4,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,y5,'Color',[0.7 0 0.7],'LineWidth',1.3,'LineStyle','-');
p6 = plot(x1,y6,'Color',[1 1 0],'LineWidth',1.3,'LineStyle','-');
p7 = plot(x1,y7,'Color',[0 0.8 0],'LineWidth',1.3,'LineStyle','-');
p8 = plot(x1,y8,'Color',[0.2 0.2 0.2],'LineWidth',1.3,'LineStyle','-');
p9 = plot(x1,y9,'Color',[0.6 0.6 0.6],'LineWidth',1.3,'LineStyle','-');
p10 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle','-');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([2007.75,2010.75])
xticks([2008 2009 2010 2011])
ylim([-3,9])
hleglines = [p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1) p8(1) p9(1)];
legend(hleglines,'1','2','3','4','5','6-10','11-15','16-20','20 plus','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',6.,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_exit_by_age_theta_z','-dpng','-r300')
%,'Orientation','horizontal'

% GDP
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = aggregates_data(:,1);
y1 = aggregates_joint(:,3);
y2 = aggregates_theta(:,3);
y3 = aggregates_z(:,3);
y4 = aggregates_data(:,2);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-8,2])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_gdp','-dpng','-r300')

% Hours
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = aggregates_data(:,1);
y1 = aggregates_joint(:,4);
y2 = aggregates_theta(:,4);
y3 = aggregates_z(:,4);
y4 = aggregates_data(:,5);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-8,2])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_hours','-dpng','-r300')

% Consumption
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = aggregates_data(:,1);
y1 = aggregates_joint(:,8);
y2 = aggregates_theta(:,8);
y3 = aggregates_z(:,8);
y4 = aggregates_data(:,3);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-8,2])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_consumption','-dpng','-r300')

% Investment
x1 = aggregates_joint(:,1);
x2 = aggregates_theta(:,1);
x3 = aggregates_z(:,1);
x4 = aggregates_data(:,1);
y1 = aggregates_joint(:,9);
y2 = aggregates_theta(:,9);
y3 = aggregates_z(:,9);
y4 = aggregates_data(:,4);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p4 = plot(x4,y4,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-30,10])
hleglines = [p1(1) p2(1) p3(1) p4(1)];
legend(hleglines,'productivity + credit','credit','productivity','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_investment','-dpng','-r300')

% Decomposition of aggregates

% Decomposition output
x1 = great_recession_decomposition(:,1);
x2 = great_recession_decomposition(:,1);
x3 = great_recession_decomposition(:,1);
x4 = great_recession_decomposition(:,1);
y1 = great_recession_decomposition(:,2);
y2 = great_recession_decomposition(:,6);
y3 = great_recession_decomposition(:,10);
y4 = great_recession_decomposition(:,14);

figure
p1 = plot(x1,y1,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle',':'); hold on
p2 = plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p4 = plot(x4,y4,'Color',[0 0.6 0],'LineWidth',1.3,'LineStyle','-.');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('quarter','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2016.75])
xticks([2008 2010 2012 2014 2016])
ylim([-4,1])
hleglines = [p3(1) p1(1) p4(1) p2(1)];
legend(hleglines,'no entry','no exit','no entry or exit','baseline','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_decomposition_output','-dpng','-r300')

% Decomposition hours
x1 = great_recession_decomposition(:,1);
x2 = great_recession_decomposition(:,1);
x3 = great_recession_decomposition(:,1);
x4 = great_recession_decomposition(:,1);
y1 = great_recession_decomposition(:,3);
y2 = great_recession_decomposition(:,7);
y3 = great_recession_decomposition(:,11);
y4 = great_recession_decomposition(:,15);

figure
p1 = plot(x1,y1,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle',':'); hold on
p2 = plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p4 = plot(x4,y4,'Color',[0 0.6 0],'LineWidth',1.3,'LineStyle','-.');
p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('quarter','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2016.75])
xticks([2008 2010 2012 2014 2016])
ylim([-4,1])
hleglines = [p3(1) p1(1) p4(1) p2(1)];
legend(hleglines,'no entry','no exit','no entry or exit','baseline','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_decomposition_hours','-dpng','-r300')


% Labor productivity
x1 = aggregates_joint(:,1);
x2 = aggregates_z(:,1);
x3 = aggregates_data(:,1);
y1 = aggregates_joint(:,11);
y2 = aggregates_z(:,11);
y3 = aggregates_data(:,7);

figure
p1 = plot(x1,y1,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-.');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([2007.75,2017])
ylim([-10,10])
hleglines = [p1(1) p3(1) p2(1)];
legend(hleglines,'productivity + credit','data','productivity','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_labor_productivity','-dpng','-r300')


