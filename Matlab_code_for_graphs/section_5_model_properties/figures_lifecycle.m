clear, clc, close all

% save('data_draft_figures.mat','data')
load ('data_figures.mat')

%life_cycle_baseline
% Column 1: age
% Column 2: exit_rate	
% Column 3: employment_size	
% Column 4: productivity	
% Column 5: capital	
% Column 6: net_debt	
% Column 7: net_debt_to_capital	
% Column 8: dividends	
% Column 9: dividends_operating	
% Column 10: debt_distribution	
% Column 11: binding	
% Column 12: average_debt_to_average_capital

%life_cycle_data
% Column 1: Age	
% Column 2: Employment size	
% Column 3: Firm exit rate


%% Lificycle properties

% specification
fontsz = 8;
paper_unit = 'inches'; % 'inches', 'normalized', 'centimeters', 'points'
paper_size = [3.2 1.55]; % [width height]
paper_position = [0 0 3.2 1.55]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'
color_graph = [0.2 0.5 1.0];
color_mark = [0.2 0.5 1.0];

% Employment by age
figure
p1 = plot(lifecycle_baseline(:,1),lifecycle_baseline(:,3),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark); hold on
p2 = plot(lifecycle_data(:,1),lifecycle_data(:,2),'LineWidth',1,'Color',[0 0 0],'Marker','o','MarkerSize',3,'markerfacecolor',[0 0 0]);
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('index, age 0 = 1','interpreter','latex','Color','k')
xlim([0,30])
ylim([0,3])
hleglines = [p1(1) p2(1)];
legend(hleglines,'model','data','Location','SouthEast')
legend boxoff
%set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_emp_age','-dpng','-r300')

% Capital by age
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,5),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark)
xlim([0,30])
ylim([0,200])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('level','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_capital_age','-dpng','-r300')

% Net debt by age
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,6),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark), hold on
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,6),'LineWidth',0.5,'Color',[0 0 0],'LineStyle',':')
xlim([0,30])
ylim([-100,50])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('level','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_debt_age','-dpng','-r300')

% Net debt-to-capital ratio by age
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,7),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark), hold on
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,7),'LineWidth',0.5,'Color',[0 0 0],'LineStyle',':')
xlim([0,30])
ylim([-.5,.4])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('ratio','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_debt_capital_age','-dpng','-r300')

% Debt share
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,10),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark)
xlim([0,30])
ylim([0,6])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('percent of total debt','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_debt_share_age','-dpng','-r300')

% Exit rate by age
figure
p1 = plot(lifecycle_baseline(2:100,1),lifecycle_baseline(2:100,2),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark); hold on
p2 = plot(lifecycle_data(2:5,1),lifecycle_data(2:5,3),'LineWidth',1,'Color',[0 0 0],'Marker','o','MarkerSize',3,'markerfacecolor',[0 0 0]);
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('percent','interpreter','latex','Color','k')
xlim([0,30])
ylim([0,40])
hleglines = [p1(1) p2(1)];
legend(hleglines,'model','data','Location','NorthEast')
legend boxoff
%set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_exit_rate_age','-dpng','-r300')

% Dividend by age
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,9),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark)
xlim([0,30])
ylim([0,6])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('level','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_dividend_age','-dpng','-r300')

% Perfect of firms for which non-negative dividend constraint is binding
figure
plot(lifecycle_baseline(:,1),lifecycle_baseline(:,11),'LineWidth',1,'Color',color_graph,'Marker','o','MarkerSize',3,'markerfacecolor',color_mark)
xlim([0,30])
ylim([70,100])
xlabel('age (years)','interpreter','latex','Color','k')
ylabel('percent of firms','interpreter','latex','Color','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_binding_age','-dpng','-r300')


%% Lificycle properties (non-targeted moments)

%life_cycle_baseline_non_targeted
% Column 1: Firm share data
% Column 2: Firm share model
% Column 3: Firm size data
% Column 4: Firm size model
% Column 5: Emp share data
% Column 6: Empshare model	
% Column 7: Job creation data
% Column 8: Job creation model
% Column 9: Job destruction data
% Column 10: Job destruction model	

% specification
fontsz = 8;
paper_unit = 'inches'; % 'inches', 'normalized', 'centimeters', 'points'
%paper_size = [5.5 1.6]; % [width height]
%paper_position = [0 0 5.5 1.6]; % [left bottom width height].
%paper_size = [3.2 2.0]; % [width height]
%paper_position = [0 0 3.2 2.0]; % [left bottom width height].
paper_size = [3.3 2.1]; % [width height]
paper_position = [0 0 3.3 2.1]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'

% Firm age distribution
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,2) lifecycle_baseline_non_targeted(:,1)]);
set(hb(1),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,30])
xlabel('age group')
ylabel('percentage of firms')
legend('Model','Data','Location','Northwest')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_age_dist','-dpng','-r300')

% firm size by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,4) lifecycle_baseline_non_targeted(:,3)]);
set(hb(1),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,4])
xlabel('age group')
ylabel('avg. size age group/avg. size all firms')
legend('Model','Data','Location','Northwest')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_size_age','-dpng','-r300')

% employment share by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,6) lifecycle_baseline_non_targeted(:,5)]);
set(hb(1),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,70])
xlabel('age group')
ylabel('percent of total employment')
legend('Model','Data','Location','Northwest')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_emp_share_age','-dpng','-r300')

% job creation by age group
labels = {'1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(2:10,8) lifecycle_baseline_non_targeted(2:10,7)]);
set(hb(1),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
xlabel('age group')
ylabel('percentage of employment')
ylim([0,100])
legend('Model','Data','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_job_creation_age','-dpng','-r300')

% job destruction by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,10) lifecycle_baseline_non_targeted(:,9)]);
set(hb(1),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
xlabel('age group')
ylabel('percentage of employment')
ylim([0,100])
legend('Model','Data','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_job_destruction_age','-dpng','-r300')
