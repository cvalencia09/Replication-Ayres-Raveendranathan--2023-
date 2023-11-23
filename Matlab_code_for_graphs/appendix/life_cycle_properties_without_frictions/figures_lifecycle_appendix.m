clear, clc, close all

% save('data_draft_figures.mat','data')
load ('data_figures.mat')

%% Lificycle properties (nontargeted moments)

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
%paper_size = [5.5 2.0]; % [width height]
%paper_position = [0 0 5.5 2.0]; % [left bottom width height].
%paper_size = [3.2 2.0]; % [width height]
%paper_position = [0 0 3.2 2.0]; % [left bottom width height].
paper_size = [3.3 2.1]; % [width height]
paper_position = [0 0 3.3 2.1]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'

% Firm age distribution
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,2) lifecycle_baseline_non_targeted(:,1) lifecycle_no_financial_friction_non_targeted(:,2) lifecycle_no_quadratic_non_targeted(:,2) lifecycle_no_irreversible_non_targeted(:,2)]);
set(hb(1),'FaceColor',[0.2 0.5 1.0],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(3),'FaceColor',[0.8 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(4),'FaceColor',[0 0.6 0],'BarWidth',1,'EdgeColor','none')
set(hb(5),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,30])
xlabel('age group')
ylabel('percentage of firms')
legend('baseline','data','no financial friction','no quadratic','no irreversible','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('appendix\fig_age_dist','-dpng','-r300')

% firm size by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,4) lifecycle_baseline_non_targeted(:,3) lifecycle_no_financial_friction_non_targeted(:,4) lifecycle_no_quadratic_non_targeted(:,4) lifecycle_no_irreversible_non_targeted(:,4)]);
set(hb(1),'FaceColor',[0.2 0.5 1.0],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(3),'FaceColor',[0.8 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(4),'FaceColor',[0 0.6 0],'BarWidth',1,'EdgeColor','none')
set(hb(5),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,30])
xlabel('age group')
ylabel('avg. size age group/avg. size all firms')
legend('baseline','data','no financial friction','no quadratic','no irreversible','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
ylim([0,4])
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('appendix\fig_size_age','-dpng','-r300')

% employment share by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,6) lifecycle_baseline_non_targeted(:,5) lifecycle_no_financial_friction_non_targeted(:,6) lifecycle_no_quadratic_non_targeted(:,6) lifecycle_no_irreversible_non_targeted(:,6)]);
set(hb(1),'FaceColor',[0.2 0.5 1.0],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(3),'FaceColor',[0.8 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(4),'FaceColor',[0 0.6 0],'BarWidth',1,'EdgeColor','none')
set(hb(5),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
ylim([0,70])
xlabel('age group')
ylabel('percent of total employment')
legend('baseline','data','no financial friction','no quadratic','no irreversible','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
ylim([0,50])
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('appendix\fig_emp_share_age','-dpng','-r300')

% job creation by age group
labels = {'1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(2:10,8) lifecycle_baseline_non_targeted(2:10,7) lifecycle_no_financial_friction_non_targeted(2:10,8) lifecycle_no_quadratic_non_targeted(2:10,8) lifecycle_no_irreversible_non_targeted(2:10,8)]);
set(hb(1),'FaceColor',[0.2 0.5 1.0],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(3),'FaceColor',[0.8 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(4),'FaceColor',[0 0.6 0],'BarWidth',1,'EdgeColor','none')
set(hb(5),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
xlabel('age group')
ylabel('percentage of employment')
ylim([0,100])
legend('baseline','data','no financial friction','no quadratic','no irreversible','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('appendix\fig_job_creation_age','-dpng','-r300')

% job destruction by age group
labels = {'0','1','2','3','4','5','6-10','11-15','16-20','21+'};
figure
hb = bar([lifecycle_baseline_non_targeted(:,10) lifecycle_baseline_non_targeted(:,9) lifecycle_no_financial_friction_non_targeted(:,10) lifecycle_no_quadratic_non_targeted(:,10) lifecycle_no_irreversible_non_targeted(:,10)]);
set(hb(1),'FaceColor',[0.2 0.5 1.0],'BarWidth',1,'EdgeColor','none')
set(hb(2),'FaceColor',[0 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(3),'FaceColor',[0.8 0 0],'BarWidth',1,'EdgeColor','none')
set(hb(4),'FaceColor',[0 0.6 0],'BarWidth',1,'EdgeColor','none')
set(hb(5),'FaceColor',[0.5 0.5 0.5],'BarWidth',1,'EdgeColor','none')
set(gca,'XTickLabels', labels);
xlabel('age group')
ylabel('percentage of employment')
ylim([0,100])
legend('baseline','data','no financial friction','no quadratic','no irreversible','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('appendix\fig_job_destruction_age','-dpng','-r300')



