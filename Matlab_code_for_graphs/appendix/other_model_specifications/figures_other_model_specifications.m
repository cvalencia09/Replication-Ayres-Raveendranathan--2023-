clear, clc, close all

load ('data_figures.mat')

%% Impulse response

% specification
fontsz = 8;
paper_unit = 'inches'; % 'inches', 'normalized', 'centimeters', 'points'
paper_size = [3.2 2.3]; % [width height]
paper_position = [0 0 3.2 2.3]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'
color_graph = [0.2 0.5 1.0];
color_mark = [0.2 0.5 1.0];

% Exit GHH
x1 = impulse_GHH(:,1);
x2 = impulse_GHH(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_GHH(:,7);
y2 = impulse_GHH(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_GHH','-dpng','-r300')

% Exit GHH habit
x1 = impulse_GHH_habit(:,1);
x2 = impulse_GHH_habit(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_GHH_habit(:,7);
y2 = impulse_GHH_habit(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_GHH_habit','-dpng','-r300')



% Exit operating cost in final good
x1 = impulse_operating_final_good(:,1);
x2 = impulse_operating_final_good(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_operating_final_good(:,7);
y2 = impulse_operating_final_good(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_operating_cost_in_final_good','-dpng','-r300')


% Exit large spread
x1 = impulse_large_spread(:,1);
x2 = impulse_large_spread(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_large_spread(:,7);
y2 = impulse_large_spread(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_larger_spread','-dpng','-r300')


% Exit random cost 1
x1 = impulse_random_cost1(:,1);
x2 = impulse_random_cost1(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_random_cost1(:,7);
y2 = impulse_random_cost1(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_random_cost1','-dpng','-r300')


% Exit random cost 2
x1 = impulse_random_cost2(:,1);
x2 = impulse_random_cost2(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_random_cost2(:,7);
y2 = impulse_random_cost2(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_random_cost2','-dpng','-r300')

