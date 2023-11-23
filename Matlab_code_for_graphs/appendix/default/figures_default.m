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

% Exit default
x1 = impulse_baseline(:,1);
x2 = impulse_default1(:,1);
x3 = impulse_default2(:,1);
y1 = impulse_baseline(:,7);
y2 = impulse_default1(:,7);
y3 = impulse_default2(:,7);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle',':');
p3 = plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p4 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-2,4])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'baseline','default v1','default v2','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit_default','-dpng','-r300')


% Spread default
x1 = impulse_default1(:,1);
x2 = impulse_default2(:,1);
y1 = impulse_default1(:,13);
y2 = impulse_default2(:,13);

figure
p1 = plot(x1,y1,'Color',[0.5 0.5 0.5],'LineWidth',1.3,'LineStyle',':'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([0,6])
hleglines = [p1(1) p2(1)];
legend(hleglines,'default v1','default v2','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_spread_default','-dpng','-r300')




