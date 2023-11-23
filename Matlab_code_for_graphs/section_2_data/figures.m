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

% recessions (NBER recession dates)
% Jan-80/Jul-80, Jul-81/Nov-82, Jul-90/Mar-91, Mar-01/Nov-01, Dec-07/Jun-09
rec1_t0 = 1980;
rec1_T = 1980.75; %rec1_T = 1980+6/12;
rec2_t0 = 1981.5; %rec2_t0 = 1981+6/12;
rec2_T = 1983; %rec2_T = 1982+10/12;
rec3_t0 = 1990.5; %rec3_t0 = 1990+6/12;
rec3_T = 1991.25; %rec3_T = 1991+2/12;
rec4_t0 = 2001; %rec4_t0 = 2001+2/12;
rec4_T = 2002; %rec4_T = 2001+10/12;
rec5_t0 = 2007.75; %rec5_t0 = 2007+11/12;
rec5_T = 2009.5; %rec5_T = 2009+5/12;

% Entry and exit time series data: trend
x1 = firm_entry_exit_time_series_data(:,1);
x2 = firm_entry_exit_time_series_data(:,1);
x3 = firm_entry_exit_time_series_data(:,1);
x4 = firm_entry_exit_time_series_data(:,1);
y1 = firm_entry_exit_time_series_data(:,2);
y2 = firm_entry_exit_time_series_data(:,3);
y3 = firm_entry_exit_time_series_data(:,4);
y4 = firm_entry_exit_time_series_data(:,5);

figure
box on;
ymax1_plot = 7; ymin1_plot = 15;
v1 = [rec1_t0 ymin1_plot; rec1_t0 ymax1_plot; rec1_T ymax1_plot; rec1_T ymin1_plot];
f1 = [1 2 3 4];
p1 = patch('Faces',f1,'Vertices',v1,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]); hold on
v2 = [rec2_t0 ymin1_plot; rec2_t0 ymax1_plot; rec2_T ymax1_plot; rec2_T ymin1_plot];
f2 = [1 2 3 4];
p2 = patch('Faces',f2,'Vertices',v2,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v3 = [rec3_t0 ymin1_plot; rec3_t0 ymax1_plot; rec3_T ymax1_plot; rec3_T ymin1_plot];
f3 = [1 2 3 4];
p3 = patch('Faces',f3,'Vertices',v3,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v4 = [rec4_t0 ymin1_plot; rec4_t0 ymax1_plot; rec4_T ymax1_plot; rec4_T ymin1_plot];
f4 = [1 2 3 4];
p4 = patch('Faces',f4,'Vertices',v4,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v5 = [rec5_t0 ymin1_plot; rec5_t0 ymax1_plot; rec5_T ymax1_plot; rec5_T ymin1_plot];
f5 = [1 2 3 4];
p5 = patch('Faces',f5,'Vertices',v5,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
p6 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-');
p7 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-');
p8 = plot(x3,y3,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle',':');
p9 = plot(x4,y4,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle',':');
%p5 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percent','Color','k')
xlim([1978,2018])
ylim([7,15])
hleglines = [p6(1) p7(1)];
legend(hleglines,'entry','exit','Location','North')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_entry_exit_time_series_data_trend','-dpng','-r300')

% Entry and exit time series data: cycle
x1 = firm_entry_exit_time_series_data(:,1);
x2 = firm_entry_exit_time_series_data(:,1);
y1 = firm_entry_exit_time_series_data(:,6);
y2 = firm_entry_exit_time_series_data(:,7);

figure
box on;
ymax1_plot = 2; ymin1_plot = -2;
v1 = [rec1_t0 ymin1_plot; rec1_t0 ymax1_plot; rec1_T ymax1_plot; rec1_T ymin1_plot];
f1 = [1 2 3 4];
p1 = patch('Faces',f1,'Vertices',v1,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]); hold on
v2 = [rec2_t0 ymin1_plot; rec2_t0 ymax1_plot; rec2_T ymax1_plot; rec2_T ymin1_plot];
f2 = [1 2 3 4];
p2 = patch('Faces',f2,'Vertices',v2,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v3 = [rec3_t0 ymin1_plot; rec3_t0 ymax1_plot; rec3_T ymax1_plot; rec3_T ymin1_plot];
f3 = [1 2 3 4];
p3 = patch('Faces',f3,'Vertices',v3,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v4 = [rec4_t0 ymin1_plot; rec4_t0 ymax1_plot; rec4_T ymax1_plot; rec4_T ymin1_plot];
f4 = [1 2 3 4];
p4 = patch('Faces',f4,'Vertices',v4,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
v5 = [rec5_t0 ymin1_plot; rec5_t0 ymax1_plot; rec5_T ymax1_plot; rec5_T ymin1_plot];
f5 = [1 2 3 4];
p5 = patch('Faces',f5,'Vertices',v5,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.6 0.6 0.6]);
p6 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-');
p7 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-');
p8 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([1978,2018])
ylim([-2,2])
hleglines = [p6(1) p7(1)];
legend(hleglines,'entry','exit','Location','South')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_entry_exit_time_series_data_cycle','-dpng','-r300')

% Exit by age data
x1 = firm_exit_by_age_data(:,1);
y1 = firm_exit_by_age_data(:,2);
y2 = firm_exit_by_age_data(:,3);
y3 = firm_exit_by_age_data(:,4);
y4 = firm_exit_by_age_data(:,5);
y5 = firm_exit_by_age_data(:,6);
y6 = firm_exit_by_age_data(:,7);
y7 = firm_exit_by_age_data(:,8);
y8 = firm_exit_by_age_data(:,9);
y9 = firm_exit_by_age_data(:,10);

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
legend(hleglines,'1','2','3','4','5','6-10','11-15','16-20','20 plus','Location','NorthWest')
legend boxoff
set(gca,'Fontsize',6.,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_great_recession_exit_by_age_data','-dpng','-r300')




