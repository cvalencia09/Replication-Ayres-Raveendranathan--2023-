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
% Column 1: quarter	
% Column 2: shock z	
% Column 3: output z	
% Column 4: employment z	
% Column 5: firm debt to output z	
% Column 6: firm entry z	
% Column 7: firm exit z	
% Column 8: consumption z	
% Column 9: investment z	
% Column 10: r z	
% Column 11: w z	

% Column 12: shock theta	
% Column 13: output theta
% Column 14: employment theta
% Column 15: firm debt to output theta
% Column 16: firm entry theta
% Column 17: firm exit theta
% Column 18: consumption theta	
% Column 19: investment theta
% Column 20: r theta
% Column 21: w theta

% Column 22: zero

% Entry
x1 = impulse_elastic_entry(:,1);
x2 = impulse_elastic_entry(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_elastic_entry(:,6);
y2 = impulse_elastic_entry(:,16);
y3 = data_average_recession_entry_exit(:,2);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p6 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-10,6])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_entry','-dpng','-r300')





%Exit

x1 = impulse_elastic_entry(:,1);
x2 = impulse_elastic_entry(:,1);
x3 = data_average_recession_entry_exit(:,1);
y1 = impulse_elastic_entry(:,7);
y2 = impulse_elastic_entry(:,17);
y3 = data_average_recession_entry_exit(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'); hold on
p2 = plot(x2,y2,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--');
p3 = plot(x3,y3,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p6 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,20])
ylim([-10,6])
hleglines = [p1(1) p2(1) p3(1)];
legend(hleglines,'productivity','collateral constraint','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_ir_exit','-dpng','-r300')




