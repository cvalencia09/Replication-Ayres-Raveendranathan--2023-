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
figure
plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-'), hold on
plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':')
%xlabel('period','Color','k')
ylabel('percentage deviation','Color','k')
xlim([1979.5,1984])
ylim([-6,2])
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_1980_recession_productivity','-dpng','-r300')

% Entry
x1 = aggregates_joint(:,1);
x2 = firm_entry_exit_1980_recession(:,1);
y1 = aggregates_joint(:,6);
y2 = firm_entry_exit_1980_recession(:,2);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([1979.5,1984])
ylim([-2,2])
hleglines = [p1(1) p2(1)];
legend(hleglines,'baseline','data','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_1980_recession_entry','-dpng','-r300')

% Exit
x1 = aggregates_joint(:,1);
x2 = firm_entry_exit_1980_recession(:,1);
y1 = aggregates_joint(:,7);
y2 = firm_entry_exit_1980_recession(:,3);

figure
p1 = plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','--'); hold on
p2 = plot(x2,y2,'Color',[0 0 0],'LineWidth',1.3,'LineStyle','-');
p3 = plot(x1,0*y1,'Color','k','LineWidth',0.5,'LineStyle',':');
%xlabel('period','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([1979.5,1984])
ylim([-2,2])
hleglines = [p1(1) p2(1)];
legend(hleglines,'baseline','data','Location','SouthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('fig_1980_recession_exit','-dpng','-r300')



