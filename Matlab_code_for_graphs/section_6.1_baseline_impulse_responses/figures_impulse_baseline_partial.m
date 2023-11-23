clear, clc, close all

load ('data_figures.mat')

%% Impulse response - partial equilibrium analysis

% specification
fontsz = 6;
paper_unit = 'inches'; % 'inches', 'normalized', 'centimeters', 'points'
paper_size = [2.2 1.9]; % [width height]
paper_position = [0 0 2.2 1.9]; % [left bottom width height].
paper_orientation = 'portrait'; % 'portrait (default)', 'landscape'
font_name = 'Times New Roman'; % 'Helvetica', 'Times New Roman'

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

% Productivity shock: entry and exit rate
x1 = impulse_baseline(:,1);
x2 = impulse_baseline(:,1);
x3 = impulse_baseline(:,1);
x4 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,6);
y2 = impulse_baseline(:,6);
y3 = impulse_baseline_w_r_fixed(:,7);
y4 = impulse_baseline(:,7);

figure
plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x4,y4,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,8])
ylim([-4,6]) 
xticks((0:1:8))
legend('entry: partial equilibrium','entry: general equilibrium','exit: partial equilibrium','exit: general equilibrium','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_prodty_pe_exit','-dpng','-r300')

% Productivity shock: r
x1 = impulse_baseline(:,1);
x2 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,10);
y2 = impulse_baseline(:,10);
figure
plot(x1,y1,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,8])
ylim([-1,1]) 
xticks((0:1:8))
legend('partial equilibrium','general equilibrium','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_prodty_pe_r','-dpng','-r300')

% Productivity shock: w
x1 = impulse_baseline(:,1);
x2 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,11);
y2 = impulse_baseline(:,11);
figure
plot(x1,y1,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage deviation','Color','k')
xlim([0,8]) 
ylim([-2,4])
xticks((0:1:15))
legend('partial equilibrium','general equilibrium','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_prodty_pe_w','-dpng','-r300')

% Credit shock: entry and exit rate
x1 = impulse_baseline_w_r_fixed(:,1);
x2 = impulse_baseline(:,1);
x3 = impulse_baseline(:,1);
x4 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,16);
y2 = impulse_baseline(:,16);
y3 = impulse_baseline_w_r_fixed(:,17);
y4 = impulse_baseline(:,17);
figure
plot(x1,y1,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x3,y3,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x4,y4,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,8])
ylim([-4,6]) 
xticks((0:1:8))
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_theta_pe_exit','-dpng','-r300')

% Credit shock: r
x1 = impulse_baseline(:,1);
x2 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,20);
y2 = impulse_baseline(:,20);
figure
plot(x1,y1,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage point deviation','Color','k')
xlim([0,8]) 
ylim([-1,1]) 
xticks((0:1:15))
legend('partial equilibrium','general equilibrium','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_theta_pe_r','-dpng','-r300')

% Credit shock: w
x1 = impulse_baseline(:,1);
x2 = impulse_baseline(:,1);
y1 = impulse_baseline_w_r_fixed(:,21);
y2 = impulse_baseline(:,21);
figure
plot(x1,y1,'Color',[0.8 0 0],'LineWidth',1.3,'LineStyle','--'), hold on
plot(x2,y2,'Color',[0.2 0.5 1.0],'LineWidth',1.3,'LineStyle','-')
plot(x1,0*y1,'Color','k','LineWidth',0.3,'LineStyle','-')
xlabel('quarter','Color','k')
ylabel('percentage deviation','Color','k')
xlim([0,8]) 
ylim([-2,4])
xticks((0:1:15))
legend('partial equilibrium','general equilibrium','Location','NorthEast')
legend boxoff
set(gca,'Fontsize',fontsz,'FontName',font_name,'xcolor','k','ycolor','k')
fig = gcf;
fig.PaperUnits = paper_unit;
fig.PaperSize = paper_size;
fig.PaperPosition = paper_position;
fig.PaperOrientation = paper_orientation;
print('partial_equilibrium\fig_ir_theta_pe_w','-dpng','-r300')





