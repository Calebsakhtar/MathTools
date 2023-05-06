clear all;
close all;

% Set Style
cal_orange = [237 125 49]/255;
cal_blue = [50 142 216]/255;
cal_grey_light = [0.95 0.95 0.95];
cal_grey_dark = [0.35 0.35 0.35];

% Read the data output from the C++ RK4 Method
x_mat = readmatrix('x.csv');
x_integ_analytic_mat = readmatrix("xAnalytic.csv");
x_integ_numeric_mat = readmatrix("xNumeric.csv");

% Extract the relevant data from the variables
x = x_mat(:,2);
x_integ_analytic = x_integ_analytic_mat(:,2);
x_integ_numeric = x_integ_numeric_mat(:,2);

% Plot y1
figure(1)
% set figure position and size:
set(gcf,'position',[160 200 500 400])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(x,x_integ_numeric,"--b",'Linewidth',1,'DisplayName','Numeric')
plot(x,x_integ_analytic,":r",'Linewidth',1,'DisplayName','Analytic')
hold off

xlabel('x');
ylabel('x^2','Interpreter','tex');

%leg = legend('location','northoutside','orientation','horizontal');
legend();

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

% % For a vectorial figure (for latex):
% print -deps2c y1.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 y1.png