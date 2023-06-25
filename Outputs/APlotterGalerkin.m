clear all;
close all;

% Set Style
cal_orange = [237 125 49]/255;
cal_blue = [50 142 216]/255;
cal_grey_light = [0.95 0.95 0.95];
cal_grey_dark = [0.35 0.35 0.35];

% Read the data output from the C++ RK4 Method
x_mat = readmatrix('xDistribution.csv');
PDF_mat = readmatrix("Exact.csv");
PDF_approx_mat = readmatrix("Approx.csv");

% Extract the relevant data from the variables
x = x_mat(:,2);
PDF = PDF_mat(:,2);
PDF_approx = PDF_approx_mat(:,2);

% Plot PDF
figure(1)
% set figure position and size:
set(gcf,'position',[160 200 500 400])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(x,PDF,"b",'Linewidth',1)
plot(x, PDF_approx,":r",'Linewidth',1 )
hold off

xlabel('x');
ylabel('PDF','Interpreter','tex');

leg = legend('Exact','Approx','location','northoutside','orientation','horizontal');

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

xlim([0,14]);
ylim([0,1]);

% % For a vectorial figure (for latex):
% print -deps2c PDF.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 PDF.png

