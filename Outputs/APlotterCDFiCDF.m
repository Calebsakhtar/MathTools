clear all;
close all;

% Set Style
cal_orange = [237 125 49]/255;
cal_blue = [50 142 216]/255;
cal_grey_light = [0.95 0.95 0.95];
cal_grey_dark = [0.35 0.35 0.35];

% Read the data output from the C++ RK4 Method
x_mat = readmatrix('xDistribution.csv');
y_mat = readmatrix('yDistribution.csv');
PDF_mat = readmatrix("PDF.csv");
CDF_mat = readmatrix("CDF.csv");
iCDF_mat = readmatrix("iCDF.csv");

% Extract the relevant data from the variables
x = x_mat(:,2);
y = y_mat(:,2);
PDF = PDF_mat(:,2);
CDF = CDF_mat(:,2);
iCDF = iCDF_mat(:,2);

% Plot PDF
figure(1)
% set figure position and size:
set(gcf,'position',[160 200 500 400])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(x,PDF,"b",'Linewidth',1)
hold off

xlabel('x');
ylabel('PDF','Interpreter','tex');

%leg = legend('location','northoutside','orientation','horizontal');

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

% % For a vectorial figure (for latex):
% print -deps2c PDF.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 PDF.png

% Plot CDF
figure(2);

% set figure position and size:
set(gcf,'position',[160 200 500 400])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(x,CDF,"b",'Linewidth',1)
hold off

xlabel('x');
ylabel('CDF','Interpreter','tex');

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

% % For a vectorial figure (for latex):
% print -deps2c samples.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 samples.png

% Plot iCDF
figure(3);

% set figure position and size:
set(gcf,'position',[160 200 500 400])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(y,iCDF,"b",'Linewidth',1)
hold off

xlabel('y');
ylabel('iCDF (x)','Interpreter','tex');

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

% % For a vectorial figure (for latex):
% print -deps2c samples.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 samples.png

