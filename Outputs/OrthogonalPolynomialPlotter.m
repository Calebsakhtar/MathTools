clear all;
close all;

% Set Style
cal_orange = [237 125 49]/255;
cal_blue = [50 142 216]/255;
cal_grey_light = [0.95 0.95 0.95];
cal_grey_dark = [0.35 0.35 0.35];

% Read the data output from the C++ RK4 Method
x_mat = readmatrix('xPolynomials.csv');
P0_mat = readmatrix("Poly00.csv");
P1_mat = readmatrix("Poly01.csv");
P2_mat = readmatrix("Poly02.csv");
P3_mat = readmatrix("Poly03.csv");
P4_mat = readmatrix("Poly04.csv");
P5_mat = readmatrix("Poly05.csv");

% Extract the relevant data from the variables
x = x_mat(:,2);
P0 = P0_mat(:,2);
P1 = P1_mat(:,2);
P2 = P2_mat(:,2);
P3 = P3_mat(:,2);
P4 = P4_mat(:,2);
P5 = P5_mat(:,2);

% Plot y1
figure(1)
% set figure position and size:
set(gcf,'position',[160 200 500 500])

% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')

hold on
plot(x,P0,"r",'Linewidth',1,'DisplayName','P0')
plot(x,P1,"g",'Linewidth',1,'DisplayName','P1')
plot(x,P2,"b",'Linewidth',1,'DisplayName','P2')
plot(x,P3,"m",'Linewidth',1,'DisplayName','P3')
plot(x,P4,"y",'Linewidth',1,'DisplayName','P4')
plot(x,P5,"k",'Linewidth',1,'DisplayName','P5')
hold off

xlabel('x');
ylabel('Orthogonal Polynomials','Interpreter','tex');

xlim([-2,2]);
ylim([-50,50]);

leg = legend('location','northoutside','orientation','horizontal');

% set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)

% % For a vectorial figure (for latex):
% print -deps2c y1.eps
% 
% % For an image figure at resolution 300 pixels/inch (for word):
% print -dpng -r300 y1.png