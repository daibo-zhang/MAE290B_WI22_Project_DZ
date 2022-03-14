%% Preliminaries

% Author: Daibo Zhang (A13591601)
% UCSD MAE290B WI22 Final Project

% Clean up
close all;
format long;

% Load parameter list
fpParam  = readtable('fpParam.csv');

% Define constants
alpha = fpParam.alpha; 
a = fpParam.a;
Omega = fpParam.Omega;
dt = fpParam.dt;
h = fpParam.h;
kappa = alpha * dt / (2 * h^2);

% Define spatial grid
xBound = [fpParam.xBound1,fpParam.xBound2];
yBound = [fpParam.yBound1,fpParam.yBound2];
[xMesh,yMesh] = meshgrid(xBound(1):h:xBound(2),yBound(1):h:yBound(2));
N = length(xMesh) - 2;

% Identify location of point of interest (0.55,0.45)
pLoc = [find(xBound(1):h:xBound(2)==0.55),...
    find(yBound(1):h:yBound(2)==0.45)];

% Import solution
bData = importdata('solutionsPartB.mat');
uOutB = bData.uOutB;
pointOutB = bData.pointOutB;
tOutB = bData.tOutB;
dData = importdata('solutionsPartD.mat');
uOutD = dData.uOutD;
pointOutD = dData.pointOutD;
tOutD = dData.tOutD;

%% Part b: visualizing results

% Plot time evolution of solution at point
figure;
myPlot(tOutB,pointOutB,'Time $t$','Temperature $T$',...
    {'Temporal Evolution of Temperature','at $x=0.55$, $y=0.45$'},20);

% Plot steady state solution as a contour plot
figure;
uSsB = uOutB(:,:,end);
contourf(xMesh,yMesh,uSsB,'ShowText','on');
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title({'Steady State Temperature Distribution',...
    strcat('Obtained at $t$ = ',string(round(tOutB(end),3)))},...
    'Interpreter','latex');

%% Part c: error analysis

% Display steady state error at the point of interest
pointExp = sin(4*pi*0.55) * sin(8*pi*0.45) / (3.2*pi*pi);
pointErrorB = (pointOutB(end) - pointExp) / pointExp;
fprintf('Steady state error on the point is %1.4e \n',pointErrorB);

% Define expected steady state solution
uSsExp = sin(4*pi*xMesh) .* sin(8*pi*yMesh) / (3.2*pi*pi);

% Plot analytical steady state solution 
figure;
contourf(xMesh,yMesh,uSsExp,'ShowText','on');
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title('Expected Steady State Temperature Distribution',...
    'Interpreter','latex');

% Calculate steady state error on each interior point
uSsError = (uSsB - uSsExp);

% Plot error side by side with analytical solution 
figure;
myColorMap = [linspace(1,1)',linspace(0,1)',linspace(0,1)';...
    linspace(1,0)',linspace(1,0)',linspace(1,1)'];
errH = heatmap(flip(uSsError,2),'colorMap',myColorMap,...
    'gridVisible','off','colorLimit',[-3E-4,3E-4]);
hmLabel = string(xBound(1):h:xBound(2));
hmLabel(mod(0:length(hmLabel)-1,8)~=0) = " ";
set(gca,'FontSize',20)
errH.FontName = 'CMU Serif';
errH.XDisplayLabels = hmLabel;
errH.YDisplayLabels = flip(hmLabel);
errH.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
errH.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
errH.NodeChildren(3).Title.Interpreter = 'latex';
errH.XLabel = '$x$';
errH.YLabel = '$y$';
errH.Title = {'Error in Computed Steady State','Temperature Distribution'};

% Plot error correlation
figure;
scatter(reshape(uSsExp,(N+2)*(N+2),1),...
    reshape(uSsError,(N+2)*(N+2),1));
set(gca,'FontSize',20,'TickLabelInterpreter','latex');
xlabel('$T_e^{Exp}$','Interpreter','latex');
ylabel('$T_e^{Num}-T_e^{Exp}$','Interpreter','latex');
title('Steady State Error vs. Expected Solution','Interpreter','latex');

%% Part b: visualizing results

% Plot time evolution of solution at point
figure;
hold on;
plot(tOutB,pointOutB,'lineWidth',2);
myPlot(tOutD,pointOutD,'Time $t$','Temperature $T$',...
    {'Temporal Evolution of Temperature','at $x=0.55$, $y=0.45$'},20,...
    'linear',{'$T_0 = 0.01\sin(\pi x)\sin(\pi y)$','$T_0 = 0$'},'best');

% Plot temperature evolution at several select timepoints
figure;
tps = [2,6,1201];
for i = 1:length(tps)
    t = string(tOutB(tps(i)));
    subplot(length(tps),2,i*2-1);
    contourf(xMesh,yMesh,uOutB(:,:,tps(i)),'ShowText','on');
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    title(strcat('Part B: $t$ = ',t),'Interpreter','latex');
    subplot(length(tps),2,i*2);
    contourf(xMesh,yMesh,uOutD(:,:,tps(i)),'ShowText','on');
    set(gca,'FontSize',20,'TickLabelInterpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    title(strcat('Part D: $t$ = ',t),'Interpreter','latex');
end
