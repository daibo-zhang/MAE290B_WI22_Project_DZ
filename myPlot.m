function myPlot(xValue,yValue,xAxisLabel,yAxisLabel,myTitle,fontSize,...
    varargin)

% A custom function for plotting
% xVal - x-data to be plotted
% yVal - y-data to be plotted
% xl - label of x-axis
% yl - label of y-axis
% ttl - title of the plot
% fontSize - fontSize of the title, axes labels, and legend
% varargin{1} - axis scaling: linear, semilog, or loglog
% varargin{2} - legend of the plot

% Select plotter
if nargin < 7
    plotter = 'linear';
else
    plotter = varargin{1};
end

% Plot with specified plotter
switch plotter
    case 'linear'
        plot(xValue,yValue,'linewidth',2);
    case 'semilogx'
        semilogx(xValue,yValue,'linewidth',2);
    case 'semilogy'
        semilogy(xValue,yValue,'linewidth',2);
    case 'loglog'
        loglog(xValue,yValue,'linewidth',2);
end

% Label the plot
set(gca,'FontSize',fontSize,'TickLabelInterpreter','latex');
xlabel(xAxisLabel,'Interpreter','latex');
ylabel(yAxisLabel,'Interpreter','latex');
title(myTitle,'Interpreter','latex');

% Make legend if specified
if nargin == 8
    legend(varargin{2},'Location','best','Interpreter','latex');
elseif nargin == 9
    legend(varargin{2},'Location',varargin{3},'Interpreter','latex');
end

end