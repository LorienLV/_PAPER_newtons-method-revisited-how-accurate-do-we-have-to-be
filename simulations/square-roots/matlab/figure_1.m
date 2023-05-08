% This script reproduces Fig 1. of the article: 
% Carl Christian Kjelgaard Mikkelsen, Lorién López-Villellas, Pablo García-Risueño
% "Newton's method revisited: How accurate do we have to be?"
% 
% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2022-04-22 Initial programming and testing
%   2023-05-08 Documentation and file path updated

% Close all graphics windows
close all

% Select the sample points
s=linspace(1,4,101);

% Select the random seed
seed=2021;

% Get a handle to a new figure
h=figure();

% Offset
t=0;

% Set the units and the position for the figure
h.Units='Pixels';
h.Position=[0, 0, 1600+t, 800];

% //////////////////////////////////////////////////////////////////
%    A large value of epsilon
% //////////////////////////////////////////////////////////////////

% Select the maximum number of iterations
maxit=7;

% Select the size of the random errors
epsilon=1e-2;

% Run the experiment
[y, rel]=newton_sqrt(s,maxit,seed,epsilon);

% Plot the error curves
sf1=subplot(1,3,1); plot(s,log10(abs(rel)),'LineWidth',2); 
sf1.Units='Pixels'
sf1.Position=[t+75 75 425 650];
lgd=legend(strcat('k=',string(0:maxit)));
lgd.Units='Pixels';
lgd.Position=[t+184 750 208 40];
lgd.NumColumns = 4;

% Finalize graphics
xlabel('\alpha'); ylabel('log10 relative error'); grid; grid minor;

% //////////////////////////////////////////////////////////////////
%    A critical value of epsilon
% //////////////////////////////////////////////////////////////////

% Select the maximum number of Newton steps
maxit=3;

% Select the size of the random errors
epsilon=1e-8;

% Run the experiment
[y, rel]=newton_sqrt(s,maxit,seed,epsilon);

% Plot the error curves
sf2=subplot(1,3,2); plot(s,log10(abs(rel)),'LineWidth',2); 
sf2.Units='Pixels'
sf2.Position=[t+600 75 425 650];
lgd=legend(strcat('k=',string(0:maxit)));
lgd.Units='Pixels';
lgd.Position=[t+709 750 208 40];
lgd.NumColumns = 4;

% Finalize graphics
xlabel('\alpha'); ylabel('log10 relative error'); grid; grid minor;

% //////////////////////////////////////////////////////////////////
%    A tiny value of epsilon
% //////////////////////////////////////////////////////////////////

% Select the maximum number of iterations
maxit=3;

% Select the size of the random errors
epsilon=1e-12;

% Run the experiment
[y, rel]=newton_sqrt(s,maxit,seed,epsilon);

% Plot the error curves
sf3=subplot(1,3,3); plot(s,log10(abs(rel)),'LineWidth',2); 
sf3.Units='Pixels'
sf3.Position=[t+1125 75 425 650];
lgd=legend(strcat('k=',string(0:maxit)));
lgd.Units='Pixels';
lgd.Position=[t+1233 750 208 40];
lgd.NumColumns = 4;

% Finalize graphics
xlabel('\alpha'); ylabel('log10 relative error'); grid; grid minor;

% Set the fontsizes
sf1.FontSize=18;
sf2.FontSize=18;
sf3.FontSize=18;

% Convert to landscape and export
h.PaperOrientation='landscape';

% Export the figure to pdf
exportgraphics(h,'../fig/figure_1.pdf');