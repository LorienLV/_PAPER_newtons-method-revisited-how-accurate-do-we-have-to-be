% This script generates Figures 2a, 2b and Figure 3 from the paper
% Carl Christian Kjelgaard Mikkelsen, Lorién López-Villellas, Pablo García-Risueño
% "Newton's method revisited: How accurate do we have to be?"

% Close all figures
close all

% Choose from the standard colors with cyclic reuse
mycolors=['r','g','b','c','m','k','r','g','b'];

% //////////////////////////////////////////////////////////////////
% Read the data directly from the relevant text files
% //////////////////////////////////////////////////////////////////

% The read the relative constraint violation
reshis=readmatrix('../data/reshis.txt');

% Read the normwise relative error between the correction t_k used for the
% symmetric quasi-Newton method and the correction s_k needed for Newton's method
normEk=readmatrix('../data/normEk.txt');

% Read the normwise relative error between the Lagrange multiplier, i.e., the zero z 
% and the current approximation x(k) generated using the symmetric quasi-Newton method
rk=readmatrix('../data/rk.txt');

% The string 'rk' is short for "normwise relative error after k steps".

% ////////////////////////////////////////////////////////////////////
% Select the data to display
% ////////////////////////////////////////////////////////////////////

% Lorién did 50000 timestep and recorded information every 2500 timesteps, so cols = 20.
% Lorién did 10 quasi-Newton steps, so rows=10;
[rows, cols]=size(normEk);

% It the strict mathematical sense we have
%
% r(k+1) <= 0.5*L*K(1+||D(k)||)*||z||*r(k)^2 + ||E(k)||*K*M*(1+||D(k)||*r(k) + ||D(k)||.
%
% but as long as we are not stagnating we therefore have reason to hope that
%
% r(k+1) <= ||E(k)||*K*M*(1+||D(k)||*r(k) 
% 
% In short, we have reason to expect lienar decay
ts=2500*(1:20);

% We eventually stagnate, so reset the number of rows manually
% Reset rows accordingly
rows=7;

% ///////////////////////////////////////////////////////////////////
% Generate figure 2a.
% ///////////////////////////////////////////////////////////////////

% Get a handle to the first figure
fig1=figure(); ax1=axes(fig1); hold(ax1);

% First plot the residual history aka. the relative constraint violation
for i=1:rows-1
  plot(ts, log10(reshis(i,:)),'Color',mycolors(i),'Linewidth',2,'Marker','*'); xlabel('time step'); ylabel('log10 relative constraint violation');
end
plot(ts, log10(reshis(rows,:)),'Color',mycolors(rows),'Linewidth',2,'Marker','s');
grid(ax1);
fig1.Units='Pixels';
fig1.Position=[75 75 825 650];
lgd1=legend(strcat('k=',string(0:rows-1)));
lgd1.Units='Pixels';
lgd1.Position=[350 610 150 30];
lgd1.NumColumns = 7;
ax1.FontSize=18;
exportgraphics(fig1,'../fig/figure_2a.pdf');

% //////////////////////////////////////////////////////////////////
% Generate figure 2b.
% //////////////////////////////////////////////////////////////////

% Get a handle to a second figure
fig2=figure(); ax2=axes(fig2); hold(ax2);
for i=1:rows-1
    plot(ts, log10(rk(i,:)),'Color',mycolors(i),'Linewidth',2,'Marker','*'); xlabel('time step'); ylabel('log10 relative error');
end
plot(ts, log10(rk(rows,:)),'Color',mycolors(rows),'Linewidth',2,'Marker','s');
grid(ax2); % axis(ax2,[0 50000 -16 0]);
fig2.Units='Pixels';
fig2.Position=[800 75 825 650];
lgd2=legend(strcat('k=',string(0:rows-1)));
lgd2.Units='Pixels';
lgd2.Position=[350 610 150 30];
lgd2.NumColumns = 7;
ax2.FontSize=18;
exportgraphics(fig2,'../fig/figure_2b.pdf');


% //////////////////////////////////////////////////////////////////////
% Generate figure 3
% //////////////////////////////////////////////////////////////////////

% Generate the data that correlates the measured values of r(k) and E(k)
data=zeros(rows-1,cols);
for j=1:cols
    for i=1:rows-1 
        data(i,j)=rk(i+1,j)/(rk(i,j)*normEk(i,j));
    end
end

% Get a handle to a third figure
fig3=figure();
for i=1:rows-1
   subplot(2,3,i); plot(ts,data(i,:),'Color',mycolors(i),'Linewidth',2,'Marker','*'); grid;
   xlabel('time step'); ylabel(strcat('\nu_',string(i-1))); 
   legend(strcat('k=',string(i-1)),'Location','Northoutside'); ax=gca; ax.FontSize=18;
end
% subplot(rows-1,1,rows-1); 
fig3.Units='Pixels';
% fig3.Position=[1525 75 825 650];
fig3.Position=[1525 75 1374 1139];
exportgraphics(fig3,'../fig/figure_3.pdf');