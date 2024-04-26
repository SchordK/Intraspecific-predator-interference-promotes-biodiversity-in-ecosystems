% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.
clear
tic
clc


% set the model parameters
 Alpha = 1.4; par.a1 = 0.1; par.a2 = 0.1;
par.k1 = 0.12; par.k2 = 0.12;
par.d1 = 0.3; par.d2 = 0.3;
par.u1 =Alpha*par.a1; par.u2 = Alpha*par.a2;
par.p1 = 0.00; par.p2 = 0.00;
par.b = par.u1; 
par.w1 = 0.15; par.w2 = 0.15;
%par.w1 = 0.15; par.w2 = 0.15;
par.R00 = 5.5; par.Ra = 5.5;
par.D1 = 0.0125; par.D2 = 0.0120;
par.v1 = 0.5;par.v2 = 0.5;
par.K0=300;
par.e = 5;

 %define time mesh
t0=1e8;
tspan = 0:1:t0;
tmesh=linspace(0,t0,t0);

%define initial condition
init=[15;10;10;0;0;0;0;0];
y0 = [0.0 0.0 0.0 0.0 0.0 10 10 15];
% It would be better to set the initial value R0 less than par.K0


%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,par),tspan,y0);

%SSA simulation
tra= SSA(par,tmesh,init);

figure;
hold on
a=tra(1,:)+tra(4,:)+tra(5,:);b=tra(2,:)+tra(4,:)+2*tra(6,:)+tra(8,:);c=tra(3,:)+tra(5,:)+2*tra(7,:)+tra(8,:);
plot(tmesh,a,'r')   %R
plot(tmesh,b,'g')   %C1
plot(tmesh,c,'b')   %C2
plot(t,y(:,8),'r','linewidth',2);
plot(t,y( :,6),'g','linewidth',2);
plot(t,y( :,7),'b','linewidth',2);
set(gca,'XScale','log');
set(gca,'YScale','log');

% ode=[t,y( :,6),y( :,7),y( :,8),];
% SSA=[tmesh',b',c',a'];
figure;
plot3(y(:,7),y(:,6),y( :,8),'k','linewidth',2);hold on  
plot3(b,c,a,'r','linewidth',2);hold on  
xlabel('\fontname{Times New Roman}\fontsize{30}\it{R}');
ylabel('\fontname{Times New Roman}\fontsize{30}\it{C1}');
zlabel('\fontname{Times New Roman}\fontsize{30}\it{C2}');
set(gca,'FontName','Times New Roman','FontSize',24,'linewidth',2);
% % axis([-0.02,0.6,0,2.1]);
% %  print(gcf,'-dpng','-r600','figle.png')



toc


