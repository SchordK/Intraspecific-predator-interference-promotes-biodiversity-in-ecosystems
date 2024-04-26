% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance  and position.
% Click the plot to show the (stochastic) trajectories for species abundances vs time and position.

clear
clc
tic
close all
%define time mesh
t0=1e5;
tmesh=linspace(0,t0,t0);
T=tmesh(1:end-1);
tspan = 0:1:t0;
%define initial condition
init=[80;60;60;0;0;0;0];
y0 = [0.0 0.0 0.0 0.0 60 60 80];


% set the model parameters
width=50;  par.r=5; %x=[-width,width]
L=2*width;
r= par.r;
vC=1;vR=0.1;
ll=200;l2=200;
par.a=ll*2*r*sqrt((vC)^2+(vR)^2)/(L*L);   %a1 d1 k1
par.a1=l2*2*r*sqrt(2)*vC/(L*L);%a1' d1'
par.k1=0;   %k1'
par.aa=ll*2*r*sqrt((vC)^2+(vR)^2 )/(L*L);    %a2 d2 k2
par.aa1=l2*2*r*sqrt(2)*vC/(L*L);%a2' d2'
par.kk1=0;     %k2'
par.R00=100;par.K0=1000;
par.DD = 0.06;par.D =0.0605;
par.w = 0.33; par.ww = 0.33;
par.k = 0.2; par.kk = 0.2;
par.d1 = 0.8; par.dd1 = 0.8; %d'
par.d=0.7;par.dd=0.7;

[t,y] = ode45(@(t,y) odefcn(t,y,par),tspan,y0);  % ODEs simulation
time=[0 1e1 1e2 1e3 1e4 1e5 1e6];
[position ibm]=simulation_model_movie(init,width,tmesh,par,time); % IBM simulation

par1=1;
a=ibm(1,:)+ibm(4,:)+ibm(5,:);b=ibm(2,:)+ibm(4,:)+2*ibm(6,:);c=ibm(3,:)+ibm(5,:)+2*ibm(7,:);
%IBM=[1:par1:t0,bb,cc,aa];


figure%('visible','off')
hold on
plot(t,y(:,7),'r','linewidth',2);   %R
plot(t,y(:,5),'g','linewidth',2);   %C1
plot(t,y(:,6),'b','linewidth',2);   %C2
plot(T,a,'r');   %R
plot(T,b,'g')   %C1
plot(T,c,'b')   %C2
% legend('Resource','Consumer1','Consumer2')
% title(['\Delta = (D_{1}-D_{2})/D_{2}=',num2str((par.D-par.DD)/par.DD)])
xlabel('Time')
ylabel('Number')
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1,inf,1e1,inf])

toc

