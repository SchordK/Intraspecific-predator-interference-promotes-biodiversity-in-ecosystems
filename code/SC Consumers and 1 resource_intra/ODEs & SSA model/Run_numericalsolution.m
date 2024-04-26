% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the (stochastic) trajectories for species abundances vs time.

clear
tic
clc
% eg.Appendix-figure 9A
N=20; % Type of consumer species
sn=1;%seed
s=randi([100,1000],1,sn);

% set the model parameters
Alpha = 1.25;
par.a = 0.1;par.u = Alpha*par.a; 
par.d= 0.5;
par.w = 0.22;par.k= 0.1;par.v = 0.6; par.p = 0.0; 
par.R00 =350; par.K0 = 1e6; 
rng(120)
par.D = normrnd(1.,0.35,1,N)*0.016;

%define initial condition

%C0=10*ones(1,N);
x0= zeros(1,N);
Y0= zeros(1,N);
C0= zeros(1,N);
col=N;
lin=1;
rng(110)
%p1=unifrnd(0,1,lin,col);
p1=abs(normrnd(0,1,lin,col));
for i=1:col
    if p1(1,i)>=0.03
        p1(1,i)=1;
    else p1(1,i)=0;
    end
end

for j=1:col
   if p1(1,j)==1
       rng(110)
         C0(1,j)=randperm(20,1);
   else C0(1,j)=0;

    end
end

R0=30;
% It would be better to set the initial value R0 less than par.K0

y0 = [x0 Y0 C0 R0];

%define time mesh
t0=1e5;
tspan = 0:0.5:t0;
tmesh=linspace(0,t0,t0);

%ODEs simulation
[t,y] = ode45(@(t,y) odefcn(t,y,N,par),tspan,y0);

%ode=y(t,1+2*N:3*N);

figure
plot(t,y(:,1+2*N:3*N));hold on
plot(t,y(:,1+3*N),'r','linewidth',3);
set(gca,'XScale','log')
set(gca,'YScale','log')
axis([1,t0,1e-2,1*1e5]);


tt2=cell(1,sn);
tt1=cell(1,sn);
tt=cell(1,sn);
TT=cell(1,sn);

%SSA simulation
    ss=1;
tra= SSA(par,tmesh,y0,N,s(ss));
tt{ss}=tra;


c=tra(1:N,:)+2*tra(1+N:2*N,:)+tra(1+2*N:3*N,:);r=sum(tra(1:N,:))+tra(1+3*N,:);
ssa=[tmesh',c',r'];
tt1{ss}=tra(1:N,t0)+2*tra(1+N:2*N,t0)+tra(1+2*N:3*N,t0);
tt2{ss}=sum(tra(1:N,t0))+tra(1+3*N,t0);
figure;
plot(tmesh,r,'r') ; hold on %R
plot(tmesh,c) ; hold on %Ci
xlabel('Time')
ylabel('Consumers and Resource')
set(gca,'XScale','log')
set(gca,'YScale','log')





toc
