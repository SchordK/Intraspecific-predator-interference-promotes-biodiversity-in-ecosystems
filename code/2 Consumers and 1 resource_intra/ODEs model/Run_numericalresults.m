% Implements the chasing pair & intraspecific interference model as described in:
% This script runs simulations while scanning the parameters for species abundance.
% Click the plot to show the trajectories for species abundances vs time.

clear
tic
clc
 % set the model parameters
Alpha = 1.25;
a1 = 0.5; a2 = 0.5; u1 = Alpha*a1; u2 = Alpha*a2; 
p1 = 0; p2 = 0; d2 = 0.5; d1 = 0.5; 
v1 = 0.02; v2 = 0.02; k1 = 0.4; k2 = 0.4;
D2 = 0.022; D1 = 0.0286;
w1 = 0.2; w2 = 0.2; 
Ra= 0.5; K0 = 10;  



%Analytical soluions
alpha1=D1/(w1*k1);alpha2=D2/(w2*k2);
beta1=u1/v1;beta2=u2/v2;
K1=(d1+k1)/a1;K2=(d2+k2)/a2;

o1=Ra/K0-k1/(2*beta1*K1)-k2/(2*beta2*K2);
o2=k1*(1-alpha1)/(2*beta1*alpha1*K1*K1)+k2*(1-alpha2)/(2*beta2*alpha2*K2*K2);
RR=(-o1+sqrt(o1*o1+4*o2*Ra))/(2*o2);
CC1=((1-alpha1)*RR*RR-K1*alpha1*RR)/(2*beta1*K1*K1*alpha1*alpha1);
CC2=((1-alpha2)*RR*RR-K2*alpha2*RR)/(2*beta2*K2*K2*alpha2*alpha2);

%Steady-state numerical solution
checkroll=10;
checkrollm=checkroll-1;
syms x1 x2 y1 y2 C1 C2 R
 S=vpasolve([a1*(R-x1-x2)*(C1-x1-2*y1)-(d1+k1)*x1==0, a2*(R-x1-x2)*(C2-x2-2*y2)-(d2+k2)*x2==0, u1*(C1-x1-2*y1)*(C1-x1-2*y1)-(v1+p1)*y1==0, u2*(C2-x2-2*y2)*(C2-x2-2*y2)-(v2+p2)*y2==0, w1*k1*x1-p1*y1-D1*C1==0, w2*k2*x2-p2*y2-D2*C2==0, Ra*(1-R/K0)-(k1*x1+k2*x2)==0],[x1,x2,y1,y2,C1,C2,R]);
       ll=length(S.x1);
        Sfinal=zeros(ll,checkroll);
        Sfinal(:,1)=S.x1;% RVC1
        Sfinal(:,2)=S.x2;% RVC2
        Sfinal(:,3)=S.y1;% C1VRVC1
        Sfinal(:,4)=S.y2;% C2VRVC2
        Sfinal(:,5)=S.C1;% C1_tot
        Sfinal(:,6)=S.C2;% C2_tot
        Sfinal(:,7)=S.R;% R_tot
        Sfinal(:,8)=S.R-S.x1-S.x2;% R_free
        Sfinal(:,9)=S.C1-S.x1-2*S.y1;% C1_free
        Sfinal(:,10)=S.C2-S.x2-2*S.y2;% C2_free
  
        
Sscore=Sfinal>0;
ssum=sum(Sscore');
tempmax=max(ssum); 
results=tempmax>checkrollm;
maxrow=find(ssum==tempmax);
C1Ta=Sfinal(maxrow,5);
C2Ta=Sfinal(maxrow,6);
RTa=Sfinal(maxrow,7);
result=[RR,RTa,CC1,C1Ta,CC2,C2Ta]
result1=[CC1,CC2,RR];
     




%ODEs simulation
t0=1e5;
tspan = 0:1:t0;
% initial species abundances
y0 = [0.0 0.0 0.0 0.0 1 1 2];

[t,y] = ode45(@(t,y) odefcnAbiotic(t,y,u1,u2,a1,d1,k1,a2,d2,k2,v1,p1,v2,p2,Ra,K0,w1,w2,D1,D2),tspan,y0);
figure;
plot(t,y(:,5),'g','linewidth',1);hold on
plot(t,y(:,6),'b','linewidth',1);hold on
plot(t,y(:,7),'r','linewidth',1);hold on
ode=[t,y(:,5),y(:,6),y(:,7)];

oo=[ode(1:500,:);ode(501:10:1000,:);ode(1001:100:10000,:);ode(10001:1000:100000,:)];

%save ode
plot(0.8*t0,CC1,'g--.',0.8*t0,CC2,'b--.',0.8*t0,RR,'r--.','LineWidth',2,'MarkerSize',10)
plot(0.85*t0,CC1,'g--.',0.85*t0,CC2,'b--.',0.85*t0,RR,'r--.','LineWidth',2,'MarkerSize',10)
plot(0.9*t0,CC1,'g--.',0.9*t0,CC2,'b--.',0.9*t0,RR,'r--.','LineWidth',2,'MarkerSize',10)
plot(0.95*t0,CC1,'g--.',0.95*t0,CC2,'b--.',0.95*t0,RR,'r--.','LineWidth',2,'MarkerSize',10)  
set(gca,'XScale','log');
set(gca,'YScale','log');
axis([1,inf,1e-1,1*1e2]);
     

