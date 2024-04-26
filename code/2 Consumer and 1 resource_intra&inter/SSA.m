function tra= SSA(par,tmesh,init)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate stochastic population dynamics
% the 'type' variable determines of species or pair abundance

rng(120)
%initialize
t=0;
current=init;
tra=zeros(8,length(tmesh));
tra(:,1)=init;
index=2;
C=zeros(8,1);
m=[par.R00 par.a1 par.d1 par.k1 par.a2 par.d2 par.k2 par.u1 par.v1 par.p1 par.u2 par.v2 par.p2 par.b par.e par.D1 par.D2];
mass1=0;mass2=0;
while t<tmesh(end)
    current(current<0)=0;
    current_R=current(1)+current(4)+current(5);
    current_C1=current(2)+current(4)+2*current(6)+current(8);
    current_C2=current(3)+current(5)+2*current(7)+current(8);
    a=[(1-current_R/par.K0) current(1)*current(2) current(4) current(4) current(1)*current(3) current(5) current(5) current(2)*(current(2)-1)/2 current(6) current(6) current(3)*(current(3)-1)/2 current(7) current(7) current(1)*current(2) current(8) current_C1 current_C2].*m;
    
    sto=[1 -1 1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0   
        0 -1 1 1 0 0 0 -2 1 1 0 0 0 -1 1 -1 0   
        0 0 0 0 -1 1 1 0 0 0 -2 1 1 -1 1 0 -1   
        0 1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0  
        0 0 0 0 1 -1 -1 0 0 0 0 0 0 0 0 0 0 
        0 0 0 0 0 0 0 1 -1 -1 0 0 0 1 -1 0 0   
        0 0 0 0 0 0 0 0 0 0 1 -1 -1 1 -1 0 0   
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0];  
    C=sto*a';
    
    %select random reaction time
    a0=sum(a);
    r1=rand();
    dt=-log(r1)/a0;
    
    %select random reaction
    M=cumsum(a);
    r2=rand();
    for i=1:length(M)
        if M(i)>r2*a0
            break
        end
    end
    
    if i==4
        mass1=mass1+par.w1;
        sto(2,4)=sto(2,4)+floor(mass1);
        mass1=mass1-floor(mass1);
    elseif i==7
        mass2=mass2+par.w2;
        sto(3,7)=sto(3,7)+floor(mass2);
        mass2=mass2-floor(mass2);
    end
    current=current+sto(:,i);
    t=t+dt;
    
    while index<=length(tmesh) &&t>tmesh(index)
        tra(:,index)=current;
        index=index+1;
    end
end