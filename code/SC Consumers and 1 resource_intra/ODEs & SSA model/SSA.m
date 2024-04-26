function tra= SSA(par,tmesh,y0,N,seed)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate stochastic population dynamics
% the 'type' variable determines of species or pair abundance

rng(110)
%initialize
t=0;
current=y0;
tra=zeros(3*N+1,length(tmesh));
tra(:,1)=y0;  %[CR CC C R]
index=2;
m1=[par.R00];
m2=[par.a par.d par.k];
m3=[par.u par.v par.p];
mass=zeros(1,N);
while t<tmesh(end)
    current(current<0)=0;
    current_R=sum(current(1:N))+current(1+3*N);
    current_C=current(1:N)+2*current(N+1:2*N)+current(1+2*N:3*N);
    a=[(1-current_R/par.K0)].*m1;
    for i=1:N
        a=[a, [current(1+3*N)*current(i+2*N) current(i) current(i)].*m2, ...
            [current(i+2*N).*(current(i+2*N)-1)/2 current(N+i) current(N+i)].*m3, ...
            (current(i)+2*current(N+i)+current(i+2*N))*par.D(i)];
    end
    
    
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
    
    if i==1
        current(1+3*N)=current(1+3*N)+1;
    else
        kk=floor((i-2)/7);
        ii=i-1-kk*7;
        switch ii
            case {1}
                current(1+3*N)=current(1+3*N)-1;
                current(kk+1+2*N)=current(kk+1+2*N)-1;
                current(kk+1)=current(kk+1)+1;
            case {2}
                current(1+3*N)=current(1+3*N)+1;
                current(kk+1+2*N)=current(kk+1+2*N)+1;
                current(kk+1)=current(kk+1)-1;
            case {3}
                mass(kk+1)=mass(kk+1)+par.w;
                current(kk+1+2*N)=current(kk+1+2*N)+1+floor(mass(kk+1));
                mass(kk+1)=mass(kk+1)-floor(mass(kk+1));
                current(kk+1)=current(kk+1)-1;
            case {4}
                current(kk+1+2*N)=current(kk+1+2*N)-2;
                current(kk+1+N)=current(kk+1+N)+1;
            case {5}
                current(kk+1+2*N)=current(kk+1+2*N)+2;
                current(kk+1+N)=current(kk+1+N)-1;
            case {6}
                current(kk+1+2*N)=current(kk+1+2*N)+1;
                current(kk+1+N)=current(kk+1+N)-1;
            case {7}
                current(kk+1+2*N)=current(kk+1+2*N)-1;
        end
    end
    
    
    t=t+dt;
    
    while index<=length(tmesh) &&t>tmesh(index)
        tra(:,index)=current;
        index=index+1;
    end
end