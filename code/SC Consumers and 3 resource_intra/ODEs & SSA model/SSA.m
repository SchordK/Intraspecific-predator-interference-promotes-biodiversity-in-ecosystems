function tra= SSA(par,tmesh,y0,N,seed)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate stochastic population dynamics
% the 'type' variable determines of species or pair abundance


rng(110)
%initialize
t=0;
current=y0;
tra=zeros(5*N+3,length(tmesh));
tra(:,1)=y0;  %[CR CC C R]
index=2;
mass=zeros(1,N);
while t<tmesh(end)
    current(current<0)=0;
    current_R1=sum(current(1:N))+current(1+5*N);
    current_R2=sum(current(1+N:2*N))+current(2+5*N);
    current_R3=sum(current(1+2*N:3*N))+current(3+5*N);
    current_C=current(1:N)+current(1+N:2*N)+current(1+2*N:3*N)+2*current(3*N+1:4*N)+current(1+4*N:5*N);
    a=[(1-current_R1/par.K01) (1-current_R2/par.K02) (1-current_R3/par.K03)].*[par.R01 par.R02 par.R03];
    a=[a current(1+4*N:5*N)*current(1+5*N).*par.a1 current(1+4*N:5*N)*current(2+5*N).*par.a2 current(1+4*N:5*N)*current(3+5*N).*par.a3 ...
        current(1:N).*par.d1 current(1+N:2*N).*par.d2 current(1+2*N:3*N).*par.d3 ...
        current(1:N).*par.k1 current(1+N:2*N).*par.k2 current(1+2*N:3*N).*par.k3 ...
        current(1+4*N:5*N).^2.*par.u current(1+3*N:4*N).*par.v current(1+3*N:4*N).*par.p ...
        current_C.*par.D];     
    
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
    
    if i<=3
        current(i+5*N)=current(i+5*N)+1;
    else
        kk=floor((i-4)/N);
        ii=i-3-kk*N;
        switch kk
            case {0,1,2}
                current(kk+1+5*N)=current(kk+1+5*N)-1;
                current(ii+4*N)=current(ii+4*N)-1;
                current(ii+kk*N)=current(ii+kk*N)+1;
            case {3,4,5}
                current(kk-2+5*N)=current(kk-2+5*N)+1;
                current(ii+4*N)=current(ii+4*N)+1;
                current(ii+(kk-3)*N)=current(ii+(kk-3)*N)-1;
            case {6}
                mass(ii)=mass(ii)+par.w1(ii);
                current(ii+4*N)=current(ii+4*N)+1+floor(mass(ii));
                mass(ii)=mass(ii)-floor(mass(ii));
                current(ii+(kk-6)*N)=current(ii+(kk-6)*N)-1;
            case {7}
                mass(ii)=mass(ii)+par.w2(ii);
                current(ii+4*N)=current(ii+4*N)+1+floor(mass(ii));
                mass(ii)=mass(ii)-floor(mass(ii));
                current(ii+(kk-6)*N)=current(ii+(kk-6)*N)-1;
            case {8}
                mass(ii)=mass(ii)+par.w3(ii);
                current(ii+4*N)=current(ii+4*N)+1+floor(mass(ii));
                mass(ii)=mass(ii)-floor(mass(ii));
                current(ii+(kk-6)*N)=current(ii+(kk-6)*N)-1;
            case {9}
                current(ii+4*N)=current(ii+4*N)-2;
                current(ii+3*N)=current(ii+3*N)+1;
            case {10}
                current(ii+4*N)=current(ii+4*N)+2;
                current(ii+3*N)=current(ii+3*N)-1;
            case {11}
                current(ii+4*N)=current(ii+4*N)+1;
                current(ii+3*N)=current(ii+3*N)-1;
            case {12}
                current(ii+4*N)=current(ii+4*N)-1;
        end
    end
    
    
    t=t+dt;
    
    while index<=length(tmesh) &&t>tmesh(index)
        tra(:,index)=current;
        index=index+1;
    end
end