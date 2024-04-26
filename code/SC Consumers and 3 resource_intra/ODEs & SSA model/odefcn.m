function dydt = odefcn(t,y,N,par)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance

dydt = zeros(5*N+3,1);

for i=1:N

dydt(i)= par.a1(i)*(y(i+4*N)-y(i)-y(i+N)-y(i+2*N)-2*y(i+3*N))*(y(5*N+1)-sum(y(1:N)))-(par.d1(i)+par.k1(i))*y(i);   %xi_chasing pair abundance

dydt(i+N)= par.a2(i)*(y(i+4*N)-y(i)-y(i+N)-y(i+2*N)-2*y(i+3*N))*(y(5*N+2)-sum(y(1+N:2*N)))-(par.d2(i)+par.k2(i))*y(i+N);   %x(i+N)_chasing pair abundance

dydt(i+2*N)= par.a3(i)*(y(i+4*N)-y(i)-y(i+N)-y(i+2*N)-2*y(i+3*N))*(y(5*N+3)-sum(y(1+2*N:3*N)))-(par.d3(i)+par.k3(i))*y(i+2*N);   %x(i+2N)_chasing pair abundance

dydt(i+3*N)= par.u(i)*(y(i+4*N)-y(i)-y(i+N)-y(i+2*N)-2*y(i+3*N))*(y(i+4*N)-y(i)-y(i+N)-y(i+2*N)-2*y(i+3*N))-(par.v(i)+par.p(i))*y(i+3*N);   %yi_intraspecific interference pair abundance

dydt(i+4*N) =(par.w1(i)*par.k1(i)*y(i)+par.w2(i)*par.k2(i)*y(i+N)+par.w3(i)*par.k3(i)*y(i+2*N))-par.D(i)*y(i+4*N);   %Ci_species abundance

end

dydt(1+5*N) =par.R01*(1-y(1+5*N)/par.K01)-par.k1*y(1:N);   %R1_species abundance

dydt(2+5*N) =par.R02*(1-y(2+5*N)/par.K02)-par.k2*y(1+N:2*N);   %R2_species abundance

dydt(3+5*N) =par.R03*(1-y(3+5*N)/par.K03)-par.k3*y(1+2*N:3*N);   %R3_species abundance

end





