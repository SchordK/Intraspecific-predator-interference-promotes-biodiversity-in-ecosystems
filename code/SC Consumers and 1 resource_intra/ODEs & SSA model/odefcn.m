function dydt = odefcn(t,y,N,par)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance

dydt = zeros(3*N+1,1);

for i=1:N

dydt(i)= par.a*(y(1+3*N)-sum(y(1:N)))*(y(i+2*N)-y(i)-2*y(i+N))-(par.d+par.k)*y(i); %xi_chasing pair abundance

dydt(i+N) =par.u*(y(i+2*N)-y(i)-2*y(i+N))*(y(i+2*N)-y(i)-2*y(i+N))-(par.v+par.p)*y(i+N); %yi_intraspecific interference pair abundance

dydt(i+2*N) =par.w*par.k*y(i)-par.p*y(i+N)-par.D(i)*y(i+2*N);  %Ci_species abundance

end

dydt(1+3*N) =par.R00*(1-y(1+3*N)/par.K0)-par.k*sum(y(1:N));  %R_species abundance
end




