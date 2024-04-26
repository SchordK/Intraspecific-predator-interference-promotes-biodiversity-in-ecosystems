function dydt = odefcn(t,y,par)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance

dydt = zeros(7,1);
dydt(1) = par.a*(y(7)-y(1)-y(2))*(y(5)-y(1)-2*y(3))-(par.d+par.k)*y(1);   %x1_chasing pair abundance
dydt(2) =par.aa*(y(7)-y(1)-y(2))*(y(6)-y(2)-2*y(4))-(par.dd+par.kk)*y(2);  %x2_chasing pair abundance
dydt(3) =par.a1*(y(5)-y(1)-2*y(3))*(y(5)-y(1)-2*y(3))-(par.d1+par.k1)*y(3);   %y1_intraspecific interference pair abundance
dydt(4) =par.aa1*(y(6)-y(2)-2*y(4))*(y(6)-y(2)-2*y(4))-(par.dd1+par.kk1)*y(4);  %y2_intraspecific interference pair abundance
dydt(5) =par.w*par.k*y(1)-par.k1*y(3)-par.DD*y(5);  %C1_species abundance
dydt(6) = par.ww*par.kk*y(2)-par.kk1*y(4)-par.D*y(6);  %C2_species abundance
dydt(7) =par.R00*(1-y(7)/par.K0)-(par.k*y(1)+par.kk*y(2));  %R_species abundance
end





