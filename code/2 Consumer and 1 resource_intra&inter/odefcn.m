function dydt = odefcn(t,y,par)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance
dydt = zeros(8,1);
dydt(1) = par.a1*(y(8)-y(1)-y(2))*(y(6)-y(1)-2*y(4)-y(3))-(par.d1+par.k1)*y(1);  %x1_chasing pair abundance
dydt(2) = par.a2*(y(8)-y(1)-y(2))*(y(7)-y(2)-2*y(5)-y(3))-(par.d2+par.k2)*y(2);  %x2_chasing pair abundance
dydt(3) = par.b*(y(6)-y(1)-2*y(4)-y(3))*(y(7)-y(2)-2*y(5)-y(3))-par.e*y(3); %z_interspecific interference pair abundance
dydt(4) = par.u1*(y(6)-y(1)-2*y(4)-y(3))*(y(6)-y(1)-2*y(4)-y(3))-par.v1*y(4); %y1_intraspecific interference pair abundance
dydt(5) = par.u2*(y(7)-y(2)-2*y(5)-y(3))*(y(7)-y(2)-2*y(5)-y(3))-par.v2*y(5); %y2_intraspecific interference pair abundance
dydt(6) = par.w1*par.k1*y(1)-par.p1*y(4)-par.D1*y(6); %C1_species abundance
dydt(7) = par.w2*par.k2*y(2)-par.p2*y(5)-par.D2*y(7); %C2_species abundance
dydt(8) = par.Ra*(1-y(8)/par.K0)-par.k1*y(1)-par.k2*y(2);  %R_species abundance
end
