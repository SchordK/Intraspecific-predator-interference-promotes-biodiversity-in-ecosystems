function dydt = odefcnAbiotic(t,y,u1,u2,a1,d1,k1,a2,d2,k2,v1,p1,v2,p2,Ra,K0,w1,w2,D1,D2)
%Chasing pair & intraspecific interference model 
% returns a function that can be used to simulate population dynamics
% the 'type' variable determines of species or pair abundance

dydt = zeros(7,1);
dydt(1) = a1*(y(7)-y(1)-y(2))*(y(5)-y(1)-2*y(3))-(d1+k1)*y(1); %x1_chasing pair abundance
dydt(2) =a2*(y(7)-y(1)-y(2))*(y(6)-y(2)-2*y(4))-(d2+k2)*y(2); %x2_chasing pair abundance
dydt(3) =u1*(y(5)-y(1)-2*y(3))*(y(5)-y(1)-2*y(3))-(v1+p1)*y(3); %y1_intraspecific interference pair abundance
dydt(4) =u2*(y(6)-y(2)-2*y(4))*(y(6)-y(2)-2*y(4))-(v2+p2)*y(4); %y2_intraspecific interference pair abundance
dydt(5) =w1*k1*y(1)-p1*y(3)-D1*y(5); %C1_species abundance
dydt(6) = w2*k2*y(2)-p2*y(4)-D2*y(6); %C2_species abundance
dydt(7) =Ra*(1-y(7)/K0)-(k1*y(1)+k2*y(2)); %R_species abundance
end





