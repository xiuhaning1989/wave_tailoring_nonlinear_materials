function [dx] = impact_equation_of_motion_asymetric(t,x,nonlinear_spring_info,material_info,F)

%{
    ==========================================================================
    Inputs
    
    - t: time variable
    - y: output variable (i.e., displacement)
    - nonlinear_spring_info: Information for the nonlinear spring of unit cell
                             i.e., coeffcients for polynomial
    - a: unit cell height
    - h: unit cell depth (= sample depth)
    - m_lattice: lattice mass
    - m_impactor: impactor mass (commonly larger than lattice's one)
    - eta: (viscous) damping parameter
    - E_lin: linear modulus
    - N_asp: numbere of asperities per unit cell length
   
    Outputs
    
    - dydt: time derivative of displacement (i.e., velocity)
    ==========================================================================
%}

% initialization of some parameters
num_mass = length(x)/2;               % number of masses
dx = zeros(length(x),1);            % vector of time derivative of displacement
c1=nonlinear_spring_info(1);c2=nonlinear_spring_info(2);c3=nonlinear_spring_info(3);
m=material_info(1);m_impact=material_info(2);C_impact=material_info(3);
zeta=material_info(4);

i=1;                                   %the 1st unit cell of 1D chain   
dx(2*i-1) = x(2*i);   
dx(2*i)=(F-impact_force(C_impact,x(2*i-1),x(2*(i+1)-1)))/m_impact;
% dx(2*i)=(F+contactforce(a,h,x(2*(i+1)-1)-x(2*i-1),E_lin,N_asp))/m_impact;


% if num_mass>2
i=2;
dx(2*i-1) = x(2*i);   
if x(2*(i+1)-1)>=x(2*i-1)
    dx(2*i)=(c1*(x(2*(i+1)-1)-x(2*i-1)) ...
        +impact_force(C_impact,x(2*(i-1)-1),x(2*i-1))-2*zeta*(-x(2*(i+1))+x(2*i)))/m;
else
    dx(2*i)=(c1*(x(2*(i+1)-1)-x(2*i-1))-c2*(x(2*(i+1)-1)-x(2*i-1))^2+c3*(x(2*(i+1)-1)-x(2*i-1))^3+ ...
        +impact_force(C_impact,x(2*(i-1)-1),x(2*i-1))-2*zeta*(-x(2*(i+1))+x(2*i)))/m;
end
% if x(2*(i+1)-1)>=x(2*i-1)
%     dx(2*i)=(c1*(x(2*(i+1)-1)-x(2*i-1))+c2*(x(2*(i+1)-1)-x(2*i-1))^2+c3*(x(2*(i+1)-1)-x(2*i-1))^3+ ...
%         -contactforce(a,h,x(2*i-1)-x(2*(i-1)-1),E_lin,N_asp))/m;
% else
%     dx(2*i)=(c1*(x(2*(i+1)-1)-x(2*i-1))-c2*(x(2*(i+1)-1)-x(2*i-1))^2+c3*(x(2*(i+1)-1)-x(2*i-1))^3+ ...
%         -contactforce(a,h,x(2*i-1)-x(2*(i-1)-1),E_lin,N_asp))/m;
% end

for i=3:num_mass-1
    dx(2*i-1) = x(2*i);  
    if x(2*(i+1)-1)>=x(2*i-1)
        Force_spring1=c1*(x(2*(i+1)-1)-x(2*i-1));
    else
        Force_spring1=c1*(x(2*(i+1)-1)-x(2*i-1))-c2*(x(2*(i+1)-1)-x(2*i-1))^2+c3*(x(2*(i+1)-1)-x(2*i-1))^3;
    end
    if x(2*i-1)>=x(2*(i-1)-1)
        Force_spring2=-c1*(x(2*i-1)-x(2*(i-1)-1));
    else
        Force_spring2=-c1*(x(2*i-1)-x(2*(i-1)-1))+c2*(x(2*i-1)-x(2*(i-1)-1))^2-c3*(x(2*i-1)-x(2*(i-1)-1))^3;
    end
    dx(2*i)=(Force_spring1+Force_spring2-2*zeta*(-x(2*(i+1))+2*x(2*i)-x(2*(i-1))))/m;
end

% end

i=num_mass;
dx(2*i-1) = x(2*i);
% if i==2
% if 0>=x(2*i-1)
%     dx(2*i)=(c1*(0-x(2*i-1))+c2*(0-x(2*i-1))^2+c3*(0-x(2*i-1))^3+ ...
%         +impact_force(C_impact,x(2*(i-1)-1),x(2*i-1))-2*zeta*(0+x(2*i)))/m;
% else
%     dx(2*i)=(c1*(0-x(2*i-1))-c2*(0-x(2*i-1))^2+c3*(0-x(2*i-1))^3+ ...
%         +impact_force(C_impact,x(2*(i-1)-1),x(2*i-1))-2*zeta*(0+x(2*i)))/m;
% end
% else
if 0>=x(2*i-1)
    Force_spring1=c1*(0-x(2*i-1));
else
    Force_spring1=c1*(0-x(2*i-1))-c2*(0-x(2*i-1))^2+c3*(0-x(2*i-1))^3;
end
if x(2*i-1)>=x(2*(i-1)-1)
    Force_spring2=-c1*(x(2*i-1)-x(2*(i-1)-1));
else
    Force_spring2=-c1*(x(2*i-1)-x(2*(i-1)-1))+c2*(x(2*i-1)-x(2*(i-1)-1))^2-c3*(x(2*i-1)-x(2*(i-1)-1))^3;
end
dx(2*i)=(Force_spring1+Force_spring2-2*zeta*(x(2*i)-x(2*(i-1)))-2*zeta*(x(2*i)))/m;
% end


