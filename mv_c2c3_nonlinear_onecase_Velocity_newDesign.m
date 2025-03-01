clc;
clear all;
close all;


% load c2_c3_data_N20.mat
%% Initialize the required input parameters (Unit: m)
% Quadrant on which the stress-strain curve is defined
a=0.125;
k1=1045.45;
k2=1059.8;
m1=398.8e-3;
m2=400.4e-3;
c_non=a*sqrt(k1/m1);
c_lin=a*sqrt(k2/m2);

mass_physical=40e-3;
M_impactor1=mass_physical/m1;
M_impactor2=mass_physical/m2;
velocity_physical=1.37;
V_impactor1=velocity_physical/c_non;
V_impactor2=(velocity_physical+0.0)/c_lin;

% M_impactor1=0.1;M_impactor2=0.1;
% mass_physical=[M_impactor1*m1 M_impactor2*m2];
% V_impactor1=0.2202;V_impactor2=0.2202;
% velocity_physical=[V_impactor1*c_non V_impactor2*c_lin];

N = 20;                      % number of layers per samplheight
m_layer=1;
M0_impactor=N/2;
C_impact =1.2e4; %C_impact1/(a1*K_layer_Linear1)*a1^1.5;

zeta_bis=0.005;       % nondimensional viscous damping parameter
zeta_lin=0.003; 

% C3=1781.748;%2238;%199.53;
% C2=-87.261;%-101.8;%-29.11;
C3=1782.701;%2238;%199.53;
C2=-86.118;%-101.8;%-29.11;
C1=1;

% Setup the options for ode45 solver
options_ode45 = [];
options2 = odeset('RelTol',1e-10,'AbsTol',1e-10.*ones(1,2*(N+1)));

k_f=1;%max(1,1+2*nonlinear_spring_info(2)+3*nonlinear_spring_info(3));
f0 = 1/(2*pi)*sqrt(1/M0_impactor);
dt_cyc0 = 1/f0;
cycles = 5;
outputpercycle = 5000;
f = 1/(2*pi)*sqrt(k_f/M0_impactor);
dt_cyc = 1/f;
dt = dt_cyc/outputpercycle;
T = dt_cyc0*cycles;
downsample = 1e+1;
time_range = [0 T];

material_info=[m_layer M_impactor1 C_impact zeta_bis];
initialvals = zeros(2*(N+1),1);
initialvals(2) = V_impactor1;
nonlinear_spring_info = [C1 C2 C3];
[t, X_non] = ode45(@(t,x) impact_equation_of_motion_asymetric(t,x,...
            nonlinear_spring_info,material_info,0),0:dt:T, initialvals, options2);
velocity_x_non=X_non(:,2:2:2*(N+1));
displacement_x_non=X_non(:,3:2:2*(N+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialvals = zeros(2*(N+1),1); 
initialvals(2) = V_impactor2;
material_info=[m_layer M_impactor2 C_impact zeta_lin];
nonlinear_spring_info = [C1 0 0];
[t, X_lin] = ode45(@(t,x) impact_equation_of_motion_asymetric(t,x,...
    nonlinear_spring_info,material_info,0),0:dt:T, initialvals, options2);
velocity_x_lin=X_lin(:,2:2:2*(N+1));
displacement_x_lin=X_lin(:,3:2:2*(N+1));
% 



Relative_disp=displacement_x_non(:,2:N)-displacement_x_non(:,1:N-1);
Relative_disp(:,N)=0-displacement_x_non(:,N);
max_relative_disp=a*max(abs(Relative_disp))

%%
Kinetic_energy_non=1/2*(c_non*velocity_x_non).^2*m1;
Kinetic_energy_lin=1/2*(c_lin*velocity_x_lin).^2*m2;

figure;plot(t/sqrt(k2/m2),Kinetic_energy_lin(:,N/2),'r-','linewidth',1.2)
ylabel('Kinetic (J)');
xlabel('Time (s)');

for ii = 1:20
    max_kinetic_energy_lin(ii)=max(Kinetic_energy_lin(:,ii+1));
    max_kinetic_energy_non(ii)=max(Kinetic_energy_non(:,ii+1));
end
figure;plot(1:N,max_kinetic_energy_lin./max_kinetic_energy_non,'r-o', LineWidth=2,MarkerSize=6)
xlabel('Number of unit cell');
ylabel('KE ratio of Linear/Bistable');
grid on;

x5=5*ones(1,length(t));
x20=20*ones(1,length(t));

figure;plot3(t(1:50:end)/sqrt(k1/m1),x5(1:50:end),c_non*velocity_x_non(1:50:end,6),'b-', LineWidth=2)
hold on;plot3(t(1:50:end)/sqrt(k1/m1),x20(1:50:end),c_non*velocity_x_non(1:50:end,21),'r-', LineWidth=2)
% plot(t(1:50:end)/sqrt(k1/m1),0.5*m_layer*(velocity_x_non(1:50:end,21)).^2/(0.5*M_impactor1*V_impactor1^2),'k-', LineWidth=2)
xlabel('Time');
ylabel('Unit cell');
zlabel('Velocity');
grid on;

figure;plot3(t(1:50:end)/sqrt(k2/m2),x5(1:50:end),c_lin*velocity_x_lin(1:50:end,6),'b-', LineWidth=2)
hold on;plot3(t(1:50:end)/sqrt(k2/m2),x20(1:50:end),c_lin*velocity_x_lin(1:50:end,21),'r-', LineWidth=2)
% plot(t(1:50:end)/sqrt(k2/m2),0.5*m_layer*(velocity_x_lin(1:50:end,21)).^2/(0.5*M_impactor2*V_impactor2^2),'k-', LineWidth=2)
xlabel('Time');
ylabel('Unit cell');
zlabel('Velocity');
grid on;





figure;pcolor(1:N,t/sqrt(k1/m1),0.5*m_layer*(velocity_x_non(:,2:end)).^2/(0.5*M_impactor1*V_impactor1^2))
shading flat
colormap parula
clim([0.0001 1])
xlabel('Number of unit cell');
ylabel('Time (s)');
title 'Kinetic Energy vs Time'
c=colorbar;
ylim([0 1.5])
set(gca,'colorscale','log')
c.Label.String = 'KE/TE';
axis square


figure;pcolor(1:N,t/sqrt(k2/m2),0.5*m_layer*(velocity_x_lin(:,2:end)).^2/(0.5*M_impactor2*V_impactor2^2))
shading flat
colormap parula
clim([0.0001 1])
xlabel('Number of unit cell');
ylabel('Time (s)');
title 'Kinetic Energy vs Time'
c=colorbar;
ylim([0 1.5])
set(gca,'colorscale','log')
c.Label.String = 'KE/TE';
axis square



%%
[t0,location]=min(abs(t-1*sqrt(k1/m1)))

figure;
startColor = [1,0.2,0.2]; endColor=[0,0,1];
for ii = 1:20
    plot(t(1:10:location)/sqrt(k1/m1), a*displacement_x_non(1:10:location,ii),'Color',startColor + (endColor-startColor)*(ii-1) /19,linewidth=2);
    hold on;
end
xlabel('Time/s');
ylabel('Displacement/m');
% legend('1st','2nd','3rd','4th')
grid on;

figure;
startColor = [1,0.2,0.2]; endColor=[0,0,1];
for ii = 1:20
    plot(t(1:10:location)/sqrt(k1/m1), a*displacement_x_lin(1:10:location,ii),'Color',startColor + (endColor-startColor)*(ii-1) /19,linewidth=2);
    hold on;
end
xlabel('Time/s');
ylabel('Displacement/m');
% legend('1st','2nd','3rd','4th')
grid on;

figure;
startColor = [1,0.2,0.2]; endColor=[0,0,1];
for ii = 1:20
    plot(t(1:40:location)/sqrt(k1/m1), c_non*velocity_x_non(1:40:location,ii+1),'Color',startColor + (endColor-startColor)*(ii-1) /19,linewidth=2);
    hold on;
end
xlabel('Time/s');
ylabel('Displacement/m');
% legend('1st','2nd','3rd','4th')
grid on;

figure;
startColor = [1,0.2,0.2]; endColor=[0,0,1];
for ii = 1:20
    plot(t(1:40:location)/sqrt(k1/m1), c_lin*velocity_x_lin(1:40:location,ii+1),'Color',startColor + (endColor-startColor)*(ii-1) /19,linewidth=2);
    hold on;
end
xlabel('Time/s');
ylabel('Displacement/m');
% legend('1st','2nd','3rd','4th')
grid on;


