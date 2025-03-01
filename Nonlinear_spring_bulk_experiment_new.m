clear all;
% close all
clc

%Target_Bistable = [396.11,-41.064, 1.0, 0.0];

N=1;
c01=1;
c03=396.9449%201.8;%190.9159
c02=-41.1872%-29.72;%-28.1210
L=1;
x=linspace(0,0.12*L,100);
x2=linspace(0,0.12*L,20);
y=c01*x+c02*x.^2+c03*x.^3;
y2=c01*x2+c02*x2.^2+c03*x2.^3;
figure;
plot(x,y,'r-',"LineWidth",2)
hold on;plot(x2,y2,'ro',MarkerSize=8,LineWidth=2)
xlabel("Nondimensionl displacement")
ylabel("Nondimensionl force")
grid on;

c01=1;
c03=396.11;
c02=-41.064;
L=1;
x=linspace(0,0.12*L,100);
x2=linspace(0,0.12*L,20);
y=c01*x+c02*x.^2+c03*x.^3;
y2=c01*x2+c02*x2.^2+c03*x2.^3;
hold on
plot(x,y,'b-',"LineWidth",2)
hold on;plot(x2,y2,'bo',MarkerSize=8,LineWidth=2)

legend('Final design','','Target')


%%

c01=1;
c03=396.11/59^2*125^2;
c02=-41.064/59*125;
L=1;
x=linspace(0,0.05*L,100);
x2=linspace(0,0.05*L,20);
y=c01*x+c02*x.^2+c03*x.^3;
y2=c01*x2+c02*x2.^2+c03*x2.^3;
a=0.001;%0.125;
k1=1045.45;
figure;
plot(x*a*1000,y*a*k1,'k--',"LineWidth",2)
% hold on;plot(x2*a*1000,y2*a*k1,'kx',MarkerSize=10,LineWidth=2)
xlabel("Displacement (mm)")
ylabel("Force (N)")
grid on;

c03=396.9449/59^2*125^2;
c02=-41.1872/59*125;
y=c01*x+c02*x.^2+c03*x.^3;
y2=c01*x2+c02*x2.^2+c03*x2.^3;
hold on;plot(x*a*1000,y*a*k1,'r-',"LineWidth",2)

c03=1782.701;%2238;%199.53;
c02=-86.118;%-101.8;%-29.11;
y=c01*x+c02*x.^2+c03*x.^3;
y2=c01*x2+c02*x2.^2+c03*x2.^3;
hold on;plot(x*a*1000,y*a*k1,'g-',"LineWidth",2)
% hold on;plot(x2*a*1000,y2*a*k1,'rx',MarkerSize=10,LineWidth=2)

legend('Target','Final design','Fitting')
