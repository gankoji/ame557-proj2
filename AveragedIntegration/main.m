tic
format long;
clear all;close all;clc;

%Constants
mu = 3.98600442*1e5;                    %km^3/s^2
R = 6378.137;                           %km
J2 = 0.0010826267;
BC = 1/12.741621/0.00013071*1000^2;     %kg/km^2

% Initial orbital elements
% a0 = 6657.455850;
% e0 = 0.002607;
% i0 = 42.747954;
% Omega0 = 345.325641;
% omega0 = 124.282654;
% f0 = 287.855161;
a0 = 6657.391879;
e0 = 0.002595;
i0 = 42.748082;
Omega0 = 345.325883;
omega0 = 124.412552;
f0 = 287.718564;

x0 = [a0;e0;Omega0;omega0];

Tmax = 200*24*3600;
dt = 10;
NsimMax = Tmax/dt;
xPlot = zeros(NsimMax,4);
Perigee_altitude = zeros(NsimMax,1);
tPlot = zeros(NsimMax,1);

t0=0;
%% Main Simulation Loop
for i=1:NsimMax
    
    tspan = [t0:(dt/10):t0+dt];
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t,x] = ode45(@averaged_oe_eom,tspan,x0,options,mu,BC,J2,R,i0);

    nPlot = i;    
    xPlot(i,:) = x(end,:);
    tPlot(i) = t(end);
    x0 = x(end,:)';
    t0 = t(end);
    
    Perigee_altitude(i) = xPlot(i,1)*(1-xPlot(i,2)) - R;
    
    [ i , Perigee_altitude(i) xPlot(i,1) xPlot(i,2) ]
    
    %Error Handling
    if xPlot(i,2) < 0 || isnan(xPlot(i,2));
        break
    end
end

xPlot_Fin = xPlot(1:nPlot,:);
tPlot_Fin = tPlot(1:nPlot);

figure(1)
subplot(4,1,1)
plot(tPlot_Fin,xPlot_Fin(:,1))
title('Semi-major axis')
xlabel('Time')
subplot(4,1,2)
plot(tPlot_Fin,xPlot_Fin(:,2))
title('Eccentricity')
xlabel('Time')
subplot(4,1,3)
plot(tPlot_Fin,xPlot_Fin(:,3))
title('Ascending Node (deg)')
subplot(4,1,4)
plot(tPlot_Fin,xPlot_Fin(:,4))
title('Argument of Periapsis (deg)')
toc