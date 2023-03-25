%% AE502 hwk 2
%% Aidan Buesing
clc

%% Question 3

%Earth Constants
J2 = 0.00108;
mu = 396800;
R = 6370;

%initializing array to store propogation
OEs = zeros(8,100);


%initial conditions
a = 26600;
i = 1.10654;
e = 0.74;
w = 5*pi/180;
omega = 90*pi/180;
M = 10*pi/180; %need to get theta from this...
days = 100;

theta = thetaFromM(M,e); %This function is just a bisection solver for theta 

%need h
h = sqrt(mu*(1+e)*a*(1-e));

for dt=1:1:days
    M = M + dt*sqrt(mu/a^3);
    E = EFromM(e,M); %Also a bisection solver. Used two different funtions
    % for this and theta as we only needed theta once but now we need 
    % E every time
    %need r
    r = a*(1-e*sin(E));
    % Curtis uses "u" to make the equations look a little neater
    u = theta + w;
    %Curtis Equations
    dh = -3/2 * (J2*mu*R^2)/r^3 *sin(i)^2*sin(2*u);
    de = 3/2 * (J2*mu*R^2)/(r^3*h) * (((h^2/(mu*r)*sin(theta)*(3*sin(i)^2*sin(u)^2-1))-(sin(2*u)*sin(i)^2)*((3+e*cos(theta))*cos(theta)+e)));
    dtheta = h/r^2 + 3/2 * (J2*mu*R^2)/(r^3*h*e) * ((h^2/(mu*r)*cos(theta)*(3*sin(i)^2*sin(u)^2-1))+((2+e*cos(theta))*sin(2*u)*sin(i)^2*sin(theta)));
    domega = -3*(J2*mu*R^2)/(r^3*h)*sin(u)^2*cos(i);
    di = -3/4 * (J2*mu*R^2)/(r^3*h)*sin(2*u)*sin(2*i);
    dw = 3/2 * (J2*mu*R^2)/(r^3*h*e)*((h^2/(mu*r)*cos(theta)*(1-3*sin(i)^2*sin(u)^2))-((2+e*cos(theta))*sin(2*u)*sin(i)^2*sin(theta))+(2*e*cos(i)^2*sin(u)^2));
    
    %Using what we found from Curtis' Equations to propogate the orbit
    h = h + dh*dt;
    e = e + de*dt;
    theta = theta + dtheta*dt;
    omega = omega + domega*dt;
    i = i + di*dt;
    w = w + dw*dt;
    a = h^2 / (mu*(1-e^2));
    
    % Storing our elements so we can plot them
    OEs(:,dt)=[h,e,theta,omega,i,w,a,M];
end

dt =linspace(1,days,days);
subplot(4,2,1)
plot(dt,OEs(1,:));
title('Spec. Ang. Mom. (h)')
xlabel('time (days)')
ylabel('(km^2/s)')
subplot(4,2,2)
plot(dt,OEs(2,:));
title('Eccentricity')
xlabel('time (days)')
subplot(4,2,3)
plot(dt,OEs(3,:));
title('True Anomaly');
xlabel('time (days)')
ylabel('(rad)')
subplot(4,2,4)
plot(dt,OEs(4,:));
title('Right of the Ascending Node')
xlabel('time (days)')
ylabel('(rad)')
subplot(4,2,5)
plot(dt,OEs(5,:));
title('Inclination')
xlabel('time (days)')
ylabel('(rad)')
subplot(4,2,6)
plot(dt,OEs(6,:));
title('Argument of Perigee')
xlabel('time (days)')
ylabel('(rad)')
subplot(4,2,7)
plot(dt,OEs(7,:));
title('Semi-Major Axis')
xlabel('time (days)')
ylabel('(km)')
% I know hwk said to not plot M, but why not? I had an extra space left 
% on my subplot anyways
subplot(4,2,8)
plot(dt,OEs(8,:));
title('Mean Anomaly')
xlabel('time (days)')
ylabel('(rad)')

