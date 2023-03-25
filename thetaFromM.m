function theta = thetaFromM(M,e)

% theta = thetaFromM(M,e)
% this function is a bisection solver to find the true anomaly
% of an orbit when given the Mean Anomaly and Eccentricity
%
% Inputs:
% M = Mean Anomaly of the orbit
% e = eccentricity
%
% Output:
% theta = true anomaly

theta = 0;
Et = 0;
diff = 1;
tol = 1e-3;

while(abs(diff)>tol)
    diff = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2))-(e*sqrt(1-e^2)*sin(theta))/(1+e*cos(theta))-M;
    theta = theta + .001;
end

