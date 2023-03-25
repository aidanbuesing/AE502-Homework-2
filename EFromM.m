function E = EFromM(e,M)

% E = EFromM(M,e)
% this function is a bisection solver to find the eccentric anomaly
% of an orbit when given the Mean Anomaly and Eccentricity
%
% Inputs:
% M = Mean Anomaly of the orbit
% e = eccentricity
%
% Output:
% E = eccentric Anomaly

E = 0;
diff = 1;
tol = 1e-3;

while(abs(diff)>tol)
    diff = E - e*sin(E) - M;
    E = E +.001;
end