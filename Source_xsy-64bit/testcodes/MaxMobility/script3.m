clear
 
% Physical properties.
rhos = 2650.0;
rho = 998.60;
s = rhos/rho;
g = 9.81;
nu = 1.0678e-6;
t = 18.0;
 
 
% Problem parameters.
x = [-1:0.1:3];
dstar = 10.^x;
theta = zeros(size(dstar));
for k = 1:length(dstar)
    theta(k) = 0.30/(1.0 + 1.2*dstar(k)) + 0.055*(1.0 - exp(-0.020*dstar(k)));
end
 
loglog(dstar,theta)
