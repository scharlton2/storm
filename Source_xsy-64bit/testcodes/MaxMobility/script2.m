clear
 
% Physical properties.
rhos = 2650.0;
rho = 998.60;
s = rhos/rho;
g = 9.81;
nu = 1.0678e-6;
t = 18.0;

d = 2.0;
d = d*0.001;
dstar = (g*(s - 1.0)/nu^2)^(1.0/3.0)*d;
tau = 0.30/(1.0 + 1.2*dstar) + 0.055*(1.0 - exp(-0.020*dstar));
tau = g*(rhos-rho)*d*tau

 
