clear
 
% Physical properties.
rhos = 2650.0;
rho = 999.7;
s = rhos/rho;
g = 9.81;
nu = 1.306e-6;
t = 10.0;
 
%nu = 1.003e-6;
%rho = 998.2;
%s = rhos/rho;
%t = 20;
 
% Problem parameters.
x = [-2:0.1:1];
d = 10.^x;
d = d*1.0e-3;
tau = zeros(size(d));
dstar = zeros(size(d));
for k = 1:length(d)
    dstar(k) = (g*(s - 1.0)/nu^2)^(1.0/3.0)*d(k);
    tau(k) = 0.30/(1.0 + 1.2*dstar(k)) + 0.055*(1.0 - exp(-0.020*dstar(k)));
    tau(k) = g*(rhos-rho)*d(k)*tau(k);
end
 
loglog(d,tau)
hold on

nu = 1.003e-6;
rho = 998.2;
s = rhos/rho;
t = 20;
 
% Problem parameters.
x = [-2:0.1:1];
d = 10.^x;
d = d*1.0e-3;
tau = zeros(size(d));
dstar = zeros(size(d));
for k = 1:length(d)
    dstar(k) = (g*(s - 1.0)/nu^2)^(1.0/3.0)*d(k);
    tau(k) = 0.30/(1.0 + 1.2*dstar(k)) + 0.055*(1.0 - exp(-0.020*dstar(k)));
    tau(k) = g*(rhos-rho)*d(k)*tau(k);
end

loglog(d,tau,'r')
hold off