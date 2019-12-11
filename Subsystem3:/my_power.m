function W = my_power(x)

%x is a 1-by-5 double;
%x(1) = beta, x(2) = P_in, x(3) = P_out, x(4) = r, x(5) = v;

rho = 1000;
g = 9.81;
eff = 0.7;
Q = pi*(x(4)^2)*x(5);
H = (x(3) - x(2))/(rho*g) + (x(4)^2)/(2*g);
W = x(1)*rho*g*Q*H/eff;

end