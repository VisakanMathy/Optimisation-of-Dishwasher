function [c,ceq] = system_level_cons2(x)

%----Parameters----
m = 13; % mass (kg)
c = 3976.7; % specific heat capacity water (J/kgK)
T1 = 335.5; % Final temp of water (K)
T2 = 328; %Temp of water fed in (K)
t = 20*60; % time (s)
V = 230; % Voltage (V)

%----Parameter calculations----
P = (m*c*(T1-T2))/t; % Power (W)
I = P/V; % Current (I)
R = P/(I^2); % Resistance required (R)

%----Material dependant variables----
r = 1.3.*10.^-6; % Resistivity of material A (rho)

%----Inequality constraints----
c = R - ((r.*x(4))/(pi.*(x(1).^2)));

%----Equality constraints----
ceq = [];

end