lb = [0.95, 100000, 100000, 0.009, 4]; %lower bound for each parameter
ub = [1, 150000, 150000, inf, inf]; %upper bound for each parameter

x_0 = [0.955, 500000 , 512345, 0.02, 0.24]; %initial value to start optimisation
A = [0,1,-1,0,0]; % Linear contraints on Pin and Pout
b = -25000; % Linear constraints

Aeq = [];
Beq = [];

nonlcon = @my_nonlcon; %nonlinear constraint representing Q = pi*(r^2)*v >= 0.001

option = optimoptions(@fmincon, 'Algorithm', 'sqp', 'MaxIterations', 1000);

%if option not activated then algorithm used for fmincon =
%inter-point method (default);
%if option activated then algorithm used = Sequential
%Quadratic Programming (SQP).

min = fmincon(@my_power,x_0,A,b,Aeq,Beq,lb,ub,nonlcon,option); % applying fmincon
pump_power = my_power(min); % putting min value 

disp('Solution: ')
disp(['Optimized power: ' num2str(pump_power)])
disp(['Optimized beta value: ' num2str(min(1,1))])
disp(['Optimized inlet pressure: ' num2str(min(1,2))])
disp(['Optimized outlet pressure: ' num2str(min(1,3))])
disp(['Optimized radius: ' num2str(min(1,4))])
disp(['Optimized fluid velocity: ' num2str(min(1,5))])

