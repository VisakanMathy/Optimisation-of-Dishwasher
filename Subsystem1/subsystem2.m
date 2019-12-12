clear all
clc

% A.x <= b
A = [0,0,1,0,0,0,0,1.5;   % eq32
    0,0,-1,0,0,0,0,-1;    % eq31
    0,0,0,0,1,0,1,0.5;    % eq30
    0,0,0,0,-1,0,-1,0.5;  % eq29
    0,0,0,0,0,1,0,1.5;    % eq23
    1,0,0,0,0,1,0,0.5;    % eq22
    -1,0,0,0,0,-1,0,-0.5; % eq21
    0,0,0,0,0,0,-1,1;     % eq20
    0,0,0,0,0,-1,0,1];     % eq19

B = [0.3;
    -0.2;  % eq31
    0.3;
    -0.2;  % eq29
    0.075; % eq23
    0.3;
    -0.2;
    0;
    0];

Aeq = [0,1,0,0,0,1,0,1.5;   % eq25
    0,0,0,1,0,0,1,1.5;      % eq26
    -1,0,0,0,1,-1,1,0;      % eq27
    0,0,-1,0,1,0,1,-1.5];    % eq28

Beq = [0.075;0.4;0;0];

lb = [0,0,0,0,0,0,0,0];     % can't have negative lengths of pipe
ub = [];                    % no upper bounds on the lengths, not needed

x0 = [0.2, 0, 0.2125, 0.2625, 0.1375, 0.0375, 0.1, 0.025];  %Initial guess, to provide some numbers

nonlcon = [];

% Uncomment below to run sqp algorithm
%
sqp_option = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

tic
sqp_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, sqp_option);
sqp_time = toc

sqp_pennies = cost_function(sqp_solution)/100
%} 


% Uncomment below to run Active_set Algorithm
%
AS_option = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'active-set');

tic
AS_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, AS_option);
AS_time = toc

AS_pennies = cost_function(AS_solution)/100
%}

% Uncomment below to run interior-point algorithm
%
interior_option = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
tic
interior_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, interior_option);
interior_time = toc

interior_pennies = cost_function(interior_solution)/100
%}


% Uncomment below to run sqp-legacy algorithm
%
legacy_option = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp-legacy');
tic
legacy_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, legacy_option);
legacy_time = toc

legacy_pennies = cost_function(legacy_solution)/100
%}

% Uncomment to compare times (must uncomment relative algorithms above)
%
AS_time
interior_time
sqp_time
legacy_time

% Below shows they all come to the same value
AS_pennies
interior_pennies
sqp_pennies
legacy_pennies
%}

values = [AS_solution; interior_solution; sqp_solution; legacy_solution];

function cost = cost_function(x)
    L2 = x(1);
    L3 = x(2);
    L4 = x(3);
    L5 = x(4);
    L6 = x(5);
    R0 = x(6);
    R1 = x(7);
    D = x(8);
    k = 0.0000015;       % surface roughness
    u0 = 0.00015/(pi*(D/2)^2);     % initial flow velocity
    u1 = u0/2;  
    Re0 = u0*D/(10^-6);
    Re1 = u1*D/(10^-6);
    f0 = 0.25/((log10((k/(3.7*D))+(5.74/(Re0^0.9))))^2);
    f1 = 0.25/((log10((k/(3.7*D))+(5.74/(Re1^0.9))))^2);
    Le0 =22.2126*(Re0*(D/R0)^2)^0.7888 * Re0^-0.71438;
    Leq0 = Le0*D + pi*R0/2;
    L0 = L2 + Leq0 + L3;
    P0 = 4*f0*(L0/D)*(1/2)*u0^2;
    Le1 = 22.2126*(Re1*(D/R1)^2)^0.7888 * Re1^-0.71438;
    Leq2 = 0.6*D + 3*D;
    Leq3 = 0.2*D + 3*D;
    Leq1 = Le1*D + pi*R1/2;
    L1 = Leq2 + L4 + Leq3 + L5 + Leq1 + L6;
    P0 = 4*f0*(L0/D)*(1/2)*997*u0^2;
    P1 = 4*f1*(L1/D)*(1/2)*997*u1^2;
    P = P0 + P1 + 2696;
    L_pipe = L2 + L3 + L4 + L5 + L6 + (R0+R1)*pi/2 + 4*D;
    watt_p = 0.01;     % 1 watt provides approximately 100Pa of pressure
    cost_p = watt_p*14.37/3600000;     % cost per Pa at 14.7p/3.6MJ
    cost_pressure = cost_p*P*(90*3600)*(52*9.5);  % cost*Pa*time in use
    volume = L_pipe*pi*(D/2)^2                  % amount of water stored in the pipes
    cost_volume = 319*volume*3*(52*9.5);          % the cost of this water
    cost_pipe = 400*L_pipe;                      % cost of the piping
    cost = cost_pressure + cost_volume + cost_pipe; % total cost function
end