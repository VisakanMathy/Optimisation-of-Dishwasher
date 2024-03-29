clear all
clc

% A.x <= b
A = [0,0,1,0,0,0,0,1.5;   % g9
    0,0,-1,0,0,0,0,-1;    % g8
    0,0,0,0,1,0,1,0.5;    % g7
    0,0,0,0,-1,0,-1,0.5;  % g6
    0,0,0,0,0,1,0,1.5;    % g5
    1,0,0,0,0,1,0,0.5;    % g4
    -1,0,0,0,0,-1,0,-0.5; % g3
    0,0,0,0,0,0,-1,1;     % g2
    0,0,0,0,0,-1,0,1];    % g1

B = [0.3;  % g9
    -0.2;  % g8
    0.3;   % g7
    -0.2;  % g6
    0.075; % g5
    0.3;   % g4
    -0.2;  % g3
    0;     % g2
    0];    % g1

Aeq = [0,1,0,0,0,1,0,1.5;   % h1
    0,0,0,1,0,0,1,1.5;      % h2
    -1,0,0,0,1,-1,1,0;      % h3
    0,0,-1,0,1,0,1,-1.5];   % h4

Beq = [0.075;0.4;0;0];

lb = [0,0,0,0,0,0,0,0];     % can't have negative lengths of pipe
ub = [];                    % no upper bounds on the lengths, not needed

x0 = [0.2, 0, 0.2125, 0.2625, 0.1375, 0.0375, 0.1, 0.025];  %Initial guess, to provide some numbers

nonlcon = [];

% Uncomment below to run sqp algorithm
%
sqp_option = optimoptions('fmincon', 'Algorithm', 'sqp');

tic
sqp_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, sqp_option);
sqp_time = toc

sqp_pennies = cost_function(sqp_solution)/100
%} 


% Uncomment below to run Active_set Algorithm
%
AS_option = optimoptions('fmincon', 'Algorithm', 'active-set');

tic
AS_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, AS_option);
AS_time = toc

AS_pennies = cost_function(AS_solution)/100
%}

% Uncomment below to run interior-point algorithm
%{
interior_option = optimoptions('fmincon', 'Algorithm', 'interior-point');
tic
interior_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, interior_option);
interior_time = toc

interior_pennies = cost_function(interior_solution)/100
%}


% Uncomment below to run sqp-legacy algorithm
%{
legacy_option = optimoptions('fmincon', 'Algorithm', 'sqp-legacy');
tic
legacy_solution = fmincon(@cost_function, x0, A, B, Aeq, Beq, lb, ub, nonlcon, legacy_option);
legacy_time = toc

legacy_pennies = cost_function(legacy_solution)/100
%}

% Uncomment to compare times (must uncomment relative algorithms above)
%{
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

%values = [AS_solution; interior_solution; sqp_solution; legacy_solution];

function cost = cost_function(x)
    L2 = x(1);
    L3 = x(2);
    L4 = x(3);
    L5 = x(4);
    L6 = x(5);
    R0 = x(6);
    R1 = x(7);
    D = x(8);
    k = 0.0000015;               % surface roughness
    u0 = 0.00015/(pi*(D/2)^2)    % initial flow velocity
    u1 = u0/2;                   % flow halves after the tee-piece
    Re0 = u0*D/(10^-6);          % Reynolds number
    Re1 = u1*D/(10^-6);          % Reynolds number after tee-piece
    f0 = 0.25/((log10((k/(3.7*D))+(5.74/(Re0^0.9))))^2);    % Swammee-Jain equation of Darcy Friction Factor
    f1 = 0.25/((log10((k/(3.7*D))+(5.74/(Re1^0.9))))^2);    % Swammee-Jain equation of Darcy Friction Factor
    Le0 =22.2126*(Re0*(D/R0)^2)^0.7888 * Re0^-0.71438;      % Reference 2 from main report
    Leq0 = Le0*D + pi*R0/2;                                 % Total equivalent length of first bend
    L0 = L2 + Leq0 + L3;                                    % Total equivalent length before tee-piece
    Le1 = 22.2126*(Re1*(D/R1)^2)^0.7888 * Re1^-0.71438;     % Reference 2 from main report
    Leq2 = 0.6*D + 3*D;                                     % Total equivalent length branching off tee-piece
    Leq3 = 0.2*D + 3*D;                                     % Total equivalent length through tee-peice
    Leq1 = Le1*D + pi*R1/2;                                 % Total equivalent length of upper bend 
    L1 = Leq2 + L4 + Leq3 + L5 + Leq1 + L6;                 % Total equivalent length after and including tee-piece
    P0 = 4*f0*(L0/D)*(1/2)*997*u0^2;                        % Pressure drop due to piping before tee-piece
    P1 = 4*f1*(L1/D)*(1/2)*997*u1^2;                        % Pressure drop due to piping during and after tee-piece
    P = P0 + P1 + 2696;                                     % Total pressure drop in the system, including due to height gain
    L_pipe = L2 + L3 + L4 + L5 + L6 + (R0+R1)*pi/2 + 4*D;   % Total length of physical pipe in the system
    watt_p = 0.01;                                          % 1 watt provides approximately 100Pa of pressure
    cost_p = watt_p*14.37/3600000;                          % cost per Pa at 14.7p/3.6MJ
    cost_pressure = cost_p*P*(90*3600)*(52*9.5);            % cost*Pa*time in use
    volume = L_pipe*pi*(D/2)^2                              % amount of water stored in the pipes
    cost_volume = 319*volume*3*(52*9.5);                    % the cost of this water
    cost_pipe = 400*L_pipe;                                 % cost of the piping
    cost = cost_pressure + cost_volume + cost_pipe;         % total cost function
end