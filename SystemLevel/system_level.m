
%----SYSTEM LEVEL COST ANALYSIS----

% SUBSYSTEM 1 & 2
% ----Initial guess---- 
x0 = [0.2, 0, 0.2125, 0.2625, 0.1375, 0.0375, 0.1, 0.025, 0.955, 500000 , 512345, 0.02, 0.24];

% ----Inequality linear constraints---- 
A = [
    0 0 0 1.5 0 0 0 0 0 0 0 0 0;  
    0 0 -1 0 0 0 0 -1 0 0 0 0 0;    
    0 0 0 0 1 0 1 0.5 0 0 0 0 0;   
    0 0 0 0 -1 0 -1 0.5 0 0 0 0 0; 
    0 0 0 0 0 1 0 1.5 0 0 0 0 0;    
    1 0 0 0 0 1 0 0.5 0 0 0 0 0;    
    -1 0 0 0 0 -1 0 -0.5 0 0 0 0 0; 
    0 0 0 0 0 0 -1 1 0 0 0 0 0;    
    0 0 0 0 0 -1 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 -1 0 0
    ];

b = [
    0.3;
    -0.2;  
    0.3;
    -0.2;  
    0.075; 
    0.3;
    -0.2;
    0;
    0;
    -25000];

% ----Equality constraints----
Aeq = [0 1 0 0 0 1 0 1.5 0 0 0 0 0;   
    0 0 0 1 0 0 1 1.5 0 0 0 0 0;     
    -1 0 0 0 1 -1 1 0 0 0 0 0 0;     
    0 0 -1 0 1 0 1 -1.5 0 0 0 0 0];
beq = [0.075;0.4;0;0];

% ----Upper and lower bounds---
lb = [0 0 0 0 0 0 0 0 0.95 100000 100000 0.009 4];
ub = [Inf Inf Inf Inf Inf Inf Inf Inf 1 150000 150000 inf inf];

costperpence = ((14.37/(10.^6)))/3.6; %pence per J

%----solver----
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

format long %small number

% ----Start timer----
tic
obj = fmincon(@cost_function1,x0,A,b,Aeq,beq,lb,ub,@system_level_cons,options);
% ----End timer----
time_for_algorithm = toc;

cost_value = cost_function1(obj); 

disp(['Final Objective: ' num2str(cost_value)])
disp(['Total price for Subsystem 1 & 3 (pence): ' num2str((cost_value))])
disp(['Time (s): ' num2str(time_for_algorithm) ])

% SUBSYSTEM 2 (needs a seperate solver because it maximises)
    % subsystem 2
    
%----Parameters---- 
h = 1525; % convective heat transfer coefficient of water (W/m^2K)
T1 = 335.5; % Final temp of water (K)
T2 = 328; %Temp of water fed in (K)
kB = [5.92, 13.46, 3.46]; % Thermal conductivities (W/mK)
eB = [1420e03, 726e03, 278e03]; % Dielectric strengths (V/m)
% ----Material dependant properties---- (use Aluminia)
k = kB(2);
e = eB(2);

% ----Initial guess----
x0 = [0.00006,0.0003,0.0003,1];

% ----Inequality linear constraints---- 
A = [1 -1 0 0;0 1 -1 0];
b = [-(230/e);-0.0000045];

% ----Equality constraints----
Aeq = [];
beq = [];

% ----Upper and lower bounds---
lb = [-Inf, 0.0002, -Inf, -Inf];
ub = [Inf,Inf,0.021,2.14];

%----solver----
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
format long %small number

% ----Start timer----
tic
obj2 = fmincon(@cost_function2,x0,A,b,Aeq,beq,lb,ub,@system_level_cons2,options);
% ----End timer----
time_for_algorithm = toc;

cost_value2 = cost_function2(obj2); 

disp(['Final Objective: ' num2str(cost_value2)])
disp(['Total price for Subsystem 2 (pence): ' num2str(-(cost_value2))])
disp(['Time (s): ' num2str(time_for_algorithm) ])
disp(['Total Price for Subsystem 1, 2 & 3 (pence): ' num2str((-(cost_value2))+ (cost_value))])

% Function combining subsystem 1 & 3
function X = cost_function1(x)
    % subsystem 1
    L2 = x(1);
    L3 = x(2);
    L4 = x(3);
    L5 = x(4);
    L6 = x(5);
    R0 = x(6);
    R1 = x(7);
    D = x(8);
    k = 0.0000015;% surface roughness
    u0 = 0.00015/(pi*(D/2)^2);% initial flow velocity
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
    watt_p = 0.01;% 1 watt provides approximately 100Pa of pressure
    cost_p = watt_p*14.37/3600000;% cost per Pa at 14.7p/3.6MJ
    cost_pressure = cost_p*P*(90*3600)*(52*9.5);% cost*Pa*time in use
    volume = L_pipe*pi*(D/2)^2;% amount of water stored in the pipes
    cost_volume = 319*volume*3*(52*9.5);% the cost of this water
    cost_pipe = 400*L_pipe;% cost of the piping
    cost = cost_pressure + cost_volume + cost_pipe; % total cost function
    % subsystem 3
    rho = 1000;
    g = 9.81;
    eff = 0.7;
    Q = pi*(x(12)^2)*x(13);
    H = (x(11) - x(10))/(rho*g) + (x(12)^2)/(2*g);
    W = (x(9)*rho*g*Q*H*((((14.37/(10.^6)))/3.6)*90*60*52*9.5))/eff;
    % combine
    X = W + cost;
end

% Function for subsystem 2
function Y = cost_function2(x)
    % subsystem 2
    
    %----Parameters---- 
    h = 1525; % convective heat transfer coefficient of water (W/m^2K)
    T1 = 335.5; % Final temp of water (K)
    T2 = 328; %Temp of water fed in (K)
    kB = [5.92, 13.46, 3.46]; % Thermal conductivities (W/mK)
    eB = [1420e03, 726e03, 278e03]; % Dielectric strengths (V/m)
    
    costperpence = ((14.37/(10.^6)))/3.6; %pence per J
    
    % ----Material dependant properties----
    k = kB(2);
    e = eB(2);
    
    % function including cost
    Y = ((-(((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))))*((2*pi*(x(3)))*x(4)*(20*60)*52*9.5))*costperpence);
end