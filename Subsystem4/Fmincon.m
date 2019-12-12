close all
clear all
objectiveFn = @objective;
x0=[20,8,0.0025,0.05,0.05,6];

A = [];
b=[];
Aeq = [];
beq = [];
lb=[0,4,0.0015,0.02,0,2];
ub = [90,10,0.003,0.2,0.1,20];
nonlincon = @nlcon;

problem = createOptimProblem('fmincon','x0',x0,'objective',objectiveFn,'nonlcon',nonlincon,'lb',lb,'ub',ub);
x = run(GlobalSearch,problem);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
objectivePlot(x);
x_optimal = x;

%% Number of holes with fixed integers
for i = 4:10
    x0=[50,8,0.0025,0.05,0.05,6];
    lb=[0,i,0.0015,0.02,0,2];
    ub = [90,i,0.003,0.2,0.1,20];
    numberOfHoles(i-3)=i;
    problem = createOptimProblem('fmincon','x0',x0,'objective',objectiveFn,'nonlcon',nonlincon,'lb',lb,'ub',ub);
    x = run(GlobalSearch,problem);
    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
    omega(i-3)=fval;
end
figure
plot(numberOfHoles,omega)
title('The effect of number of holes on omega')
xlabel('Number of holes')
ylabel('Angular velocity (rad/s)')


%% Testing Non Gradient methods
ip = [0,0];
sqp = [0,0];
as = [0,0];
runs=10;
for i =1:runs
x0 = [randi([0 90]),4,randi([15,30])*0.0001,randi([2,20])*0.001,0.1,randi([2,20])];
lb=[0,4,0.0015,0.02,0,2];
ub = [90,4,0.003,0.2,0.1,20];
tic;
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x, fval]=fmincon(objectiveFn,x0,A,b,Aeq,beq,lb,ub,nonlincon,options)
a = toc;
ip(1) = ip(1) + a;
ip(2) = ip(2) + fval;
tic
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x, fval]=fmincon(objectiveFn,x0,A,b,Aeq,beq,lb,ub,nonlincon,options)
a = toc;
sqp(1) = sqp(1) + a;
sqp(2) = sqp(2) + fval;
tic
options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
[x, fval]=fmincon(objectiveFn,x0,A,b,Aeq,beq,lb,ub,nonlincon,options);
a = toc;
as(1) = as(1) + a;
as(2) = as(2) + fval;
end
ip/runs;
sqp/runs;
as/runs;
%% Post global search algorithm values
lb=[0,4,0.0015,0.02,0,2];
ub = [90,4,0.003,0.2,0.1,20];
ip = [0,0];
sqp = [0,0];
as = [0,0];
runs=10;
for i =1:runs
x0 = [randi([0 90]),4,randi([15,30])*0.0001,randi([2,20])*0.001,0.1,randi([2,20])];
problem = createOptimProblem('fmincon','x0',x0,'objective',objectiveFn,'nonlcon',nonlincon,'lb',lb,'ub',ub);
x = run(GlobalSearch,problem);
tic
options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
a = toc;
ip(1) = ip(1) + a;
ip(2) = ip(2) + fval;
tic
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
a = toc;
sqp(1) = sqp(1) + a;
sqp(2) = sqp(2) + fval;
tic
options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
a = toc;
as(1) = as(1) + a;
as(2) = as(2) + fval;
end
ip/runs;
sqp/runs;
as/runs;
%% Sensitivity analysis of speed
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
for j = 2:19
    x = x_optimal;
    lb=[0,4,0.0015,0.02,0,j];
    ub = [90,4,0.003,0.2,0.1,j+1];
    [x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
    speed(j-1) = j+1;
    omega_speed(j-1)=fval;
end
figure
plot(speed,omega_speed)
title("Jet velocity's effect on omega")
xlabel('Jet velocity (m/s)')
ylabel('Angular velocity (rad/s)')

%% Sensitivity analysis of radius of hole
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
for j = 1:15
    x=x_optimal;
    lb=[0,4,(14+j)*0.0001,0.02,0,2];
    ub = [90,4,(14+j)*0.0001,0.2,0.1,20];
    [x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
    radius(j) = (14+j)*0.0001;
    omega_radius(j)=fval;
end
figure
plot(radius,omega_radius)
title("Radius of hole's effect on omega")
xlabel('Radius (m)')
ylabel('Angular velocity (rad/s)')


%% Sensitivity analysis of angles
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
for j = 1:10
    x=x_optimal;
    lb=[(j-1)*10,4,0.0015,0.02,0,2];
    ub = [(j-1)*10,4,0.003,0.2,0.1,20];
    [x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
    angle(j) = (j-1)*10;
    omega_angle(j)=fval;
end
figure
plot(angle,omega_angle)
title("Angle of hole's effect on omega")
xlabel('Angle (deg)')
ylabel('Angular velocity (rad/s)')