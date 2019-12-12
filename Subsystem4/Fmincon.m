close all
clear all
objectiveFn = @objective;
x0=[50,8,0.0025,0.05,0.05,6];


A = [];
b=[];
Aeq = [];
beq = [];
lb=[0,4,0.0015,0.02,0,2];
ub = [90,10,0.003,0.2,0.1,20];
nonlincon = @nlcon;

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');

[x, fval]=fmincon(objectiveFn,x0,A,b,Aeq,beq,lb,ub,nonlincon,options)

%problem = createOptimProblem('fmincon','x0',x0,'objective',objectiveFn,'nonlcon',nonlincon,'lb',lb,'ub',ub);
%x = run(GlobalSearch,problem);
%options = optimoptions('fmincon','Display','iter','Algorithm','sqp')
%[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options)
%objectivePlot(x);

%%

 x_initial = [];
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
    if i == 4
        x_initial = x;
    end
end
figure
plot(numberOfHoles,omega)
title('Number of holes effect on omega')
xLabel('Number of holes')
yLabel('Angular velocity (rad/s)')


%% Sensitivity analysis of speed

for j = 2:19
    x=x_initial;
lb=[0,4,0.0015,0.02,0,j];
ub = [90,4,0.003,0.2,0.1,j+1];
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
speed(j-1) = j+1;
omega_speed(j-1)=fval;
end
figure
plot(speed,omega_speed)
title("Jet velocity's effect on omega")
xLabel('Jet velocity (m/s)')
yLabel('Angular velocity (rad/s)')

%% Sensitivity analysis of radius of hole

for j = 1:15
    x=x_initial;
lb=[0,4,(14+j)*0.0001,0.02,0,2];
ub = [90,4,(14+j)*0.0001,0.2,0.1,20];
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
radius(j) = (14+j)*0.0001;
omega_radius(j)=fval;
end
figure
plot(radius,omega_radius)
title("Radius of hole's effect on omega")
xLabel('Radius (m)')
yLabel('Angular velocity (rad/s)')


%% Sensitivity analysis of angles

for j = 1:8
    x=x_initial;
lb=[j*10,4,0.0015,0.02,0,2];
ub = [j*10,4,0.003,0.2,0.1,20];
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options);
angle(j) = j*10;
omega_angle(j)=fval;
end
figure
plot(angle,omega_angle)
title("Angle of hole's effect on omega")
xLabel('Angle (deg)')
yLabel('Angular velocity (rad/s)')


    
