close all
objectiveFn = @objective;
x0=[20,7,0.0025,0.05,0.05,6];


A = [];
b=[];
Aeq = [];
beq = [];
lb=[0,3,0.0015,0.02,0,2];
ub = [90,10,0.003,0.2,0.1,20];
nonlincon = @nlcon;

options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');

[x, fval]=fmincon(objectiveFn,x0,A,b,Aeq,beq,lb,ub,nonlincon,options)
problem = createOptimProblem('fmincon','x0',x0,'objective',objectiveFn,'nonlcon',nonlincon,'lb',lb,'ub',ub);
x = run(GlobalSearch,problem);
[x, fval]=fmincon(objectiveFn,x,A,b,Aeq,beq,lb,ub,nonlincon,options)
objectivePlot(x);


