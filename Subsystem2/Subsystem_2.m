clc
clear


%% Initial parametric analysis

disp([ ' -------------------- Initial parametric analysis -------------------- '])

% ---- rA changing ----
rA = 0:0.0005:0.005;
% parameters - average vales based on upper and lower bounds 
rB = ((0.005+0.02)/2);
rC = ((0.02+0.021)/2);
L = ((1.57+2.14)/2); 
kB = 5.92; % Using material Mica
kC = 18.68; % Using material Stainless steel
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(1)
plot(rA,y)
xlabel('rA (m)')
ylabel('Heat flux (W/m^2)')
title('Effect of rA changing on objective')

clc
clear

% ---- rB changing ----
rB = 0.005:0.001:0.02;
% parameters - average vales based on upper and lower bounds 
rA = ((0+0.005)/2);
rC = ((0.02+0.021)/2);
L = ((1.57+2.14)/2); 
kB = 5.92; % Using material Mica
kC = 18.68; % Using material Stainless steel
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(2)
plot(rB,y)
xlabel('rB (m)')
ylabel('Heat flux (W/m^2)')
title('Effect of rB changing')

clc
clear

% ---- rC changing ----
rC = 0.02:0.001:0.021;
% parameters - average vales based on upper and lower bounds 
rA = ((0+0.005)/2);
rB = ((0.005+0.02)/2);
L = ((1.57+2.14)/2); 
kB = 5.92; % Using material Mica
kC = 18.68; % Using material Stainless steel
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(3)
plot(rC,y)
xlabel('rC (m)')
ylabel('Heat flux (W/m^2)')
title('Effect of rC changing')

clc
clear

% ---- L changing ----
L = 1.57:0.1:2.14;
% parameters - average vales based on upper and lower bounds 
rA = ((0+0.005)/2);
rB = ((0.005+0.02)/2);
rC = ((0.02+0.021)/2);
kB = 5.92; % Using material Mica
kC = 18.68; % Using material Stainless steel
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(4)
plot(L,y)
xlabel('L (m)')
ylabel('Heat flux (W/m^2)')
title('Effect of L changing')

clc
clear

% ---- L changing ----
kB = 0:0.1:10;
% parameters - average vales based on upper and lower bounds 
rA = ((0+0.005)/2);
rB = ((0.005+0.02)/2);
rC = ((0.02+0.021)/2);
L = ((1.57+2.14)/2);
kC = 18.68; % Using material Stainless steel
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(5)
plot(kB,y)
xlabel('kB (W/mK)')
ylabel('Heat flux (W/m^2)')
title('Effect of kB changing')

clc
clear

% ---- kC changing ----
kC = 0:0.1:20;
% parameters - average vales based on upper and lower bounds 
rA = ((0+0.005)/2);
rB = ((0.005+0.02)/2);
rC = ((0.02+0.021)/2);
L = ((1.57+2.14)/2);
kB = 5.92; % Using material Mica
% heat flux  
y = (-((7.5)./(((log(rB./rA))./(2.*pi.*L.*kB))+((log(rC./rB))./(2.*pi.*L.*kC))+(1./(1525.*2.*pi.*L.*rC)))));
% plot
figure(6)
plot(kC,y)
xlabel('kC (W/mK)')
ylabel('Heat flux (W/m^2)')
title('Effect of kC changing')

clc
clear


%% Material Data - for 3 materials

% Each list is in order of: Mica, Aluminia, Mullite

% Materials
materials = ["Mica", "Aluminia", "Mullite"];
% Thermal conductivities (W/mK)
kB = [5.92, 13.46, 3.46]; 
% Dielectric strengths (V/m)
eB = [1420e03, 726e03, 278e03]; 
% Densities (kg/m^3) (averaged from min and max values)
density = [((2599+3210)/2),((3432+3515)/2),((2698+2989)/2)];
% CO2 impacts of materials (J/kg) (averaged from min and max values)
CO2_impact = [(((21189e04)+(2349e05))/2), (((4954e04)+(5466e04))/2), (((5280e04)+(5838e04))/2)]; 



%% Fmincon: SQP, Active-set, Interior-point - for 3 materials

% Empty lists add values into
% List for optimsed objective functions
fval_list = [];
% List for optimised variables
variable_list = [];

disp([ ' -------------------- Fmincon for 2 algorithms using 3 different materials -------------------- '])

for i = 1:3 % 3 times for each algorithm
    
    algorithms = ["sqp", "active-set", "interior-point"];
    
    for j = 1:length(materials) % 3 times for each material

        %----Parameters----
        h = 1525; %convective heat transfer coefficient of water (W/m^2K)
        T1 = 335.5; % Final temp of water (Temp of outide layer A) (K)
        T2 = 328; %Temp of water fed in (Temp of outide layer C)(K)
        % ----Material dependant properties----
        k = kB(j);
        e = eB(j);
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
        % ----Non linear constraints----
        nonlcon = @cons2;
        % ----Objective function----
        % VARIABLES - x(1) = rA, x(2) = rB, x(3) = rC, x(4) = L
        fun = @(x) (-((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))));
        %----solver----
        options = optimoptions('fmincon','Display','iter','Algorithm',algorithms(i));
        format long %small number
        % ----Start timer----
        tic
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        % ----End timer----
        time_for_algorithm = toc;
        % ----Add values to list----
        fval_list = [fval_list, fval];
        variable_list = [variable_list, [x(1),x(2),x(3),x(4)]];
        % ----Display values----
        disp([materials(j)])
        disp([algorithms(i)])
        disp(['Initial Objective: ' num2str(fun(x0))])
        disp(table(x(1),x(2),x(3),x(4),'VariableNames',{'rA', 'rB', 'rC', 'Length'}))
        disp(['Final Objective: ' num2str(fun(x))])
        disp(['Heat flux (W/m^2): ' num2str(-(fun(x)))])
        disp(['Time (s): ' num2str(time_for_algorithm) ])
        disp(['% improvement of objective: ' num2str(((-fun(x)-(-fun(x0)))/(-fun(x0)))*100)])
    end
end

%% Global search - for 3 materials

disp([ ' -------------------- Global search for 3 materials -------------------- '])

%----List for optimised variables----
variable_list_global_search = [];

for j = 1:length(materials) % 3 times for each material
    
    % ----Material dependant properties----
    k = kB(j);
    e = eB(j);    
    %----Parameters----
    h = 1525; %convective heat transfer coefficient of water (W/m^2K)
    T1 = 335.5; % Final temp of water (Temp of outide layer A)(K)
    T2 = 328; %Temp of water fed in (Temp of outide layer C)(K)
    % ----Objective function----
    % VARIABLES - x(1) = rA, x(2) = rB, x(3) = rC, x(4) = L
    fun = @(x) (-((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))));
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
    % ----Non linear constraints----
    nonlcon = @cons2;
    rng default
    % ----Start timer----
    tic
    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',nonlcon,'lb',lb,'ub',ub);
    x = run(GlobalSearch,problem);
    % ----End timer----
    time_for_algorithm = toc;
    %----Output variables from each material---
    variable_list_global_search = [variable_list_global_search, [x(1),x(2),x(3),x(4)]];
    disp([materials(j)])
    disp(table(x(1),x(2),x(3),x(4),'VariableNames',{'rA', 'rB', 'rC', 'Length'}))
    disp(['Final Objective: ' num2str(fun(x))])
    disp(['Heat flux (W/m^2): ' num2str(-(fun(x)))])
    disp(['Time (s): ' num2str(time_for_algorithm) ])
    disp(['% improvement of objective: ' num2str(((-fun(x)-(-fun(x0)))/(-fun(x0)))*100)])

end

%% Post optimality parametric analysis
% Using MICA values

disp([ ' -------------------- Post optimality parametric analysis -------------------- '])

% ---- Change in TEMPERATURE ----

 % ----Material dependant properties----
 k = kB(1);
 e = eB(1);    
 %----Parameters----
 h = 1525; % convective heat transfer coefficient of water (W/m^2K)
 T1 = 335.5; % Final temp of water (K)
 T2 = 280; %Temp of water fed in (K)
 % ----Objective function----
 % VARIABLES - x(1) = rA, x(2) = rB, x(3) = rC, x(4) = L
 fun = @(x) (-((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))));
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
 % ----Non linear constraints----
 nonlcon = @cons3;
 rng default
 problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',nonlcon,'lb',lb,'ub',ub);
 x = run(GlobalSearch,problem);
 % output variables from each material
 variable_list_global_search = [variable_list_global_search, [x(1),x(2),x(3),x(4)]];
 disp(table(x(1),x(2),x(3),x(4),'VariableNames',{'rA', 'rB', 'rC', 'Length'}))
 disp(['Final Objective: ' num2str(fun(x))])
 disp(['Heat flux (W/m^2): ' num2str(-(fun(x)))])
    
% ---- TIME FOR DISHWASHER TO RUN changing ----

 % ----Material dependant properties----
 k = kB(1);
 e = eB(1);    
 %----Parameters----
 h = 1525; %convective heat transfer coefficient of water %UNIT
 T1 = 335.5; % Temp of water final (K)
 T2 = 328; %Temp of water fed in (K)
 % ----Objective function----
 % VARIABLES - x(1) = rA, x(2) = rB, x(3) = rC, x(4) = L
 fun = @(x) (-((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))));
 % ----Initial guess----
 x0 = [0.00006,0.0003,0.0003,1];
 % ----Inequality linear constraints---- = g1, g4, g5
 A = [1 -1 0 0;0 1 -1 0];
 b = [-(230/e);-0.0000045];
 % ----Equality constraints----
 Aeq = [];
 beq = [];
 % ----Upper and lower bounds---
 lb = [-Inf, 0.0002, -Inf, -Inf];
 ub = [Inf,Inf,0.021,2.14];
 % ----Non linear constraints----
 nonlcon = @cons4;
 rng default
 problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',nonlcon,'lb',lb,'ub',ub);
 x = run(GlobalSearch,problem);
 % output variables from each material
 variable_list_global_search = [variable_list_global_search, [x(1),x(2),x(3),x(4)]];
 disp(table(x(1),x(2),x(3),x(4),'VariableNames',{'rA', 'rB', 'rC', 'Length'}))
 disp(['Final Objective: ' num2str(fun(x))])
 disp(['Heat flux (W/m^2): ' num2str(-(fun(x)))])
 
%---- MASS OF WATER changing ----

 % ----Material dependant properties----
 k = kB(1);
 e = eB(1);    
 %----Parameters----
 h = 1525; %convective heat transfer coefficient of water %UNIT
 T1 = 335.5; % Final temp of water (K)
 T2 = 328; %Temp of water fed in (K)
 % ----Objective function----
 % VARIABLES - x(1) = rA, x(2) = rB, x(3) = rC, x(4) = L
 fun = @(x) (-((T1-T2)/(((log(x(2)/x(1)))/(2*pi*x(4)*k))+((log(x(3)/x(2)))/(2*pi*x(4)*18.68))+(1/(h*2*pi*x(4)*x(3))))));
 % ----Initial guess----
 x0 = [0.00006,0.0003,0.0003,1];
 % ----Inequality linear constraints---- = g1, g4, g5
 A = [1 -1 0 0;0 1 -1 0];
 b = [-(230/e);-0.0000045];
 % ----Equality constraints----
 Aeq = [];
 beq = [];
 % ----Upper and lower bounds---
 lb = [-Inf, 0.0002, -Inf, -Inf];
 ub = [Inf,Inf,0.021,2.14];
 % ----Non linear constraints----
 nonlcon = @cons5;
 rng default
 problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'nonlcon',nonlcon,'lb',lb,'ub',ub);
 x = run(GlobalSearch,problem);
 % output variables from each material
 variable_list_global_search = [variable_list_global_search, [x(1),x(2),x(3),x(4)]];
 disp(table(x(1),x(2),x(3),x(4),'VariableNames',{'rA', 'rB', 'rC', 'Length'}))
 disp(['Final Objective: ' num2str(fun(x))])
 disp(['Heat flux (W/m^2): ' num2str(-(fun(x)))])

%% Multiobjective for Heat flux v embodied CO2 - for 3 materials

disp([ ' -------------------- Multiobjective for 3 different materials -------------------- '])

% Using Global search values

% Empty arrays for pass and fail values
pass_values = table();
fail_values = table();

thickness_range =  0:0.0001:0.001; % thickness range for rB, gives us 10 thicknesses

for y = 1:length(thickness_range) %10 times for 10 thicknesses
    for i = 1:length(materials) %3 times for each material
        
        %----Current thickness----
        current_thickness = thickness_range(y);
        
        %----Material dependant properties----
        k = kB(i);
        e = eB(i);
        
        materials(i);
        
        if i == 1
            a = 1;
        end
        
        if i == 2
            a = 5;
        end
            
        if i == 3
            a = 9;
        end
        
        b = a+1;
        c = b+1;
        d = c+1;
        
        %----Extract optimised variables from Global Search----
        rA = variable_list_global_search(:,a); %extract rA
        x_2 = variable_list_global_search(:,b); %extract rB
        x_3 = variable_list_global_search(:,c); %extract rC
        L = variable_list_global_search(:,d); %extract L
        
        %----Change objective according to new thickness value for rB-----
        rB = current_thickness + rA; %new rB is rA + current thickness
        rC = rB + (x_3-x_2); %new rC is new rB + (Old rC - rB)
        
        %----Find heat flux for each thickness----
        heatflux = (-((7.5)./(((log(rB/rA))./(2.*pi.*L.*k))+((log(rC./rB))/(2.*pi.*L.*18.68))+(1./(1525.*2.*pi.*L.*rC)))));
        
        %----Find CO2 impact of layer B----
        volume = (((rB.^2)*pi)-((rA.^2).*pi)).*L;
        mass = density(i).*volume;
        CO2impact = mass.*CO2_impact(i);
        
        %----Current material----
        current_material = materials(i);  
        
        %----Append data----
        newrow = {current_material, current_thickness, CO2impact, heatflux};
        
        % Check dielectric strength, if current thickness > thickness
        % for sufficient dielectric strength = PASS 
        if  thickness_range(y) > (230/eB(i));
            pass_values = [pass_values; newrow];

        else
            fail_values = [fail_values; newrow];
            
        end  
    end
end


fail_values;

% Uncomment to display pass values 
% pass_values


%----Plot multiobjective problem----
figure(7)
scatter(pass_values{:,4}, pass_values{:,3})  % PASS VALUES - heat flux value, CO2 value
hold on
scatter(fail_values{:,4}, fail_values{:,3}) % FAIL VALUES - heat flux value, CO2 value
grid on
ylabel('CO2 Impact (J/kg)')
xlabel('Heat Flux (W/m^2)')
title('Multi-Objective Optimisation - Pareto set')
legend('Pass', 'Fail')



