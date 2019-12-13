# DE-4 Optimisation: Optimsing the environmental impact of a dishwasher

Christopher Turner (Subsystem 1), Francesca Suer (Subsystem 2), Seung Huh (Subsystem 3), Visakan Mathivannan (Subsystem 4)

A dishwasher is a complicated kitchen appliance used to wash plates and utensils, commonly found in homes under the kitchen countertop. The aim of this project was to reduce the environmental impact of a dishwasher. This was achieved on a system level by minimising the cost throughout the lifetime of the dishwasher through reducing the energy consumption of the pump, friction losses in the piping, maximizing the efficiency of the heating element and maximizing the cleaning ability of the spray arms. The problem was split into four subsystems:

- Piping
- Heating element
- Pump
- Spray mechanism

This enabled each author to create models of their assigned subsystem before optimizing each part.

Afterwards the models were combined to create a single unified optimization problem. This was solved using SQP algorithm yielding the result of £197.18 total cost.

## Dependancies
The project was completed entirely using MATLAB with some addition toolboxes which are required to run these programs.

- [Deep Learning Toolbox](https://uk.mathworks.com/products/deep-learning.html?s_tid=AO_PR_info)
- [Global Optimization Toolbox](https://uk.mathworks.com/products/global-optimization.html)

## System Level

This combines subsystem 1,2 and 3 to find the monetary cost value over the dishwashers lifetime. Subsystem 4 has been ommited because its contributions are negliable.

### How to run

Run the [system_level_main.m](System/system_level_main.m) file.

### Summary

The system is optimised in terms of cost.

## Sub system 1 - Piping

This system is to optimise the pipe size, diameter and layout for the water distribution within the dishwasher. 
Contained in the Subsystem1 folder is another MATLAB document [Best_pipe_Cross_section.m](Subsystem1/Best_pipe_Cross_section.m) which shows and proves that a circle is the optimum shape for a pipe cross section.

### How to run

Run [subsystem1.m](Subsystem1/subsystem1.m) file, found in Subsystem1 folder

### Summary

Pipes are optimised

## Sub system 2 - Heating element

This subsystem optimises heat flux (a measure of watts per square meter) for the heating element.

### How to run

Run the [Subsystem_2.m](Subsystem2/Subsytem_2.m) file.

### Summary

The heating element design is optimised

## Sub system 3 - Pump

This subsystem optimises power output (wattage) of centrifugal pump's motor.

### How to run

Run the [optimization_subsystem_3.m](Subsystem3/optimization_subsystem_3.m) file.

### Summary

Pump is optimised

## Sub system 4 - Spray mechanism

The spray mechanism is the subsystems which designates how the water is distributed across the plates and the goal of this problem was to find the optimal configuration to optimise the rotation of the spray arm in order to increase the cleaning efficiency

### How to run

Run the [spray_mechanism_main.m](Subsystem4/spray_mechanism_main.m) file

### Summary

Ultimately a configuration is defined and the rotation is optimised to have an agular velocity 18 rad/s

## Report
For further indepth analysis take a look at the [report](FinalReport_Group2.pdf)
