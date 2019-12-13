# DE-4 Optimisation: Optimsing the environmental impact of a dishwasher
Christopher Turner (Subsystem 1), Francesca Suer (Subsystem 2), Seung Huh (Subsystem 3), Visakan Mathivannan (Subsystem 4)

This repository contains the teams project work for solving an optimisation problem. You can find a copy of our final report alongside code for each individual subsystem and for the system level problem

The dishwasher could be broken down into a series of components each subsystem chosen was designed to optimise the transfer of either heat or pressure and they were as follows:
* Piping
* Heating element
* Pump
* Spray mechanism

The  project was completed entirely using MATLAB with some addition toolboxes which are required to run these programs.

* [Deep Learning Toolbox](https://uk.mathworks.com/products/deep-learning.html?s_tid=AO_PR_info)
* [Global Optimization Toolbox](https://uk.mathworks.com/products/global-optimization.html)
## System Level
This combines subsystem 1,2 and 3 to find the monetary cost value over the dishwashers lifetime. Subsystem 4 has been ommited because its contributions are negliable.
### How to run
Run the system_level_main.m file.
### Summary
The system is optimised in terms of cost.
## Sub system 1 - Piping
This system is to optimise the pipe size, diameter and layout for teh water distribution within the dishwasher.
Contained in the Subsystem1 folder is another MATLAB document [Best_pipe_Cross_section.m](Subsystem1/Best_pipe_Cross_section.m) which shows and proves that a circle is the optimum shape for a pipe cross section.
### How to run
Run [subsystem1.m](Subsystem1/subsystem1.m) file, found in Subsystem1 folder
### Summary
Pipes are optimised
## Sub system 2 - Heating element
This subsystem optimises heat flux (a measure of watts per square meter) for the heating element.
### How to run
Run the Subsystem_2.m file.
### Summary
The heating element design is optimised
## Sub system 3 - Pump
This subsystem optimises centrifugal pump's motor power (wattage) output.
### How to run
Run the Main.m file 
### Summary
Pipes are optimised
## Sub system 4 - Spray mechanism
The spray mechanism is the subsystems which designates how the water is distributed across the plates and the goal of this problem was to find the optimal configuration to optimise the rotation of the spray arm in order to increase the cleaning efficiency
### How to run
Run the [spray_mechanism_main.m](Subsystem4/spray_mechanism_main.m) file 
### Summary
Ultimately a configuration is defined and the rotation is optimised to have an agular velocity 18 rad/s
