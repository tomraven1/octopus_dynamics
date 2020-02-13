# octopus_dynamics
opti_4sec can be used for running control experiments with objetive functions in evalu_4sec.
Load learnedmodel.mat for a prelearned dynamic model

For using the dynamic model of the octopus arm :
Use explorenew.m
define the time of simulation and the forces to be applied in the matrix 'ini' (size of 3 x time*100)
To change number of sections in the arm change 'npie' variable in piecewise_driver. Note that only the base section is actuated
Single section manipulator is much faster to simulate. 

