README file - A hybrid framework for the simulation of stochastic reaction-diffusion processes

--------------------------------------------------------------------------------------------

Contact information:
Name: 
Email: 

--------------------------------------------------------------------------------------------

This is the README file for the various parts of code that can be found in the supplimentary material. The directory contains several folders, each of which contains the code for different parts of the paper. The example numbers refer to the examples in the main text of the manuscript.

The code is written using MATLAB.

--------------------------------------------------------------------------------------------

Examples 1,2: Pure Diffusion
----------------------------

The initial conditions and parameter values are stored in Initial_Conditions_Pure_Diffusion.m. This should be edited to alter all parameter values and the initial condition. Initial condition "unif" refers to the first example and initial condition "left" refers to the second example. 

The script Meso_Micro_Pure_Diffusion_v2_4.m runs the compartment-to-Brownian hybrid method with the parameters and initial condition specified in Initial_Conditions_Pure_Diffusion.m. The analytical solution can also be run by using the script analytical.m.

Both Meso_Micro_Pure_Diffusion_v2_4.m and analytical.m save the final solution into the Plotting folder, under the appropriate designation (Data_hybrid_compartment_brownian and then TEST_PROBLEM_1 or 2 depending on the initial condition). These can be used by the plotting tool.

Example 3: Morphogen Gradient
-----------------------------

The initial conditions and parameter values are stored in Initial_Conditions_MG.m. This should be edited to alter all parameter values.

The script Meso_Micro_MG_v1_5.m runs the compartment-to-Brownian hybrid method with the parameters and initial condition specified in Initial_Conditions_MG.m. The analytical solution can also be run by using the script Analytical_MG_v1_2.m.

Both Meso_Micro_MG_v1_5.m and Analytical_MG_v1_2.m save the final solution into the Plotting folder, under the appropriate designation (Data_hybrid_compartment_brownian and then TEST_PROBLEM_3). These can be used by the plotting tool.

Example 4: Bimolecular Reactions
--------------------------------

The initial conditions and parameter values are stored in Initial_Conditions_2ND.m. This should be edited to alter all parameter values.

The script Meso_Micro_2ND_v1_6.m runs the compartment-to-Brownian hybrid method with the parameters and initial condition specified in Initial_Conditions_2ND.m. The fully Brownian solution can also be run by using the script Particle_2ND_v1_5.m.

Both Meso_Micro_2ND_v1_6.m and Particle_2ND_v1_5.m save the final solution into the Plotting folder, under the appropriate designation (Data_hybrid_compartment_brownian and then TEST_PROBLEM_4). These can be used by the plotting tool.

All other functions are helper functions. See their code for more details.