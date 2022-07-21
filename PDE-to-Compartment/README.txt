# PDE - Compartment

pure_diffusion.py contains the code to simulate test problem 1 and 2. The class PureDiffusion has an argument start that can be equal to 'left' or 'equilibrium'. And morphogen_gradient.py is the code for the test problem 3.

### DATA

The data for the three first test problems are in the experiment folder.

### Run simulations

To run a simulation, select the chemical_reaction equal to 'TEST_PROBLEM_1', 'TEST_PROBLEM_2' or 'TEST_PROBLEM_3' in run_sims.py. Be careful, run_sims.py will overwrite the data in the experiment folder.

### Plot the DATA

To plot the data, run plot_sims.py with the appropriate chemical_reaction.

### Packages
numpy
tqdm