TEST 1 - MAINTAINING UNIFORM GRADIENT

Pure diffusion system with a uniform start, particles are uniformly distributed accross the domain:

Parameter Values and Initial Conditions:

n0=1000;          %initial number of particles

X=[0,1];          %One dimentional domain

repeats=1000;     %number of repeats

I=[1/3,2/3];      %Interval Positions

D = 1;            %Diffusion coefficient

t_final = 1;      %Total time

k=30;             %Number of Compartments

H=1/30;           %Compartment height

dt = 10^(-4);     %Time step (for particle diffusion)

time_int = 0.01   %Time step for recording

