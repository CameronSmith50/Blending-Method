TEST 3 - Morphogen gradient, uniform start

Morphogen gradient sytem with particle entry through left hand side of domain and unifrom degradation rate across 
the domain, particles start uniformly distributed across the domain.

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

lambda=10000      %lambda value for particle entry though left hand boundary

mu=10             %particle degradation rate
