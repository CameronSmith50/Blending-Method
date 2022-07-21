TEST 4 - 2nd order slanted start

Second order system with reactions: 2A->A  empty->A

Parameter Values and Initial Conditions:

n0=465;          %initial number of particles

X=[0,10];          %domains
Y=[0,1];          %domains
Z=[0,1];          %domains

repeats=1000;     %number of repeats

I=[10/3,20/3];      %Interval Positions X-axis

D = 1;            %Diffusion coefficient

k1=0.1;           %Reaction rate

VOL=10;

k2=89.7           %production rate: k1*(300*(299))/(VOL^2);

rho=0.06          %reaction radius

t_final = 1;      %Total time

k=20;             %Number of Compartments

H=1/30;           %Compartment height

dt = 25^(-6);     %Time step (for particle diffusion)

time_int = 0.01   %Time step for recording


