%%Initial Conditions%%

n0=1000; %initial number of particles
X=[0,1]; %One dimentional domain
repeats=1000; %number of repeats

%Set the interval positions:
I=[1/3,2/3];

%Diffusion variables
D = 1;

%Time variables
t_final = 1; %Total time

%compartment number
k=30;

%Calculate comaprtment height
H=((X(2)-X(1))/k);

dt = 10^(-4); %Time step
%dt=(H^2)/(2*D);


T = (0:dt:t_final); %Vector of time intervals
t_steps = length(T); %no of time steps
time_int=0.01; %set the recording time interval
time_int_vec=(0:time_int:t_final); %vector of recording time intervals
time_ints=length(time_int_vec); 


%Boolean to determine initial positions
left=1;
uniform=0;
right=0;
    
if (left+uniform+right)~=1
    disp('ERROR: Inital particle positions boolean values must sum to 1')
    return
end
  