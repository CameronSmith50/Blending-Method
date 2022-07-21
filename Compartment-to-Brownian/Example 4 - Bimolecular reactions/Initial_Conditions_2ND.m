%%Initial Conditions%%

%clear

repeats=5; %number of repeats
n0=465; %initial number of particles

X=[0,10]; %One dimentional domain

%compartment number
k=30;


%Calculate comaprtment height
H=((X(2)-X(1))/k);

Y=[0,1]; 
Z=[0,1];

VOL=(X(2)-X(1))*(Y(2)-Y(1))*(Z(2)-Z(1));

ComVOL=VOL/k;

%Set the interval positions:
I=[10/3,20/3];

%Diffusion variables
D = 1;

%Second order reaction rate
k1=0.1;
k1_saved=k1;


%zeroth order production rate
k2=k1_saved*(300*(300-1))/(VOL^2); %*4 to account for changed size

%Time variables
t_final = 1; %Total time

dt = 0.000025; %Time step0


T = (0:dt:t_final); %Vector of time intervals
t_steps = length(T); %no of time steps
time_int=0.01; %set the recording time interval
time_int_vec=(0:time_int:t_final); %vector of recording time intervals
time_ints=length(time_int_vec);