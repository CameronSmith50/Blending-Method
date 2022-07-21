%Solves the diffusion equation analytically based on the initial conditions
%where the equation has been made using fourier analysis and 

clear
close all

Initial_Conditions_MG %load initial conditions

F=2000; %set the number of sums for the fourier series;

C=3000;       %number of compartments for the mesh

%work out the coefficients for the sum based on the initial condition

%discretise and calculate the analtical value at t_final
x = X(1):(1/(10*k)):X(2);
t=t_final;

sigma=sqrt(mu/D);

u=zeros(length(time_int_vec),length(x));
for q=1:time_ints
    
    t=time_int_vec(q);
    for j=1:length(x)
        temp=lambda*(cosh(sigma*(x(j)-1))/(sigma*sinh(sigma)));
        for i=1:F
            temp=temp-((2*lambda*D)*...
                (exp(-(((i*pi)^2)*D+mu)*t)*cos(i*pi*x(j)))/(D*(i^2)*(pi^2)+mu));
        end
        u(q,j)=temp;
        
        
    end
    
end

save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_3/TEST3_ANALYTICAL','u')
