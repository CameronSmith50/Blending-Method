%Solves the diffusion equation analytically based on the initial conditions
%where the equation has been made using fourier analysis and 

clear

Initial_Conditions_Pure_Diffusion %load initial conditions

F=1000; %set the number of sums for the fourier series;

%work out the coefficients for the sum based on the initial condition
a=zeros(1,F);
if uniform ==1
for j=1:F

    if j==1
        a(1) = n0/(X(2)-X(1));
    else
        a(j) = 0;
    end
end

elseif left==1
for j=1:F

    if j==1
        a(1) = n0/(X(2)-X(1));
    else
        a(j) = 2*n0/((j-1)*pi*(I(1)-X(1)))*sin((j-1)*pi*(I(1)-X(1))/(X(2)-X(1)));
    end
end

elseif right==1
for j=1:F

    if j==1
        a(1)=(n0/(k));
    else
        a(j)=(-(8*n0)/(k*(j-1)*pi))*sin((3*pi*(j-1))/4);
    end
end
end


u=zeros(length(time_int_vec),10*k);
for q=1:time_ints
    %discretise and calculate the analtical value at t_final
    x = X(1):(1/(10*k)):X(2);
    
    t=time_int_vec(q);
    for j=1:length(x)
        temp=a(1);
        
        for ii=1:(F-1)
            temp = temp + a(ii+1)*exp(-D*ii^2*pi^2/((X(2)-X(1))^2)*t)*cos(ii*pi*x(j)/(X(2)-X(1)));
        end
        u(q,j)=temp;
    end
end

if uniform == 1
    save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_1/TEST1_ANALYTICAL','u')
elseif left == 1
    save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_2/TEST2_ANALYTICAL','u')
end

