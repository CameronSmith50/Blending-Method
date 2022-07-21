%clear

Initial_Conditions_2ND_v1_1

RHO=0.06;%reaction radius rho
GAM = sqrt(4*D*dt)/RHO;  % Dimensionless parameter
KAPPA = k1*dt/(RHO^3);  % Dimensionless parameter


%CALCULATING UNKNOWN PARAMETERS - TAKEN  FROM CAMERON
S = 50;  % Upper integral bound
h1 = 0.01;  % Mesh size in [0,1]
n2 = 3000;  % Number of mesh points in (1,S]
h2 = (S-1)/n2;  % Second mesh width
r_VEC = [0:h1:1,1+h2:h2:S];  % Radial mesh vector
r_VEC_LEN = length(r_VEC);
K_MAT = zeros(r_VEC_LEN);  % Diffusion kernel matrix
for ii = 1:r_VEC_LEN
    for jj = 1:r_VEC_LEN
        K_MAT(ii,jj) = K_fun(r_VEC(ii),r_VEC(jj),GAM);
    end
end
fprintf('\nFinding Plambda...')
P_LAMBDA = Find_Params(GAM,KAPPA,K_MAT,r_VEC);  % Calculate parameters

%Make a recording matrix that records particle positions in boxes at chosen
%time intervals while summing over repeats
rec_com_sum=zeros(length(time_int_vec),k);

%make a recording vector which records the number of particles at each
%recording step
rec_part=zeros(length(time_int_vec),1);

% Progress bar
n_count = 0;
if repeats > 1
    count = 0;
    fprintf('\nProgress of mesoscopic method:\n\n   0                        50                      100 (Percentage)\n   |                        |                        |\n   |')
end

phi_part=zeros(1,1000);

tic

for l=1:repeats
    % Progress counter
    if repeats > 1
        if mod(floor(l/repeats*100),2) == 0 && n_count ~= floor(l/repeats*100)
            fprintf('|')
            n_count = n_count + 2;
        end
    end
    
    
  %  com=(n0/k)*ones(1,k);
    
    

    %
    NPART=n0;
    
 
   com=31-[1:30];    
    
    %Calculate initial number of particles
         part_x=[];
         for i=1:k
             part_x = [part_x,H*(i-1) + rand(1,com(i))*(X(2)-X(1))/(k)];
    
         end
    %part_x = rand(1,n0)*(X(2)-X(1));
    part_y = rand(1,sum(com(1:k)))*(Y(2)-Y(1));
    part_z = rand(1,sum(com(1:k)))*(Z(2)-Z(1));
    part=[part_x',part_y',part_z'];        
    
    %update the compartment numbers in the blended region and part
    %region
    for i=1:k
        com(i)=sum((ceil(part(:,1)*(k/(X(2)-X(1)))))==i);
    end

    n_part=sum(com);
    
    
    %initialise time
    t=0;
    

    %initial sum recording
    rec_com_sum(1,:)=rec_com_sum(1,:)+com;

    %initial particle number recording
    rec_part(1)=rec_part(1)+n0;
    
    %find pairwise distances between particles
    part=transpose([part_x;part_y;part_z]);
    
    %to be used later in the sum recording section
    z=2; %z is current recording level
    
    %Do compartment method using algorithm 9
    while t<t_final-dt/2
        
        t=t+dt;
        
        % Zeroth order reaction
        % Calculate the probability that a particle is introduced
        PROB = k2*dt*VOL;
        
        % Draw a random number and check against this. Add particle
        % appropriately
        if rand < PROB
            part = [part;X(1)+X(2)*rand,Y(1)+Y(2)*rand,Z(1)+Z(2)*rand]; %#ok
            NPART = NPART + 1;
        end
        
        
        dists=Pairwise_Dist(part);
        
        %find out if particles are in a reaction radius
        logical = dists > 0 & dists < RHO;
        parts_to_react=Find_Parts(logical);
        
        if sum(sum(logical)) > 0
            REACT = P_LAMBDA > rand(sum(sum(logical)),1);
        else
            REACT = [];
        end
        
        % Find the rows in PARTS that correspond to the reacting pairs
        PAIRS = parts_to_react(REACT == 1,:);
        PAIRS_REM = Remove_Duplicates(PAIRS,NPART);  % OK
        
        % Find the number of pairs to react
        [TO_REACT_2,~] = size(PAIRS_REM);  % OK
        
        % Turn PAIRS_REM into a vector
        PAIRS_VEC = PAIRS_REM(:);  % OK
        
        logical2=rand(length(PAIRS_VEC),1)<0.5;
        PAIRS_VEC=PAIRS_VEC(find(PAIRS_VEC.*(logical2)));
        
        part(PAIRS_VEC,:) = [];  % OK

        % Update N
        NPART = NPART - length(PAIRS_VEC);

        
        %if PROB>0.1
        %    error('PROB too big, choose smaller time step')
        %end
        
        
        %Now update particle positions in 3 dimensions
        
        %update particle positions according to brownian motion
        part(:,1)=part(:,1)+((sqrt(2*D*dt)).*randn(length(part(:,1)),1));
        part(:,2)=part(:,2)+((sqrt(2*D*dt)).*randn(length(part(:,2)),1));
        part(:,3)=part(:,3)+((sqrt(2*D*dt)).*randn(length(part(:,3)),1));
        
        
        %apply right hand boundary condition
        part(:,1)=part(:,1)+2*(part(:,1)>X(2)).*(X(2)-part(:,1));
        part(:,2)=part(:,2)+2*(part(:,2)>Y(2)).*(Y(2)-part(:,2));
        part(:,3)=part(:,3)+2*(part(:,3)>Z(2)).*(Z(2)-part(:,3));
        
        %apply left hand boundary condition
        part(:,1)=part(:,1)+2*(part(:,1)<X(1)).*(X(1)-part(:,1));
        part(:,2)=part(:,2)+2*(part(:,2)<Y(1)).*(Y(1)-part(:,2));
        part(:,3)=part(:,3)+2*(part(:,3)<Z(1)).*(Z(1)-part(:,3));
        
        for y=1:k
            com(y)=sum((ceil(part(:,1)*(k/(X(2)-X(1))))==y));
        end
        
        if t>25
        phi_part(sum(com))=phi_part(sum(com))+1;
        end
        
        
        
        
        while t+dt/10>=(z-1)*time_int %for recording
            
            %find number of "compartments" in particle region
            
            %rec_com_sum(z)=rec_com_sum(z)+ com;
            rec_com_sum(z,:)=rec_com_sum(z,:)+ com;
            rec_part(z)=rec_part(z) + sum(com);
            z=z+1;
            
        end
        
    end
end

time_part=toc;

second_order_particle=rec_com_sum/repeats;

second_order_sum=sum(second_order_particle,2);

NUM_PART_PART=rec_part/repeats;

phi_part=phi_part/(sum(phi_part));

save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_4/TEST4_particle','second_order_particle')







