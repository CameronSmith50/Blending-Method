%Solve a pure diffusion system using a meso-micro blending method between
%a compartment based model and a brownian model.
% clear

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

%Calculate comaprtment height (x axis)
H=((X(2)-X(1))/k);

%calculate compartment volume
V=H*(Y(2)-Y(1))*(Z(2)-Z(1));

%calculate the number of compartments in each region
NCOMP=round((I(1)-X(1))/H);
NBLEND=((I(2)-I(1))/H);
NPART=((X(2)-I(2))/H);

%error check for compartments across boundaries
if NCOMP~=ceil(NCOMP)
    disp 'ERROR: Compartments cannot intersect interval boundaries'
    return
end

prod=0;
reac=0;

%Find Dc the diffusion coefficients for left jumps
Dcl=zeros(1,k-NPART);
xc=(X(1)+(H/2):H:I(2)); %midpoints of compartments

for i=1:k-NPART
    if xc(i)<I(1)
        Dcl(i)=D;
    else
        Dcl(i)=D*((I(2)-(xc(i)-H/2))/(I(2)-I(1)));
        
    end
end

%Find Dc the diffusion coefficients for right jumps
Dcr=zeros(1,k-NPART);
xc=(X(1)+(H/2):H:I(2)); %midpoints of compartments

for i=1:k-NPART
    if xc(i)<I(1)
        Dcr(i)=D;
    else
        Dcr(i)=D*((I(2)-(xc(i)+H/2))/(I(2)-I(1)));
        
    end
end

%find d for each compartment
dl=Dcl(1:(k-NPART))/((H)^2);
dr=Dcr(1:(k-NPART))/((H)^2);

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

tic

PART=[];

for l=1:repeats
    % Progress counter
    if repeats > 1
        if mod(floor(l/repeats*100),2) == 0 && n_count ~= floor(l/repeats*100)
            fprintf('|')
            n_count = n_count + 2;
        end
    end
    
    
    %com=21-[1:20];
    
    com=31-[1:30];
    %com=ones(1,k)*(n0/k);
    
    
    %Calculate initial number of particles
         part_x=[];
         for i=NCOMP+1:k
             part_x = [part_x,H*(i-1) + rand(1,com(i))*(X(2)-X(1))/(k)];
    
         end
    part_y = rand(1,sum(com(NCOMP+1:k)))*(Y(2)-Y(1));
    part_z = rand(1,sum(com(NCOMP+1:k)))*(Z(2)-Z(1));
    part=[part_x',part_y',part_z'];

    alpha_coeff=[dl,dr,ones(1,NCOMP+NBLEND)*k1*(k/VOL),k2*VOL];

    n_part=sum(com);
    
    
    %initialise time
    t=0;
    
    %initial sum recording
    rec_com_sum(1,:)=rec_com_sum(1,:)+com;
    
    %initial particle number recording
    rec_part(1)=rec_part(1)+n0;
    
    %to be used later in the sum recording section
    z=2; %z is current recording level
    
    tb=t+dt; %initial next brownian timestep
    
    part_rec=[];
    
    %Do compartment method using algorithm 9
    while t<=t_final
        
        
        %generate a random number
        r1=rand;
        
        %find alpha
        alpha=[com(1:(NCOMP+NBLEND)),com(1:(NCOMP+NBLEND)),com(1:(NCOMP+NBLEND)).*((com(1:(NCOMP+NBLEND)))-1),1].*alpha_coeff;
                
        %Find Alpha0
        alpha0=sum(alpha);
        
        %find the time until the next reaction
        tau=(1/(alpha0))*log(1/r1);
        
        tc=t+tau; %next reaction
        
        %find out wheter a compartment jump occurs in a given timestep
        if tc<tb %compartment reaction
            t=tc;
            
            r2a=rand*alpha0; %define r2a as a random number [0,1] times alpha0
            
            %initialise the cumalitive sum which will be built iteratively
            cumsumalpha=alpha(1);
            j=1;
            while cumsumalpha<r2a
                j=j+1;
                cumsumalpha=cumsumalpha+alpha(j);
            end
            
            if j<=NCOMP+NBLEND && j>1%left jump
                
                if j>1 &&j<=NCOMP %left jump in pure compartment region
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    
                elseif j==NCOMP+1 %left jump over I1 boundary
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    CHOOSE_PART=find(part(:,1).*(((j-1)*H)<part(:,1)).*(part(:,1)<((j)*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART,:) = [];
                else
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    CHOOSE_PART=find(part(:,1).*(((j-1)*H)<part(:,1)).*(part(:,1)<((j)*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART,1) = part(RAND_PART,1)-H;
                end
                
                
            elseif j>NCOMP+NBLEND && j<=2*(NCOMP+NBLEND) %right jump
                
                if j<(2*NCOMP)+NBLEND %right jump in pure compartment region
                    com(j-(NCOMP+NBLEND))=com(j-(NCOMP+NBLEND))-1;
                    com(j+1-(NCOMP+NBLEND))=com(j+1-(NCOMP+NBLEND))+1;
                    
                elseif j==(2*NCOMP)+NBLEND %right jump over I1 boundary
                    com(j-(NCOMP+NBLEND))=com(j-(NCOMP+NBLEND))-1;
                    com(j+1-(NCOMP+NBLEND))=com(j+1-(NCOMP+NBLEND))+1;
                    
                    %calculate y coordinates of trapezium for jumping
                    y1=(com(NCOMP-1)*(3*(H/2))+com(NCOMP+1)*((H/2)))/(2*H);
                    y2=(com(NCOMP-1)*((H/2))+com(NCOMP+1)*(3*(H/2)))/(2*H);
                    
                    y1_bar=y1/max(y1,y2);
                    y2_bar=y2/max(y1,y2);
                    
                    check=0;
                    while check==0
                        r=rand;
                        part_new=I(1)+(rand*H);
                        if r<(y1_bar*(I(1)+H-part_new)+y2_bar*(part_new-I(1)))/H
                            part=[part;part_new,rand*(Y(2)-Y(1)),rand*(Z(2)-Z(1))];
                            check=1;
                        end
                    end
                    
                elseif (2*NCOMP)+NBLEND<j && j~= (2*(NCOMP+NBLEND)) %right jump in blended region
                    com(j-(NCOMP+NBLEND))=com(j-(NCOMP+NBLEND))-1;
                    com(j+1-(NCOMP+NBLEND))=com(j+1-(NCOMP+NBLEND))+1;
                    CHOOSE_PART=find(part(:,1).*(((j-1-(NCOMP+NBLEND))*H)<part(:,1)).*(part(:,1)<((j-(NCOMP+NBLEND))*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART,1) = part(RAND_PART,1)+H;
                end
                
                
            elseif j>2*(NCOMP+NBLEND) && j<=2*(NCOMP+NBLEND)+(NCOMP+NBLEND) %2nd order reaction
                
                reac=reac+1;
                
                if j>2*(NCOMP+NBLEND) && j<=2*(NCOMP+NBLEND)+NCOMP
                    com(j-(2*(NCOMP+NBLEND)))=com(j-(2*(NCOMP+NBLEND)))-1;
                else
                    com(j-(2*(NCOMP+NBLEND)))=com(j-(2*(NCOMP+NBLEND)))-1;
                    CHOOSE_PART=find(part(:,1).*(((j-1-(2*(NCOMP+NBLEND)))*H)<part(:,1)).*(part(:,1)<((j-(2*(NCOMP+NBLEND)))*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART,:) = [];
                    n_part=n_part-1;
                end
      
            end
            
        elseif tb<tc
            t=tb;
            
            %apply particle reactions for particles in pure brownian
            %region
            
            brown_part=(part.*(I(2)<part(:,1)));%find particles in this compartment
            brown_part(brown_part==0)=[];
            brown_part=reshape(brown_part,[length(brown_part)/3,3]);
            

            dists=Pairwise_Dist(brown_part);
            
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
            PAIRS_REM = Remove_Duplicates(PAIRS,length(brown_part));  % OK
            
            % Find the number of pairs to react
            [TO_REACT_2,~] = size(PAIRS_REM);  % OK
            
            % Turn PAIRS_REM into a vector
            PAIRS_VEC = PAIRS_REM(:);  % OK
            
            % Remove the particles in PAIRS_REM, and add the new particles in
            % NEW_PARTS_2
       
            logical2=rand(length(PAIRS_VEC),1)<0.5;
            PAIRS_VEC=PAIRS_VEC(find(PAIRS_VEC.*(logical2)));
%         
%             
            brown_part(PAIRS_VEC,:) = [];  % OK
            
            n_part = n_part - length(PAIRS_VEC);
            
            reac=reac+length(PAIRS_VEC)/2; %reac counter
            
            part1=(part.*(I(2)>part(:,1)));
            part1(part1==0)=[];
            part1=reshape(part1,[length(part1)/3,3]);
            
            part=[part1;brown_part];
            
            % Zeroth order reaction
            % Calculate the probability that a particle is introduced
            PROB = k2*dt*VOL;
            
            % Draw a random number and check against this. Add particle
            % appropriately
            if rand < PROB
                
                prod = prod+1;
                
                new_part=[(X(1)+X(2))*rand,Y(1)+Y(2)*rand,Z(1)+Z(2)*rand];
                
                if new_part(1)<I(1)
                    new_com=ceil(new_part(1)*(k)/(X(2)-X(1)));
                    com(new_com)=com(new_com)+1;
                    
                else
                    
                    part = [part;new_part];
                end
                n_part=n_part+1;
                
            end
            
            Db=zeros(1,length(part));
            for i=1:length(part)
                if part(i)>I(2)
                    Db(i)=D;
                elseif I(1)<part(i) && part(i)<I(2)
                    Db(i)=D*(1-((I(2)-part(i))/(I(2)-I(1))));
                else
                    Db(i)=0;
                    
                    
                end
            end
            
            %update particle positions according to brownian motion
            
            temp=part(:,1)+((sqrt(2*Db*dt))'.*randn(length(part(:,1)),1))+(part(:,1)>I(1)).*(part(:,1)<I(2))*((D/(I(2)-I(1)))*dt);
            part(:,1)=temp;
            %part(:,1)=part(:,1)+((sqrt(2*D*dt)).*randn(length(part(:,1)),1));
            
            temp=part(:,2)+((sqrt(2*Db*dt))'.*randn(length(part(:,2)),1));
            part(:,2)=temp;
            temp=part(:,3)+((sqrt(2*Db*dt))'.*randn(length(part(:,3)),1));
            part(:,3)=temp;
            
            %apply right hand boundary condition
            part(:,1)=part(:,1)+2*(part(:,1)>X(2)).*(X(2)-part(:,1));
            part(:,2)=part(:,2)+2*(part(:,2)>Y(2)).*(Y(2)-part(:,2));
            part(:,3)=part(:,3)+2*(part(:,3)>Z(2)).*(Z(2)-part(:,3));
            
            %apply left hand boundary condition
            part(:,1)=part(:,1)+2*(part(:,1)<I(1)).*(I(1)-part(:,1));
            part(:,2)=part(:,2)+2*(part(:,2)<Y(1)).*(Y(1)-part(:,2));
            part(:,3)=part(:,3)+2*(part(:,3)<Z(1)).*(Z(1)-part(:,3));
            
                        
            %update the compartment numbers in the blended region and part
            %region
            for y=NCOMP+1:k
                com(y)=sum((ceil(part(:,1)*(k/(X(2)-X(1)))))==y);
            end
            
            %reset tb
            tb=tb+dt;
            
 
        end
        
        
        while t+dt/10>=(z-1)*time_int %for recording
                      
            %find number of "compartments" in particle region
            rec_com_sum(z,:)=rec_com_sum(z,:)+ com;
            rec_part(z)=rec_part(z) + sum(com);
            z=z+1;           
            
        end
        if z>time_ints
            break
        end
        
        
    end
    PART=[PART;part];
end

NUM_PART_HYBRID=rec_part/repeats;

second_order_blend=rec_com_sum/repeats;

Time_second=toc;

save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_4/TEST4_BLENDING','second_order_blend')
