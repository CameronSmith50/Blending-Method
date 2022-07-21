%Solve a pure diffusion system using a meso-micro blending method between
%a compartment based model and a brownian model.

clear

%Run initial conditions
Initial_Conditions_MG

%Calculate comaprtment height
H=((X(2)-X(1))/k);

%calculate the number of compartments in each region
NCOMP=round((I(1)-X(1))/H);
NBLEND=round((I(2)-I(1))/H);
NPART=round((X(2)-I(2))/H);



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
    
    %Calculate initial number of particles
    determ_prop=floor(n0/k)*k;
    rand_init_pos=rand(1,n0-determ_prop);
    
    
    %Calculate initial number of particles per compartment (including recording
    %compartments in pure particle region)
    rand_com=zeros(1,k);
    for y=1:k
        rand_com(y)=sum(ceil(rand_init_pos*(k))==y);
    end
    
    com=(determ_prop/k)*ones(1,k)+rand_com;
    
    %combine coefficients to be used to calculate alpha
    alpha_coeff=[dl,dr,ones(1,k)*mu,D*lambda];

    %initialise time
    t=0;

    %Calculate initial number of particles
    part=[];
    for i=NPART+1:k
        part= [part,H*(i-1) + rand(1,com(i))/(k)]; %#ok
    end

    %initial sum recording
    rec_com_sum(1,:)=rec_com_sum(1,:)+com;
    
    %initial particle number recording
    rec_part(1)=rec_part(1)+n0;
    
    %to be used later in the sum recording section
    z=2; %z is current recording level
    
    tb=t+dt; %initial next brownian timestep
    
    %Do compartment method using algorithm 9
    while t<=t_final
        
        %generate a random number
        r1=rand;
        
        %find alpha
        alpha=[com(1:(NCOMP+NBLEND)),com(1:(NCOMP+NBLEND)),com,1].*alpha_coeff;
        
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
            
            if j<=NCOMP+NBLEND %left jump
                
                if j>1 &&j<=NCOMP %left jump in pure compartment region
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    
                elseif j==NCOMP+1 %left jump over I1 boundary
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    CHOOSE_PART=find(part.*(((j-1)*H)<part).*(part<((j)*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART) = [];
                    
                elseif j>NCOMP+1
                    com(j)=com(j)-1;
                    com(j-1)=com(j-1)+1;
                    CHOOSE_PART=find(part.*(((j-1)*H)<part).*(part<((j)*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART) = part(RAND_PART)-H;
                else
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
                            part=[part,part_new];
                            check=1;
                        end
                    end
                                        
                elseif (2*NCOMP)+NBLEND<j && j~= (2*(NCOMP+NBLEND)) %right jump in blended region
                    com(j-(NCOMP+NBLEND))=com(j-(NCOMP+NBLEND))-1;
                    com(j+1-(NCOMP+NBLEND))=com(j+1-(NCOMP+NBLEND))+1;
                    CHOOSE_PART=find(part.*(((j-1-(NCOMP+NBLEND))*H)<part).*(part<((j-(NCOMP+NBLEND))*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART) = part(RAND_PART)+H;
                    
                end
            elseif j>2*(NCOMP+NBLEND) && j<=2*k+NCOMP %degredation
                if j<=2*(NCOMP+NBLEND)+(NCOMP)
                    com(j-(2*(NCOMP+NBLEND)))=com(j-(2*(NCOMP+NBLEND)))-1;
                else 
                    com(j-(2*(NCOMP+NBLEND)))=com(j-(2*(NCOMP+NBLEND)))-1;
                    CHOOSE_PART=find(part.*(((j-1-(2*(NCOMP+NBLEND)))*H)<part).*(part<((j-(2*(NCOMP+NBLEND)))*H)));%find particles in this compartment
                    RAND_PART = CHOOSE_PART(ceil(length(CHOOSE_PART)*rand));
                    part(RAND_PART) = [];

                end
            else
                com(1)=com(1)+1;
                
            end
                        
        elseif tb<tc
            t=tb;
            
    
            
            %calculate Db for each particle
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
            part=part+((sqrt(2*Db*dt)).*randn(1,length(part)))+(part<I(2))*((D/(I(2)-I(1)))*dt);
            
            %apply right hand boundary condition        
            part=part+2*(part>X(2)).*(X(2)-part);
           
            %apply left hand boundary condition
            part=part+2*(part<I(1)).*(I(1)-part);
            
            
            %error check for particles passing left hand boundary
            if sum(part<I(1))~=0
                disp 'Rogue particle in compartment region'
                return
            end
            
            %update the compartment numbers in the blended region and part
            %region
            for y=(NCOMP+1):k
                com(y)=sum((ceil(part*k))==y);
            end
            
            %reset tb
            tb=tb+dt;
        end
        
        while t>=(z-1)*time_int %for recording
            
            %find number of "compartments" in particle region
            
            rec_com_sum(z,:)=rec_com_sum(z,:)+ com;
            rec_part(z)=rec_part(z) + sum(com);
            z=z+1;
            
        end
        
    end
    PART=[PART,part]; %#ok
end

MG=rec_com_sum/repeats;
NUM_PART=rec_part/repeats;


Time_MG=toc;

save('./../../Plotting/Data_hybrid_compartment_brownian/TEST_PROBLEM_3/TEST3_BLENDING','MG')











