function [arcx, arcy, arcz, kBestInterp, omegaBestInterp, deltaSBestInterp, BestCosts, swarmVecPre, swarmVecCurr, swarmCost] = SpeciesPSO (gamma1, gamma2, swarmIndAvoiding, posIniIndiv, posIniSwarmRel, posIniSwarmRelNext, x0, y0, z0, zForbid, k, omega, deltaS, arcxPre, arcyPre, arczPre, smoothField, EextField, Maxiter, N, init_curv, init_omega, rnd_init) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Species based PSO 
% 
% gamma1, gamma2: coeficients for internal and external energies 
% 
% swarmIndAvoiding: index of swarm members in avoidance 
% 
% posIniIndiv: initial positions of swarm members 
% posSwarmRel: previous position of the swarm center 
% posSwarmRelNext: current position of the swarm center 
% 
% x0, y0: current positions of swarm members (starting points for trajectory planning)
% 
% k, omega, deltaS: three dimensions of PSO, curvature, slope angle and arc
% length (fix in this version) 
% 
% arcxPre, arcyPre, arczPre: previous trajectory 
% 
% smoothField: working space of UAVs 
% 
% EextField: external energy field 
% 
% Maxiter: iterations of PSO 
% 
% N: number of interpolation points of each arc 

[~, ~, swarmSize]=size(x0); 
searchBound1=size(k, 2); 
searchBound2=size(omega, 2); 
searchBound3=size(deltaS, 2); 
searchBound=[searchBound1, searchBound2, searchBound3]; 
%==========================================================================
% Parameters of PSO
%==========================================================================
PopulationSize = 200; % Enter Population Size (Swarm size)
N_variables = 3; % Enter number of Design Variables
searchBoundSet=searchBound(ones(PopulationSize, 1),:,ones(1, swarmSize)); 
%==========================================================================
% Constriction Coefficients
%==========================================================================
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

% w = chi;                    % Intertia Coefficient
% wdamp = 1;                  % Damping Ratio of Inertia Coefficient
c1 = chi*phi1;              % Personal Acceleration Coefficient
c2 = chi*phi2;              % Social Acceleration Coefficient

% Generate Random Solution for Particle Position 
speciesNum = swarmSize; 
if rnd_init
    temp1=ones(PopulationSize, 1)*(searchBound - 1); 
    rand_seed = rand(1, N_variables, speciesNum); 
    Particle_Position = 1 + temp1(:,:,ones(speciesNum, 1)).*rand_seed(ones(PopulationSize, 1), :, :); 
else
    Particle_Position = zeros(PopulationSize, N_variables, speciesNum); 
    for agent = 1: speciesNum 
        init_omega_temp = init_omega(agent); 
        init_curv_temp = init_curv(agent); 
        Particle_Position(:,:,agent) = [normrnd(init_curv_temp, 1, [PopulationSize, 1]), ...
            normrnd(init_omega_temp, 1, [PopulationSize, 1]), ...
            normrnd(0, 1, [PopulationSize, 1])]; 
    end 
end 

Particle_Position = max(Particle_Position, 1);
% Particle_Position = min(Particle_Position, searchBound1);
Particle_Position = min(Particle_Position, searchBoundSet);
    
% Generate Random Solution for Particle Position 
Particle_Velocity = rand(PopulationSize, N_variables, speciesNum);

Particle_Cost=zeros(PopulationSize, 1, speciesNum) + inf; 
% Particle_Cost=sum(Particle_Cost, 3); 
% Particle_Cost=Particle_Cost(:,:,ones(1, speciesNum)); 

% Update the Personal Best 
Particle_BestPosition = Particle_Position;
Particle_BestCost = Particle_Cost;

% Finding Global Best 
[M,I] = min(Particle_BestCost,[],1); % M is GlobalBest_Cost 
GlobalBest = zeros(1, N_variables, speciesNum); 
GlobalBestPlt = zeros(1, N_variables, speciesNum); 
for ind_i=1: speciesNum 
    GlobalBest(:,:,ind_i) = Particle_BestPosition(I(ind_i), :, ind_i); 
    GlobalBestPlt(:,:,ind_i) = Particle_BestPosition(I(ind_i), :, ind_i); 
end 

% Array to Hold Best Cost Value on Each Iteration 
BestCosts = zeros(Maxiter, 1, speciesNum);

dg=zeros(Maxiter,1,speciesNum); 
dmin=zeros(Maxiter,1,speciesNum); 
dmax=zeros(Maxiter,1,speciesNum); 
evolFactor=zeros(Maxiter,1,speciesNum); 
w=zeros(Maxiter,1,speciesNum); 
w(1,:,:)=0.9; 
swarmCost=zeros(1,1, speciesNum); 
swarmVecPre=zeros(3,1, speciesNum); 
swarmVecCurr=zeros(3,1, speciesNum); 
%==========================================================================
% Main Loop of PSO 
%==========================================================================
for it=1: Maxiter 
        
    Particle_Velocity = w(it,:,:).*Particle_Velocity...
        + c1*rand(PopulationSize, N_variables, speciesNum).*(Particle_BestPosition - Particle_Position) ...
        + c2*rand(PopulationSize, N_variables, speciesNum).*(GlobalBest(ones(PopulationSize,1),:,:) - Particle_Position); 
    Particle_Position = Particle_Position + Particle_Velocity; 
    
    Particle_Position = max(Particle_Position, 1);
%     Particle_Position = min(Particle_Position, searchBound1);
    Particle_Position = min(Particle_Position, searchBoundSet);

    ParticleXLeft = Particle_Position(:,1,:); 
    ParticleYLeft = Particle_Position(:,2,:); 
    ParticleXRight = permute(Particle_Position(:,1,:), [2 1 3]); 
    ParticleYRight = permute(Particle_Position(:,2,:), [2 1 3]); 
    dParticle = sum(sqrt((ParticleXLeft(:, ones(1, PopulationSize), :)-ParticleXRight(ones(1, PopulationSize), :, :)).^2 + (ParticleYLeft(:, ones(1, PopulationSize), :)-ParticleYRight(ones(1, PopulationSize), :, :)).^2), 2)./(PopulationSize); 
    
    kInterp=zeros(1,PopulationSize,speciesNum);
    omegaInterp=zeros(1,PopulationSize,speciesNum);
    deltaSInterp=zeros(1,PopulationSize,speciesNum);
    for ind_i=1: speciesNum
        kInterp(:,:,ind_i)=interp1(k(:,:,ind_i),Particle_Position(:,2,ind_i));
        omegaInterp(:,:,ind_i)=interp1(omega(:,:,ind_i),Particle_Position(:,1,ind_i));
        deltaSInterp(:,:,ind_i)=interp1(deltaS(:,:,ind_i),Particle_Position(:,3,ind_i));
%         deltaSInterp(:,:,ind_i)=5*ones(size(Particle_Position(:,2,ind_i)));
    end
    [arcxTemp, arcyTemp, x1, y1]=arcParameterization (x0, y0, kInterp, omegaInterp, deltaSInterp, N); 
%     arczTemp = interp2(smoothField, arcxTemp, arcyTemp); 
    z0=z0(ones(1, N), ones(1, PopulationSize), :); 
%     arczTemp(arczTemp>z0) = z0; 
    arczTemp = z0; 
    
    [E, Eint, Eext, Eswm, vecPre, vecCurr] = snakeCostFun (gamma1, gamma2, swarmIndAvoiding, posIniSwarmRel, posIniSwarmRelNext, posIniIndiv, arcxTemp, arcyTemp, arczTemp, zForbid, arcxPre, arcyPre, arczPre, EextField, smoothField);
    E=permute(E, [2 1 3]); 
    
    Particle_Cost=E; 
%     Particle_Cost=sum(Particle_Cost, 3); 
%     Particle_Cost=Particle_Cost(:,:,ones(1, speciesNum)); 
    
    for ind_i=1: speciesNum
        Particle_BestPosition(Particle_Cost(:,:,ind_i) < Particle_BestCost(:,:,ind_i),:,ind_i) = Particle_Position(Particle_Cost(:,:,ind_i) < Particle_BestCost(:,:,ind_i),:,ind_i);
        Particle_BestCost(Particle_Cost(:,:,ind_i) < Particle_BestCost(:,:,ind_i), :, ind_i) = Particle_Cost(Particle_Cost(:,:,ind_i) < Particle_BestCost(:,:,ind_i),:,ind_i);
    end 
    
    % Update Global Best
    [M,I] = min(Particle_BestCost,[],1); 
    for ind_i=1: speciesNum
        GlobalBest(:,:,ind_i) = Particle_BestPosition(I(ind_i), :, ind_i);
        GlobalBestPlt(:,:,ind_i) = Particle_BestPosition(I(ind_i), :, ind_i); 
        swarmCost(:,:,ind_i)=Eswm(:,I(ind_i),ind_i); 
        swarmVecPre(:,:,ind_i)=vecPre(:,I(ind_i),ind_i); 
        swarmVecCurr(:,:,ind_i)=vecCurr(:,I(ind_i),ind_i); 
    end 
    
    mu=0; 
    
    sigma=1; 
    Gaussrnd=normrnd(mu,sigma^2); 
    GlobalBest=GlobalBest+searchBound(:,:,ones(speciesNum, 1))/10.*Gaussrnd; 
    
    GlobalBest = max(GlobalBest, 1);
    GlobalBest = min(GlobalBest, searchBound);
    
    dg(it,:,:)=sum(sqrt(sum((GlobalBest(ones(PopulationSize,1),:,:)-Particle_Position).^2,2)), 1)./(PopulationSize);
    dmin(it,:,:)=min(dParticle); 
    dmax(it,:,:)=max(dParticle); 
    evolFactor(it,:,:)=(dg(it,:,:)-dmin(it,:,:))./(dmax(it,:,:)-dmin(it,:,:)); 
    
    w(it,:,:)=1./(1+1.5*(10.^(-2.6.*evolFactor(it,:,:))));
    
    % Store the Best Cost Value
    BestCosts(it,:,:) = M(:,1,:);
    
    kBestInterp=zeros(1,1,speciesNum);
    omegaBestInterp=zeros(1,1,speciesNum);
    deltaSBestInterp=zeros(1,1,speciesNum);
    for ind_i=1: speciesNum
        kBestInterp(:,:,ind_i)=interp1(k(:,:,ind_i),GlobalBestPlt(:,2,ind_i));
        omegaBestInterp(:,:,ind_i)=interp1(omega(:,:,ind_i),GlobalBestPlt(:,1,ind_i));
        deltaSBestInterp(:,:,ind_i)=interp1(deltaS(:,:,ind_i),GlobalBestPlt(:,3,ind_i));
%         deltaSBestInterp(:,:,ind_i)=5*ones(size(GlobalBestPlt(:,2,ind_i))); 
    end 
    [arcx, arcy, x1, y1]=arcParameterization (x0, y0, kBestInterp, omegaBestInterp, deltaSBestInterp, N);
    arcz = interp2(smoothField, arcx, arcy); 
end 

kNew=zeros(1,1,speciesNum); 
omegaNew=zeros(1,1,speciesNum); 
deltaSNew=zeros(1,1,speciesNum); 
for ind_i=1: speciesNum 
    kNew(:,:,ind_i)=interp1(k(:,:,ind_i),GlobalBestPlt(:,2,ind_i)); 
    omegaNew(:,:,ind_i)=interp1(omega(:,:,ind_i),GlobalBestPlt(:,1,ind_i)); 
    deltaSNew(:,:,ind_i)=interp1(deltaS(:,:,ind_i),GlobalBestPlt(:,3,ind_i)); 
end 