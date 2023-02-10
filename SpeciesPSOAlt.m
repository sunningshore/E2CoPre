function [tBestInterp] = SpeciesPSOAlt (xEnd, yEnd, zEnd, AltEnd, t, guardDistv2v, Maxiter) 
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
[searchBound1, ~, N_variables]=size(t); 
searchBound=searchBound1; 
%==========================================================================
% Parameters of PSO
%==========================================================================
PopulationSize = 100; % Enter Population Size (Swarm size)
% N_variables = 1; % Enter number of Design Variables
searchBoundSet=searchBound(ones(PopulationSize, 1),ones(N_variables, 1)); 
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
% temp1=ones(PopulationSize, 1)*(searchBound - 1); 
Particle_Position = searchBound/2 + rand(PopulationSize, N_variables); 
Particle_Position = max(Particle_Position, 1);
Particle_Position = min(Particle_Position, searchBoundSet);
    
% Generate Random Solution for Particle Position 
Particle_Velocity = rand(PopulationSize,N_variables);

Particle_Cost = zeros(PopulationSize, 1)+Inf; 

% Update the Personal Best 
Particle_BestPosition = Particle_Position;
Particle_BestCost = Particle_Cost;

% Finding Global Best 
[M,I] = min(Particle_BestCost,[],1); % M is GlobalBest_Cost 

GlobalBest = Particle_BestPosition(I, :); 
GlobalBestPlt = Particle_BestPosition(I, :); 

% Array to Hold Best Cost Value on Each Iteration 
BestCosts = zeros(Maxiter, 1); 

% dg=zeros(Maxiter,1); 
% dmin=zeros(Maxiter,1); 
% dmax=zeros(Maxiter,1); 
% evolFactor=zeros(Maxiter,1); 
w=zeros(Maxiter,1); 
% w(1,:,:)=0.9; 
wdamp=0.99; 

% figure(100); 
% h=plot3(Particle_Position(:,1), Particle_Position(:,2), Particle_Position(:,3), '.'); 
% axis([0, searchBound, 0, searchBound, 0, searchBound]); 
%==========================================================================
% Main Loop of PSO 
%==========================================================================
for it=1: Maxiter 
    
    if it == 1 
        w(it,:,:) = 0.9; 
    else 
        w(it,:,:) = w(it-1,:,:)*wdamp;
    end 
        
    Particle_Velocity = w(it,:,:).*Particle_Velocity...
        + c1*rand(PopulationSize, N_variables).*(Particle_BestPosition - Particle_Position) ...
        + c2*rand(PopulationSize, N_variables).*(GlobalBest(ones(1, PopulationSize),:) - Particle_Position); 
    Particle_Position = Particle_Position + Particle_Velocity; 
    
    Particle_Position = max(Particle_Position, 1);
    Particle_Position = min(Particle_Position, searchBoundSet);
    
%     dParticle = 0; 
%     for ind_i=1: N_variables 
%         ParticleXLeft = Particle_Position(:,ind_i);
%         ParticleXRight = Particle_Position(:,ind_i)'; 
%         dParticle=dParticle+(ParticleXLeft(:, ones(1, PopulationSize))-ParticleXRight(ones(1, PopulationSize), :)).^2;
%     end 
%     dParticle = sum(sqrt(dParticle), 2)./(PopulationSize); 
    
    tInterp=zeros(PopulationSize,N_variables);
    for ind_i=1: N_variables 
        tInterp(:,ind_i)=interp1(t(:,ind_i),Particle_Position(:,ind_i)); 
    end 
%% Swarm Energy (V2V Distance) %% 
    combInd=nchoosek([1: N_variables],2); 
    AltTemp=tInterp+AltEnd(ones(PopulationSize, 1), :); 
    xTemp=xEnd(ones(PopulationSize, 1), :); 
    yTemp=yEnd(ones(PopulationSize, 1), :); 
    dSwmXY=(xTemp(:,combInd(:,1))-xTemp(:,combInd(:,2))).^2 + ... 
        (yTemp(:,combInd(:,1))-yTemp(:,combInd(:,2))).^2; 
    dSwm=sqrt(dSwmXY+(AltTemp(:,combInd(:,1))-AltTemp(:,combInd(:,2))).^2); 
%     dSwm=sqrt((xTemp(:,combInd(:,1))-xTemp(:,combInd(:,2))).^2 + ... 
%         (yTemp(:,combInd(:,1))-yTemp(:,combInd(:,2))).^2 + ... 
%         (AltTemp(:,combInd(:,1))-AltTemp(:,combInd(:,2))).^2); 
    EswmMin=min(dSwm, [], 2); 
    Eswm=EswmMin; 
    Eswm(Eswm <= guardDistv2v)=Inf; 
    E = 1*Eswm; 
%% Gravity Energy (Gravity) %% 
    Eg=sum(tInterp.^2, 2); 
    E = E + Eg; 
    
%     E = sum(abs(tInterp), 2); 
    
    Particle_Cost=E; 
    
    Particle_BestPosition(Particle_Cost < Particle_BestCost, :) = Particle_Position(Particle_Cost < Particle_BestCost, :);
    Particle_BestCost(Particle_Cost < Particle_BestCost) = Particle_Cost(Particle_Cost < Particle_BestCost);
    
    % Update Global Best
    [M,I] = min(Particle_BestCost,[],1); 
    GlobalBest = Particle_BestPosition(I,:);
    GlobalBestPlt = Particle_BestPosition(I,:); 
    
    mu=0; 
    
%     sigma=1; 
%     Gaussrnd=normrnd(mu,sigma^2); 
%     GlobalBest=GlobalBest+searchBound/10.*Gaussrnd; 
    
    GlobalBest = max(GlobalBest, 1);
    GlobalBest = min(GlobalBest, searchBound);
    
%     dg(it,:)=sum(sqrt(sum((GlobalBest(ones(PopulationSize,1),:)-Particle_Position).^2,2)), 1)./(PopulationSize);
%     dmin(it,:)=min(dParticle); 
%     dmax(it,:)=max(dParticle); 
%     evolFactor(it,:)=(dg(it,:)-dmin(it,:))./(dmax(it,:)-dmin(it,:)); 
%     w(it,:,:)=1./(1+1.5*(10.^(-2.6.*evolFactor(it,:)))); 
    
    % Store the Best Cost Value
    BestCosts(it,:) = M; 
    
    tBestInterp=zeros(1,N_variables);
    for ind_i=1: N_variables
        tBestInterp(:,ind_i)=interp1(t(:,ind_i), GlobalBestPlt(ind_i));
    end 
%     [arcx, arcy, x1, y1]=arcParameterization (x0, y0, kBestInterp, omegaBestInterp, deltaSBestInterp, N);
%     arcz = interp2(smoothField, arcx, arcy); 
    
%     pause(.1); 
%     h.XData= Particle_Position(:,1); 
%     h.YData= Particle_Position(:,2); 
%     h.ZData= Particle_Position(:,3); 
end 