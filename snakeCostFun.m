function [E, Eint, Eext, Evel, vecPre, vecCurr] = snakeCostFun (gamma1, gamma2, swarmIndAvoiding, posSwarmRel, posSwarmRelNext, posIniIndiv, x, y, z, zForbid, xPre, yPre, zPre, EextField, smoothField) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The cost function according to Snake energy. 
% E = Eint + Eext.
% 
% gamma1, gamma2: coeficients for internal and external energies 
% 
% swarmIndAvoiding: index of swarm members in avoidance 
% 
% posSwarmRel: previous position of the swarm center 
% posSwarmRelNext: current position of the swarm center 
% posIniIndiv: initial positions of swarm members 
% 
% xPre, yPre, zPre: previous trajectories 
% x, y, z: current trajectories 
% EextField: external energy field 
% 
% E, Eint, Eext, Eswm: Snake energy, internal energy, external energy and
% swarm formation energy (not used in this version)
% 
% vecPre, vecCurr: not used in this version 
[N,PopulationSize,swarmSize]=size(x); 
x0=ones(1,1,swarmSize); 
x0(1,1,:)=posIniIndiv(1,:); 
x0=x0(:,ones(PopulationSize,1),:); 
y0=ones(1,1,swarmSize); 
y0(1,1,:)=posIniIndiv(2,:); 
y0=y0(:,ones(PopulationSize,1),:); 
%% Previous traject %% 
[N,PopulationSize,swarmSize]=size(x); 

xPre=xPre(:,ones(PopulationSize,1),:); 
yPre=yPre(:,ones(PopulationSize,1),:); 
zPre=zPre(:,ones(PopulationSize,1),:); 

%% Overall traject = Previous traject + Current traject %% 
% x1=cat(1,xPre([1:end-10,end],:,:),x([10:end],:,:)); 
% y1=cat(1,yPre([1:end-10,end],:,:),y([10:end],:,:)); 
x1=cat(1,xPre([1:10:end],:,:),x([10:10:end],:,:)); 
y1=cat(1,yPre([1:10:end],:,:),y([10:10:end],:,:)); 
%% Internal energy = first order deriv (continuouty) + second order deriv (curvature) %% 
Econt=0; 
% Econt = sum((x1(1:end-1,:,:)-x1(2:end,:,:)).^2+(y1(1:end-1,:,:)-y1(2:end,:,:)).^2, 1); 
Ecurv = sum(sqrt((x1(3:end,:,:)-2*x1(2:end-1,:,:)+x1(1:end-2,:,:)).^2+(y1(3:end,:,:)-2*y1(2:end-1,:,:)+y1(1:end-2,:,:)).^2), 1); 

% Econt=0; 
% Econt = sum(abs((x1(1:end-1,:,:)-x1(2:end,:,:)))+abs((y1(1:end-1,:,:)-y1(2:end,:,:))), 1); 
% Ecurv = sum(abs((x1(3:end,:,:)-2*x1(2:end-1,:,:)+x1(1:end-2,:,:)))+abs((y1(3:end,:,:)-2*y1(2:end-1,:,:)+y1(1:end-2,:,:))), 1); 

Eint = Econt + Ecurv; 
% Eint = 0; 
%% External energy (Gradient) %% 
[S1,~,~]=size(x1); 
if size(EextField, 3) == 1 
    EextField=cat(3, EextField, EextField); 
end 
z1=ones(S1*PopulationSize, 1)*swarmIndAvoiding;
Eext1 = interp3(EextField,x1(:),y1(:),z1(:));
Eext1(isnan(Eext1))=0; 
Eext1 = reshape(Eext1, [S1, PopulationSize, swarmSize]); 
Eext = sum(abs(Eext1).^1, 1); 
E = gamma1*Eint - gamma2*Eext; 
%% External energy (Difference of Intensity, DOI) %% 
% [S1,~,~]=size(x1); 
% xBench=x(1,:,:); 
% yBench=y(1,:,:); 
% Eext0 = interp2(smoothField,xBench(:),yBench(:)); 
% Eext0 = reshape(Eext0, [1, PopulationSize, swarmSize]); 
% Eext0 = Eext0(ones(1, N), :, :); 
% Eext1 = interp2(smoothField,x(:),y(:)); 
% Eext1(isnan(Eext1))=0; 
% Eext1 = reshape(Eext1, [N, PopulationSize, swarmSize]); 
% % % Eext1(Eext1 > Eext0)=Inf; 
% % % Eext1(Eext1 > min(zForbid(:)))=Inf; 
% Eext = sum(abs(Eext0-Eext1).^1, 1); 
% % Eext = sum(abs(diff(Eext1)).^1, 1); 
% % Eext(Eext1(end,:,:) > Eext0(end,:,:))=Inf; 
% % % Eext(Eext1(end,:,:) > min(zForbid(:)))=Inf; 
% E = gamma1*Eint + gamma2*Eext; 

% Eswm=zeros(size(E)); 
% E = E + Eswm; 
%% Velocity Consistence Force %% 
Evel=sum(xPre.*yPre+x.*y, 1); 
% E = 1*E - .00001*Evel;
%% Swarm Formation Force - Vector Dot Product (Not in use in this version, kept for debugging) %% 
centPre=posSwarmRel; 
vecPre=[(x0-centPre(1)*ones(1,PopulationSize, swarmSize)); (y0-centPre(2)*ones(1,PopulationSize, swarmSize))]; 
vecPre=[vecPre; zeros(1,PopulationSize, swarmSize)]; 
centCurr=posSwarmRelNext; 
vecCurr=[(x(end,:,:)-centCurr(1)*ones(1,PopulationSize, swarmSize)); (y(end,:,:)-centCurr(2)*ones(1,PopulationSize, swarmSize))]; 
vecCurr=[vecCurr; zeros(1,PopulationSize, swarmSize)]; 