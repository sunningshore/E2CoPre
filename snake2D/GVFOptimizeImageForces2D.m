function Fext=GVFOptimizeImageForces2D(Fext, Mu, Iterations, Sigma)
% This function "GVFOptimizeImageForces" does gradient vector flow (GVF)
% on a vector field. GVF gives the edge-forces a larger capature range,
% to make the snake also reach concave regions
%
% Fext = GVFOptimizeImageForces2D(Fext, Mu, Iterations, Sigma) 
% 
% inputs,
%   Fext : The image force vector field N x M x 2
%   Mu : Is a trade of scalar between noise and real edge forces
%   Iterations : The number of GVF itterations
%   Sigma : Used when calculating the Laplacian
% 
% outputs,
%   Fext : The GVF optimized image force vector field
%
% Function is written by D.Kroon University of Twente (July 2010)

% Iterations=100; 
% Mu=.2; 
% Fext=[]; 
% Sigma=2; 
% I=im2double(imread('test.png'));
% I=I(:,:,1); 
% Eext = ExternalForceImage2D(I,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
% Fx=ImageDerivatives2D(Eext,Sigma,'x'); 
% Fy=ImageDerivatives2D(Eext,Sigma,'y'); 
% Fext(:,:,1)=-Fx*2*Sigma^2; 
% Fext(:,:,2)=-Fy*2*Sigma^2; 

% Fext=FextOrig; 
% Iterations = Options.GIterations ; 
% Mu=Options.Mu;
% Sigma=Options.Sigma3;
% Squared magnitude of force field
Fx= Fext(:,:,1);
Fy= Fext(:,:,2);

% Calculate magnitude
sMag = Fx.^2+ Fy.^2;

% Set new vector-field to initial field
u=Fx;  v=Fy;
  
% Iteratively perform the Gradient Vector Flow (GVF)
for i=1:Iterations 
  % Boundary Condition 
%   u = BoundMirrorEnsure(u);
%   v = BoundMirrorEnsure(v);

  % Calculate Laplacian 
  Uxx=ImageDerivatives2D(u,Sigma,'xx');
  Uyy=ImageDerivatives2D(u,Sigma,'yy');
  
  Vxx=ImageDerivatives2D(v,Sigma,'xx');
  Vyy=ImageDerivatives2D(v,Sigma,'yy');

  % Update the vector field
  u = u + Mu*(Uxx+Uyy) - sMag.*(u-Fx);
  v = v + Mu*(Vxx+Vyy) - sMag.*(v-Fy);
end

Fext(:,:,1) = u;
Fext(:,:,2) = v;
%% Plotting %% 
% figure(10);
% imshow(I); 
% hold; 
% [x,y]=ndgrid(1:10:size(Fext,1),1:10:size(Fext,2));
% quiver(y,x,Fext(1:10:end,1:10:end,2),Fext(1:10:end,1:10:end,1)); 
% hold; 