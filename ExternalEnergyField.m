function [binField, EextField, px, py, u, v] = ExternalEnergyField (worldField, z0) 
%% External energy field (Gradient) %% 

binField=worldField; % fitst, we get the binary field for gradients 
binField(binField<=z0)=0; 
binField(binField>z0)=1; 

Sigma = sqrt(2); 
Ix=ImageDerivatives2D(binField,Sigma,'x');
Iy=ImageDerivatives2D(binField,Sigma,'y');
Ixx=ImageDerivatives2D(binField,Sigma,'xx');
Ixy=ImageDerivatives2D(binField,Sigma,'xy');
Iyy=ImageDerivatives2D(binField,Sigma,'yy');

Eline = imgaussian(binField,Sigma);
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((1+Ix.^2 + Iy.^2).^(3/2));
Eedge = sqrt(Ix.^2 + Iy.^2); % Gradients for external energy 

EextField = Eedge; 

px=0; 
py=0; 
u=0; 
v=0; 