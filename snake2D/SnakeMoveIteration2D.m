function P=SnakeMoveIteration2D(B,P,Fext,gamma,kappa,delta,ExpdContFlag) 
% This function will calculate one iteration of contour Snake movement
%
% P=SnakeMoveIteration2D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)
%   gamma : Time step
%   kappa : External (image) field weight
%   delta : Balloon Force weight
%
% outputs,
%   P : The (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)
% B=S; 
% P=P(:,:,flagEvolving); 
% Fext=Fext; 
% gamma=Options.Gamma; 
% kappa=Options.Kappa; 
% delta=Options.Delta; 
% ExpdContFlag=ExpdContFlag; 

[~,~,n]=size(P); 
% Fext=max(Fext(:))-Fext; 
% Clamp contour to boundary
P(:,1,:)=min(max(P(:,1,:),1),size(Fext,1));
P(:,2,:)=min(max(P(:,2,:),1),size(Fext,2));

% Get image force on the contour points 
for ind_i = 1: n 
    Fext1(:,:,ind_i,1)=kappa*interp2(Fext(:,:,ind_i,1),P(:,2,ind_i),P(:,1,ind_i));
    Fext1(:,:,ind_i,2)=kappa*interp2(Fext(:,:,ind_i,2),P(:,2,ind_i),P(:,1,ind_i));
end 
% Interp2, can give nan's if contour close to border
Fext1(isnan(Fext1))=0;

% Calculate the baloonforce on the contour points
N=GetContourNormals2D(P);
Fext2=delta*N;

% Update contour positions
ssx = gamma*P(:,1,:) + ExpdContFlag* Fext1(:,:,:,1) + ExpdContFlag* Fext2(:,1,:); 
ssy = gamma*P(:,2,:) + ExpdContFlag* Fext1(:,:,:,2) + ExpdContFlag* Fext2(:,2,:); 
% ssx = gamma*P(:,1); 
% ssy = gamma*P(:,2); 
for ind_i = 1: n
    P(:,1,ind_i) = B * ssx(:,:,ind_i);
    P(:,2,ind_i) = B * ssy(:,:,ind_i);
end 

% Clamp contour to boundary
% P(:,1)=min(max(P(:,1),1),size(Fext,1));
% P(:,2)=min(max(P(:,2),1),size(Fext,2));

    
