function [px,py,u,v] = GVF(f, mu, ITER) 
%GVF Compute gradient vector flow.
%   [u,v] = GVF(f, mu, ITER) computes the
%   GVF of an edge map f.  mu is the GVF regularization coefficient
%   and ITER is the number of iterations that will be computed.  

%   Chenyang Xu and Jerry L. Prince 6/17/97
%   Copyright (c) 1996-99 by Chenyang Xu and Jerry L. Prince
%   Image Analysis and Communications Lab, Johns Hopkins University

%   modified on 9/9/99 by Chenyang Xu
%   MATLAB do not deal their boundary condition for gradient and del2 
%   consistently between MATLAB 4.2 and MATLAB 5. Hence I modify
%   the function to take care of this issue by the code itself.
%   Also, in the previous version, the input "f" is assumed to have been
%   normalized to the range [0,1] before the function is called. 
%   In this version, "f" is normalized inside the function to avoid 
%   potential error of inputing an unnormalized "f".

% clear;
% clc; 
% mu=0.2; 
% ITER=100;  
% f=imread('example.bmp'); 
% f=rgb2gray(f); 
% f=im2double(imread('test.png'));
% f=f(:,:,1); 
% load('weight_EC_cube.mat','weight_EC_cube'); 
% field=weight_EC_cube; 

[frmRow0,frmCol0,frmNum]=size(f); 
pxSet=zeros(frmRow0,frmCol0,frmNum); 
pySet=zeros(frmRow0,frmCol0,frmNum); 

% for ind_frm=1: frmNum 
    
% f0=imresize(field(:,:,ind_frm),1); 
% f=f0; 

% [pxG,pyG] = gradient(-double(f)); 
% magG = sqrt(pxG.*pxG+pyG.*pyG);
% pxGU = pxG./(magG+1e-10); pyGU = pyG./(magG+1e-10);

[m,n] = size(f);
fmin  = min(f(:));
fmax  = max(f(:));
f = (f-fmin)/(fmax-fmin);  % Normalize f to the range [0,1]

f = BoundMirrorExpand(f);  % Take care of boundary condition
[fx,fy] = gradient(f);     % Calculate the gradient of the edge map
u = fx; v = fy;            % Initialize GVF to the gradient
SqrMagf = fx.*fx + fy.*fy; % Squared magnitude of the gradient field

% Iteratively solve for the GVF u,v
for i=1:ITER 
  u = BoundMirrorEnsure(u);
  v = BoundMirrorEnsure(v);
  u = u + mu*4*del2(u) - SqrMagf.*(u-fx);
  v = v + mu*4*del2(v) - SqrMagf.*(v-fy);
%   fprintf(1, '%3d', i);
%   if (rem(i,20) == 0)
%      fprintf(1, '\n');
%   end 
end
% fprintf(1, '\n');

u = BoundMirrorShrink(u); 
v = BoundMirrorShrink(v); 
% disp(' Nomalizing the GVF external force ...');
mag = sqrt(u.*u+v.*v);
px = u./(mag+1e-10); py = v./(mag+1e-10);
% px=u; 
% py=v; 

% pxSet(:,:,ind_frm)=px; 
% pySet(:,:,ind_frm)=py; 
%% Plotting %% 
% figure(10);
% imshow(I); 
% hold; 
% [x,y]=ndgrid(1:10:size(px,1),1:10:size(px,2));
% quiver(y,x,px(1:10:end,1:10:end),px(1:10:end,1:10:end)); 
% hold; 
% % figure (1); 
% % imshow(f); 
% % hold; 
% % set(gca,'ydir','normal'); 
% h=figure (1); 
% mesh(f0); hold; 
% alpha .1 
% quiver(gca, px,py,'b'); 
% hold; 
% axis off; axis equal; axis 'ij';     % fix the axis
% eval(['title(''Normalized GVF Field frm ',num2str(ind_frm),''');' ]); 
% 
% g=figure (2); 
% mesh(f0); hold; 
% alpha .1 
% quiver(gca, pxGU,pyGU,'r'); 
% hold; 
% axis off; axis equal; axis 'ij';     % fix the axis
% eval(['title(''Normalized Gradient Field frm ',num2str(ind_frm),''');' ]); 
% 
% f=figure(3); 
% mesh(f0);
% axis off; axis equal; axis 'ij';     % fix the axis
% eval(['title(''Potential Field frm ',num2str(ind_frm),''');' ]); 
% pause(.1); 
% % hold; 
% end 
% save('GVF.mat','pxSet','pySet'); 