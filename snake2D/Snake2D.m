function [PreorderLow, PreorderUp, P, Fext]=Snake2D(I, P, x0Temp, y0Temp, posIniDest, iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function SNAKE implements the basic snake segmentation. A snake is an 
% active (moving) contour, in which the points are attracted by edges and
% other boundaries. To keep the contour smooth, an membrame and thin plate
% energy is used as regularization.
%
% [O,J]=Snake2D(I,P,Options)
%  
% Literature: 
% - Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes Active Contour Models", 1987 
% - Jim Invins and John Porril, "Everything you always wanted to know about snakes (but were afraid to ask) 
% - Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New external force for Snakes 
% - Christoph Lurig, Leif Kobbelt, Thomas Ertl, "Hierachical solutions for the Deformable Surface Problem in Visualization"
% 
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%   
% outputs,
%   O : List with coordinates of the final contour M x 2
%   J : Binary image with the segmented region
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default 100
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
%  Options.Sigma1 : Sigma used to calculate image derivatives, default 10
%  Options.Wline : Attraction to lines, if negative to black lines otherwise white
%                    lines , default 0.04
%  Options.Wedge : Attraction to edges, default 2.0
%  Options.Wterm : Attraction to terminations of lines (end points) and
%                    corners, default 0.01
%  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
%                    image (which gives the image force), default 20
%
% options (Gradient Vector Flow)
%  Options.Mu : Trade of between real edge vectors, and noise vectors,
%                default 0.2. (Warning setting this to high >0.5 gives
%                an instable Vector Flow)
%  Options.GIterations : Number of GVF iterations, default 0
%  Options.Sigma3 : Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.2
%  Options.Delta : Baloon force, default 0.1
%  Options.Kappa : Weight of external image force, default 2
%
%
% Literature:
%   - Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes : Active
%       Contour Models", 1987
%   - Jim Ivins amd John Porrill, "Everything you always wanted to know
%       about snakes (but wer afraid to ask)
%   - Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New
%       external force for Snakes
%
% Example, Basic:
%

% writerObj=VideoWriter('SnakeMotion.avi');
% writerObj.FrameRate = 10;
% open(writerObj);

% % %   I = imread('test.png'); 
%   ind_u = 1:9; 
%   ind_p = 1; 
%   x0Temp=x0(:,:,ind_u); 
%   y0Temp=y0(:,:,ind_u); 
%   I = binField(:,:,ind_u);  
%   iteration = 100; 
%   P = trajectStart(:,:,ind_u); 
%   figure(100), imshow(I(:,:,ind_p)); view(0, 90); 
% % %   [y,x] = getpts;
% % %   P=[x(:) y(:)];
%   hold; 
%   plot3(x0Temp(ind_p), y0Temp(ind_p), ones(size(x0Temp(ind_p))), 'rs'); 
%   h=plot(P(:, 2, ind_p),P(:, 1, ind_p),'r.'); 
  
  PIni=P; 
  Options=struct;
  Options.Verbose=true;
  Options.Iterations=iteration;
  Options.Wedge=2;
  Options.Wline=0;
  Options.Wterm=0;
  Options.Kappa=4;
  Options.Sigma1=8;
  Options.Sigma2=8;
  Options.Alpha=0;
  Options.Beta=1;
  Options.Mu=0.2;
  Options.Delta=-0.1;
  Options.GIterations=600; 
  Options.Gamma = 1; 
  
%   [O,J]=Snake2D(I,P,Options);
%   
% Function is written by D.Kroon University of Twente (July 2010)

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',100,'Wline',0.04,'Wedge',2,'Wterm',0.01,'Sigma1',10,'Sigma2',20,'Alpha',0.2,'Beta',0.2,'Delta',0.1,'Gamma',1,'Kappa',2,'Iterations',100,'GIterations',0,'Mu',0.2,'Sigma3',1);
if(~exist('Options','var')) 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))) 
        warning('snake:unknownoption','unknown options found');
    end
end

% Convert input to double
I = double(I); 

% If color image convert to grayscale
% if(size(I,3)==3), I=rgb2gray(I); end

% The contour must always be clockwise (because of the balloon force)
P = MakeContourClockwise2D(P); 

% Make an uniform sampled contour description
P = InterpolateContourPoints2D(P, Options.nPoints); 

%% Field Inside Contour 
[Irow, Icol, n]=size(I); 
% [gridX, gridY]=meshgrid([1:Irow],[1:Icol]);
% [idx,odx] = inpolygon(gridX,gridY,PIni(:,1),PIni(:,2));
% odx=~idx; 
% odx=1; 
% IInside=I.*idx';
% field=(1+zeros(size(I))).*idx'; 

% origin=I;
% sigma0=3*sqrt(2);
% scale=sigma0*sqrt(2);
% f=fspecial('gaussian',[3,floor(3*scale)],scale);
% L1=double(origin);
% L2=conv2(origin,f,'same');
% L3=conv2(L2,f','same'); 
% DOG=L3; 

% Transform the Image into an External Energy Image
field=I; 
Eext = ExternalForceImage2D(field, Options.Wline, Options.Wedge, Options.Wterm, Options.Sigma1);

% Make the external force (flow) field.
Fx = ImageDerivatives2D(Eext,Options.Sigma2,'x'); 
Fy = ImageDerivatives2D(Eext,Options.Sigma2,'y'); 

% Fx=Fx.*odx'; 
% Fy=Fy.*odx'; 

Fext(:,:,:,1)=-Fx*2*Options.Sigma2^2; 
Fext(:,:,:,2)=-Fy*2*Options.Sigma2^2; 

FextOrig=Fext; 
% Fext=Fext.*idx'; 
% Fext=Fext.*(~I); 

% Do Gradient vector flow, optimalization 
% Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, Options.GIterations, Options.Sigma3); 
% Fext=Fext.*idx'; 
% Fext=GVFOptimizeImageForces2D(Fext, Options.Mu, 1, Options.Sigma3); 
% mu=.2;
% ITER=50;
% [px,py,u,v] = GVF(field, mu, ITER); 
% px=px.*odx'; 
% py=py.*odx'; 
% Fext=cat(3,px,py); 

% Make the interal force matrix, which constrains the moving points to a
% smooth contour
[b,A,S]=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);
% h=[]; 
ExpdContFlag = 1; % -1 for Expand; +1 for Contract 
dMinSet=[]; 

flagEvolving = [1: n]'; 
flagEvolved = []; 
swarmInd = [1: n]'; 
it=0; 
% for i=1: 1 %Options.Iterations 
while ~isempty(flagEvolving) 
    it=it+1; 
    itSet=it(:,:,ones(1, n));  
    
    P(:,:,flagEvolving) = SnakeMoveIteration2D(S,P(:,:,flagEvolving),Fext(:,:,flagEvolving,:), ... 
        Options.Gamma,Options.Kappa,Options.Delta,ExpdContFlag);

    % Show current contour 
%     eval(['title(''it=',num2str(it),''')']); 
%     h.XData = P(:, 2, ind_p);
%     h.YData = P(:, 1, ind_p);
%     pause(.1); 
    
%     if(Options.Verbose) 
%         if(ishandle(h)), delete(h), end 
%         h=plot(P(:, 2, ind_p),P(:, 1, ind_p),'r.'); 
%         c=i/Options.Iterations; 
%         plot([P(:, 2, ind_p);P(1, 2, ind_p)],[P(:, 1, ind_p);P(1, 1, ind_p)],'-','Color',[c 1-c 0]); 
%         drawnow
%     end 

%     [dMax, indMax] = max(sqrt((x0Temp-P(:,2,:)).^2+(y0Temp-P(:,1,:)).^2));
    [dMax, indMax] = min(sqrt((posIniDest(1)-P(:,2,:)).^2+(posIniDest(2)-P(:,1,:)).^2)); 
    [dMin, indMin] = min(sqrt((x0Temp-P(:,2,:)).^2+(y0Temp-P(:,1,:)).^2));
    dMinSet = cat(1, dMinSet, dMin); 
    
    thr = 1; 
%     flagEvolved = union(flagEvolved, swarmInd(dMin <= thr)); 
%     flagEvolving = setdiff(flagEvolving, swarmInd(dMin <= thr)); 
    
    if it > 1 
        stopFlag = dMinSet(end, :, :)-dMinSet(end-1, :, :); 
    else 
        stopFlag = -1; 
    end 
    flagEvolved = union(flagEvolved, swarmInd(stopFlag > .01 | dMin <= thr | itSet >= Options.Iterations));
    flagEvolving = setdiff(flagEvolving, flagEvolved); 
    
%         frame = getframe(gca);
%         writeVideo(writerObj,frame);
end 
PreorderLow=cell(1,1,n); 
PreorderUp=cell(1,1,n); 
for ind_i = 1: n 
    if indMin(ind_i) >= indMax(ind_i) 
        PreorderLow{ind_i} = [P(indMin(ind_i):end, [2 1], ind_i); P(1:indMax(ind_i), [2 1], ind_i)]; 
        PreorderUp{ind_i} = [P(indMin(ind_i):-1:indMax(ind_i), [2 1], ind_i)]; 
    else 
        PreorderLow{ind_i} = P(indMin(ind_i):indMax(ind_i), [2 1], ind_i); 
        PreorderUp{ind_i} = [P(indMin(ind_i):-1:1, [2 1], ind_i); P(end:-1:indMax(ind_i), [2 1], ind_i)]; 
    end 
end
Prevert = P(:, [2 1], :); 
% hold; 
% close(writerObj); 

% figure(10); 
% imshow(field); 
% hold; 
% plot([P(:,2);P(1,2)],[P(:,1);P(1,1)],'r-');  
% plot([PIni(:,2);PIni(1,2)],[PIni(:,1);PIni(1,1)],'b-'); 
% quiver(y,x,Fext(1:10:end,1:10:end,2),Fext(1:10:end,1:10:end,1)); 
% hold; 

% if(nargout>1)
%      J=DrawSegmentedArea2D(P,size(I));
% end

