function [x, y, x1, y1]=arcParameterization (x0, y0, k, omega, deltaS, N) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw arcs according to start point, curvature, omega angle, arc length
% and interpolating points. 
% 
% x0, y0: starting points 
% k, omega, deltaS: curvature, slope angle, arc length 
% N: number of interpolating points 
% 
% x, y: arc 
% x1, y1: tangent vector 
[~, PopulationSize, swarmSize]=size(k); 
x0=x0(:,ones(PopulationSize,1),:); 
y0=y0(:,ones(PopulationSize,1),:); 
%% Calculating %% 
% N=100; % number of interpolation points. 
x1=deltaS.*cosd(omega)+x0;
y1=deltaS.*sind(omega)+y0; 

x=zeros(N, PopulationSize* swarmSize); 
y=zeros(N, PopulationSize* swarmSize); 

linSpcInd = linspace(0,1,N); 
alpha1=-(90-omega); 
theta10=180+alpha1;
theta11=theta10(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta1=ones(N,1)*theta11(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta10(:)-theta11(:))'); 

alpha2=90+omega; 
theta20=-(180-alpha2);
theta21=theta20(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta2=ones(N,1)*theta20(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta21(:)-theta20(:))'); 

alpha3=-(90-omega); 
theta30=-(180-alpha3);
theta31=theta30(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta3=ones(N,1)*theta31(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta30(:)-theta31(:))'); 

alpha4=-(270-omega); 
theta40=180+alpha4;
theta41=theta40(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta4=ones(N,1)*theta40(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta41(:)-theta40(:))'); 

alpha5=-(90-omega); 
theta50=-(180-alpha5);
theta51=theta50(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta5=ones(N,1)*theta51(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta50(:)-theta51(:))'); 

alpha6=-(270-omega); 
theta60=(180+alpha6);
theta61=theta60(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta6=ones(N,1)*theta60(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta61(:)-theta60(:))'); 

alpha7=-(450-omega);  
theta70=(180+alpha7);
theta71=theta70(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta7=ones(N,1)*theta71(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta70(:)-theta71(:))'); 

alpha8=-(270-omega);  
theta80=-(180-alpha8);
theta81=theta80(:)-k(:).*deltaS(:)/pi*180; 
deltaTheta8=ones(N,1)*theta80(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(theta81(:)-theta80(:))'); 

xTemp=ones(N,1)*x0(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(x1(:)-x0(:))'); 
yTemp=ones(N,1)*y0(:)'+linSpcInd'*ones(1,length(k(:))).*(ones(N,1)*(y1(:)-y0(:))'); 
if ~isempty (x(:, k==0)) 
    x(:, k==0)=xTemp(:,k==0);
    y(:, k==0)=yTemp(:,k==0);
end 
if ~isempty (x(:, k>0 & omega<=90)) 
    x(:, k>0 & omega<=90)= ones(N,1)*abs(1./k(:, k>0 & omega<=90)).*cosd(deltaTheta1(:, k>0 & omega<=90))+ones(N,1)*(abs(1./k(:, k>0 & omega<=90)).*cosd(alpha1(:, k>0 & omega<=90)))+ones(N,1)*x0(:, k>0 & omega<=90);
    y(:, k>0 & omega<=90)= ones(N,1)*abs(1./k(:, k>0 & omega<=90)).*sind(deltaTheta1(:, k>0 & omega<=90))+ones(N,1)*(abs(1./k(:, k>0 & omega<=90)).*sind(alpha1(:, k>0 & omega<=90)))+ones(N,1)*y0(:, k>0 & omega<=90);
end
if ~isempty (x(:, k<0 & omega<=90))
    x(:, k<0 & omega<=90)= ones(N,1)*abs(1./k(:, k<0 & omega<=90)).*cosd(deltaTheta2(:, k<0 & omega<=90))+ones(N,1)*(abs(1./k(:, k<0 & omega<=90)).*cosd(alpha2(:, k<0 & omega<=90)))+ones(N,1)*x0(:, k<0 & omega<=90);
    y(:, k<0 & omega<=90)= ones(N,1)*abs(1./k(:, k<0 & omega<=90)).*sind(deltaTheta2(:, k<0 & omega<=90))+ones(N,1)*(abs(1./k(:, k<0 & omega<=90)).*sind(alpha2(:, k<0 & omega<=90)))+ones(N,1)*y0(:, k<0 & omega<=90);
end
if ~isempty (x(:, omega <= 180 & omega > 90 & k > 0))
    x(:, omega <= 180 & omega > 90 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 180 & omega > 90 & k > 0)).*cosd(deltaTheta3(:, omega <= 180 & omega > 90 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 180 & omega > 90 & k > 0)).*cosd(alpha3(:, omega <= 180 & omega > 90 & k > 0)))+ones(N,1)*x0(:, omega <= 180 & omega > 90 & k > 0);
    y(:, omega <= 180 & omega > 90 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 180 & omega > 90 & k > 0)).*sind(deltaTheta3(:, omega <= 180 & omega > 90 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 180 & omega > 90 & k > 0)).*sind(alpha3(:, omega <= 180 & omega > 90 & k > 0)))+ones(N,1)*y0(:, omega <= 180 & omega > 90 & k > 0);
end
if ~isempty (x(:, omega <= 180 & omega > 90 & k < 0))
    x(:, omega <= 180 & omega > 90 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 180 & omega > 90 & k < 0)).*cosd(deltaTheta4(:, omega <= 180 & omega > 90 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 180 & omega > 90 & k < 0)).*cosd(alpha4(:, omega <= 180 & omega > 90 & k < 0)))+ones(N,1)*x0(:, omega <= 180 & omega > 90 & k < 0);
    y(:, omega <= 180 & omega > 90 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 180 & omega > 90 & k < 0)).*sind(deltaTheta4(:, omega <= 180 & omega > 90 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 180 & omega > 90 & k < 0)).*sind(alpha4(:, omega <= 180 & omega > 90 & k < 0)))+ones(N,1)*y0(:, omega <= 180 & omega > 90 & k < 0);
end
if ~isempty (x(:, omega <= 270 & omega > 180 & k > 0))
    x(:, omega <= 270 & omega > 180 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 270 & omega > 180 & k > 0)).*cosd(deltaTheta5(:, omega <= 270 & omega > 180 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 270 & omega > 180 & k > 0)).*cosd(alpha5(:, omega <= 270 & omega > 180 & k > 0)))+ones(N,1)*x0(:, omega <= 270 & omega > 180 & k > 0);
    y(:, omega <= 270 & omega > 180 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 270 & omega > 180 & k > 0)).*sind(deltaTheta5(:, omega <= 270 & omega > 180 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 270 & omega > 180 & k > 0)).*sind(alpha5(:, omega <= 270 & omega > 180 & k > 0)))+ones(N,1)*y0(:, omega <= 270 & omega > 180 & k > 0);
end
if ~isempty (x(:, omega <= 270 & omega > 180 & k < 0))
    x(:, omega <= 270 & omega > 180 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 270 & omega > 180 & k < 0)).*cosd(deltaTheta6(:, omega <= 270 & omega > 180 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 270 & omega > 180 & k < 0)).*cosd(alpha6(:, omega <= 270 & omega > 180 & k < 0)))+ones(N,1)*x0(:, omega <= 270 & omega > 180 & k < 0);
    y(:, omega <= 270 & omega > 180 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 270 & omega > 180 & k < 0)).*sind(deltaTheta6(:, omega <= 270 & omega > 180 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 270 & omega > 180 & k < 0)).*sind(alpha6(:, omega <= 270 & omega > 180 & k < 0)))+ones(N,1)*y0(:, omega <= 270 & omega > 180 & k < 0);
end
if ~isempty (x(:, omega <= 359 & omega > 270 & k > 0))
    x(:, omega <= 359 & omega > 270 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 359 & omega > 270 & k > 0)).*cosd(deltaTheta7(:, omega <= 359 & omega > 270 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 359 & omega > 270 & k > 0)).*cosd(alpha7(:, omega <= 359 & omega > 270 & k > 0)))+ones(N,1)*x0(:, omega <= 359 & omega > 270 & k > 0);
    y(:, omega <= 359 & omega > 270 & k > 0)= ones(N,1)*abs(1./k(:, omega <= 359 & omega > 270 & k > 0)).*sind(deltaTheta7(:, omega <= 359 & omega > 270 & k > 0))+ones(N,1)*(abs(1./k(:, omega <= 359 & omega > 270 & k > 0)).*sind(alpha7(:, omega <= 359 & omega > 270 & k > 0)))+ones(N,1)*y0(:, omega <= 359 & omega > 270 & k > 0);
end
if ~isempty (x(:, omega <= 359 & omega > 270 & k < 0)) 
    x(:, omega <= 359 & omega > 270 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 359 & omega > 270 & k < 0)).*cosd(deltaTheta8(:, omega <= 359 & omega > 270 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 359 & omega > 270 & k < 0)).*cosd(alpha8(:, omega <= 359 & omega > 270 & k < 0)))+ones(N,1)*x0(:, omega <= 359 & omega > 270 & k < 0);
    y(:, omega <= 359 & omega > 270 & k < 0)= ones(N,1)*abs(1./k(:, omega <= 359 & omega > 270 & k < 0)).*sind(deltaTheta8(:, omega <= 359 & omega > 270 & k < 0))+ones(N,1)*(abs(1./k(:, omega <= 359 & omega > 270 & k < 0)).*sind(alpha8(:, omega <= 359 & omega > 270 & k < 0)))+ones(N,1)*y0(:, omega <= 359 & omega > 270 & k < 0);
end

x(:,x(1,:)~=x0(:)' | y(1,:)~=y0(:)')=flipud(x(:,x(1,:)~=x0(:)' | y(1,:)~=y0(:)')); 
y(:,x(1,:)~=x0(:)' | y(1,:)~=y0(:)')=flipud(y(:,x(1,:)~=x0(:)' | y(1,:)~=y0(:)')); 

x=reshape(x, [N, PopulationSize, swarmSize]); 
y=reshape(y, [N, PopulationSize, swarmSize]); 
