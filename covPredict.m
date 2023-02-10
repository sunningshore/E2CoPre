function [xPred, yPred, P] = snakePredict (iteration, binField, trajectStart, ...
    x0, y0, posIniDest, arcxPre, arcyPre, arczPre, trajectSet) 
%% Prediction %% 
[~,~,n] = size(x0); 
N = 100; 
[PreorderLow, PreorderUp, P, Fext] = Snake2D(binField, trajectStart, x0, y0, posIniDest, iteration);

% [xPoly,yPoly] = curvInterp (boundContReg(:,2)',boundContReg(:,1)', N); 
stepPred = 20;
PolyLowPred = []; 
PolyUpPred = []; 
for ind_i = 1: n 
    PreorderLowTemp = PreorderLow{ind_i}; 
    PreorderUpTemp = PreorderUp{ind_i}; 
    [xPolyLow,yPolyLow] = curvInterp (PreorderLowTemp(:,2)',PreorderLowTemp(:,1)', N);
    [xPolyUp,yPolyUp] = curvInterp (PreorderUpTemp(:,2)',PreorderUpTemp(:,1)', N);
    
    if length(xPolyLow) >= stepPred*N 
        xPolyLowPred = xPolyLow(1 : stepPred*N)'; 
        yPolyLowPred = yPolyLow(1 : stepPred*N)'; 
    elseif (~isempty(xPolyLow)) && (length(xPolyLow) < stepPred*N )
        xPolyLowPred = [xPolyLow, xPolyLow(end)*ones(1, stepPred*N-length(xPolyLow))]'; 
        yPolyLowPred = [yPolyLow, yPolyLow(end)*ones(1, stepPred*N-length(xPolyLow))]'; 
    else 
        xPolyLowPred = x0(ind_i)*ones(1, stepPred*N)'; 
        yPolyLowPred = y0(ind_i)*ones(1, stepPred*N)'; 
    end 
    if length(xPolyUp) >= stepPred*N
        xPolyUpPred = xPolyUp(1 : stepPred*N)'; 
        yPolyUpPred = yPolyUp(1 : stepPred*N)'; 
    elseif (~isempty(xPolyUp)) && (length(xPolyUp) < stepPred*N )
        xPolyUpPred = [xPolyUp, xPolyUp(end)*ones(1, stepPred*N-length(xPolyUp))]'; 
        yPolyUpPred = [yPolyUp, yPolyUp(end)*ones(1, stepPred*N-length(xPolyUp))]'; 
    else
        xPolyUpPred = x0(ind_i)*ones(1, stepPred*N)'; 
        yPolyUpPred = y0(ind_i)*ones(1, stepPred*N)'; 
    end 
    
    PolyLowPred = cat(3, PolyLowPred, [xPolyLowPred, yPolyLowPred]); 
    PolyUpPred = cat(3, PolyUpPred, [xPolyUpPred, yPolyUpPred]); 
end 
%% Snake Energy %% 
xPre = trajectSet(:,1,:); 
yPre = trajectSet(:,2,:); 
% xPre = arcxPre; 
% yPre = arcyPre; 
xLowAll=cat(1, xPre(100:100:end,:,:), PolyLowPred(100:100:end, 2, :)); 
yLowAll=cat(1, yPre(100:100:end,:,:), PolyLowPred(100:100:end, 1, :)); 
xUpAll=cat(1, xPre(100:100:end,:,:), PolyUpPred(100:100:end, 2, :)); 
yUpAll=cat(1, yPre(100:100:end,:,:), PolyUpPred(100:100:end, 1, :)); 
% Econt=0; 

% EcontLow = sum(sqrt((xLowAll(1:end-1,:,:)-xLowAll(2:end,:,:)).^2+(yLowAll(1:end-1,:,:)-yLowAll(2:end,:,:)).^2), 1); 
% EcontUp = sum(sqrt((xUpAll(1:end-1,:,:)-xUpAll(2:end,:,:)).^2+(yUpAll(1:end-1,:,:)-yUpAll(2:end,:,:)).^2), 1); 
EcontLow = 0; 
EcontUp = 0; 
EcurvLow = sum(sqrt((xLowAll(3:end,:,:)-2*xLowAll(2:end-1,:,:)+xLowAll(1:end-2,:,:)).^2+(yLowAll(3:end,:,:)-2*yLowAll(2:end-1,:,:)+yLowAll(1:end-2,:,:)).^2), 1); 
EcurvUp = sum(sqrt((xUpAll(3:end,:,:)-2*xUpAll(2:end-1,:,:)+xUpAll(1:end-2,:,:)).^2+(yUpAll(3:end,:,:)-2*yUpAll(2:end-1,:,:)+yUpAll(1:end-2,:,:)).^2), 1); 

EintLow = EcontLow + EcurvLow; 
EintUp = EcontUp + EcurvUp; 

xAll = [PolyLowPred(:,2,:), PolyUpPred(:,2,:)]; 
yAll = [PolyLowPred(:,1,:), PolyUpPred(:,1,:)]; 
EAll = [EintLow, EintUp]; 
[EminVal, EminInd] = min(EAll, [], 2); 

xPred = []; 
yPred = []; 
for ind_i = 1: n 
    xPredTemp = xAll(:, EminInd(ind_i), ind_i); 
    yPredTemp = yAll(:, EminInd(ind_i), ind_i); 
    xPred = cat(3, xPred, xPredTemp); 
    yPred = cat(3, yPred, yPredTemp); 
end 