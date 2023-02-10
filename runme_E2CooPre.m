% function [trajectAvoiding, trajectSet] = runmeE2Coop (velAmpIntru0, posIniIndiv)
%% test %%
clear;
clc;
plotFlag=1;
addpath('GVF');
addpath('snake2D');
%% Define intruders %%
velAmpIntru1=0:10;
n1=2:10;
gamma1Set=.1:.1:.9;
varLen=length(velAmpIntru1);
timecost = zeros(1, 1); 
step = 0; 
for ind_l = 1 % Number of run 
    for ind_ll = 1 % Index of obstacle velocity 
        eval(['file_swarm_x = fopen(''swarm_x_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_swarm_y = fopen(''swarm_y_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_obs_x = fopen(''obs_x_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_obs_y = fopen(''obs_y_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_swarm_kappa = fopen(''swarm_kappa_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_swarm_omega = fopen(''swarm_omega_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_obs_omega = fopen(''obs_omega_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        eval(['file_swarm_start_end = fopen(''swarm_start_end_',num2str(ind_l),'_',num2str(ind_ll),'.txt'',''w'');']);
        
        eval(['disp(''repeat=',num2str(ind_l),''');' ]);
        eval(['disp(''variable=',num2str(ind_ll),''');' ]);
        cd('D:\Scripts\E2CoPre');

        gamma1=.5;
        gamma2 = 1-gamma1;
        velAngIntru=[-135];
        posIniIntru=[180; 180]; % Position of the intruders

        [~, intruNum]=size(posIniIntru);
        velAmpSwarm = 10; % Velocity of swarm (center)
        velAmpIntru0=velAmpIntru1(ind_ll); % Velocity of obstacles
%         velAmpIntru0 = 10;
        velAmpIntru = max(velAmpIntru0, velAmpSwarm);
        velAmpIntru = velAmpIntru * ones(intruNum, 1); 
        velVecIntru=[velAmpIntru.*cosd(velAngIntru), velAmpIntru.*sind(velAngIntru)];
        %% Define destination %%
        posIniDest=[250; 250]; % Position of the destination
        fprintf(file_swarm_start_end,'%s\r\n', 'end: x, y');
        %% Define swarm %%
        n=n1(4); % Swarm size N 
        theta=linspace(-pi,pi,n+1);
        x0=95;
        y0=95;
        r=20; % Swarm radius \tau 

        xi=r*cos(theta)+x0;
        yi=r*sin(theta)+y0;
        posIniIndivAlt=zeros(1, n);
        posIniIndiv=[xi(1:end-1); yi(1:end-1)];

        fprintf(file_swarm_start_end,'%s\r\n', 'start: 1st row x, 2nd row y');
        fprintf(file_swarm_start_end,'%12.8f %12.8f\r\n',xi(1:end-1));
        fprintf(file_swarm_start_end,'%12.8f %12.8f\r\n',yi(1:end-1));
        
        [~, swarmSize]=size(posIniIndiv);
        posIniSwarm=sum(posIniIndiv,2)./swarmSize;
        distSwarmIndiv = sqrt(sum((posIniIndiv-posIniSwarm).^2,1)); % Geometrical center of the swarm
        vecSwarm2Dest=posIniDest-posIniSwarm; % Offset vector towards the destination
        vecSwarm2DestUni=vecSwarm2Dest./norm(vecSwarm2Dest);
        posIniSwarm=posIniSwarm+vecSwarm2DestUni.*2*max(max(distSwarmIndiv), r); % Swarm center with offset
        %% Potential Fields %%
        fieldSize = 300; % Define size of the field;
        posIniSwarmRel = round(posIniSwarm);
        posIniIntruRel = round(posIniIntru);
        groundIni = zeros(fieldSize);
        [groundRow, groundCol]=size(groundIni);
        fieldMotionIni = groundIni;
        %% Parameter Settings %%
        simulTime = 50;
        reslS=30; % Granularity of deltaS (length of arc)
        deltaSMax = max(velAmpIntru0, 5)*ones(1,1,swarmSize); % Maximum deltaS
        deltaSMin = max(velAmpIntru0, 5)*ones(1,1,swarmSize); % Manimum deltaS
        linSpcInd = linspace(0,1,reslS);
        deltaS = deltaSMin+linSpcInd'.*(deltaSMax-deltaSMin); % Dimension of arc length in PSO
                
        omega=0:1:359;
        omega=omega(:,:,ones(swarmSize, 1)); % Dimension of omega in PSO
        
        linSpcInd = linspace(0,1,350);
        kappa=ones(350, 1)*[-pi./(1*deltaSMax(:)')]+linSpcInd'.*(ones(350, 1)*(pi./(1*deltaSMax(:)')+pi./(1*deltaSMax(:)')));
        kappa=reshape(kappa, 1, 350, swarmSize); % Dimension of curvature in PSO
        
        vetVelMax=30*ones(1,1,swarmSize);
        linSpcInd = linspace(0,1,2*5+1);
        tau=-vetVelMax+linSpcInd'.*(vetVelMax+vetVelMax);
        
        Maxiter=100; % Iteratons of optimization in PSO
        N=100; % number of interpolation points of each arc
        %% Pre-trajectory Settings %%
        smoothField = groundIni;
        xPre = reshape(posIniIndiv(1,:),1,1,swarmSize);
        yPre = reshape(posIniIndiv(2,:),1,1,swarmSize);
        AltPre = reshape(posIniIndivAlt,1,1,swarmSize);
        kPre = zeros(1, 1, swarmSize);
        tPre = zeros(1, 1, swarmSize);
        omegaPre = zeros(1, 1, swarmSize);
        contourPlotCell=cell(simulTime, 1);
        [arcxPre, arcyPre, arcAltPre, x1Pre, y1Pre]=arcParameterization3D (xPre, yPre, AltPre, kPre, tPre, omegaPre, 5*5, N);
        arcxPre=arcxPre-5*5;
        arczPre=interp2(smoothField, arcxPre, arcyPre); % Previous trajectories before avoidance
        x0=xPre;
        y0=yPre;
        Alt0=AltPre;
        z0=interp2(smoothField, x0, y0); % Starting points for trajectory planning
        
        arcPreOrigin=cat(2,arcxPre,arcyPre,arcAltPre);
        contourPlot=[0 0 0];
        trajectSet=[x0, y0, AltPre, zeros(size(x0)), z0];
        x0_obs = reshape(posIniIntru(1,:),1,1,intruNum);
        y0_obs = reshape(posIniIntru(2,:),1,1,intruNum);
        trajectSetObs = [x0_obs, y0_obs, zeros(size(x0_obs)), zeros(size(x0_obs)), zeros(size(x0_obs))]; 
        posSwarmSet = posIniSwarmRel;
        smoothFieldSet = smoothField;
        %% Plotting %%
        if plotFlag
            writerObj1=VideoWriter('AnimiFinal.avi');
            writerObj1.FrameRate = 1;
            open(writerObj1);
            
            fig2=figure(ind_ll);
            set(gcf,'color','w');
            g3=mesh(gca, smoothField);
            axe2=get(fig2,'CurrentAxes');
            view (0,90);
            hold (axe2);
            g1=plot3(axe2, x0(:), y0(:), AltPre(:), 'r*');
            for ind_h=1: swarmSize
                eval(['h',num2str(ind_h),'=plot3(axe2, trajectSet(:,1,ind_h), trajectSet(:,2,ind_h), ones(size(trajectSet(:,3,ind_h)))*max(smoothField(:)), ''r'',''linewidth'',2); ']);
            end
            g4=plot3(contourPlot(:,1), contourPlot(:,2), contourPlot(:,3), 'w.','markersize',3);
            for ind_k=1: swarmSize
                eval(['k',num2str(ind_k),'=plot3(axe2, arcPreOrigin(:,1,ind_k), arcPreOrigin(:,1,ind_k), ones(size(arcPreOrigin(:,1,ind_k)))*max(smoothField(:)), ''r:'',''linewidth'',2); ']);
            end
            g2=plot3(posIniIntru(1,:), posIniIntru(2,:), interp2(smoothField, posIniIntru(1,:), posIniIntru(2,:)), 'bs', 'MarkerFaceColor', 'b');
            % g5=plot3(posIniDest(1,:), posIniDest(2,:), max(smoothField(:)), 'rs');
            g6=plot3(posIniSwarm(1), posIniSwarm(2), interp2(smoothField, posIniSwarm(1), posIniSwarm(2)), 'ro', 'MarkerFaceColor', 'r');
            for ind_p=1: swarmSize
                eval(['p',num2str(ind_p),'=plot3(axe2, arcPreOrigin(:,1,ind_p), arcPreOrigin(:,1,ind_p), ones(size(arcPreOrigin(:,1,ind_p)))*max(smoothField(:)), ''b:'',''linewidth'',2); ']);
            end
        end
        %% Optimization %%
        posSwarmRel=posIniSwarmRel; % Position of swarm
        posIntruRel=posIniIntruRel; % Position of intruder
        vecSwarm2Dest=posIniDest-posSwarmRel; % Vector from swarm center to destination
        vecSwarm2DestUni=vecSwarm2Dest./norm(vecSwarm2Dest);
        distSwarm2Dest=norm(vecSwarm2Dest);
        
        distIntruIndiv = sqrt(sum((posIniIndiv(:,:,ones(intruNum, 1))-reshape(posIniIntru, 2, 1, intruNum)).^2,1));
        distIntruIndiv = min(distIntruIndiv, [], 3); % Distance from swarm member to intruder
        
        vecIntruIndiv = reshape(posIniIntru, 2, 1, intruNum)-posIniIndiv(:,:,ones(intruNum, 1));
        vecIntruIndivUni = vecIntruIndiv./ vecnorm(vecIntruIndiv); % Vector from swarm member to intruder
        
        % vecDestIndiv = posIniDest - posIniIndiv;
        vecDestIndiv = [350; 350] - posIniIndiv;
        vecDestIndivUni = vecDestIndiv./ vecnorm(vecDestIndiv); % Vector from swarm member to destination
        
        stateIndicator = dot(vecIntruIndivUni, vecDestIndivUni(:,:,ones(1, intruNum)));
        % stateIndicator is the dot product of vector UAV-intruder and vector
        % UAV-destination. It indicates when the UAVs finish avoiding. If the dot
        % product is positive, posibilities of collision exist. If the dot product
        % is nagative, the UAVs have passed the intruders, avoidance is deemed accomplished.
        stateIndicator (stateIndicator >= 0) = 0;
        stateIndicator (stateIndicator < 0) = -1;
        stateIndicator = sum(stateIndicator, 3);
        stateIndicator (stateIndicator == 0) = 1;
        
        vecIniSwarmIndiv =  posIniIndiv - posIniSwarmRel; % Records the initial relative
        % positions of swarm members. It's for formation restoration after
        % avoidance.
        
        posIndiv=posIniIndiv;
        posIndivAlt=posIniIndivAlt;
        swarmInd=1:swarmSize;
        guardDistSwarm=5;
        guardDist=20; % d_{safe} 
        
        startDist=30+guardDist+velAmpIntru0; % Distance from the intruder, where UAVs start to avoid.
        startDist=startDist(1, ones(1, swarmSize));
        deltaSTemp=deltaS(1,:,:);
        deltaSTemp=deltaSTemp(:);
        % First sort by velocity;
        [uniqS, uniqSInd]=unique(deltaSTemp);
        if length(uniqSInd) > 1
            sortMethod='ascend';
            [uniqSSort1, uniqSSortInd1]=sort(uniqS, sortMethod); % Slow ones avoid first.
            uniqSSortInd1=uniqSInd(uniqSSortInd1);
            uniqSdiff=abs(diff(uniqSSort1));
            if strcmp(sortMethod, 'descend')
                uniqSSortInd1=uniqSSortInd1(1: end-1);
            elseif strcmp(sortMethod, 'ascend')
                uniqSSortInd1=uniqSSortInd1(2: end);
            end
            uniqSSortInd1=uniqSSortInd1';
            % Then sort by dist2Obs;
            ind2sort=setdiff(swarmInd, uniqSSortInd1);
            [uniqSSort2, uniqSSortInd2]=sort(distIntruIndiv(ind2sort), 'ascend');
            uniqSSortInd2=ind2sort(uniqSSortInd2);
            if strcmp(sortMethod, 'ascend')
                uniqSSortInd=[uniqSSortInd2, uniqSSortInd1];
            elseif strcmp(sortMethod, 'descend')
                uniqSSortInd=[uniqSSortInd1, uniqSSortInd2];
            end
        else
            % Then sort by dist2Obs;
            [uniqSSort2, uniqSSortInd2]=sort(distIntruIndiv, 'ascend');
            uniqSSortInd=uniqSSortInd2;
        end
        
        uniqSLen=length(uniqSSortInd);
        guardDistv2v = 5; % d_{v2v} 
        
        buff=40;
        detectDist=startDist+buff; % Distance from the intruder, where UAVs detect the intruder.
        swarmIndStop = [];
        swarmIndAvoiding = [];
        swarmInd2Avoid = [];
        swarmCostSet = [];
        internCostSet = [];
        trajectAvoiding = [];
        velCenter = 0;
        velIntru = velAmpIntru0;
        fieldMotion = fieldMotionIni;
        vecIndivFormNorm = ones(1, swarmSize)*5+1;
        dist2intru = [];
        xPred=trajectSet(end,1,:);
        yPred=trajectSet(end,2,:);
        PSet = cell(swarmSize, 1);
%         for ind_t=1:2 %1: 3 %simulTime
        ind_t=0;
        while length(find(vecIndivFormNorm<5)) ~= swarmSize && ind_t < 40 
            ind_t=ind_t+1;
            %% Define Swarm Dir %%
            vecSwarm2Dest=posIniDest-posSwarmRel;
            distSwarm2Dest=norm(vecSwarm2Dest);
            vecSwarm2DestUni=vecSwarm2Dest./distSwarm2Dest; % Vector from swarm center to the destination
            
            if isempty(stateIndicator (stateIndicator>0))
                velCenter = 0; % If avoidance accomplished, stop moving swarm center and start formation restoration.
                %     elseif ~isempty(swarmInd2Avoid)
                %         velCenter = velCenter-velAmpIntru;
            else
                velCenter = mean(deltaSMax);
                %         velCenter = abs(velAmpIntru- velCenter);
            end
            
            posSwarmRel=posSwarmRel+vecSwarm2DestUni.*abs(velCenter); % Update swarm center movement
            %     posSwarmRel=posSwarmRel+vecSwarm2DestUni.*(velCenter);
            
            posSwarmForm = posSwarmRel;
            posIndivForm = posSwarmForm + vecIniSwarmIndiv;
            vecIndivForm = posIndivForm - posIndiv;
            vecIndivFormNorm = vecnorm(vecIndivForm);
            vecIndivFormUni = vecIndivForm./vecIndivFormNorm; % Formation restoration
            
            vecIntru=[cosd(velAngIntru); sind(velAngIntru)];
            posIntruRel = posIntruRel + vecIntru.*abs(velIntru); % Update intruder movement
            %% Construct potential field and contours %%
            comm_noise = 0; % Interference in communication channel 
            if comm_noise  
                mu = 0;
                delta = 1; 
                posIntruRel = normrnd(mu + posIntruRel, delta, size(posIntruRel, 1), size(posIntruRel, 2)); 
                posSwarmRel = normrnd(mu + posSwarmRel, delta, size(posSwarmRel, 1), size(posSwarmRel, 2));  
            end 

            obsFlag=1;
            [posCurrIntru, smoothField1] = testUAVMotion ...
                (fieldMotion, posIntruRel, velAmpIntru*20, guardDist, obsFlag); % Construct potential field
                        
            swarmInd2Avoid = swarmInd(distIntruIndiv <= detectDist);
            if isempty(swarmInd2Avoid)
                swarmInd2Avoid = zeros(1,0);
            elseif size(swarmInd2Avoid, 1) > size(swarmInd2Avoid, 2)
                swarmInd2Avoid = swarmInd2Avoid';
            end
            
            ind_tMark = 1;
            if ~isempty(swarmInd2Avoid)
                % Accumulation field without interpolation
                obsFlag=1;
                [posCurrIntru, smoothField0] = testUAVMotion (fieldMotion, ...
                    [posSwarmSet(:,ind_tMark:end)], [velAmpIntru(1,ones(ind_t-ind_tMark+1,1))*10], guardDistSwarm, obsFlag);
                smoothField = smoothField1.*(ind_t-ind_tMark+1) + smoothField0 * 1;
                % Memory of potential field, from when the obstacle is deteced.
            else
                obsFlag=1;
                [posCurrIntru, smoothField0] = testUAVMotion(fieldMotion, ...
                    [posSwarmRel], [velAmpIntru*10], guardDistSwarm, obsFlag);
                smoothField = smoothField1 + smoothField0 * 1;
                ind_tMark = ind_t; % Mark the time obstacle is detected.
            end
            smoothFieldSet=cat(3, smoothFieldSet, smoothField);
            %% PSO %%
            if plotFlag
                tic;
            end
            
            z0=interp2(smoothField, x0, y0); % Update z0 immediately after updating smoothField.
                        
            intruMember=[];
            zForbid = zeros(intruNum, 1);
            for ind_i=1: intruNum % No contours drawn within guardDist.
                zForbidMat=smoothField([max(round(posIntruRel(2,ind_i)-guardDist), 1)...
                    :round(posIntruRel(2,ind_i)+guardDist)], ...
                    [max(round(posIntruRel(1,ind_i)-guardDist), 1)...
                    :round(posIntruRel(1,ind_i)+guardDist)]);
                if ~isempty(zForbidMat)
                    zForbid(ind_i)=min(zForbidMat(:));
                end
                distIntruIndiv3D=reshape(distIntruIndiv, 1, 1, swarmSize);
            end
            
            binField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            EextField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            paraField=0;
            pxField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            pyField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            uField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            vField=zeros(size(smoothField, 1), size(smoothField, 2), swarmSize);
            for ind_s=1: swarmSize
                [binField(:,:,ind_s), EextField(:,:,ind_s), pxField(:,:,ind_s), ...
                    pyField(:,:,ind_s), uField(:,:,ind_s), vField(:,:,ind_s)] = ...
                    ExternalEnergyField (smoothField, z0(ind_s));
            end
            % Generate binary field and external energy field.
            %% PSO Here %%
            arcx=zeros(N, 1, swarmSize);
            arcy=zeros(N, 1, swarmSize);
            arcz=zeros(N, 1, swarmSize);
            arcAlt=zeros(N, 1, swarmSize);
            swarmIndAvoiding = swarmInd(distIntruIndiv <= startDist & stateIndicator >= 0);
            % Find UAVs avoiding obstacles
            if isempty(swarmIndAvoiding)
                swarmIndAvoiding=zeros(1,0);
            elseif size(swarmIndAvoiding, 1) > size(swarmIndAvoiding, 2)
                swarmIndAvoiding=swarmIndAvoiding';
            end
            
            swarmIndNonAvoiding = swarmInd(distIntruIndiv > startDist & stateIndicator >= 0);
            % Find UAVs not start avoidance yet
            swarmIndStop = swarmInd(stateIndicator < 0);
            
            swarmIndStop=union(swarmIndStop, swarmIndNonAvoiding);
            % Find UAVs finished avoidance
            
            if isempty(stateIndicator (stateIndicator>0))
                swarmIndStop = swarmInd;
                swarmIndAvoiding = [];
                swarmIndNonAvoiding = [];
            end
            
            if isempty(swarmIndAvoiding)
                % UAVs not start avoidance yet just go straight to destination.
                kPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                tPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                omegaPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                omegaPre2(1,1,:)=atand(vecSwarm2DestUni(2)./vecSwarm2DestUni(1));
                [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                    (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), Alt0(:,:,swarmIndNonAvoiding), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                
                arcx(:,1,swarmIndNonAvoiding)=arcx2;
                arcy(:,1,swarmIndNonAvoiding)=arcy2;
                arcAlt(:,1,swarmIndNonAvoiding)=arcAlt2;
                omegaPre (:,:,swarmIndNonAvoiding) = omegaPre2;
                
                % UAVs finihsed avoidance proceed to formation restoration.
                deltaSMax(:,:,swarmIndStop)=5;
                deltaSMin(:,:,swarmIndStop)=5;
                kPre2 = zeros(1, 1, length(swarmIndStop));
                tPre2 = zeros(1, 1, length(swarmIndStop));
                omegaPre2 = zeros(1, 1, length(swarmIndStop));
                angTemp=atand(vecIndivFormUni(2,swarmIndStop)./vecIndivFormUni(1,swarmIndStop));
                angTemp(vecIndivFormUni(1,swarmIndStop)<0)=180+angTemp(vecIndivFormUni(1,swarmIndStop)<0);
                omegaPre2(1,1,:)=angTemp;

                [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                    (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), Alt0(:,:,swarmIndStop), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                
                arcx(:,1,swarmIndStop)=arcx2;
                arcy(:,1,swarmIndStop)=arcy2;
                arcAlt(:,1,swarmIndStop)=arcAlt2;
                omegaPre (:,:,swarmIndStop) = omegaPre2;
                
                x0temp=x0(:,:,vecnorm(vecIndivForm)<=5);
                y0temp=y0(:,:,vecnorm(vecIndivForm)<=5);
                Alt0temp=Alt0(:,:,vecnorm(vecIndivForm)<=5);
                arcx(:,1, vecnorm(vecIndivForm)<=5)=x0temp(ones(N,1),:,:);
                arcy(:,1, vecnorm(vecIndivForm)<=5)=y0temp(ones(N,1),:,:);
                arcz(:,1, vecnorm(vecIndivForm)<=5)=ones(size(arcy(:,1, vecnorm(vecIndivForm)<=5)))*100;
                arcAlt(:,1, vecnorm(vecIndivForm)<=5)=Alt0temp(ones(N,1),:,:); 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                kappa2write = zeros(1, n);
                omega2write = zeros(1, n);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                if ~isempty(intruMember)
                    arcx=arcxPre;
                    arcy=arcyPre;
                    arcAlt=arcAltPre;
                else
                    %% Prediction %%
                    iteration = 100;
                    
                    theta1 = linspace(1*pi/2, 7*pi/4, 10);
                    theta2 = linspace(-1*pi/2, 3*pi/4, 10);
                    d = 2*r;
                    xStart1 = d*cos(theta1)+sum(x0(:))/n;
                    xStart2 = d*cos(theta2)+posIntruRel(1);
                    yStart1 = d*sin(theta1)+sum(y0(:))/n;
                    yStart2 = d*sin(theta2)+posIntruRel(2);
                    xStart = [xStart1, xStart2, xStart1(1)];
                    yStart = [yStart1, yStart2, yStart1(1)];
                    
                    trajectStart = [yStart', xStart'];
                    trajectStart = InterpolateContourPoints2D(trajectStart, N);
                    trajectStart = trajectStart(:, :, ones(1,n));
                    
                    for ind_y = 1: swarmSize
                        trajectStartTemp = PSet{ind_y};
                        if ~isempty(trajectStartTemp)
                            trajectStart(:,:,ind_y) = trajectStartTemp;
                        end
                    end
                    
                    stepPred = 20;
                    xPred = zeros(stepPred*N, 1, swarmSize);
                    yPred = zeros(stepPred*N, 1, swarmSize);
                    PSet = cell(swarmSize, 1);
                    [yPred(:,:,swarmIndAvoiding), xPred(:,:,swarmIndAvoiding), P] = covPredict (iteration, ...
                        binField(:,:,swarmIndAvoiding), trajectStart(:,:,swarmIndAvoiding), ...
                        x0(:,:,swarmIndAvoiding), y0(:,:,swarmIndAvoiding), posIniDest, ...
                        arcxPre(:,:,swarmIndAvoiding), arcyPre(:,:,swarmIndAvoiding), arczPre(:,:,swarmIndAvoiding), ...
                        trajectSet(:,:,swarmIndAvoiding));

                    step_len = 500;
                    xPredStep = xPred(N+1:N+1+step_len,:,:);
                    yPredStep = yPred(N+1:N+1+step_len,:,:); 
                    init_curv = mean(sqrt((xPredStep(3:end,:,:)-2*xPredStep(2:end-1,:,:)+xPredStep(1:end-2,:,:)).^2 ...
                        +(yPredStep(3:end,:,:)-2*yPredStep(2:end-1,:,:)+yPredStep(1:end-2,:,:)).^2), 1);
                    omega_len = step_len; 
                    init_omega = atand((yPredStep(1:omega_len,:,:) - y0)./(xPredStep(1:omega_len,:,:) - x0)); 
                    init_omega(xPredStep(1:omega_len,:,:) < x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) > y0(ones(1, omega_len), :, :)) = ...
                        init_omega(xPredStep(1:omega_len,:,:) < x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) > y0(ones(1, omega_len), :, :)) + 180; 
                    init_omega(xPredStep(1:omega_len,:,:) < x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) < y0(ones(1, omega_len), :, :)) = ...
                        init_omega(xPredStep(1:omega_len,:,:) < x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) < y0(ones(1, omega_len), :, :)) + 180;
                    init_omega(xPredStep(1:omega_len,:,:) > x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) < y0(ones(1, omega_len), :, :)) = ...
                        init_omega(xPredStep(1:omega_len,:,:) > x0(ones(1, omega_len), :, :) & yPredStep(1:omega_len,:,:) < y0(ones(1, omega_len), :, :)) + 360;
                    init_omega = mean(init_omega); 

                    kappa_diff = abs(kappa - init_curv); 
                    omega_diff = abs(omega - init_omega); 
                    [~, init_kappa_ind] = min(kappa_diff, [], 2); 
                    [~, init_omega_ind] = min(omega_diff, [], 2);
                    
                    for ind_u = 1: length(swarmIndAvoiding)
                        PSet{swarmIndAvoiding(ind_u)} = P(:,:,ind_u);
                    end
                    %% De-conflict %%
                    AltAdjust = Alt0*0;
                    if length(swarmIndAvoiding) > 1
                        dist2Obs = min(sqrt((xPred-posIntruRel(1)).^2 + (yPred-posIntruRel(2)).^2), [], 1);
                        combo = nchoosek(swarmIndAvoiding, 2);
                        if ismember(ind_l, [501: 520])
                            distV2V = sqrt((x0(end,:,combo(:,1))-x0(end,:,combo(:,2))).^2 ...
                                + (y0(end,:,combo(:,1))-y0(end,:,combo(:,2))).^2 ...
                                + (Alt0(end,:,combo(:,1))-Alt0(end,:,combo(:,2))).^2);
                        else
                            distV2V = sqrt((xPred(end,:,combo(:,1))-xPred(end,:,combo(:,2))).^2 ...
                                + (yPred(end,:,combo(:,1))-yPred(end,:,combo(:,2))).^2 ...
                                + (Alt0(end,:,combo(:,1))-Alt0(end,:,combo(:,2))).^2);
                        end 
                        
%                         eval(['disp(''t = ',num2str(ind_t),', d_{v2v}: '')']);
%                         squeeze(distV2V)
                        
                        thrV2V = guardDistv2v;
                        dangerInd = combo(distV2V < thrV2V, :);
                        %% Alt PSO %%
                        if ~isempty(dangerInd)
                            t = tau(:,:,swarmIndAvoiding);
                            xEnd = xPred(end,:,swarmIndAvoiding);
                            yEnd = yPred(end,:,swarmIndAvoiding);
                            zEnd = z0(:,:,swarmIndAvoiding);
                            AltEnd = Alt0(:,:,swarmIndAvoiding);
                            if ismember(ind_l, [271: 290])
                                mu = 0;
                                delta = 1;
                                AltEnd = normrnd(mu + AltEnd, delta, size(AltEnd, 1), size(AltEnd, 2), size(AltEnd, 3));
                            end
                            Maxiter = 100;
                            [tBestInterp] = SpeciesPSOAlt (xEnd, yEnd, zEnd, ...
                                AltEnd, t, thrV2V, Maxiter);
                        else
                            tBestInterp = 0;
                        end
                    else
                        tBestInterp = 0;
                    end
                    AltAdjust(:,:,swarmIndAvoiding) = tBestInterp;
                    altTran=Alt0(ones(N,1),:,:);
                    for ind_j = swarmIndAvoiding
                        altStart = Alt0(ind_j);
                        altEnd = AltAdjust(ind_j)+Alt0(ind_j);
                        altTran(:,:,ind_j) = linspace(altStart, altEnd, N);
                    end
                    %% Horizontal PSO %%
                    if ismember(ind_l, [501: 520])
                        rnd_init = 1;
                    else
                        rnd_init = 0;
                    end
                    Maxiter = 100;
                    [arcx1, arcy1, arcz1, kNew, omegaNew, deltaSNew, BestCosts, swarmVecPre, swarmVecCurr, swarmCost] = ...
                        SpeciesPSO (gamma1, gamma2, swarmIndAvoiding, posIniIndiv(:,  swarmIndAvoiding), posIniSwarmRel, ...
                        posSwarmRel, x0(:,:, swarmIndAvoiding), y0(:,:, swarmIndAvoiding), z0(:,:, swarmIndAvoiding), ...
                        zForbid, kappa(:,:, swarmIndAvoiding), omega(:,:, swarmIndAvoiding), deltaS(:,:, swarmIndAvoiding), ...
                        arcxPre(:,:, swarmIndAvoiding), arcyPre(:,:, swarmIndAvoiding), arczPre(:,:, swarmIndAvoiding), ...
                        smoothField, EextField, Maxiter, N, init_kappa_ind, init_omega_ind, rnd_init);
                    arcAlt1=Alt0(ones(N,1),:,swarmIndAvoiding);
                    
                    arcx(:,1, swarmIndAvoiding)=arcx1;
                    arcy(:,1, swarmIndAvoiding)=arcy1;
                    arcz(:,1, swarmIndAvoiding)=arcz1;
                    arcAlt(:,1, swarmIndAvoiding)=altTran(:,:,swarmIndAvoiding);
                    
                    deltaSMax(:,:,swarmIndStop)=5;
                    deltaSMin(:,:,swarmIndStop)=5;
                    kPre2 = zeros(1, 1, length(swarmIndStop));
                    tPre2 = zeros(1, 1, length(swarmIndStop));
                    omegaPre2 = zeros(1, 1, length(swarmIndStop));
                    angTemp=atand(vecIndivFormUni(2,swarmIndStop)./vecIndivFormUni(1,swarmIndStop));
                    angTemp(vecIndivFormUni(1,swarmIndStop)<0)=180+angTemp(vecIndivFormUni(1,swarmIndStop)<0);
                    omegaPre2(1,1,:)=angTemp;

                    [arcx2, arcy2, arcAlt2, x1Pre2, y1Pre2]=arcParameterization3D ...
                        (x0(:,:,swarmIndStop), y0(:,:,swarmIndStop), Alt0(:,:,swarmIndStop), kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndStop), N);
                    
                    arcx(:,1,swarmIndStop)=arcx2;
                    arcy(:,1,swarmIndStop)=arcy2;
                    arcz(:,1,swarmIndStop)=ones(size(arcy2))*100;
                    arcAlt(:,1,swarmIndStop)=arcAlt2;
                    omegaPre (:,:,swarmIndStop) = omegaPre2;
                    
                    kPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    tPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    omegaPre2 = zeros(1, 1, length(swarmIndNonAvoiding));
                    omegaPre2(1,1,:)=atand(vecSwarm2DestUni(2)./vecSwarm2DestUni(1));

                    [arcx2, arcy2, arcAlt2, x1Pre, y1Pre]=arcParameterization3D ...
                        (x0(:,:,swarmIndNonAvoiding), y0(:,:,swarmIndNonAvoiding), Alt0(:,:,swarmIndNonAvoiding), ...
                        kPre2, tPre2, omegaPre2, deltaSMax(:,:,swarmIndNonAvoiding), N);
                    
                    arcx(:,1,swarmIndNonAvoiding)=arcx2;
                    arcy(:,1,swarmIndNonAvoiding)=arcy2;
                    arcz(:,1,swarmIndNonAvoiding)=ones(size(arcy2))*500;
                    arcAlt(:,1,swarmIndNonAvoiding)=arcAlt2;
                    omegaPre (:,:,swarmIndNonAvoiding) = omegaPre2;
                    
                    x0temp=x0(:,:,vecnorm(vecIndivForm)<=5);
                    y0temp=y0(:,:,vecnorm(vecIndivForm)<=5);
                    Alt0temp=Alt0(:,:,vecnorm(vecIndivForm)<=5);
                    arcx(:,1, vecnorm(vecIndivForm)<=5)=x0temp(ones(N,1),:,:);
                    arcy(:,1, vecnorm(vecIndivForm)<=5)=y0temp(ones(N,1),:,:);
                    arcz(:,1, vecnorm(vecIndivForm)<=5)=ones(size(arcy(:,1, vecnorm(vecIndivForm)<=5)))*100;
                    arcAlt(:,1, vecnorm(vecIndivForm)<=5)=Alt0temp(ones(N,1),:,:);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                kappa2write(swarmIndAvoiding) = kNew;
                omega2write(swarmIndAvoiding) = omegaNew;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(file_swarm_x,'%12.8f %12.8f\r\n',x0(:)');
            fprintf(file_swarm_y,'%12.8f %12.8f\r\n',y0(:)');
            fprintf(file_obs_x,'%12.8f\r\n',posIntruRel(1));
            fprintf(file_obs_y,'%12.8f\r\n',posIntruRel(2));
            fprintf(file_swarm_kappa,'%12.8f %12.8f\r\n',kappa2write);
            fprintf(file_swarm_omega,'%12.8f %12.8f\r\n',omega2write);
            fprintf(file_obs_omega,'%12.8f\r\n',velAngIntru);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            timeCost=toc;
            step = step + 1; 
            timecost(step) = timeCost; 
            trajectTemp=cat(2, arcx, arcy, arcAlt, timeCost(ones(size(arcz))), arcz);
            trajectSet=cat(1,trajectSet,trajectTemp);

            xObsPre = trajectSetObs(end, 1); 
            yObsPre = trajectSetObs(end, 2); 
            xObs = posIntruRel(1, :); 
            yObs = posIntruRel(2, :); 
            arcxObs = zeros(N, 1, intruNum); 
            arcyObs = zeros(N, 1, intruNum); 
            for obs_ind = 1: intruNum
                arcxObs(:, :, obs_ind) = linspace(xObsPre, xObs(obs_ind), N)'; 
                arcyObs(:, :, obs_ind) = linspace(yObsPre, yObs(obs_ind), N)'; 
            end 
            trajectTempObs = cat(2, arcxObs, arcyObs, zeros(size(arcxObs)), zeros(size(arcxObs)), zeros(size(arcxObs)));
            trajectSetObs = cat(1, trajectSetObs, trajectTempObs);

            if size(trajectSetObs, 1) ~= size(trajectSet, 1)
                break; 
            end 
            trajectAvoiding=cat(1, trajectAvoiding, [trajectTemp, ind_t(ones(1,N),1,ones(1,swarmSize))]);
            
            arcxPre=trajectSet(end-N+1:end, 1, :);
            arcyPre=trajectSet(end-N+1:end, 2, :);
            arcAltPre=trajectSet(end-N+1:end, 3, :);
            x0=arcxPre(end,:,:);
            y0=arcyPre(end,:,:);
            Alt0=arcAltPre(end,:,:);
            
            posIndiv = [x0(:)'; y0(:)'];
            posIndivAlt = Alt0(:)';
            
            dist2intruTemp = zeros(intruNum, n); 
            for obs_ind = 1: intruNum
                dist2intruTemp(obs_ind, :) = vecnorm(posIntruRel(:, obs_ind)-posIndiv, 2);
            end 
            dist2intruTemp = min(dist2intruTemp, [], 1); 
            dist2intru=cat(1, dist2intru, dist2intruTemp);
            
            distIntruIndiv = sqrt(sum((posIndiv(:,:,ones(intruNum, 1))-reshape(posIntruRel, 2, 1, intruNum)).^2,1));
            distIntruIndiv = min(distIntruIndiv, [], 3);
            distIndivDest = sqrt(sum((posIndiv-posIniDest).^2,1));
            
            vecIntruIndiv = reshape(posIntruRel, 2, 1, intruNum) - posIndiv(:,:,ones(intruNum, 1));
            vecIntruIndivUni = vecIntruIndiv./ vecnorm(vecIntruIndiv);
            vecDestIndiv = [350; 350] - posIniIndiv;
            vecDestIndivUni = vecDestIndiv./ vecnorm(vecDestIndiv);
            
            stateIndicator=dot(vecIntruIndivUni, vecDestIndivUni(:,:,ones(1, intruNum)));
            stateIndicator (stateIndicator >= 0) = 0;
            stateIndicator (stateIndicator < 0) = -1;
            stateIndicator = sum(stateIndicator, 3);
            stateIndicator (stateIndicator == 0) = 1;
            
            if isempty(stateIndicator (stateIndicator>0))
                stateIndicator = -ones(1, swarmSize);
            end
            if plotFlag
                %% Plotting %%
                pause(.01);
                
                [contourPlot, sectorNum] = ContourPlot (smoothField, z0);
                contourPlotCell{ind_t}=contourPlot;
                if ~isempty(contourPlot) 
                    g4.XData=contourPlot(:,1);
                    g4.YData=contourPlot(:,2);
                    g4.ZData=contourPlot(:,3);
                end 
                    g3.ZData=smoothField;
                g1.XData=x0;
                g1.YData=y0;
                g1.ZData=Alt0+max(smoothField(:));
                g2.XData=posIntruRel(1,:);
                g2.YData=posIntruRel(2,:);
                g2.ZData=interp2(smoothField, posIntruRel(1,:), posIntruRel(2,:));

                g6.XData=posSwarmRel(1);
                g6.YData=posSwarmRel(2);
                g6.ZData=interp2(smoothField, posSwarmRel(1), posSwarmRel(2));

                for ind_h=1: swarmSize
                    eval(['h',num2str(ind_h),'.XData=trajectSet(:,1,ind_h);']);
                    eval(['h',num2str(ind_h),'.YData=trajectSet(:,2,ind_h);']);
                    eval(['h',num2str(ind_h),'.ZData=trajectSet(:,3,ind_h);']);
                end
                for ind_k=1: swarmSize
                    eval(['k',num2str(ind_k),'.XData=arcPreOrigin(:,1,ind_k);']);
                    eval(['k',num2str(ind_k),'.YData=arcPreOrigin(:,2,ind_k);']);
                    eval(['k',num2str(ind_k),'.ZData=arcPreOrigin(:,3,ind_k);']);
                end
                for ind_p=1: swarmSize
                    eval(['p',num2str(ind_p),'.XData=yPred(:,1,ind_p);']);
                    eval(['p',num2str(ind_p),'.YData=xPred(:,1,ind_p);']);
                    eval(['p',num2str(ind_p),'.ZData=ones(size(xPred(:,1,ind_p)))*max(smoothField(:));']);
                end
%                 
                eval(['title(axe2,''\fontsize{16} step ',num2str(ind_t),''');' ]);
                
                posSwarmSet=cat(2, posSwarmSet, posSwarmRel);
                
                frame = getframe(fig2);
                writeVideo(writerObj1,frame);
            end
        end
        fprintf(file_swarm_start_end,'%12.8f %12.8f\r\n', mean(x0(:)), mean(y0(:)));
        
        endX=arcx(end,:,:);
        endY=arcy(end,:,:);
        endZ=arcz(end,:,:);
        arcAlt(end,:,:)=posIniIndivAlt;
        endAlt=arcAlt(end,:,:);
        trajectTemp=cat(2, arcx(end,:,:), arcy(end,:,:), arcAlt(end,:,:), timeCost(ones(size(arcAlt(end,:,:)))), arcz(end,:,:));
        trajectSet=cat(1,trajectSet,trajectTemp);

        xObs = posIntruRel(1, :); 
        yObs = posIntruRel(2, :); 
        trajectTempObs = cat(3, xObs, yObs, zeros(size(xObs)), zeros(size(xObs)), zeros(size(xObs)));
        trajectTempObs = permute(trajectTempObs, [1 3 2]); 
        trajectSetObs = cat(1, trajectSetObs, trajectTempObs);

        for ind_h=1: swarmSize
            eval(['h',num2str(ind_h),'.XData=trajectSet(:,1,ind_h);']);
            eval(['h',num2str(ind_h),'.YData=trajectSet(:,2,ind_h);']);
            eval(['h',num2str(ind_h),'.ZData=trajectSet(:,3,ind_h);']);
        end
        pause (1);
        if plotFlag
            hold (axe2);
            close(writerObj1);
        end
        fclose(file_swarm_x);
        fclose(file_swarm_y);
        fclose(file_obs_x);
        fclose(file_obs_y);
        fclose(file_swarm_kappa);
        fclose(file_swarm_omega);
        fclose(file_obs_omega);
        fclose(file_swarm_start_end); 
    end
end
