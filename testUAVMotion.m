%% Little Lovely UAV %%
function [posCurr, fieldMotion] = testUAVMotion (fieldIni, motionPos, velAmp, guardDist, obsFlag) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw potential fields according to realtime motions of the swarm and the intruder. 
%
% fieldIni: the plane fields are drawn on. Can be seen as working
% space of swarms in practical applications, like a farm, a play ground, a
% park or any area UAVs woking on. 
%
% motionPos: center point of the field, where the field is drawn on
% the plane. Usually the realtime positions of swarms or obstacles. 
%
% velAmp: velocity of swarms or obstacles, it decided the intensity of
% fields. 
%
% fieldMotion: potential fields. 

[~, intruNum]=size(motionPos);
posIniX=motionPos (1,:);
posIniY=motionPos (2,:);
posCurrX=round(posIniX);
posCurrY=round(posIniY);

posCurr=[posCurrX; posCurrY];

allPOI=[posCurrX; posCurrY];
%% Edge & Corner Weighting %%
intensity=velAmp;
[fieldMotion]= OvalMapWeights (fieldIni, allPOI, intensity, .95, guardDist, obsFlag); 