function [contourPlot, sectorNum] = ContourPlot (fieldMotion, cont_level) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw real contours according to realtime positions of swarm members
% (using MATLAB buildin function). 
% 
% fieldMotion: the field contours are to be drawn on
% cont_level: level of contours 
% 
% 
% sectorNum: number of sectors of a contour. Contours sometimes can be
% discontinuous, imagine contours on two humps of a camel usually have two sectors or parts. 
% 
% contourPlot: coordinates of contours 

cont_level=cont_level(:); 
[cont_val, sector_set] = Contour_Detection (fieldMotion, cont_level);
contour_x=cont_val(1,:);
contour_y=cont_val(2,:);
[~,sector_len]=size(sector_set);
contourPlot=[]; 
sectorNum=zeros(size(cont_level)); 
for ind_i=1: sector_len
    sta_p=sector_set(1,ind_i);
    end_p=sector_set(2,ind_i);
    level2=sector_set(3,ind_i);
    cont_x=contour_x(sta_p:end_p);
    cont_y=contour_y(sta_p:end_p);
    cont_z=ones(size(cont_x))*level2;
    
    contTemp=[cont_x; cont_y; cont_z];
    contourPlot=cat(1, contourPlot, contTemp'); 
    
    sectorNum(cont_level==level2)=sectorNum(cont_level==level2)+1; 
end
