function [fieldPotential] = OvalMapWeights (input, allPOI, intensity, alpha, guardDist, obsFlag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw potential fields according to Electrostatic theories. 
% field = velocity/(d-d0)2 
% 
% input: background
% allPOI: field centers, where fields are drawn 
% intensity: field intensity
% alpha: offset, making the field wide 
% 
% fieldPotential: output 
% input=zeros(300); 
% allPOI=[10; 10]; 
% intensity=[1]; 
% alpha=.9; 

[map_row,map_col]=size(input); 
allPOI=round(allPOI);
edge_x=allPOI(1,:);
edge_y=allPOI(2,:);
[~, edge_len]=size(allPOI);
weight=zeros(map_row, map_col, edge_len);
weight_len = max(map_row, map_col);
for ind_i=1: edge_len 
    edge_val=intensity(ind_i);
    x1=edge_x(ind_i);
    y1=edge_y(ind_i);
    weight_range_y=[y1-weight_len:y1+weight_len];
    weight_range_x=[x1-weight_len:x1+weight_len];
    weight_range_y (weight_range_y<1)=[];
    weight_range_y (weight_range_y>map_row)=[];
    weight_range_x (weight_range_x<1)=[];
    weight_range_x (weight_range_x>map_col)=[];
    [weight_range_Y, weight_range_X]=meshgrid(weight_range_y, weight_range_x);
    
    edge_X1=x1*ones(size(weight_range_X));
    edge_Y1=y1*ones(size(weight_range_Y));
    
    d1=sqrt((weight_range_Y-edge_Y1).^2+(weight_range_X-edge_X1).^2);
    
    d2=d1-guardDist;
    if obsFlag 
        d1(d1<=guardDist)=guardDist; 
    end 
    


    weightTemp = (d1+alpha*10).^(-2); 
    
    weight_uni=weightTemp./max(max(weightTemp));
    weight_scal=weight_uni.*edge_val;
    weight(1:size(weight_scal, 2), 1:size(weight_scal, 1), ind_i)=weight_scal'; 
end
fieldPotential=input+sum(weight, 3);