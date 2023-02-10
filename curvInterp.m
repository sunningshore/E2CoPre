function [xPoly, yPoly]= curvInterp (input_x, input_y, N) 
% input_x = xStart; 
% input_y = yStart; 

% N=10; 
bound_len=length(input_x);
xPoly=[];
yPoly=[];
for ind_i = 1: bound_len-1 
    
    x_temp=[input_x(ind_i),input_x(ind_i+1)];
    y_temp=[input_y(ind_i),input_y(ind_i+1)];
    
    slope_temp=(y_temp(2)-y_temp(1))/(x_temp(2)-x_temp(1)); 
    if ismember(slope_temp, [Inf, -Inf]) 
        y_set=linspace(y_temp(1), y_temp(2), N); 
        x_set=x_temp(1)*ones(size(y_set)); 
    else 
        x_set=linspace(x_temp(1), x_temp(2), N);
        y_set=slope_temp*(x_set-x_temp(1))+y_temp(1); 
    end 
    xPoly=cat(2,xPoly,x_set);
    yPoly=cat(2,yPoly,y_set);
end