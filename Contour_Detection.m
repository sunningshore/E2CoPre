function [cont_val, sector_set] = Contour_Detection (input,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB buildin function for drawing contours. 
figure (100)
[cont_val,~]=contour(input,n);
contour_x=cont_val(1,:);
contour_y=cont_val(2,:);
contour_len=length(contour_x);
ind_ini=1;
sector_set=[];
while (ind_ini <= contour_len)
    x_ini=contour_x(ind_ini);
    y_ini=contour_y(ind_ini);
    sec_start=ind_ini+1;
    sec_end=y_ini+ind_ini;
    sector_ind=[sec_start;sec_end];
    sector_val=x_ini;
    sector_set=cat(2,sector_set,[sector_ind;sector_val]);
    ind_ini=sec_end+1;
end
close 100;
end 