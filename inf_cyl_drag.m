function [Cd] = inf_cyl_drag(Re)
% Lookup table for the drag coefficient of an infinite cylinder

% Read in plot of drag on an infinite cylinder
global loadflag
if(loadflag==1)
    load('re_v_cd.mat')
else
    re_cyl = xlsread('cylinder drag.xlsx','C3:C50');
    cd_cyl = xlsread('cylinder drag.xlsx','B3:B50');
    save('re_v_cd.mat','re_cyl','cd_cyl')
    loadflag = 1;
end

    if(Re<=re_cyl(1))
        re_vals = [re_cyl(1) re_cyl(2)];
        cd_vals = [cd_cyl(1) cd_cyl(2)];
        Cd = interp1(re_vals,cd_vals,Re);
    else
        index = find(Re<=re_cyl);
        re_vals = [re_cyl(index(1)-1) re_cyl(index(1))];
        cd_vals = [cd_cyl(index(1)-1) cd_cyl(index(1))];
        Cd = interp1(re_vals,cd_vals,Re);
    end
end

