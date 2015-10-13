function [power_reqd_parasitic,drag_total] = heli_drag_power(v_inf,M,rho_alt,num_rotors,visc,num_wings,S,c,Cd0_struts,Cd0_gear,A_gear,b)
% This function computes the parasitic power requirement of the baseline
% quadrotor-biplane tailsitter concept (AirEZ) when in its quadcopter
% configuration given a relative wind velocity greater than its stall speed

% Obtain equivalent flat plate drag area in crosswind
% Flat plate area of fuselage as viewed from the side
A_fuse = (650/144); % [ft^2]
A_struts = .36/12*(14/12)*num_rotors; % [ft^2]
Re_struts = rho_alt*v_inf*3/12/visc;
Re_gear = rho_alt*v_inf*1/visc;
Cf_struts = 1.328/sqrt(Re_struts);
Cf_gear = 1.328/sqrt(Re_gear);

d_nac = 3.5/12;
l_nac = 2*d_nac;
S_nacelle = d_nac^2*pi/4;
S_wet_nacelle = 4*pi*(d_nac/2)^2; % Approx as sphere
Re_nacelle = rho_alt*v_inf*3.5/24/visc;
Cf_nacelle = 1.328/sqrt(Re_nacelle);
CD_nacelle = inf_cyl_drag(Re_nacelle);


d_fuse = sqrt((41*17)/pi)/12; % equivalent diameter
l_fuse = 17/12;
S_fuse = 650/144;
S_wet_fuse = 2198/144;
Re_fuse = rho_alt*v_inf*d_fuse/visc;
Cf_fuse = 1.328/sqrt(Re_fuse);
Cd0_fuse = inf_cyl_drag(Re_fuse);

drag_fuselage = 0.5*rho_alt*v_inf^2*A_fuse*(Cd0_fuse+Cf_fuse)+...
    0.5*rho_alt*v_inf^2*A_struts*(Cd0_struts+Cf_struts)+0.5*rho_alt*v_inf^2*A_gear*...
    (Cd0_gear+Cf_gear)+0.5*rho_alt*v_inf^2*(CD_nacelle+Cf_nacelle)*S_nacelle;

% Approximate wing as slim rectangular prism
Re_wing = rho_alt*v_inf*b/visc; % wing length - assume a/c will align wings parallel wind
Cf_wing = 1.328/sqrt(Re_wing); % Assume laminar flow
d_wing = sqrt((c*c*0.137)/pi); % Equivalent diameter
l_wing = b;
f_ld_wing = 1+60/(l_wing/d_wing)+.0025*(l_wing/d_wing);
f_m_wing = 1-.08*M^1.45;
S_wet_wing = 2*(1+0.5*0.137)*b*c; % Approximation of wing wetted area, 13.7% thick airfoil
Cd0_wing = Cf_wing*f_ld_wing*f_m_wing*S_wet_wing/S*1.2;

drag_wings = 0.5*rho_alt*v_inf^2*Cd0_wing*S*num_wings;

% Total parasitic drag
drag_total = drag_wings+drag_fuselage;

% Parasite power in hover, with crosswind, [hp]
power_reqd_parasitic = (drag_total*v_inf)/550;

end

