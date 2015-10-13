function [power_reqd_parasitic] = biplane_drag_power(v_inf,M,rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e)
% This function computes the parasitic power requirement of the baseline
% quadrotor-biplane tailsitter concept (AirEZ) when in its biplane
% configuration given a relative wind velocity greater than its stall speed

A_fuse = (260/144); % [ft^2]
A_struts = .36/12*(20/12)*num_rotors; % [ft^2]
Re_fuse = rho_alt*v_inf*3.17/visc;
Re_wing = rho_alt*v_inf*c/visc;
Re_struts = rho_alt*v_inf*3/12/visc;
Re_gear = rho_alt*v_inf*1/visc;
Cf_fuse = 1.328/sqrt(Re_fuse);
Cf_struts = 1.328/sqrt(Re_struts);
Cf_gear = 1.328/sqrt(Re_gear);

d_nac = 3.5/12;
l_nac = 2*d_nac;
S_nacelle = d_nac^2*pi/4;
S_wet_nacelle = 4*pi*(d_nac/2)^2;
Re_nacelle = rho_alt*v_inf*3.5/24/visc;
f_ld_nacelle = 1+60/(l_nac/d_nac)+.0025*(l_nac/d_nac);
f_m_nacelle = 1-.08*M^1.45;
Cf_nacelle = 1.328/sqrt(Re_nacelle);
CD_nacelle = Cf_nacelle*f_ld_nacelle*f_m_nacelle*S_wet_nacelle/S_nacelle*1.2;

d_fuse = 17/12;
l_fuse = 41/12;
S_fuse = 260/144;
S_wet_fuse = 2198/144;
Re_fuse = rho_alt*v_inf*l_fuse/visc;
f_ld_fuse = 1+60/(l_fuse/d_fuse)+.0025*(l_fuse/d_fuse);
f_m_fuse = 1-.08*M^1.45;
Cf_fuse = 1.328/sqrt(Re_fuse);
Cd0_fuse = Cf_fuse*f_ld_fuse*f_m_fuse*S_wet_fuse/S_fuse*1.2;

% flat plate area and drag of fuselage, SI (m2,N)
drag_helicopter = 0.5*rho_alt*v_inf^2*A_fuse*(Cd0_fuse+Cf_fuse)+...
    0.5*rho_alt*v_inf^2*A_struts*(Cd0_struts+Cf_struts)+...
    0.5*rho_alt*v_inf^2*A_gear*(Cd0_gear+Cf_gear)+...
    0.5*rho_alt*v_inf^2*CD_nacelle*S_nacelle;

% compute Cl_cruise
CL = 2*tow/(rho_alt*v_inf^2*S*num_wings);

% Drag of the wing, [lb]
drag_wing = (0.5*(CL)^2/(pi*(b/c)*e)*...
    rho_alt*v_inf^2*S+...
    (0.5*rho_alt*v_inf^2*S*Cd0_wing))*num_wings;

% PARASITIC POWER, [hp from lb-ft/s]
power1 = drag_helicopter*v_inf/550;
power2 = drag_wing*v_inf/550;

% TOTAL POWER, , [hp]
power_reqd_parasitic = (power1+power2)/prop_efficiency;

end

