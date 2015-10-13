function [energy_reqd,time_reqd] = cruiseenergy(TAS,speed_sound_SL,...
                            rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                            c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,...
                            dx,dy,wind_direction,xwind)
% This determines the power required to cruise in a given direction at a
% given airspeed given a crosswind
cruise_dist = sqrt(dx^2+dy^2); % [ft]
if(dx<0)
    cruise_direction = atan(dy/dx)+pi; % [rad]
else
    cruise_direction = atan(dy/dx); % [rad]
end

ground_speedx = TAS*cos(cruise_direction)+0.5*xwind*cos(wind_direction);
ground_speedy = TAS*sin(cruise_direction)+0.5*xwind*sin(wind_direction);
ground_speed = sqrt(ground_speedx^2+ground_speedy^2);
TAS_x = TAS*cos(cruise_direction)-0.5*xwind*cos(wind_direction);
TAS_y = TAS*sin(cruise_direction)-0.5*xwind*sin(wind_direction);
TAS = sqrt(TAS_x^2+TAS_y^2);
M = TAS/speed_sound_SL;

v_cruisex = TAS*cos(cruise_direction)-xwind*cos(wind_direction);
v_cruisey = TAS*sin(cruise_direction)-xwind*sin(wind_direction);
TAS = sqrt(v_cruisex^2+v_cruisey^2); % [ft/s]

power_reqd = biplane_drag_power(TAS,M,...
        rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
        c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e);  % [hp]
energy_reqd = power_reqd*cruise_dist/ground_speed/3600; % [hp-hr]
time_reqd = cruise_dist/ground_speed; % [s]
end

