function [energy_inspection] = airezwt_inspection(v_inf,rho_alt,visc,...
    num_wings,S,c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,inspection_speed,...
    vel_stall,tip_speed,tow,ct_single_rotor,rotor_area,sig,fm,num_rotors,...
    prop_efficiency,e,speed_sound_SL,ipf,blade_length)

    % Initialize global variables for conversion factors
    global W2hp kts2ftps 
    
    %
    % Wind turbine inspection routine
    %
    
    % rotor angle, 1st blade
    blade_angle = rand*2*pi; % [rad]

    power_inspection = [0 0 0];
    for h = 1:3
        % ============================================================
        %             DRAG APPROXIMATIONS IN CROSSWIND
        % ============================================================

        % Flow characterization for specified blade angles - assume roc is
        % small, hence momentum theory is still valid
        roc = inspection_speed*sin(blade_angle+(h-1)*2*pi/3);
        M = abs((v_inf+inspection_speed*cos(blade_angle+(h-1)*2*pi/3))/speed_sound_SL);
        mu_xwind = abs((v_inf+inspection_speed*cos(blade_angle+(h-1)*...w
            2*pi/3))/tip_speed);

        % For hovering case
        if((v_inf+inspection_speed*cos(blade_angle+(h-1)*2*pi/3))<(vel_stall*kts2ftps))
            
            % PARASITE POWER, [hp]
            [power_reqd_parasitic,drag_total] = heli_drag_power(v_inf,M,rho_alt,...
                num_rotors,visc,num_wings,S,c,Cd0_struts,Cd0_gear,A_gear,b);

            climb_inflow = roc/tip_speed;

            % INFLOW COMPUTATION
            alpha_rotor = drag_total/tow;

            % Newton-Rhapson solver for inflow
            % Define first derivative of the equation of interest for use in the
            % solver (i.e., 0 = d(ueq)/d(li))
            inflow = sqrt(ct_single_rotor/2);
            ueqn = @(inflow) (climb_inflow*cos(alpha_rotor)+mu_xwind*tan(alpha_rotor)+...
                ct_single_rotor/(2*sqrt(mu_xwind^2 + inflow^2)));
            dueqn = @(inflow) (-1-0.5*ct_single_rotor*inflow*(mu_xwind^2+inflow^2)^-1.5);
            inflow = newton(ueqn,dueqn,inflow,1E-7); % Tol: 1E-7

            % INDUCED POWER PER ROTOR, [ft-lb/s]
            cp_induced = ipf*(ct_single_rotor^2)/(2*sqrt(inflow^2+mu_xwind^2));
            power_induced = cp_induced*rho_alt*rotor_area*tip_speed^3;

            % PROFILE POWER PER ROTOR, [ft-lb/s]
            cp_profile    = sig*.01/8*(1+5*mu_xwind^2);
            power_profile = cp_profile*rho_alt*rotor_area*tip_speed^3;

            % TOTAL POWER REQ'D, [hp]
            power_rotors = (power_induced + power_profile)...
                *num_rotors/fm/550+power_reqd_parasitic;

        else % For when aircraft can fly in biplane mode into the crosswind
            
            % PARASITE POWER, [hp]
            power_required_parasitic = biplane_drag_power(v_inf,M,...
                rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e);

            % CLIMB POWER, [hp]
            power_climb = roc*tow/550;

            % TOTAL POWER REQ'D, [hp]
            power_rotors = (power_required_parasitic+power_climb)/...
                prop_efficiency;
        end

        % TOTAL POWER REQ'D FOR 2-PASS INSPECTION OF 1 BLADE, [hp]
        camera_power = 16.7*W2hp;
        power_inspection(h) = (power_rotors + camera_power)*2;

    end %h
    energy_inspection=sum(power_inspection*blade_length/inspection_speed/3600);
end