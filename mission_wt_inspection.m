function [n_refuel,time_total] = mission_wt_inspection(weight,spec_energy,...
     Battery_lb,hover_dwld_factor,num_rotors,visc,num_wings,S,prop_efficiency,...
      c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,vel_stall,tip_speed,rotor_area,sig,fm,...
      speed_sound_SL,ipf,k_wndspd,c_wndspd,mu_wnddir,kappa_wnddir,k_wave,c_wave,...
      k_hwavetime,c_hwavetime,k_lwavetime,c_lwavetime,month,daylight);

    % wtinspection.m
    
    global lb2kg fps2mph kw2hp m2f kts2ftps % conversions

    % ========================================================================
    %
    %                       TURBINE INSPECTION PERFORMANCE
    %
    % ========================================================================
    % This section of the code computes the energy and time required to perform the
    % thermal and visual inspection of one set of turbine blades.  It accounts
    % for crosswinds from 0-25 mps (representing from still air to a typical
    % cut-out speed of ~55 mph) and a random rotor orientation.  The blade
    % length of the turbine inspected here is the NREL 5MW baseline.  The wind farm is
    % Robin Rigg, a 180MW OSW farm (Vestas V90 turbines) on the Scottish/British border

    blade_length = (126-3)/2*m2f; % (rotor diam. - hub diam.) [ft]
    farm_map = xlsread('Robin Rigg layout.xlsx','B4:C63')*m2f; % [ft]
    substation_map = xlsread('Robin Rigg layout.xlsx','E4:F5')*m2f; % [ft]

    h = figure(3);
    plot(farm_map(:,1),farm_map(:,2),'rx'); hold all;
    plot(substation_map(:,1),substation_map(:,2),'bx')
    % PL weight for Pumair gimbal+cameras = 1.822 kg (~4 lbf)
    % Design PL is 5 lb, assume that mission PL is equivalent with required
    % mounting hardware, wiring, and other adjustments

    colors = {'yx' 'mx' 'cx' 'gx' 'bx' 'kx'};
    
    xwind = wndspd; % [ft/s]

    T_above = [-20 0 20 40];
    len_T   = length(T_above);

    % initialization
    pres_alt = zeros(len_T,1);
    n_inspections = 0;
    n_turbines = length(farm_map);
    n_landing = 0;
    time = 0;
    disk_radius = sqrt(rotor_area*4/pi);
    global loadflag
    loadflag = 0;

    % TO condition
    tow = weight;
    wind_direction = rand*2*pi; % This describes the "to" direction of the wind vector
    bat_energy_full = spec_energy*Battery_lb*lb2kg*kw2hp; % [hp-hr]

    % Vehicle transit speed
    inspection_speed = 0.654; % [ft/s]
    inspection_time=blade_length*6/inspection_speed; % [s]
    v_br = 33.1*kts2ftps; % [fps], from "baseline" a/c in AHS report
    v_be = 24.3*kts2ftps; % [fps], from "baseline" a/c in AHS report
    v_cruise = 55*kts2ftps; % [fps], from "baseline a/c in AHS report

    for k = 1:length(xwind)

        for i = 1:length(T_above)
            farm_map_use = farm_map;

            T_above_ISA = T_above(i);

            % Reset battery capacity to full
            bat_energy = bat_energy_full;

            rho = .00238;
            % obtain true altitude in ft
            altitude = 250;
            % compute T at altitude (according to ISA conditions)
            T_alt_ISA = 15-0.001981*altitude; % alt must be in ft
            % compute p/p0 (delta)
            delta = (1-6.876e-6*altitude)^5.265; % alt must be in ft
            % compute pressure altitude (m)
            altitude_pres = (518.4/0.00357*(1-delta^0.1903));
            % compute rho/rho0 (sigma)
            sigma = 288.16/(T_alt_ISA+T_above_ISA+273.16)*...
            (1-0.001981*altitude_pres/288.16)^5.256;
            % compute height for this density ratio, aka, dens altitude
            altitude_dens = (518.4/0.00357*(1-sigma^0.235));
            % compute equivalent air density based on density altitude
            rho_alt = rho*sigma;
            % true temperature at altitude
            T_alt = T_alt_ISA + T_above_ISA;
            % record the pressure altitude
            pres_alt(i) = altitude_pres;

            ct_single_rotor = (tow*(1+hover_dwld_factor)/num_rotors)/...
                (rho_alt*rotor_area*tip_speed^2);

            % BEGIN INSPECTION MISSION - TAKEOFF, TRANSITION, AND CRUISE
            % From AHS paper: 
            %    0.8% of max battery capacity to complete heli->biplane transition 
            %       (+200 ft altitude gain, 55 kts cruise, 350 ft forward displacement)
            %    0.3% of max battery to complete biplane->heli transition 
            %       (350 ft forward displacement, no altitude change)
            %    2.0% of max battery to complete heli->biplane "tumble"
            %       transition (200 ft axial altitude gain)
            while(n_inspections<n_turbines)

                % If a/c is at substation
                if(bat_energy==bat_energy_full)
                    %
                    % Determine first turbine
                    %
                    x=farm_map_use(:,1)-substation_map(1,1);
                    y=farm_map_use(:,2)-substation_map(1,2);
                    proximity=sqrt(x.^2+y.^2);
                    index=find(min(proximity)==proximity);
                    cruise_dist = proximity(index); % [ft]
                    if(x(index)<0)
                        cruise_direction = atan(y(index)/x(index))+pi; % [rad]
                    else
                        cruise_direction = atan(y(index)/x(index)); % [rad]
                    end
                    figure(h)
                    plot(farm_map_use(index,1),farm_map_use(index,2),colors{mod(n_landing,6)+1})

                    % TAKEOFF FROM SUBSTATION (0 MSL) AND TRANSITION TO CRUISE
                    if(cruise_dist<350)
                        bat_energy = 0.992*bat_energy_full; % [hp-hr]
                        time = time+17; % From AHS report
                    else 
                        cruise_dist = cruise_dist-350;

                        % Target airspeed is adjusted to keep cruising ground
                        % speed at 55 kts +/- 1/2 the wind
                        ground_speedx = v_cruise*cos(cruise_direction)+0.5*xwind(k)*cos(wind_direction);
                        ground_speedy = v_cruise*sin(cruise_direction)+0.5*xwind(k)*sin(wind_direction);
                        ground_speed = sqrt(ground_speedx^2+ground_speedy^2);
                        TAS_x = v_cruise*cos(cruise_direction)-0.5*xwind(k)*cos(wind_direction);
                        TAS_y = v_cruise*sin(cruise_direction)-0.5*xwind(k)*sin(wind_direction);
                        TAS = sqrt(TAS_x^2+TAS_y^2);
                        M = TAS/speed_sound_SL;
                        
                        power_required = biplane_drag_power(TAS,M,...
                                rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                                c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e); % [hp]
                        bat_energy = power_required*cruise_dist/ground_speed/3600; % [hp-hr]
                        bat_energy = 0.992*bat_energy_full-bat_energy; % [hp-hr]
                        time = time+cruise_dist/ground_speed+17; % [s]
                    end

                    % TURBINE INSPECTION
                    [energy_inspection] = wt_inspection(xwind(k),rho_alt,visc,...
                        num_wings,S,c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,inspection_speed,...
                        vel_stall,tip_speed,tow,ct_single_rotor,rotor_area,sig,fm,num_rotors,...
                        prop_efficiency,e,speed_sound_SL,ipf,blade_length); % [hp-hr]
                    bat_energy = bat_energy-energy_inspection; % [hp-hr]

                    time = time+inspection_time; % [s]
                    n_inspections = n_inspections + 1;
                end

                % Continue turbine inspections
                if(n_inspections~=n_turbines)
                    %
                    % Check the amount of energy required to go home
                    %
                    x = substation_map(1,1)-farm_map_use(index,1);
                    y = substation_map(1,2)-farm_map_use(index,2);
                    [land_energy,land_time]=cruiseenergy(v_cruise,speed_sound_SL,...
                                rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                                c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,...
                                x,y,wind_direction,xwind(k));
                    %
                    % Check for ability to inspect next turbine
                    % Assume that the aircraft travels to the nearest turbine
                    %
                    x=farm_map_use(index,1);
                    y=farm_map_use(index,2);
                    farm_map_use([index],:)=[]; % Eliminate inspected turbine from index of targets
                    x=farm_map_use(:,1)-x;
                    y=farm_map_use(:,2)-y;

                    proximity=sqrt(x.^2+y.^2);
                    index=find(proximity);
                    index=find(min(proximity(index))==proximity);

                    [cruise_energy,cruise_time]=cruiseenergy(v_cruise,speed_sound_SL,...
                                rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                                c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,...
                                x(index),y(index),wind_direction,xwind(k)); 
                    bat_reqd = cruise_energy+...
                            wt_inspection(xwind(k),rho_alt,visc,...
                                num_wings,S,c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,inspection_speed,...
                                vel_stall,tip_speed,tow,ct_single_rotor,rotor_area,sig,fm,num_rotors,...
                                prop_efficiency,e,speed_sound_SL,ipf,blade_length); % [hp-hr]

                    % Reserve power reqd (5 min @ Vbe)
                    power_required = biplane_drag_power(xwind(k)+v_be,(xwind(k)+v_be)/speed_sound_SL,...
                                rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                                c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e); % [hp]
                    bat_reqd = bat_reqd+power_required*5*60/3600; % [hp-hr]
                end

                % Distance to substation from next turbine [ft]
                x = substation_map(1,1)-farm_map_use(index,1);
                y = substation_map(1,2)-farm_map_use(index,2);
                [cruise_energy,cruise_time] = cruiseenergy(v_cruise,speed_sound_SL,...
                        rho_alt,num_rotors,visc,num_wings,S,tow,prop_efficiency,...
                        c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,...
                        x,y,wind_direction,xwind(k));
                bat_reqd=cruise_energy+bat_reqd+0.0235*bat_energy_full; % [hp-hr]  

                % This is an evaluatory statement
                if(n_inspections==2)
                    energy_percent_battery(i,k) = energy_inspection/bat_energy_full*100; % [%]
                end
                if(bat_reqd>=bat_energy) % Continue inspection - land & recharge
                    n_landing = n_landing+1;
                    bat_energy = bat_energy_full;
                    time = time+land_time;
                    figure(h); hold all;
                elseif(n_inspections==n_turbines) % Inspections complete - land
                    bat_energy = bat_energy-land_energy; % [hp-hr]
                    time = time+land_time;
                else % Continue inspection - don't land
                    figure(h)
                    plot(farm_map_use(index,1),farm_map_use(index,2),colors{mod(n_landing,6)+1})
                    bat_energy = bat_energy-bat_reqd+power_required*5*60/3600; % [hp-hr]
                    n_inspections = n_inspections + 1;
                    time = time+cruise_time+inspection_time;
                end  
            end % h - inspection loop
            n_refuel(i,k) = n_landing;
            time_total(i,k)=time;
            n_inspections = 0;
            n_landing = 0;
        end %i

    end %k

    figure(4)
    set(gcf,'Color',[1,1,1])
    set(0,'DefaultAxesFontSize',27.5)
    set(0,'DefaultLineLineWidth',5)
    set(0,'DefaultLineMarkerSize',10)
    grid off;
    xlabel('Crosswind Speed (mph)')
    ylabel('Battery Used (%)')
    title('Battery Usage for Single Turbine Inspection')
    for i = 1:length(T_above)
        hold all
        plot(xwind*fps2mph,energy_percent_battery(i,1:length(xwind)))
        ltext{i} = sprintf('%g K from ISA',T_above(i));
    end
    legend(ltext)

%     figure(5)
%     set(gcf,'Color',[1,1,1])
%     set(0,'DefaultAxesFontSize',27.5)
%     set(0,'DefaultLineLineWidth',5)
%     set(0,'DefaultLineMarkerSize',10)
%     grid off;
%     xlabel('Crosswind Speed (mph)')
%     ylabel('Battery Used (%)')
%     title('Battery Usage for travel to next turb+return+land')
%     for i = 1:length(T_above)
%         hold all
%         plot(xwind*fps2mph,energy_percent_battery(i,1:length(xwind)))
%         ltext{i} = sprintf('%g K from ISA',T_above(i));
%     end
%     legend(ltext)
end