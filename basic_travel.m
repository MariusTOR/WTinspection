function [tt,tb,index] = basic_travel(v_cruise,dest_list,base_x,base_y,...
    current_x,current_y)
% This function takes in a speed, destination list, current location, and a
% base location.  It determines the closest destination to the present
% location and the amount of time required to travel to that point and then
% to the base location

% Cruise to nearest destination from present position
x=dest_list(:,1)-current_x;
y=dest_list(:,2)-current_y;
proximity=sqrt(x.^2+y.^2);
index=find(min(proximity)==proximity);
cruise_dist = proximity(index); % [ft]
tt=cruise_dist/v_cruise; % [s]

% Cruise from nearest destination to base
x=base_x-dest_list(index,1);
y=base_y-dest_list(index,2);
proximity=sqrt(x.^2+y.^2); % [ft]
tb=proximity/v_cruise; % [s]
end

