% inspection.m
%
% This code initializes data for an runs the wind turbine inspection
% missions for AirEZ, CTV, and CTV with UAS support
%

% Written by Br. Marius Strom, TOR

% clean it up
clear
close all

%%%%%%%%%%%%%%%%%%% Conversions %%%%%%%%%%%%%%%%%%%%

global W2hp mph2mps mps2mph kg2lb lb2kg mph2kts kts2mph mph2fps fps2mph ...
    kts2mps mps2kts hp2kw kw2hp m2f f2m kts2ftps ftps2kts

% Watts to hp ---> 1 Watt = .001341 hp
W2hp = .001341; hp2W = 1/W2hp;
% mph to m/s ---> 1 mph = .44704 m/s, 1 m/s = 2.2689 mph
mph2mps = .44704; mps2mph = 2.2369;
kg2lb = 2.2046; lb2kg = 1/kg2lb;
mph2kts = 0.868976242; kts2mph = 1/mph2kts;
mph2fps = 1.46667; fps2mph = 1/mph2fps;
kts2mps = 0.51444; mps2kts = 1/kts2mps;
hp2kw = hp2W/1000; kw2hp = 1/hp2kw;
m2f = 3.28; f2m = 1/m2f;
kts2ftps = kts2mph*mph2fps; ftps2kts = 1/kts2ftps;

% Load AirEZ data
load('airEZ_baseline.mat');

% CTV wave restrictions
ctv_max_wave=2.5*m2f;

% Load buoy data and elinate out of range data points
buoydata=xlsread('2012buoydata.xlsx','B3:J8074');
for i=1:length(buoydata)
  if(buoydata(i,5)<90||buoydata(i,6)<90||buoydata(i,4)>360) 
    buoymonth(i)=buoydata(i,1);
    buoyday(i)=buoydata(i,2);
    buoyhour(i)=buoydata(i,3);      % [hr]
    wnddir(i)=buoydata(i,4)*pi/180; % [rad]
    wndspd(i)=buoydata(i,5)*m2f;    % [ft]
    waveht(i)=buoydata(i,6)*m2f;    % [ft]
  end
end
daylight = xlsread('CapeMay2012Daylight.xlsx',2,'B2:M33'); % Daily qty of daylight in [s]
days=[31 29 31 30 31 30 31 31 30 31 30 31];
monthnames={'January' 'February' 'March' 'April' 'May' 'June' 'July'...
            'August' 'September' 'October' 'November' 'December'};

% Bin data by month, generate distributions for wave heights, wave
% persistence, wind speed, wind persistence, wind and direction
% Wave/wind data: Weibull distribution
%  - k (shape--b in MATLAB), c (scaling--a in MATLAB)
% Wind direction: von Mises distribution

n=0;
for i=1:12
    index=find(buoymonth==i);
    month_mean = mean(wndspd(index));

    %
    % Wind speed distribution (Weibull)
    %
    
    k_wndspd(i) = (std(wndspd(index)/month_mean)^-1.086);
    c_wndspd(i) = month_mean/gamma(1+1/k_wndspd(i));
    
    %
    % Wind direction distribution (von Mises)
    %
    
    [mu(i),~,~]=circ_mean(wnddir(index)');
    kappa(i)=circ_kappa(wnddir(index)');
    [pdf,bins]=circ_vmpdf(wnddir(index)',mu(i),kappa(i));
    wnddir_pdf(i,1:length(wnddir(index)))=pdf;

    %
    % Wave height distribution (Weibull)
    %
        
    month_mean = mean(waveht(index));
    k_wave(i) = (std(waveht(index)/month_mean)^-1.086);
    c_wave(i) = month_mean/gamma(1+1/k_wave(i));
    index=find(waveht(index)<ctv_max_wave);
    
    %
    % Wave persistence distribution (Weibull)
    % Find distribution for persistence times of waves above/below critical
    %
    
    high_counter=0;
    low_counter=0;
    for ii=1:length(index)
        if waveht(index(ii))>=ctv_max_wave
            high_counter=high_counter+1; % [hr]
        elseif high_counter~=0
            high_collector(ii)=high_counter+rand; % [hr]
            high_counter=0;
        end
        
        if waveht(index(ii))<ctv_max_wave
            low_counter=low_counter+1; % [hr]
        elseif low_counter~=0
            low_collector(ii)=low_counter+rand; % [hr]
            low_counter=0;
        end
    end
    high_collector=high_collector*3600; % [s]
    low_collector=low_collector*3600; % [s]
    
    index=find(high_collector);
    high_wave_time(i,1:length(index))=high_collector(index);
    month_mean=mean(high_collector(index));
    k_hwavetime(i) = (std(high_collector(index))/month_mean)^-1.086;
    c_hwavetime(i) = month_mean/gamma(1+1/k_hwavetime(i));
    
    index=find(low_collector);
    low_wave_time(i,1:length(index))=low_collector(index);
    month_mean=mean(low_collector(index));
    k_lwavetime(i) = (std(low_collector(index))/month_mean)^-1.086;
    c_lwavetime(i) = month_mean/gamma(1+1/k_lwavetime(i));
end

% Select date
month=0;
while(month==0)
    month=round(rand*12);
end
month_index=find(buoymonth==month);
day_index=buoyday(month_index);


% [n_landings,time]=mission_wt_inspection(weight,spec_energy,...
%     Battery_lb,hover_dwld_factor,num_rotors,visc,num_wings,S,prop_efficiency,...
%      c,Cd0_struts,Cd0_gear,Cd0_wing,A_gear,b,e,vel_stall,tip_speed,rotor_area,sig,fm,...
%      speed_sound_SL,ipf,wnddir,wndspd,daylight,wnddir_pdf,wndspd_pdf,days,month,...
%      day);