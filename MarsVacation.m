% AEM 4301 Project - Mars Vacation With Minimum dV - Spring 2020
%
% Group members: Dan Bombeck, Yinjie Liu, Jacob Meyer, RJ Nelson
%
clear all
% define constants
muE = 0.39860e6; % km^3/s^2 (gravitational parameter of Earth)
muS = 132712e6; % km^3/s^2 (gravitational parameter of Sun)
muM = 0.042828e6; % km^3/s^2 (gravitational parameter of mars)
soiE = 0.929e6; % km (sphere of influence radius for Earth)
soiM = 0.578e6; % km (sphere of influence radius for Mars)

% check for "ideal" Hohmann EM/ME transfer dates within the next ten years
[DatesEM, DatesME, TransferTime] = opportunities(10); % E->M and M->E

% generate V's, times of flight for +-50 days from ideal Earth-to-Mars and
% Mars-to-Earth transfers using a Lambert solver
tofEM = []; % seconds (attempted transfer time from Earth to Mars
tofME = [];% seconds (attempted transfer time from Mars to Earth)
HohmannWindowEM = [];
HohmannWindowME = [];
for i=1:length(DatesEM)
    HohmannWindowEM = [HohmannWindowEM;datetime(year(DatesEM(i)),month(DatesEM(i)),(day(DatesEM(i))+[-50:1:50]))];
    tofEM = [tofEM;TransferTime+(24*60*60).*[50:-1:-50]]; % seconds (Maybe modify later!!)
end
for i=1:length(DatesME)
    HohmannWindowME = [HohmannWindowME;datetime(year(DatesME(i)),month(DatesME(i)),(day(DatesME(i))+[-50:1:50]))];
    tofME = [tofME;TransferTime+(24*60*60).*[50:-1:-50]]; % seconds (Check later!!!)
end
% <debugging only>
disp(HohmannWindowEM)
disp(HohmannWindowME)
disp(TransferTime)
size(tofEM)
size(tofME)
length(HohmannWindowME(:,1))
% ^^debigging/testing only^^

LambertV1EM = [];
LambertV2EM = [];

for i=1:length(HohmannWindowEM(:,1))
    for j=1:length(HohmannWindowEM(1,:))
        % get current date
        DepYear = year(HohmannWindowEM(i,j)); % year of departure
        DepMonth = month(HohmannWindowEM(i,j)); % month of departure
        DepDay = day(HohmannWindowEM(i,j)); % day of departure
        MarsArrivalDay = ceil(DepDay+TransferTime/(60*60*24));
        arrive_at_mars = datetime(DepYear,DepMonth,MarsArrivalDay);
        ArrYear = year(arrive_at_mars); % year of Mars arrival
        ArrMonth = month(arrive_at_mars); % month of Mars arrival
        ArrDay = day(arrive_at_mars); % day of Mars arrival
        hour = 0;
        minute = 0;
        second = 0;
        TOF = tofEM(i,j);
        dir = 'pro';
        % obtain heliocentric radii (3x1 state vector) for the specified date
        % PlanetCoe = SolarSystem([thisyear thismonth thisday],0); % can't figure out 
        [coeE,rE, vE, jdE] = planet_elements_and_sv(3, DepYear, DepMonth, DepDay,0,0,0); % Question: is planetdata.m even reliable???
        [coeM,rM, vM, jdM] = planet_elements_and_sv(4, ArrYear, ArrMonth, ArrDay,0,0,0); % Question: is planetdata.m even reliable???
        rE = rE';
        rM = rM';
        [LambertV1,LambertV2,dum] = LambertSolver( rE, rM, TOF, muS, dir );
        %disp(LambertV1)
        LambertV1EM(i,j) = norm(LambertV1);
        LambertV2EM(i,j) = norm(LambertV2);
    end
end
% disp(LambertV1EM)
% disp(LambertV2EM)
plot(LambertV1EM(1,:)), hold on
plot(LambertV1EM(2,:)), hold on
plot(LambertV1EM(3,:)), hold on
plot(LambertV1EM(4,:)), hold on
plot(LambertV1EM(5,:)), hold on
legend(datestr(DatesEM));
axis([0 100 30 35])
min(LambertV1EM)
max(LambertV1EM)
% obtain the best dv, slowest time of flight from each E-M opportunity

% obtain the best dv, slowest time of flight from each M-E opportunity

% obtain "pure" Hohmann tansfer case total DV and [DVs 1-4]

% obtain DV and DVs [1-4] for fastest optimal DV round-trip case

% plots and shit
    % plot how mission delta-v varies with TOF for best-case round-trip
    
    % plot best Hohmann transfer 
    
    % plot (seperately) two best round-trip transfers
    
    % tables for the Hohmann and fastest round-trip case
        % "ideal" Hohmann case
            % DV to go from LEO to a hyperbolic orbit at Earth
            
            % DV (if any) upon leaving Earth's SOI
            
            % DV (if any) upon arriving at Mars' SOI
            
            % DV to acheive our parking orbit at Mars
        % more realistic fastest round-trip case
            % DV to go from LEO to a hyperbolic orbit at Earth
            
            % DV (if any) upon leaving Earth's SOI
            
            % DV (if any) upon arriving at Mars' SOI
            
            % DV to acheive our parking orbit at Mars

%% Sources:
% -----------------------------------------------------------------------
% muE, muS, muM
    % https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    % https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
    % https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
% soiE, soiM
    % https://en.wikipedia.org/wiki/Sphere_of_influence_(astrodynamics)
