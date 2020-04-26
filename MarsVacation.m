% AEM 4301 Project - Mars Vacation With Minimum dV - Spring 2020
%
% Group members: Dan Bombeck, Lengji Liu, Jacob Meyer, RJ Nelson
%
clear all
tic
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
    HohmannWindowEM = [HohmannWindowEM;datetime(year(DatesEM(i)),month(DatesEM(i)),(day(DatesEM(i))+[-200:10:200]))];
    tofEM = [tofEM;TransferTime+(24*60*60).*[200:-10:-200]]; % seconds (Maybe modify later!!)
end
for i=1:length(DatesME)
    HohmannWindowME = [HohmannWindowME;datetime(year(DatesME(i)),month(DatesME(i)),(day(DatesME(i))+[-200:10:200]))];
    tofME = [tofME;TransferTime+(24*60*60).*[200:-10:-200]]; % seconds (Check later!!!)
end
% <debugging only>
disp(HohmannWindowEM);
disp(HohmannWindowME);
disp(TransferTime);
size(tofEM);
size(tofME);
length(HohmannWindowME(:,1));
% ^^debigging/testing only^^

% Cont: Earth2Mars V1's, V2's
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
%         [coeE,rE, vE, jdE] = planet_elements_and_sv(3, DepYear, DepMonth, DepDay,0,0,0); % Question: is planetdata.m even reliable???
%         [coeM,rM, vM, jdM] = planet_elements_and_sv(4, ArrYear, ArrMonth, ArrDay,0,0,0); % Question: is planetdata.m even reliable???
        % obtain current true anomaly of planets ("current" = 12am, 5/1/2020)
        [rE, v, jd, coem coe2] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
        [rM, v, jd, coem coe2] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
        rE = rE;
        rM = rM;
        [LambertV1,LambertV2,dum] = LambertSolverND( rE, rM, TOF, muS, dir );
        %disp(LambertV1)
        LambertV1EM(i,j) = norm(LambertV1);
        LambertV2EM(i,j) = norm(LambertV2);
    end
end
% disp(LambertV1EM)
% disp(LambertV2EM)
figure(1)
plot(LambertV1EM(1,:)), hold on
plot(LambertV1EM(2,:)), hold on
plot(LambertV1EM(3,:)), hold on
plot(LambertV1EM(4,:)), hold on
plot(LambertV1EM(5,:)), hold on
legend(datestr(DatesEM));
axis([0 40 30 35]);
min(LambertV1EM);
max(LambertV1EM);
hold off

% Cont: Mars2Earth V1's, V2's
LambertV1ME = [];
LambertV2ME = [];
for i=1:length(HohmannWindowME(:,1))
    for j=1:length(HohmannWindowME(1,:))
        % get current date
        DepYear = year(HohmannWindowME(i,j)); % year of departure
        DepMonth = month(HohmannWindowME(i,j)); % month of departure
        DepDay = day(HohmannWindowME(i,j)); % day of departure
        EarthArrivalDay = ceil(DepDay+TransferTime/(60*60*24));
        arrive_at_earth = datetime(DepYear,DepMonth,EarthArrivalDay);
        ArrYear = year(arrive_at_earth); % year of Earth arrival
        ArrMonth = month(arrive_at_earth); % month of Earth arrival
        ArrDay = day(arrive_at_earth); % day of Earth arrival
        hour = 0;
        minute = 0;
        second = 0;
        TOF = tofME(i,j);
        dir = 'pro';
        % obtain heliocentric radii (3x1 state vector) for the specified date
        % PlanetCoe = SolarSystem([thisyear thismonth thisday],0); % can't figure out 
%         [coeE,rM, vM, jdM] = planet_elements_and_sv(4, DepYear, DepMonth, DepDay,0,0,0); % Question: is planetdata.m even reliable???
%         [coeE,rE, vE, jdE] = planet_elements_and_sv(3, ArrYear, ArrMonth, ArrDay,0,0,0); % Question: is planetdata.m even reliable???
        [rE, v, jd, coem coe2] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
        [rM, v, jd, coem coe2] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
        rE = rE;
        rM = rM;
        [LambertV1,LambertV2,dum] = LambertSolverND( rM, rE, TOF, muS, dir );
        %disp(LambertV1)
        LambertV1ME(i,j) = norm(LambertV1);
        LambertV2ME(i,j) = norm(LambertV2);
    end
end
% disp(LambertV1EM)
% disp(LambertV2EM)
figure(2)
plot(LambertV1ME(1,:)), hold on
plot(LambertV1ME(2,:)), hold on
plot(LambertV1ME(3,:)), hold on
plot(LambertV1ME(4,:)), hold on
plot(LambertV1ME(5,:)), hold on
legend(datestr(DatesME));
axis([0 40 0 35]);
min(LambertV1EM);
max(LambertV1EM);
hold off



% obtain delta-v's for every Earth2Mars opportunity
TotalDV_EM = [];
for i=1:length(HohmannWindowEM(:,1))
    for j=1:length(HohmannWindowEM(1,:))
        V1 = LambertV1EM(i,j);
        V2 = LambertV2EM(i,j);
        MarsArrDate = DatesEM(i);
        MarsArrYear = year(MarsArrDate);
        MarsArrMonth = month(MarsArrDate);
        MarsArrDay = month(MarsArrDate);
        Time2MarsDays = tofEM(i,j)/(24*60*60);
        MarsArrDate = datetime(MarsArrYear,MarsArrMonth,ceil(MarsArrDay+Time2MarsDays));
        MarsArrYear = year(MarsArrDate);
        MarsArrMonth = month(MarsArrDate);
        MarsArrDay = month(MarsArrDate);
        [rM, VMars, jd, coem coe2] = PlanetData(4, MarsArrYear, MarsArrMonth, MarsArrDay,0,0,0);
        VMars = norm(VMars);
        [DV1,Beta,DELTA] =  Dv_Departure(V1,350,'earth2mars'); % departure 
        [DV2,Beta,DELTA] = Dv_Arrive(V2,VMars,500,'earth2mars');% arrival 
        TotalDV_EM(i,j) = DV1 + DV2;
    end
end
figure(3)
plot(TotalDV_EM(1,:)), hold on
plot(TotalDV_EM(2,:)), hold on
plot(TotalDV_EM(3,:)), hold on
plot(TotalDV_EM(4,:)), hold on
plot(TotalDV_EM(5,:)), hold on
legend(datestr(DatesEM));
axis([0 40 15 35])
hold off

% obtain delta-v's for every Mars2Earth opportunity
TotalDV_ME = [];
for i=1:length(HohmannWindowME(:,1))
    for j=1:length(HohmannWindowME(1,:))
        V1 = LambertV1ME(i,j);
        V2 = LambertV2ME(i,j);
        EarthArrDate = DatesME(i);
        EarthArrYear = year(EarthArrDate);
        EarthArrMonth = month(EarthArrDate);
        EarthArrDay = month(EarthArrDate);
        Time2EarthDays = tofME(i,j)/(24*60*60);
        EarthArrDate = datetime(EarthArrYear,EarthArrMonth,ceil(EarthArrDay+Time2EarthDays));
        EarthArrYear = year(EarthArrDate);
        EarthArrMonth = month(EarthArrDate);
        EarthArrDay = month(EarthArrDate);
        [rM, VEarth, jd, coem coe2] = PlanetData(4, EarthArrYear, EarthArrMonth, EarthArrDay,0,0,0);
        VEarth = norm(VEarth);
        [DV1,Beta,DELTA] =  Dv_Departure(V1,500,'mars2earth'); % departure
        [DV2,Beta,DELTA] = Dv_Arrive(V2,VEarth,350,'mars2earth');% arrival 
        TotalDV_ME(i,j) = abs(DV1)+ abs(DV2);
    end
end
figure(4)
plot(TotalDV_ME(1,:)), hold on
plot(TotalDV_ME(2,:)), hold on
plot(TotalDV_ME(3,:)), hold on
plot(TotalDV_ME(4,:)), hold on
plot(TotalDV_ME(5,:)), hold on
legend(datestr(DatesME));
axis([0 40 -30 30])
hold off

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
toc