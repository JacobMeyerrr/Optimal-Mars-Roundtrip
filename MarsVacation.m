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

% check for "ideal" Hohmann E->M transfer dates within the next ten years
[DatesEM,~, TransferTime] = opportunities(10); % E->M
% check for "ideal" Hohmann M->E return transfer dates within the next 15 years
[~, DatesME,~] = opportunities(15); % M->E

% generate V's, times of flight for +-50 days from ideal Earth-to-Mars and
% Mars-to-Earth transfers using a Lambert solver
tofEM = []; % seconds (attempted transfer time from Earth to Mars
tofME = [];% seconds (attempted transfer time from Mars to Earth)
HohmannWindowEM = [];
HohmannWindowME = [];
for i=1:length(DatesEM)
    HohmannWindowEM = [HohmannWindowEM;datetime(year(DatesEM(i)),month(DatesEM(i)),(day(DatesEM(i))+[-200:5:200]))];
    tofEM = [tofEM;TransferTime+(24*60*60).*[200:-5:-200]]; % seconds (Maybe modify later!!)
end
for i=1:length(DatesME)
    HohmannWindowME = [HohmannWindowME;datetime(year(DatesME(i)),month(DatesME(i)),(day(DatesME(i))+[-200:5:200]))];
    tofME = [tofME;TransferTime+(24*60*60).*[200:-5:-200]]; % seconds (Check later!!!)
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

LambertV1EMx = [];
LambertV1EMy = [];
LambertV1EMz = [];
LambertV2EMx = [];
LambertV2EMy = [];
LambertV2EMz = [];
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
        
        LambertV1EMx(i,j) = LambertV1(1);
        LambertV1EMy(i,j) = LambertV1(2);
        LambertV1EMz(i,j) = LambertV1(3);
        LambertV2EMx(i,j) = LambertV2(1);
        LambertV2EMy(i,j) = LambertV2(2);
        LambertV2EMz(i,j) = LambertV2(3);
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
axis([0 80 30 35]);
min(LambertV1EM);
max(LambertV1EM);
hold off

% Cont: Mars2Earth V1's, V2's
LambertV1ME = [];
LambertV2ME = [];

LambertV1MEx = [];
LambertV1MEy = [];
LambertV1MEz = [];
LambertV2MEx = [];
LambertV2MEy = [];
LambertV2MEz = [];
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
        
        LambertV1MEx(i,j) = LambertV1(1);
        LambertV1MEy(i,j) = LambertV1(2);
        LambertV1MEz(i,j) = LambertV1(3);
        LambertV2MEx(i,j) = LambertV2(1);
        LambertV2MEy(i,j) = LambertV2(2);
        LambertV2MEz(i,j) = LambertV2(3);
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
axis([0 80 0 35]);
min(LambertV1EM);
max(LambertV1EM);
hold off



% obtain delta-v's for every Earth2Mars opportunity
TotalDV_EM = [];
for i=1:length(HohmannWindowEM(:,1))
    for j=1:length(HohmannWindowEM(1,:))
        V1 = LambertV1EM(i,j);
        V2 = LambertV2EM(i,j);
        
        V1vec = [LambertV1EMx(i,j);LambertV1EMy(i,j);LambertV1EMz(i,j)];
        V2vec = [LambertV2EMx(i,j);LambertV2EMy(i,j);LambertV2EMz(i,j)];
        
       EarthDepDate = DatesEM(i);
        EarthDepYear = year(EarthDepDate);
        EarthDepMonth = month(EarthDepDate);
        EarthDepDay = month(EarthDepDate);
        Time2MarsDays = tofEM(i,j)/(24*60*60);
        MarsArrDate = datetime(EarthDepYear,EarthDepMonth,ceil(EarthDepDay+Time2MarsDays));
        MarsArrYear = year(MarsArrDate);
        MarsArrMonth = month(MarsArrDate);
        MarsArrDay = month(MarsArrDate);
        [rM, VMars, jd, coem coe2] = PlanetData(4, MarsArrYear, MarsArrMonth, MarsArrDay,0,0,0);
        [rM, VEarth, jd, coem coe2] = PlanetData(3, EarthDepYear, EarthDepMonth, EarthDepDay,0,0,0);
        VMars = norm(VMars);
        [DV1,Beta,DELTA] =  Dv_Departure(V1vec,VEarth,350,'earth2mars'); % departure 
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
axis([0 80 4 10])
hold off

% obtain delta-v's for every Mars2Earth opportunity
TotalDV_ME = [];
for i=1:length(HohmannWindowME(:,1))
    for j=1:length(HohmannWindowME(1,:))
        V1 = LambertV1ME(i,j);
        V2 = LambertV2ME(i,j);
        
        V1vec = [LambertV1MEx(i,j);LambertV1MEy(i,j);LambertV1MEz(i,j)];
        V2vec = [LambertV2MEx(i,j);LambertV2MEy(i,j);LambertV2MEz(i,j)];
        
        MarsDepDate = DatesME(i);
        MarsDepYear = year(MarsDepDate);
        MarsDepMonth = month(MarsDepDate);
        MarsDepDay = month(MarsDepDate);
        Time2EarthDays = tofME(i,j)/(24*60*60);
        EarthArrDate = datetime(MarsDepYear,MarsDepMonth,ceil(MarsDepDay+Time2EarthDays));
        EarthArrYear = year(EarthArrDate);
        EarthArrMonth = month(EarthArrDate);
        EarthArrDay = month(EarthArrDate);
        [rM, VEarth, jd, coem coe2] = PlanetData(3, EarthArrYear, EarthArrMonth, EarthArrDay,0,0,0);
        [rM, VMars, jd, coem coe2] = PlanetData(4, MarsDepYear, MarsDepMonth, MarsDepDay,0,0,0);
        VEarth = norm(VEarth);
        [DV1,Beta,DELTA] =  Dv_Departure(V1vec,VMars,500,'mars2earth'); % departure
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
axis([0 80 4 10])
hold off



% obtain the best dv, slowest time of flight from each E-M opportunity
for i=1:length(TotalDV_EM(:,1))
    min = inf;
    for j=1:length(TotalDV_EM(1,:))
        if(TotalDV_EM(i,j) < min)
            min = TotalDV_EM(i,j);
            index = j;
            EMindex(i) = j;
        end
    end
    minDV_EM(i) = min;
    minDV_Date_EM(i) = HohmannWindowEM(i,index);
    index = inf;
end

% obtain the best dv, slowest time of flight from each M-E opportunity
for i=1:length(TotalDV_ME(:,1))
    min = inf;
    index = inf;
    for j=1:length(TotalDV_ME(1,:))
        if(TotalDV_ME(i,j) < min)
            min = TotalDV_ME(i,j);
            index = j;
            MEindex(i) = j;
        end
    end
    minDV_ME(i) = min;
    minDV_Date_ME(i) = HohmannWindowME(i,index);
end

% obtain best "roundtrip" (lowest DV) case and dates
minEM = inf; 
minME = inf;
DVmin = inf;
tofEM = tofEM./(24*3600);
for i=1:length(minDV_EM)
    for j=1:length(minDV_ME)
        if(minDV_Date_EM(i)+tofEM(i,EMindex(i))<minDV_Date_ME(j))
            DV = minDV_EM(i) + minDV_ME(j);
            if(DV < DVmin)
                minEM = minDV_Date_EM(i);
                minME = minDV_Date_ME(j);
                DVmin = DV
                minEM_index = i;
                minME_index = j;
            end
        end
    end
end

disp("Minimum DV Earth to Mars date: " + datestr(minDV_Date_EM(minEM_index)))
disp("Minimum DV Mars to Earth date: " + datestr(minDV_Date_ME(minME_index)))
disp("Minimum Round-trip DV for this trip is: "+num2str(DVmin))


% obtain "pure" Hohmann tansfer case total DV and [DVs 1-4]

% obtain DV and DVs [1-4] for fastest optimal DV round-trip case

% plots and shit
    % plot how mission delta-v varies with TOF for best-case round-trip
    
    % plot the "ideal" Hohmann transfer trajectory
    
    % plot the best round-trip transfer trajectory
    
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