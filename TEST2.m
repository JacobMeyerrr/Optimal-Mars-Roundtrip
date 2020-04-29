%
% Reconstruct Jacob's "MarsVacation script
%
% 
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
%disp(HohmannWindowEM);
%disp(HohmannWindowME);
%disp(TransferTime);
%size(tofEM);
%size(tofME);
%length(HohmannWindowME(:,1));
% ^^debigging/testing only^^
xindx = -200:5:200;
%%
% Cont: Earth2Mars V1's, V2's
LambertV1EM = zeros(5,81);
LambertV2EM = zeros(5,81);


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
        [rE, vE, ~, ~ , ~] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
        [rM, vM, ~, ~, ~] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
        %rE = rE;
        %rM = rM;
        [LambertV1,LambertV2,~] = LambertSolverND( rE, rM, TOF, muS, dir );
        %disp(LambertV1)
        LambertV1EM(i,j) = norm(LambertV1-vE)+Dv_Departure(LambertV1,vE,350,'e2m');
        LambertV2EM(i,j) = norm(LambertV2-vM)+Dv_Arrive(LambertV2,vM,500,'e2m');
        
        
    end
end
LambertEM_Tot = LambertV1EM+LambertV2EM;
% disp(LambertV1EM)
% disp(LambertV2EM)
figure(1)
xindx = -200:5:200;
plot(xindx,LambertEM_Tot)
%plot(xindx,LambertV1EM(1,:)), hold on
%plot(xindx,LambertV1EM(2,:))
%plot(xindx,LambertV1EM(3,:))
%plot(xindx,LambertV1EM(4,:))
%plot(xindx,LambertV1EM(5,:))
legend(datestr(DatesEM));
%axis([0 80 30 35]);
min(LambertV1EM);
max(LambertV1EM);
title('E2M')
hold off
%%
% Cont: Mars2Earth V1's, V2's
LambertV1ME = zeros(5,81);
LambertV2ME = zeros(5,81);


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
        [rE, vE, ~, ~, ~] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
        [rM, vM, ~, ~, ~] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
        rE = rE;
        rM = rM;
        [LambertV1,LambertV2,~] = LambertSolverND( rM, rE, TOF, muS, dir );
        %disp(LambertV1)
        LambertV1ME(i,j) = norm(LambertV1-vM)+Dv_Departure(LambertV1,vM,500,'m2e');
        LambertV2ME(i,j) = norm(LambertV2-vE)+Dv_Arrive(LambertV2,vE,350,'m2e');
        
        
    end
end
LambertME_Tot =LambertV1ME+LambertV2ME;
% disp(LambertV1EM)
% disp(LambertV2EM)
figure(2)
plot(xindx,LambertME_Tot)
hold on 
%plot(LambertV1ME(1,:)), hold on
%plot(LambertV1ME(2,:)), hold on
%plot(LambertV1ME(3,:)), hold on
%plot(LambertV1ME(4,:)), hold on
%plot(LambertV1ME(5,:)), hold on
legend(datestr(DatesME));
%axis([0 80 0 35]);
min(LambertV1ME);
max(LambertV1ME);
title('M2E')
hold off
%%

Tot = LambertEM_Tot+LambertME_Tot(1:5,:);

figure(3)
plot(xindx,Tot)
legend(datestr(DatesME))
title('TOT')
toc