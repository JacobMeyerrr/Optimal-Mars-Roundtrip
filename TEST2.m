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
% check for "ideal" Hohmann M->E return transfer dates within the next 15 years
[DatesEM,DatesME, TransferTime] = HohmannWindow(10,15); % E->M and M->E dates/TT

% % generate V's, times of flight for +-50 days from ideal Earth-to-Mars and
% % Mars-to-Earth transfers using a Lambert solver
% tofEM = []; % seconds (attempted transfer time from Earth to Mars
% tofME = [];% seconds (attempted transfer time from Mars to Earth)
% HohmannWindowEM = [];
% HohmannWindowME = [];
% for i=1:length(DatesEM)
%     HohmannWindowEM = [HohmannWindowEM;datetime(year(DatesEM(i)),month(DatesEM(i)),(day(DatesEM(i))+[-500:5:200]))];
%     tofEM = [tofEM;TransferTime+(24*60*60).*[500:-5:-200]]; % seconds (Maybe modify later!!)
% end
% for i=1:length(DatesME)
%     HohmannWindowME = [HohmannWindowME;datetime(year(DatesME(i)),month(DatesME(i)),(day(DatesME(i))+[-500:5:200]))];
%     tofME = [tofME;TransferTime+(24*60*60).*[500:-5:-200]]; % seconds (Check later!!!)
% end
% % <debugging only>
% %disp(HohmannWindowEM);
% %disp(HohmannWindowME);
% %disp(TransferTime);
% %size(tofEM);
% %size(tofME);
% %length(HohmannWindowME(:,1));
% % ^^debigging/testing only^^
% xindx = -500:5:200;
% %%
% % Cont: Earth2Mars V1's, V2's
% LambertV1EM = zeros(5,length(xindx));
% LambertV2EM = zeros(5,length(xindx));
% 
% 
% for i=1:length(HohmannWindowEM(:,1))
%     for j=1:length(HohmannWindowEM(1,:))
%         % get current date
%         DepYear = year(HohmannWindowEM(i,j)); % year of departure
%         DepMonth = month(HohmannWindowEM(i,j)); % month of departure
%         DepDay = day(HohmannWindowEM(i,j)); % day of departure
%         MarsArrivalDay = ceil(DepDay+TransferTime/(60*60*24));
%         arrive_at_mars = datetime(DepYear,DepMonth,MarsArrivalDay);
%         ArrYear = year(arrive_at_mars); % year of Mars arrival
%         ArrMonth = month(arrive_at_mars); % month of Mars arrival
%         ArrDay = day(arrive_at_mars); % day of Mars arrival
%         hour = 0;
%         minute = 0;
%         second = 0;
%         TOF = tofEM(i,j);
%         dir = 'pro';
%         % obtain heliocentric radii (3x1 state vector) for the specified date
%         % PlanetCoe = SolarSystem([thisyear thismonth thisday],0); % can't figure out 
% %         [coeE,rE, vE, jdE] = planet_elements_and_sv(3, DepYear, DepMonth, DepDay,0,0,0); % Question: is planetdata.m even reliable???
% %         [coeM,rM, vM, jdM] = planet_elements_and_sv(4, ArrYear, ArrMonth, ArrDay,0,0,0); % Question: is planetdata.m even reliable???
%         % obtain current true anomaly of planets ("current" = 12am, 5/1/2020)
%         [rE, vE, ~, ~ , ~] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
%         [rM, vM, ~, ~, ~] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
%         %rE = rE;
%         %rM = rM;
%         [LambertV1,LambertV2,~] = LambertSolverND( rE, rM, TOF, muS, dir );
%         %disp(LambertV1)
%         LambertV1EM(i,j) = norm(LambertV1-vE)+Dv_Departure(LambertV1,vE,350,'e2m');
%         LambertV2EM(i,j) = norm(LambertV2-vM)+Dv_Arrive(LambertV2,vM,500,'e2m');
%         
%         
%     end
% end
% LambertEM_Tot = LambertV1EM+LambertV2EM;
% % disp(LambertV1EM)
% % disp(LambertV2EM)
% figure(1)
% %xindx = -300:5:200;
% plot(xindx,LambertEM_Tot)
% %plot(xindx,LambertV1EM(1,:)), hold on
% %plot(xindx,LambertV1EM(2,:))
% %plot(xindx,LambertV1EM(3,:))
% %plot(xindx,LambertV1EM(4,:))
% %plot(xindx,LambertV1EM(5,:))
% legend(datestr(DatesEM));
% %axis([0 80 30 35]);
% min(LambertV1EM);
% max(LambertV1EM);
% title('E2M')
% hold off
% %%
% % Cont: Mars2Earth V1's, V2's
% LambertV1ME = zeros(5,length(xindx));
% LambertV2ME = zeros(5,length(xindx));
% 
% 
% for i=1:length(HohmannWindowME(:,1))
%     for j=1:length(HohmannWindowME(1,:))
%         % get current date
%         DepYear = year(HohmannWindowME(i,j)); % year of departure
%         DepMonth = month(HohmannWindowME(i,j)); % month of departure
%         DepDay = day(HohmannWindowME(i,j)); % day of departure
%         EarthArrivalDay = ceil(DepDay+TransferTime/(60*60*24));
%         arrive_at_earth = datetime(DepYear,DepMonth,EarthArrivalDay);
%         ArrYear = year(arrive_at_earth); % year of Earth arrival
%         ArrMonth = month(arrive_at_earth); % month of Earth arrival
%         ArrDay = day(arrive_at_earth); % day of Earth arrival
%         hour = 0;
%         minute = 0;
%         second = 0;
%         TOF = tofME(i,j);
%         dir = 'pro';
%         % obtain heliocentric radii (3x1 state vector) for the specified date
%         % PlanetCoe = SolarSystem([thisyear thismonth thisday],0); % can't figure out 
% %         [coeE,rM, vM, jdM] = planet_elements_and_sv(4, DepYear, DepMonth, DepDay,0,0,0); % Question: is planetdata.m even reliable???
% %         [coeE,rE, vE, jdE] = planet_elements_and_sv(3, ArrYear, ArrMonth, ArrDay,0,0,0); % Question: is planetdata.m even reliable???
%         [rE, vE, ~, ~, ~] = PlanetData(3, DepYear, DepMonth, DepDay,0,0,0);
%         [rM, vM, ~, ~, ~] = PlanetData(4, ArrYear, ArrMonth, ArrDay,0,0,0);
%         
%         [LambertV1,LambertV2,~] = LambertSolverND( rM, rE, TOF, muS, dir );
%         %disp(LambertV1)
%         LambertV1ME(i,j) = norm(LambertV1-vM)+Dv_Departure(LambertV1,vM,500,'m2e');
%         LambertV2ME(i,j) = norm(LambertV2-vE)+Dv_Arrive(LambertV2,vE,350,'m2e');
%         
%         
%     end
% end
% LambertME_Tot =LambertV1ME+LambertV2ME;
% % disp(LambertV1EM)
% % disp(LambertV2EM)
% figure(2)
% plot(xindx,LambertME_Tot)
% hold on 
% %plot(LambertV1ME(1,:)), hold on
% %plot(LambertV1ME(2,:)), hold on
% %plot(LambertV1ME(3,:)), hold on
% %plot(LambertV1ME(4,:)), hold on
% %plot(LambertV1ME(5,:)), hold on
% legend(datestr(DatesME));
% %axis([0 80 0 35]);
% min(LambertV1ME);
% max(LambertV1ME);
% title('M2E')
% hold off
% %%
% 
% Tot = LambertEM_Tot+LambertME_Tot(3:7,:);
% 
% figure(3)
% plot(xindx,Tot)
% legend(datestr(DatesME))
% title('TOT')
% 
% 
% %%
% muS   = 1.327e11; 
% rE    = 149.6e6; 
% rM    = 1.524*149.6e6;
% 
% [DE,DM] = HohmannWindow(10);
% %%
% Dep = datetime(DE);
% Ret = datetime(DM);
% %zeros(length(Dep),length(-20:20))
% %J0D = J0(year(Dep),month(Dep),day(Dep));
% E2M_start = {zeros(length(Dep),1),zeros(1,length(-20:20))};
% M2E_return ={zeros(length(Ret),1),zeros(1,length(-20:20))};
% %%
% for i = 1:length(Dep)
%     
%     E2M_start{i,1} = Dep(i);
%     E2M_start{i,2} = (Dep(i)+days(-20:20));
%         
% end
% 
% for i = 1:length(Ret)
%     
%     M2E_return{i,1} = Ret(i);
%     M2E_return{i,2} = (Ret(i)+days(-20:20));
%     
% end
% 
% 
% [E2M,M2E] = Mission_start_day(DE,DM,20,20,1);
% 
% %%
% muS   = 1.327e11; 
% rE    = 149.6e6; 
% rM    = 1.524*149.6e6;
% [DE,DM] = HohmannWindow(10);
% D1 = datevec(DE(2,:));
% D2 = datevec(DM(2,:));
% [psE, ~, ~, ~ ,~] = PlanetData(3,D1(1),D1(2),D1(3),D1(4),D1(5),round(D1(6)));
% [psM, ~, ~, ~, ~] = PlanetData(4,D1(1),D1(2),D1(3),D1(4),D1(5),round(D1(6)));
% thEi = atan2(psE(2),psE(1));
% thMi = atan2(psM(2),psM(1));
% [tW,tT,tSyn] = HohmannTransferAnimation( rE, rM, thEi, thMi, muS );
% 
% [psE2, ~, ~, ~ ,~] = PlanetData(3,D2(1),D2(2),D2(3),D2(4),D2(5),round(D2(6)));
% [psM2, ~, ~, ~, ~] = PlanetData(4,D2(1),D2(2),D2(3),D2(4),D2(5),round(D2(6)));
% thEi2 = atan2(psE2(2),psE2(1));
% thMi2 = atan2(psM2(2),psM2(1));
% [tW2,tT2,tSyn2] = HohmannTransferAnimation( rM, rE, thMi2, thEi2, muS );
% Lines 22-213 intentionally commented out!!!!!

%%
close all 
 
% S = datetime('now') + days(1:5);
% D = 1:10;
% TOF = 1:20;
S = DatesEM
D = [-200:20:200] % D represents the columns in the DV vector
TOF = TransferTime+(24*60*60).*[200:-20:-200] % TOF is the row in DV vector
DV_EM = DV_matrix(S,D,TOF);

% For Earth-to-Mars
for i=1:length(DV_EM) % iterate over all dates
    disp(i)
     for j=1:length(TOF)
        for k=1:length(D)
            Edep = DV_EM{i,1}+D(k);
            Marr = DV_EM{i,1}+TOF(j)/(24*3600);
            [rE, vE, ~, ~, ~] = PlanetData(3, year(Edep), month(Edep), day(Edep),hour(Edep),minute(Edep),second(Edep));
            [rM, vM, ~, ~, ~] = PlanetData(4, year(Marr), month(Marr), day(Marr), hour(Marr), minute(Marr), second(Marr));
            
            [LambertV1,LambertV2,~] = LambertSolverND( rE, rM, TOF(j), muS, 'pro' );

            LambertV1EM = norm(LambertV1-vE)+Dv_Departure(LambertV1,vE,350,'e2m');
            LambertV2EM = norm(LambertV2-vM)+Dv_Arrive(LambertV2,vM,500,'e2m');
            dv = LambertV1EM + LambertV2EM; %
            DV_EM{i,2}(j,k) = dv;
        end
     end
end
%dv = DV_matrix;

S = DatesME
D = [-200:20:200] % D represents the columns in the DV vector
TOF = TransferTime+(24*60*60).*[200:-20:-200] % TOF is the row in DV vector
DV_ME = DV_matrix(S,D,TOF);
%%
% For Mars-to-Earth 
for i=1:length(DV_ME) % iterate over all dates
    disp(i)
     for j=1:length(TOF)
        for k=1:length(D)
            Mdep = DV_ME{i,1}+D(k);
            Earr = DV_ME{i,1}+TOF(j)/(24*3600);
            [rM, vM, ~, ~, ~] = PlanetData(4, year(Mdep), month(Mdep), day(Mdep),hour(Mdep),minute(Mdep),second(Mdep));
            [rE, vE, ~, ~, ~] = PlanetData(3, year(Earr), month(Earr), day(Earr), hour(Earr), minute(Earr), second(Earr));
            
            [LambertV1,LambertV2,~] = LambertSolverND( rM, rE, TOF(j), muS, 'pro' );

            LambertV1ME = norm(LambertV1-vM)+Dv_Departure(LambertV1,vM,500,'m2e');
            LambertV2ME = norm(LambertV2-vE)+Dv_Arrive(LambertV2,vE,350,'m2e');
            dv = LambertV1ME + LambertV2ME; %
            DV_ME{i,2}(j,k) = dv;
        end
     end
end
%%
V_EM = Plot_DV(DV_EM,D,TOF,0,"e2m")

%%
V_ME = Plot_DV(DV_ME,D,TOF,5,"m2e")

%%
RoundtripDVs = ["Hello"];
for i=1:length(DatesEM)
    for j=1:length(DatesME)
        if(DatesEM(i)+TOF(V_EM{i,3}(1))/(24*3600)<DatesME(j))
             temp = num2str(V_EM{i,2}(1) + V_ME{i,2}(1));
             RoundtripDVs(i,j) = string(temp);
        else
            RoundtripDVs(i,j) = "Not Possible";
        end
    end
end
%%
%disp(RoundtripDVs);
RoundtripDVs_Size = size(RoundtripDVs);
RoundtripDVsTable = ["" "" "" "" "" "" "" "";"" "" "" "" "" "" "" "";...
                    "" "" "" "" "" "" "" "";"" "" "" "" "" "" "" "";...
                    "" "" "" "" "" "" "" ""];

RoundtripDVsTable(2:RoundtripDVs_Size(1)+1,2:RoundtripDVs_Size(2)+1) = RoundtripDVs;
DatesMEstr = [];

for i=1:length(DatesME)
    RoundtripDVsTable(1,i+1) = string(year(DatesME(i)))+"-"+string(month(DatesME(i)))+"-"+string(day(DatesME(i)));
end
DatesEMstr = [];

for i=1:length(DatesEM)
    RoundtripDVsTable(i+1,1) = string(year(DatesEM(i)))+"-"+string(month(DatesEM(i)))+"-"+string(day(DatesEM(i)));
end
RoundtripDVsTable(1,1) = "E->M/M->E Dates";

disp(RoundtripDVsTable)
disp(".......................DONE.....................")
        
%%
A = RoundtripDVs(:,1);
B = RoundtripDVs(:,2);
C = RoundtripDVs(:,3);
D = RoundtripDVs(:,4);
E = RoundtripDVs(:,5);
F = RoundtripDVs(:,6);
G = RoundtripDVs(:,7);

Table = timetable(A,B,C,D,E,F,G,'RowTimes',DatesEM)

toc