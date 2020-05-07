function  EM_speed_round(D_EM,D_ME,TofEM,TofME,i,LambertV1_EM,LambertV1_ME,titl)

% INPUT 
% "D" THE Date start the mission (datetime)
% "T" TOF per trip (days)
% "i" figure number 
muS = 1.3271e11;
AU = 1.496e8; % 1 Astronomical Unit in km
% D = DatesEM(1);
% R = DatesME(1);
% 
% Dvec = (D+days(-200:20:200));
% Rvec = (R+days(-200:20:200));
% 
% n = length(Dvec);
% m = length(Rvec);
% 
% Dmin = Dvec(end);
% Rmin = Rvec(end);
% Diff = abs(Rvec(end)-Dmin(1));
% %%
% for i = 1:n
%     for j = 1:m
%         X(i,j) = 1;
%        DiffT = abs(Rvec(i) - Dvec(j));
%        if(DiffT<Diff)
%            Dmin = Dvec(i);
%            Rmin = Rvec(j);
%            Diff = DiffT;
%        end
%         
%     end
% end
% Dmin
% Rmin
%  TOF = seconds(Diff);    
if(nargin<3||nargin<2||nargin<1)
i = 100;
D = datetime('now');
TOF = 50*86400;
end

TOF = TofEM*86400;
[rE, vE, ~, ~, ~] = PlanetData(3, year(D_EM), month(D_EM), day(D_EM),hour(D_EM),minute(D_EM),second(D_EM));


[rM, vM, ~, ~, ~] = PlanetData(4, year(D_ME), month(D_ME), day(D_ME), hour(D_ME), minute(D_ME), second(D_ME));

DaE = D_ME + TofME;

[rEn,vEn,~,~,~] = PlanetData(3, year(DaE), month(DaE), day(DaE), hour(DaE), minute(DaE), second(DaE));

% [LambertV1,~,~] = LambertSolverND( rE, rM, TOF, muS, 'pro' );
% 
% [LambertV2,~,~] = LambertSolverND( rM,rEn,TOF,muS,'pro');

[a,inc,W,w,e,thE,~,~,~] = OrbitalElementsFromRV(rE,vE,muS);
%coeE = [a,inc,W,w,e,th];
[th1,th2] = TrueAnomalyRange( e, 1e-1 );
th_all = linspace(th1,th2);
r = RVFromCOE( a,inc,W,w,e,th_all, muS )./AU;



% Plot Earth 
figure(i)
plot(r(1,:),r(2,:))
grid on ,axis equal, hold on 
plot(0,0,'r.','MarkerSize',40)
title(titl)
xlabel("X-Axis Distance From Sun (Astronomical Units)")
ylabel("Y-Axis Distance From Sun (Astronomical Units)")




% Plot Mars

[a,inc,W,w,e,thM,~,~,~] = OrbitalElementsFromRV(rM,vM,muS);
coeM = [a,inc,W,w,e,thM];
[th1,th2] = TrueAnomalyRange( e, 1e-1 );
th_all = linspace(th1,th2);
r = RVFromCOE( a,inc,W,w,e,th_all, muS )./AU;
plot(r(1,:),r(2,:))

% Plot Trojectory 

[a,inc,W,w,e,thL,~,~,~] = OrbitalElementsFromRV(rE,LambertV1_EM,muS);
coeL = [a,inc,W,w,e,thL];

tx = linspace(0,TOF);
thx = TrueAnomFromTime(tx,a,e,muS,thL);

% compute the position and velocity from the COE at each true anomaly "thx"
rx = RVFromCOE( a,inc,W,w,e,thx, muS )./AU;

plot(rx(1,:),rx(2,:),'g','linewidth',3)
plot(rx(1,1),rx(2,2),'go','linewidth',2,'markersize',14)
plot(rx(1,end),rx(2,end),'rs','linewidth',2,'markersize',14)

TOF = TofME*86400;
[a,inc,W,w,e,thL,~,~,~] = OrbitalElementsFromRV(rM,LambertV1_ME,muS);
coeL = [a,inc,W,w,e,thL];

tx = linspace(0,TOF);
thx = TrueAnomFromTime(tx,a,e,muS,thL);

% compute the position and velocity from the COE at each true anomaly "thx"
rx = RVFromCOE( a,inc,W,w,e,thx, muS )./AU;

plot(rx(1,:),rx(2,:),'b','linewidth',3)
plot(rx(1,1),rx(2,2),'go','linewidth',2,'markersize',14)
plot(rx(1,end),rx(2,end),'rs','linewidth',2,'markersize',14)

end