function  EM_speed_round(D,T,i)

% INPUT 
% "D" THE Date start the mission (datetime)
% "T" TOF per trip (days)
% "i" figure number 
muS = 1.3271e11;

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

TOF = T*86400;
[rE, vE, ~, ~, ~] = PlanetData(3, year(D), month(D), day(D),hour(D),minute(D),second(D));

Da = D+seconds(TOF);

[rM, vM, ~, ~, ~] = PlanetData(4, year(Da), month(Da), day(Da), hour(Da), minute(Da), second(Da));

DaE = Da +seconds(TOF);

[rEn,vEn,~,~,~] = PlanetData(3, year(DaE), month(DaE), day(DaE), hour(DaE), minute(DaE), second(DaE));

[LambertV1,~,~] = LambertSolverND( rE, rM, TOF, muS, 'pro' );

[LambertV2,~,~] = LambertSolverND( rM,rEn,TOF,muS,'pro');

[a,inc,W,w,e,thE,~,~,~] = OrbitalElementsFromRV(rE,vE,muS);
%coeE = [a,inc,W,w,e,th];
[th1,th2] = TrueAnomalyRange( e, 1e-1 );
th_all = linspace(th1,th2);
r = RVFromCOE( a,inc,W,w,e,th_all, muS );



% Plot Earth 
figure(i)
plot(r(1,:),r(2,:))
grid on ,axis equal, hold on 
plot(0,0,'r.','MarkerSize',40)
title('Fastest Lamber Earth Mars round trip ')


% Plot Marth 

[a,inc,W,w,e,thM,~,~,~] = OrbitalElementsFromRV(rM,vM,muS);
coeM = [a,inc,W,w,e,thM];
[th1,th2] = TrueAnomalyRange( e, 1e-1 );
th_all = linspace(th1,th2);
r = RVFromCOE( a,inc,W,w,e,th_all, muS );
plot(r(1,:),r(2,:))

% Plot Trojectory 

[a,inc,W,w,e,thL,~,~,~] = OrbitalElementsFromRV(rE,LambertV1,muS);
coeL = [a,inc,W,w,e,thL];

tx = linspace(0,TOF);
thx = TrueAnomFromTime(tx,a,e,muS,thL);

% compute the position and velocity from the COE at each true anomaly "thx"
rx = RVFromCOE( a,inc,W,w,e,thx, muS );

plot(rx(1,:),rx(2,:),'g','linewidth',3)
plot(rx(1,1),rx(2,2),'go','linewidth',2,'markersize',14)
plot(rx(1,end),rx(2,end),'rs','linewidth',2,'markersize',14)

[a,inc,W,w,e,thL,~,~,~] = OrbitalElementsFromRV(rM,LambertV2,muS);
coeL = [a,inc,W,w,e,thL];

tx = linspace(0,TOF);
thx = TrueAnomFromTime(tx,a,e,muS,thL);

% compute the position and velocity from the COE at each true anomaly "thx"
rx = RVFromCOE( a,inc,W,w,e,thx, muS );

plot(rx(1,:),rx(2,:),'b','linewidth',3)
plot(rx(1,1),rx(2,2),'go','linewidth',2,'markersize',14)
plot(rx(1,end),rx(2,end),'rs','linewidth',2,'markersize',14)


end