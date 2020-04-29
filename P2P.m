function d = P2P

%% Start date
startYMD=[2020 5 1];
p1 = 3;
p2 = 4;
year = startYMD(1);
month = startYMD(2);
day = startYMD(3);
jD0 = J0(year,month,day);   % Julian date for our earliest start date

%% Julian day to calendar date function...
JDay2Date = @(x) datevec(x-1721058.5);

%% Methodology
%
% We will try to find a good start date by first planning a Hohmann transfer. 
% This approach will assume circular orbits for Earth and Saturn. 
% 
% We will then consider different start dates and different durations, using
% the ACTUAL eccentric orbits for both planets.
% 
% So... we will have 2 independent variables:
% 
%     * jD1: Start Date
%     * tT:  Transfer time (duration)

mu = 1.327e11; 
% First get the radius and initial phase angle for both planets...
[pos1,~,~,coe1] = PlanetData(p1, year, month, day, 0, 0, 0, mu);
[pos2,~,~,coe2] = PlanetData(p2, year, month, day, 0, 0, 0, mu);
% use the semi-major axis for each planet's circular orbit radius
r1 = coe1(1); 
r2 = coe2(1);
% look at the position in the ecliptic plane to compute the initial phase
% angle
th10 = atan2(pos1(2),pos1(1));
th20 = atan2(pos2(2),pos2(1));
% Compute the delta-vs
hData = OrbMvrHohmann(r1,r2,r1,r2,mu);
% Compute the wait time, transfer time, and synodic period
[tW,tT,tSyn] = HohmannTransferAnimation( r1, r2, th10, th20, mu );

% If we had circular orbits, we would want to do a Hohmann transfer
% starting at day: jD0 + tW/86400
jDStartH = jD0 + tW/86400;

% But we have slightly eccentric orbits, so let's consider a range of start
% dates NEAR this date, and a range of transfer times NEAR the Hohmann
% transfer time...
jD1Vec = jDStartH + (-200 : 10 : 200); % e.g. from 50 days before to 50 days after
tTVec  = tT + (-200 : 10 : 200)*86400; % e.g. from 50 days shorter to 50 days longer

% make big empty matrices for first/second delta-v's
nST = length(jD1Vec); % number of different start times
nTT = length(tTVec); % number of different transfer times
dV1 = zeros(nST,nTT);
dV2 = zeros(nST,nTT);
dV3 = zeros(nST,nTT);
dV4 = zeros(nST,nTT);
% For each combination of start date and transfer time ... use Lambert to find
% the transfer orbit.
for i=1:length(jD1Vec)
  
  % start date
  date1 = JDay2Date(jD1Vec(i));
  y1 = date1(1);
  m1 = date1(2);
  d1 = date1(3);
  
  for j=1:length(tTVec)
    % position of Earth at start date
    [r1,v1] = PlanetData(p1,y1,m1,d1,0,0,0,mu);
    
    % TOF, end date
    TOF = tTVec(j);
    date2 = JDay2Date( jD1Vec(i) + tTVec(j)/86400 );
    y2 = date2(1);
    m2 = date2(2);
    d2 = date2(3);

    % position of Saturn at start date
    [r2,v2] = PlanetData(p2,y2,m2,d2,0,0,0,mu);
    
    % compute transfer maneuver
    [v1T,v2T] = LambertSolverND(r1,r2,TOF,mu,'pro');
    
    
    % compute and record delta-v's
   % DD = Dv_Departure(v1T,v1,350,'e2m')
   % DA = Dv_Arrive(v2T,v2,500,'e2m')
   
   date3 = JDay2Date(date2+tTVec(j)/86400);
    y3 = date3(1);
    m3 = date3(2);
    d3 = date3(3);
    
   [r1R,v1R] = PlanetData(p1,y3,m3,d3,0,0,0,mu);
   [v2TR,v1TR] = LambertSolverND(r2,r1R,TOF,mu,'pro');
   D1 = Dv_Departure(v1T,v1,350,'e2m');
   A1 = Dv_Arrive(v2T,v2,500,'e2m');
   D2 = Dv_Departure(v2TR,v2,500,'m2e');
   A2 = Dv_Arrive(v1TR,v1,350,'m2e');
    dV1(i,j) = norm(v1T-v1)+D1;
    dV2(i,j) = norm(v2T-v2)+A1;
    dV3(i,j) = norm(v2TR-v2)+D2;
    dV4(i,j) = norm(v1TR-v1R)+A2;
  end
  
end
dVTot = dV1+dV2+dV3+dV4;

% Plot the results of the Direct planet to planet analysis...
figure, plot(jD1Vec-jDStartH,dVTot), grid on, xlabel('Start day variation'), ylabel('DV')
figure, surf(jD1Vec-jDStartH,(tTVec-tT)/86400,dVTot','edgecolor','none'), grid on, rotate3d on
xlabel('Start day variation (days)'), ylabel('Transfer time variation (days)'), zlabel('DV')
hold on, plot3(0,0,hData.dVTot,'r.','markersize',30)
shading interp
[v,row]=min(dVTot); [v,col]=min(v); row=row(col);
plot3(jD1Vec(row)-jDStartH,(tTVec(col)-tT)/86400,dVTot(row,col),'g.','markersize',30)
legend('Lambert Solver DVs','Hohmann Solution','Min DV from Lambert','location','northeast')

d.dVTot_min = min(dVTot(:));
d.dV2_min = min(dV2(:));
d.dV1_min = min(dV1(:));

d.days_minDV  = tTVec(col)/86400;

d.jD1 = jD1Vec(row);
d.jD2 = d.jD1 + d.days_minDV;