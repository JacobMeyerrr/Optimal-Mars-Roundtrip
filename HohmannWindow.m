function  [DatesEM,DatesME,TransferTime] = HohmannWindow(yearsEM,yearsME)
% Find HohmannTransfer oppounity in following years starting from 2020-5-1
%
% INPUT
% "year" Number of years into the future which you check for 
%                     Hohmann transfer opportunities
%
% OUTPUT
% "EMdate" [n,6] The Earth to Mars Hohmann Transfer oppounity  date array
%          Output Format [YYYY, MM, DD, HH, min, ss]
% "MEdate" [n,6] The Mars to Earth Hohmann Transfer oppounity  date array
%          Output Format [YYYY, MM, DD, HH, min, ss]
%
%

rE    = 149.6e6; 
rM    = 1.524*149.6e6;

muS   = 1.327e11; 
%year = 10;
% obtain position of planets ("current" = 12am, 5/1/2020)
[psE, ~, ~, ~ ,~] = PlanetData(3,2020,5,1,0,0,0);

[psM, ~, ~, ~, ~] = PlanetData(4,2020,5,1,0,0,0);

thEi = atan2(psE(2),psE(1)); 
thMi = atan2(psM(2),psM(1));

nE = sqrt(muS/rE^3);
nM = sqrt(muS/rM^3);
tT = pi/sqrt(muS)*((rE+rM)/2)^(3/2);
tSyn = 2*pi/abs(nM-nE); % 

% Di = datetime([2020 5 1]);% The mission starting date
% DF = datetime(Di + years(year),'InputFormat','yyyy-MM-dd'); % Duration end date

% lots of this stuff is simply reused from my final project -- yay!!
% define a simplified, circular orbit model of Earth-Mars-Sun system.
muS = 132712e6; % km^3/s^2 (gravitational parameter of Sun)
Re = 1.49597870e8; % km (average distance of Earth from sun)
aE = Re; % km (semi-major axis of Earth - circular orbit assumption)
Rm = 1.524*Re; % km (average distance of Mars from sol - assuming circular)
aM = Rm; % km (semi-major axis of Mars' orbit - circular orbit assumption)
vE = sqrt(muS/Re); % km/s (orbital velocity of Earth - circular assumption)
vM = sqrt(muS/Rm); % km/s (orbital velocity of Mars - circular assumption)
TE = 2*pi/sqrt(muS/aE^3); % seconds (period of circular Earth orbit)
TM = 2*pi/sqrt(muS/aM^3); % seconds (period of circular Mars' orbit)
nE = sqrt(muS/aE^3)*24*3600; % rad/day (mean orbit rate of Earth)
nM = sqrt(muS/aM^3)*24*3600; % rad/day (mean orbit rate of Mars)

% Earth-to-Mars transfer time in seconds
TransferTime = pi/sqrt(muS)*((rE+rM)/2)^(3/2)

% synotic period of Earth relative to Mars
Tsyn = 2*pi/(nE-nM);
n_syn = 2*pi/Tsyn; % rad/day

% Earth/Mars "phase angle" on 5/1/2020
phi = 57*pi/180;

% time of flight from Earth-to-Mars, vice-versa
t12 = ((pi/sqrt(muS))*((Re+Rm)/2)^(3/2))/(24*3600); % days

% Earth-to-Mars departure transfer angle (spacecraft leaves Earth)
phi0 = pi - nM*t12; % radians 

% Earth-to-Mars arrival transfer angle (spacecraft reaches Mars)
phif = pi - nE*t12; % radians 

% find next Earth-to-Mars opportunity
twaitEM = (phi-phi0)/n_syn;

% find next Mars-to-Earth opportunity
twaitME = (-2*phif-2*pi)/(nM-nE); % days until the next opportunity

% figure out what the ideal E-M transfer angle (difference of TA's) is
TransferAngleEM = pi - TransferTime*2*pi/TM; % radians (TAM - TAE)

% figure out what the ideal M-E transfer angle (difference of TA's) is
TransferAngleME = TransferTime/TE*2*pi-pi; % rad (TAM - TAE)

% Earth-to-Mars earth departure date
dateEMdeparture = datetime(2020,5,ceil(1+twaitEM));

% Mars-to-Earth Hohmann transfer Mars departure date 
dateMEdeparture = datetime(2020,5,ceil(1+twaitEM+t12+twaitME));

firstEMdate = dateEMdeparture
firstMEdate = dateMEdeparture
Tsyn = tSyn/(24*3600); % synotic period in days
disp(dateEMdeparture)
DatesEM = [dateEMdeparture]
DatesME = [dateMEdeparture]

COUNT = 1; % loop counter
while(DatesEM(COUNT)+Tsyn < datetime(2020+yearsEM,5,1))
    DatesEM(COUNT+1) = DatesEM(COUNT)+Tsyn
    COUNT = COUNT + 1;
end
COUNT = 1; % reset loop counter
while(DatesME(COUNT)+Tsyn < datetime(2020+yearsME,5,1))
    DatesME(COUNT+1) = DatesME(COUNT)+Tsyn
    COUNT = COUNT + 1;
end
return 