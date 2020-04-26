function [DatesEM, DatesME, TransferTime] = opportunities(Years)
%
% Function to obtain a list of all possible Earth-to-Mars and Mars-to-Earth 
% Hohmann transfer windows within the number of years specified from today.
%
% For simplicity, this function assumes circular orbits (e = 0) for both 
% Earth and Mars. We make the proper corrections later on in our script 
% with our Lambert solver and real orbit data obtained from NASA JPL.
%
% Synopsis:
%     [DatesEM, DatesME, TransferTime] = opportunities(Years)
%
% Input:
%     Years        =  Number of years into the future which you check for 
%                     Hohmann transfer opportunities
%
% Output:
%     DatesEM      = "Ideal" opportunity Earth-to-Mars Hohmann transfer window
%                     dates within the amount time (years) specified
%     DatesME      = "Ideal" opportunity Mars-to-Earth window dates within
%                     the time specified
%     TransferTime = Hohmann transfer half-orbit time (circular assumption)
%
%
% By: Jacob J. Meyer -- Apr. 2020

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

% obtain current true anomaly of planets ("current" = 12am, 5/1/2020)
[r, v, jd, coem coe2] = PlanetData(3,2020,5,1,0,0,0);
TAE = mod(coem(6),2*pi); % rad (starting true anomaly of Earth - Heliocentric Frame)
[r, v, jd, coem coe2] = PlanetData(4,2020,5,1,0,0,0);
TAM = mod(coem(6),2*pi)-pi*1.5/3; % rad (starting true anomaly of Mars - Heliocentric Frame)

% obtain current true anomaly of planets ("current" = 12am, 5/1/2020)
% to switch to the actual current date, don't input anything to SolarSys.m
% PlanetCoe = SolarSystem([2020 5 1],0); % obtain heliocentric COEs (all planets)
% TAE = PlanetCoe(3).trueAnom; % rad (starting true anomaly of Earth) - Heliocentric
% TAM = PlanetCoe(4).trueAnom; % rad (starting true anomaly of Mars) - Heliocentric

% obtain period of Hohmann transfer orbit from Earth-to-Mars (vice versa)
aT = (Re+Rm)/2; % km (semi-major axis of heliocentric transfer orbit)
Tt = 2*pi/sqrt(muS/aT^3); % seconds (period of Hohmann transfer orbit)
TransferTime = Tt/2; % seconds (actual transfer time from E-M, M-E, etc.)

% figure out what the ideal E-M transfer angle (difference of TA's) is
TransferAngleEM = pi - TransferTime*2*pi/TM; % radians (TAM - TAE)

% figure out what the ideal M-E transfer angle (difference of TA's) is
TransferAngleME = TransferTime/TE*2*pi-pi; % rad (TAM - TAE)

% ^kinda based off of the "symnonic period" stuff that we learned in class^

% define our range of potential days (and dates) for our E-M transfer
% e.g. day 1 = 5/1/2020, day 2 = 5/2/2020, last day = 5/1/(2020+Years)
allDates = datetime(2020,5,1:ceil(365.25*Years)); % Start on May 1st, 2020
Day = [1:1:ceil(365.25*Years)]; % correspond to May 1st, 2020 to allDates(length(allDates))

% obtain true anomalies corresponding to ~12am on each day
EarthTrueAnoms = mod((TAE + (2*pi).*(24*3600).*(Day-1)./TE),(2*pi));
MarsTrueAnoms = mod((TAM + (2*pi).*(24*3600).*(Day-1)./TM),(2*pi));

% obtain all symtotic angles for E-M transfer
SymEM = [];
for i=1:length(allDates)
    Earth = EarthTrueAnoms(i);
    Mars = MarsTrueAnoms(i);
    if((0<=Earth && Earth<=pi) && (0<=Mars && Mars<=pi))
        SymEM(i) = MarsTrueAnoms(i)-EarthTrueAnoms(i);
    elseif((pi<Earth && Earth<2*pi) && (pi<Mars && Mars<2*pi))
        SymEM(i) = MarsTrueAnoms(i)-EarthTrueAnoms(i);
    elseif((0<=Earth && Earth<=pi) && (pi<Mars && Mars<2*pi))
        SymEM(i) = MarsTrueAnoms(i)-EarthTrueAnoms(i);
    elseif((pi<Earth && Earth<2*pi) && (0<=Mars && Mars<=pi))
        SymEM(i) = (MarsTrueAnoms(i)+2*pi)-EarthTrueAnoms(i);
    end
end
SymEM;

% figure out best transfer times in the next Y number of years specified
DatesEM = allDates(find(abs(SymEM./TransferAngleEM-1)<0.01));
DatesME = allDates(find(abs(SymEM./TransferAngleME-1)<0.01));
BestEMopps = [];
BestMEopps = [];
prev_year = year(DatesEM(1));
j = 1;
for i=1:length(DatesEM)
    current = abs(SymEM(find(allDates==DatesEM(i)))/TransferAngleEM-1);
    best = abs(SymEM(find(allDates==DatesEM(j)))/TransferAngleEM-1);
    if(year(DatesEM(i)) == prev_year)
        if(current<best)
            j = i;
        end
        if(i==length(DatesEM))
            BestEMopps = [BestEMopps DatesEM(j)];
        end            
    else
        BestEMopps = [BestEMopps DatesEM(j)];        
        if(i==length(DatesEM))
            j = i;
            BestEMopps = [BestEMopps DatesEM(j)]; 
        end
        j = i;
        prev_year = year(DatesEM(i));
    end
end
DatesEM = BestEMopps;  % complete extraction of each single best yearly date
prev_year = year(DatesME(1));
j = 1;
for i=1:length(DatesME)
    current = abs(SymEM(find(allDates==DatesME(i)))/TransferAngleME-1);
    best = abs(SymEM(find(allDates==DatesME(j)))/TransferAngleME-1);
    if(year(DatesME(i)) == prev_year)
        if(current<best)
            j = i;
        end
        if(i==length(DatesME))
            BestMEopps = [BestMEopps DatesME(j)];
        end            
    else
        BestMEopps = [BestMEopps DatesME(j)];        
        if(i==length(DatesME))
            j = i;
            BestMEopps = [BestMEopps DatesME(j)]; 
        end
        j = i;
        prev_year = year(DatesME(i));
    end
end
DatesME = BestMEopps; % complete extraction of each single best yearly date

% debugging only
% plot(SymEM-TransferAngleEM), hold on
% plot(EarthTrueAnoms), hold on
% plot(MarsTrueAnoms), hold on
% legend("SymEM","Earth","Mars")
% tic toc is an interesting MATLAB feature...
%% Sources:
% -----------------------------------------------------------------------
% General constants, formulas, etc.
    % SolarSystem.m (Joseph Mueller, NASA "X-files" contributor)
    % Various other NASA websites
    % Orbital Mechanics For Engineering Students (Howard D. Curtiss)
% Check to see if my Hohmann windows are roughly correct
    % https://spectrum.ieee.org/tech-talk/aerospace/robotic-exploration/china-says-its-mars-landing-technology-is-ready-for-2020
% Why choose a standard Hohmann transfer and not a bi-elliptic transfer?
    % https://space.stackexchange.com/questions/14618/how-does-bi-elliptic-transfer-compare-to-low-energy-transfers-interplanetary-t