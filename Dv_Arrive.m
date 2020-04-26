function [dvR,Beta,DELTA] = Dv_Arrive(vA,vP,h,dir)
%
% Calculate the delta v require to end up the parking orbit ( Planetary
% Rendezous) 
%
% IMPORTANT
% All input and output are not vector form
%
% INPUT
% "vA"   Helicentric approach velocity (KM/s)
% "vP"   The planet velocity (KM/s)  
% "h"    The (minmmun) altitude of parking orbit (KM) 
% "mu_p" The Gravtional constnat of the planet arrival to (KM^3/s^2)
% "dir"  Direction of the planetary transfer
%      'earth2mars' or 'e2m' Earth to Mars transfer (inner to outer)
%      'mars2earth' or 'm2e' Mars to Earth transfer (outer ot inner)
%
% OUTPUT
% "dvR"   Required delta v to arrivel 
% "Beta"  The angle between apse line and asymptote line for hyperbola
% "DELTA" Aiming ardius 'жд'
%
% NOTE : This function only returen the magnitude of deltia V relative to
% the Planet. NOT HELICENTRIC !!! 
%% Constants 

muE = 398600.44;
muM = 42828;
radE  = 6378.14;
radM  = 3389.90;
mu = 0;
Rp = 0;
%% Input Check

if(~ischar(dir))
    error('Must indicate the planetary transfer direction "e2m" or "m2e"')
end

switch lower(dir)
    case{'earth2mars','e2m'}
        vinf = norm(vP-vA);
        mu = muM;
        Rp = radM;
    case{'mars2earth','m2e'}
        vinf = norm(vA - vP);
        mu = muE;
        Rp = radE;
end

if(mu ==0 || Rp == 0)
    error('Something when wrong in function "Dv_Arrive" Check the input!')
end

%%
rp = Rp + h;

e = 1+(rp*vinf^2)/mu;

H = rp*sqrt(vinf^2+2*mu/rp);

DELTA = rp*sqrt(1+2*mu/(rp*vinf^2)); % aiming radius (eqn 8.57)

vp_hyp = sqrt(vinf^2 + 2*mu/rp); % eqn 8.58

%vp_capture = sqrt(mu*(1+e)/rp); % eqn 8.59

vp_capture = sqrt(mu/rp); % ciculor capature orbit assumption  

dvR = vp_hyp-vp_capture; % eqn 8.60

Beta = acos(1/e);




return