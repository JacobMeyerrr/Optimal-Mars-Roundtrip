function [dvR,Beta,DELTA] =  Dv_Departure(vT,vP,h,dir)    
% 
% Calculate the delta v require to departure from a planet (hyperbola
% terojectory) 
%
% IMPORTANT
% All input and output are not vector form
%
% INPUT
% "vT,   The velocity aimed to reach (KM/s)[3,1]
% 'vP'   The velocity of the Planet (km/s) [3,1]
% "h"    The altitude of parking orbit (KM) 
% "dir"  Direction of the planetary transfer
%      'earth2mars' or 'e2m' Earth to Mars transfer (inner to outer)
%      'mars2earth' or 'm2e' Mars to Earth transfer (outer ot inner)
%
% OUTPUT
% "dvR"   Required delta v to departure 
% "Beta"  The angle between apse line and asymptote line for hyperbola
% "DELTA" Aiming ardius 'жд'
%
% NOTE : This function only returen the magnitude of deltia V relative to
% the Planet. NOT HELICENTRIC !!! 
% 

%% Input Check 
if(~ischar(dir))
    error('Must indicate the planetary transfer direction "e2m" or "m2e"')
end

%% Constants
muE = 398600.44;
muM = 42828;
radE  = 6378.14;
radM  = 3389.90;
%% 
    switch lower(dir)
        case{'earth2mars','e2m'}

        % departure fomr earth. circular parking orbit
        mu_p = muE; 
        R = radE+h;
        vinf = norm(vT) - norm(vP);
        case{'mars2earth','m2e'}
        mu_p = muM;
        R = radM+h;
        vinf = norm(vP) - norm(vT);
    end
    
    
    
        Vc = sqrt(mu_p/R);              % parking orbit velocity(Assuming they are all circular)
        %vinf = norm(vT);                % hyperbolic excess velocity
        h = R*sqrt(vinf^2+2*mu_p/R);    % angular momentum of hyperbola
        e = 1+R*vinf^2/mu_p;            % eccentricity of hyperbola
        vP = h/R;                       % velocity of hyperobola at perigee 
        dvR = abs(vP - Vc);             % required delta-v from parking orbit
        Beta = acos(1/e);               % angle between apse line and asymptote line for hyperbola
        sma = h^2/mu_p/(1-e^2);         % Smei-major axis 
        DELTA = (R+abs(sma))*sin(Beta); % aiming radius


        

         
    
    
return 