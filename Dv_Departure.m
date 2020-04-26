function [dvR,Beta,DELTA] =  Dv_Departure(vT,h,dir)    
% 
% Calculate the delta v require to departure from a planet (hyperbola
% terojectory) 
%
% IMPORTANT
% All input and output are not vector form
%
% INPUT
% "vT,   The velocity aimed to reach (KM/s)
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

        Vc = sqrt(mu_p/R);              % parking orbit velocity(Assuming they are all circular)
        vinf = vT;                      % hyperbolic excess velocity
        h = R*sqrt(vinf^2+2*mu_p/R);    % angular momentum of hyperbola
        e = 1+R*vinf^2/mu_p;            % eccentricity of hyperbola
        vP = h/R;                       % velocity of hyperobola at perigee 
        dvR = vP - Vc;                  % required delta-v from parking orbit
        Beta = acos(1/e);               % angle between apse line and asymptote line for hyperbola
        sma = h^2/mu_p/(1-e^2);         % Smei-major axis 
        DELTA = (R+abs(sma))*sin(Beta); % aiming radius


        case{'mars2earth','m2e'}

         % Departure from Mars. From Planetary Rendezous. The capature orbit
         % is not circular.
         mu_p = muM;
         R = radM + h;
         vinf = vT;
         H = R*sqrt(vinf^2+2*mu_p/R);    % angular momentum of hyperbola
         e = 1+R*vinf^2/mu_p;            % eccentricity of hyperbola
         vPc = sqrt(muM*(1+e)/R);         % velocity of capature orbit at perigee 
         vPh = sqrt(vinf^2+2*muM0R);
         dvR = vPh - vPc;                % required delta-v from capature orbit
         Beta = acos(1/e);               % angle between apse line and asymptote line for hyperbola  
         sma = H^2/mu_p/(1-e^2);         % Smei-major axis 
         DELTA = (R+abs(sma))*sin(Beta); % aiming radius

     % PS: THe capature orbit is a ellipical orbit. However, We don have
     % enought information to calculate the orbit element for capature orbit.
     % But we know that at perigee v_capiture and v_hyperbola is the same so
     % delta_v require is simply vinf-vP 

    end
    
    
return 