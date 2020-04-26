function [dV] = DVfromSOI(planet,r,v)

% Calculates the change in velocity of a satellite due to the sphere of
% influence based on its entry position to its parking orbit. 
%
% Input data:
%
% Planet data
% 3 - Earth
% 4 - Mars
%
% r = radius
% v = velocity

muS = 1.327e11;

if planet = 3
    rSOI = 925000;
    mu = 398600;
    R = 149.6e6;
    po = 350; %parking orbit
elseif planet = 4
    rSOI = 577000;
    mu = 52484;
    R = 277.9e6;
    po = 500;
else
    disp('Warning: function only built for Earth(3) and Mars(4)')
    break
end

a = -mu/norm(r)^3*r - muS/R^3*r;

dx = rSOI - po;

v2 = sqrt(v^2+2*a*dx);

dV = v2-v;

