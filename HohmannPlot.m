function[f] = HohmannPlot(num)

rE    = 149.6e6; 
rM    = 1.524*149.6e6;
muE   = 398600.44;
muM   = 42828;
muS   = 1.327e11; 
radE  = 6378.14;
radM  = 3389.90;
AU = 1.496e8; % 1 Astronomical Unit in km

thx = linspace(0,2*pi,300);
figure(num)
plot(0,0,'ro','markersize',10,'MarkerFaceColor',[1,0,0]), hold on, axis equal, grid on
plot(rE/AU,0,'bo','markersize',10,'MarkerFaceColor',[0,0,1])
plot(-rM/AU,0,'yo','markersize',10,'MarkerFaceColor',[1,1,0])
plot(rE*cos(thx)/AU,rE*sin(thx)/AU,'b');
plot(rM*cos(thx)/AU,rM*sin(thx)/AU,'b','LineStyle','- -');
title("Trajectory of Earth-to-Mars Hohmann Transfer, Circular Orbits")
xlabel("X-Axis Distance From Sun (Astronomical Units)")
ylabel("Y-Axis Distance From Sun (Astronomical Units)")


aT = 0.5*(rE+rM);
eT = abs(rE-rM)/(rE+rM);
nT = sqrt(muS/aT^3);
TT = 2*pi/nT;
t = 0:100:TT;
thT = TrueAnomFromTime(t,aT,eT,muS,0);
%M = [cos(th1s), -sin(th1s), 0; sin(th1s), cos(th1s), 0; 0 0 1];
[rT,xT,yT] = PerifocalOrbit( aT, eT, thT );
rrT = [xT; yT; 0*rT];
plot(rrT(1,:)/AU,rrT(2,:)/AU,'m')
legend('Sun','Earth','Mars','Earth orbit','Mars orbit','Transfer orbit')