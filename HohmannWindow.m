function  [EMdate,MEdate] = HohmannWindow(year)
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

Di = datetime([2020 5 1]);% The mission starting date
DF = datetime(Di + years(year),'InputFormat','yyyy-MM-dd'); % Duration end date


tW = (thMi - thEi + tT*nM - pi)/(nE-nM);
tWE = (thEi-thMi+tT*nE-pi)/(nM-nE);

while( tW < 0 )
  tW = tW + tSyn;
end
while( tW > tSyn )
  tW = tW - tSyn;
end

while( tWE < 0 )
  tWE = tWE + tSyn;
end
while( tWE > tSyn )
  tWE = tWE - tSyn;
end

D1 =  Di + seconds(tW); % First Earth to Mars HohmannTransfer date 
D1M = Di+ seconds(tWE); % First Mars to Earth HohmannTransfer date

%[tW,tT,tSyn] = HohmannTransferAnimation( rE, rM, thEi, thMi, muS )
Dc = D1;
DcM = D1M;
n = 0;
j = 0;
% Find how may oppounity we have in time period E2M
while(Dc < DF)
    Dc = Dc + seconds(tSyn);
    n = n+1;
end
%Find how may oppounity we have in time period M2E
while(DcM < DF)
    DcM = DcM + seconds(tSyn);
    j = j+1;
end

%D2 = D1+66*seconds(tW);

TransferDateEM = zeros(n-1,6);
TransferDateME = zeros(j-1,6);


% Find all the oppounity date in the duration interested. 
for i = 1:n
    TransferDateEM(i,:) = datevec(D1+(i-1)*seconds(tSyn));
end

for i = 1:j
    TransferDateME(i,:) = datevec(D1M + (i-1)*seconds(tSyn));
end



EMdate = TransferDateEM;
MEdate = TransferDateME;

return 