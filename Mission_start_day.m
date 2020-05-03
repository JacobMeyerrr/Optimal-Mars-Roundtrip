function [E2M,M2E] = Mission_start_day(D_set,R_set,LB,UB,dT)
%
% Return E2M, M2E Mission starting day in formate of cell array with 
% time veriance spcificed
% 
% INPUT 
% "D_set" [n,1] datetime vector for E2M Hohmann transfer oppounity.
% "R_set" [n,1] datetime vector for M2E Hohmann transfer oppounity.
% "LB"          Lower bound for time veriance interested.(postive int)
% "UB"          Upper bound for time veriance interested.(postive int)
% "dT'          Change in time( day )
%
% OUTPUT 
% "E2M"   Earth to Mars Mission starting date(cell array formate) 
% "M2E"   Mars to Earth mission startign date(ce;; arrau formate)

%{

OutPut cell array will in formate { Hohmann transfer date(n,1), time
veriety coorsponding to the Hohmann transfer date at that row(n,m)}


%}

Dep = datetime(D_set);
Ret = datetime(R_set);

delta_T = -LB:dT:UB;
DT = length(delta_T);

E2M = {zeros(length(Dep),1),zeros(1,DT)};
M2E = {zeros(length(Ret),1),zeros(1,DT)};

for i = 1: length(Dep)
    E2M{i,1} = Dep(i);
    E2M{i,2}= (Dep(i)+days(delta_T));
end

for i =1: length(Ret)
    M2E{i,1} = Ret(i);
    M2E{i,2} =(Ret(i)+days(delta_T));
end

return 