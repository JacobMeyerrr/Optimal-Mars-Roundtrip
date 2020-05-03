function DV = DV_matrix (S,D,TOF)

% Create a empety cell array for storing DV
% 
% INPUT 
% "S"     [n,1] datetime vector 
% "D"     [n,1] Transfer starting date vector 
% "TOF"   [1,n]  TOF vector 
%  
% OUTPUT 
% "DV"    a empet cell array to store DV

% IMPORTANT NOTE 
% THIS function only creat EMPTY cell array for storing DV 
% The matrix is EMPTY!!!!!!




if(nargin<1)
    warning('NO input detected!')
    I = 5; J = 5; K = 5;
    for i = 1:I

        DV{i,1} = datetime([1111 11 1 1 1 1]);
        for j = 1:J
             dv(j,1:I) = j;
        end
        DV{i,2} = dv;
    end


else
    
    I = length(S);
    J = length(D);
    K = length(TOF);
    dv = zeros(K,J);
    for i = 1:I
        DV{i,1} = S(i);
        DV{i,2} = zeros(K,J);
    end
end
warning('This function return a EMPTY cell array to store the value!')
warning('The Cell array return from function "DV_matrix" is EMPTY')
warning('This matrix is EMPTY for now !')
return