function DV = DV_matrix (S,D,TOF)

% Create a empety cell array for storing DV
% 
% INPUT 
% "S"     [i,1] datetime vector 
% "D"     [m,1] Transfer starting date vector 
% "TOF"   [1,n]  TOF vector 
%  
% OUTPUT 
% "DV"    a empet cell array to store DV 
%         DV{:,1} is the Hohmann transfer date
%         DV{:,2} is the starting date and transfer time variation
%         combination matrix for storing the dv value. [n,m]

% IMPORTANT NOTE 
% THIS function only creat EMPTY cell array for storing DV 
% The matrix is EMPTY!!!!!!




if(nargin<1 || nargin<2 || nargin<3)
    warning('NO input detected!')
    warning('date = 1:10')
    warning('TOF = 1:20')
    I = 5; J = 10; K = 20;
    for i = 1:I

        DV{i,1} = datetime([i*1111 i i i i i*10]);
        for j = 1:K
            for k = 1:J
                dv(j,k) = rand(1);
            end
        end
        DV{i,2} = dv;
    end


else
    
    I = length(S);
    J = length(D);
    K = length(TOF);
    %dv = zeros(K,J);
    for i = 1:I
        DV{i,1} = S(i);
        DV{i,2} = zeros(K,J);
    end
end
warning('This function return a EMPTY cell array to store the value!')
warning('The Cell array return from function "DV_matrix" is EMPTY')
warning('This matrix is EMPTY for now !')
return