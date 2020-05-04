function DV = Plot_DV(D,date,TOF)

%DV_plot plot DV 
% Plot the surface plot for eace DV 
% 
% INPUT 
% "D"     A {k,2} cell arrray with first cell array column in formate of 
%         datetime Second cell colums in formate of [n,m] double matrixs
%         coorspinding to the datetime in cell array column 1.
% "date"  A [m,1] Start day variation(days) vector 
% "TOF"   A [n,1] Transfer time variation(days) vector 
%
% OUTPUT 
% "DV"    A {k,3} cell array with first cell arry column in formate of
%         datetime. second cell arry column return the min DV coorsponding
%         to the datetime in colume 1. Thired cell array return the
%         row,colum index number coorsponding to the min DV 
%         IN format [row,col]



%D = dv;

% Preallocating the array 
[i,j] = size(D);
string = strings(i,1);
 
for k = 1:i
    
   string(k) = D{k,1};
     
end


% Ploting the surface plot for each DV 
for k = 1:i
    figure(k)
    d_v = D{k,2};
    %[n,m] = size(d_v);
    xind = date;
    yind = TOF ;
    surf(xind,yind,d_v,'edgecolor','none');
    grid on,rotate3d on, xlabel('n'),ylabel('m'),colorbar
    hold on
    shading interp
    title(string(k));
    [v,row]=min(d_v); [v,col]=min(v); row=row(col);
    plot3(xind(col),yind(row),d_v(row,col),'g.','markersize',30)
    legend('Total delta v ','Min DV','location','northeast')
    hold off
    DV{k,1} = string(k);
    DV{k,2} = d_v(row,col);
    DV{k,3} = [row,col];
end

return