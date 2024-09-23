

function [C]=slope(B)
ndx=size(B,2);

A =[];
for i=1:ndx-1
   X=B(:,(i+1))-B(:,i);
   A =[A X];
   clear X
end
C=[B A];