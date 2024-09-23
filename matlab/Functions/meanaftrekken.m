function [B,Y]=meanaftrekken(A)
B=[];
ndx=size(A,2);
Y=mean(A,2);
for i=1:ndx
   B(:,i)=A(:,i)-Y;
end