function [B,Y]=maxdelen(A)
B=[];
ndx=size(A,2);
[Y]=max(A,[],2);
for i=1:ndx
   B(:,i)=A(:,i)./Y;
   end