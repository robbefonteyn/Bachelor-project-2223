function [P,U,MU,S] = pcagene(A)

Genes = size(A,1);

MU = mean(A);
A=A-(ones(Genes,1)*MU);


[U1,S1,V] = svd(A,0);

S1=diag(S1);
S=S1.^2;

U=V';

P=A*V;