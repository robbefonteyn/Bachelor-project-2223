function [P,U,MU,S] = pcasam(A)

% [P,U,MU,S] = pcasam(A)
% This function performs Principal Component Analysis for the column vectors in A 
% (No preprocessing, except substraction of the mean column vector, is done)
% The output is given by:
% - P contains the new coordinates (projection onto principal components, after substraction 
%   of the mean column vector) for the column vectors in A
% - U contains the principal components (Number of PC = Number of samples - 1) in the original coordinates
% - MU contains the mean of the samples in A in the original coordinates
% - S contains the eigenvalues associated with the principal components (column vector)

Genes = size(A,1);
Samples = size(A,2);

MU = mean(A')';

A=A-(MU*ones(1,Samples));

[U1,S1,V] = svd(A,0);

if Genes >= Samples
	U=U1(:,1:(Samples-1));
	S1=diag(S1);
   S=S1(1:Samples-1).^2;
else
   U=U1;
   S1=diag(S1);
   S=S1.^2;
end

P=U'*A;
