function [A,MURES,SIGMRES] = normalise(A,FLAG1,FLAG2,MU,SIGM)

% [A,MURES,SIGMRES] = normalise(A,FLAG1,FLAG2,MU,SIGM)
% This function normalises the expression levels in A (each row is normalised separately)
% Each expression level X on row j is replaced by:
% (X-(FLAG1>0)*MURES(j))/((FLAG2>0)*SIGMRES(j)+(FLAG2==0))
% MURES/SIGMRES is set equal to MU/SIGM if FLAG1/FLAG2 is different from 1
% MURES/SIGMRES is re-calculated if FLAG1/FLAG2 is equal to 1 (column vectors,
% containing the means/std's for each row)
       
Samples = size(A,2);

if FLAG1==1
   MU=mean(A')';
end

if FLAG2==1
   SIGM = std(A')';
   SIGM = SIGM + (SIGM==0);
end

if FLAG1>0
   A=A-(MU*ones(1,Samples));
end

if FLAG2>0
   A=A./(SIGM*ones(1,Samples));
end

MURES=MU;
SIGMRES=SIGM;



