function [A,MURES,SIGMRES] = normalise_misval(A,FLAG1,FLAG2,MU,SIGM)

% [A,MURES,SIGMRES] = normalise_misval(A,FLAG1,FLAG2,MU,SIGM)
% This function normalises the expression levels in A (each row is normalised separately)
% Each expression level X on row j is replaced by:
% (X-(FLAG1>0)*MURES(j))/((FLAG2>0)*SIGMRES(j)+(FLAG2==0))
% MURES/SIGMRES is set equal to MU/SIGM if FLAG1/FLAG2 is different from 1
% MURES/SIGMRES is re-calculated if FLAG1/FLAG2 is equal to 1 (column vectors,
% containing the means/std's for each row)
       
Samples = size(A,2);
nr_mis_val = sum(sum(~finite(A)));

if FLAG1==1
   if nr_mis_val==0
      MU=mean(A')';
   else      
      Q=sum(finite(A),2);
      Q(find(Q==0))=NaN;
      B=A;
      B(find(~finite(A)))=0;
      MU=sum(B,2)./Q;      
   end   
end

if FLAG2==1
   if nr_mis_val==0
      SIGM = std(A')';
      SIGM = SIGM + (SIGM==0);
   else
      Nrrow=size(A,1);
		for row=1:Nrrow
   		M=std(A(row,find(finite(A(row,:)))));
   		if isempty(M)
      		M=NaN;
   		end
         SIGM(row,1)=M;
      end 
      SIGM = SIGM + (SIGM==0);
   end
end

if FLAG1>0
   A=A-(MU*ones(1,Samples));
end

if FLAG2>0
   A=A./(SIGM*ones(1,Samples));
end

MURES=MU;
SIGMRES=SIGM;



