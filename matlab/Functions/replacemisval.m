function AOUT = replacemisval(A)

nr_mis_val = sum(sum(~finite(A)));

if nr_mis_val==0
   AOUT=A;
else
   Nrrow=size(A,1);
	Q=sum(finite(A),2);
	for row=1:Nrrow
   	M=sum(A(row,find(finite(A(row,:)))));
   	if Q(row)==0
      	M==0;
   	else
      	M=M/Q(row);
   	end
   	A(row,~finite(A(row,:)))=M;
	end
	AOUT=A;
end



