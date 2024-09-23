function [C,CM] = kmeans3(A,N,I,DIM, PLOT)

% [C,CM] = kmeans2(A,N,I,DIM)
% This function performs k-means analysis on the column vectors (samples) of A (No preprocessing is done)
% It is advised that PCA is performed (see function pcasam) before using this function (only the first (=DIM) PC's are used!)
% N is the number of clusters that has to be found
% I is the number of maximum iterations of the algorithm
% DIM is the number of components/dimensions that is used
% The outputs are:
% C is the class (or in this case cluster) distinction (row-)vector
% CM contains the cluster means (column vectors)
% aangepate versie van Frank 1230601
samples = size(A,2);
if nargin < 4 | isempty(DIM)
    DIM=size(A,1);
end
A=A(1:DIM,:);
CONV=0;

CI=(1:N)';
CI=CI*ones(1,ceil(samples/N));
CI=CI(1:samples);
Q=CI;

for loop=1:I
   CMI=[];
   for cluster=1:N
      E = find(CI==cluster);
      if isempty(E)
         CMI(:,cluster)=QM(:,cluster);
      else
         CMI(:,cluster)=mean_misval(A(:,E),2);
      end
   end
   CI=[];
   for s=1:samples
      sam=A(:,s);
      dis=dist_misval(CMI,sam,2);
      [m,p]=min(dis);
      CI=[CI p];
   end
   QM=CMI;
   if sum(Q==CI)==samples
      disp(['		Convergence after ',num2str(loop),' iterations'])
      CONV=1;
      break
   end
   Q=CI;
end
if CONV==0
    CMI=[];
    for cluster=1:N
   	E = find(CI==cluster);
   	    if isempty(E)
      	    CMI(:,cluster)=QM(:,cluster);
   	    else
      	    CMI(:,cluster)=mean_misval(A(:,E),2);
        end
        disp('		No convergence after maximum iterations')
    end
end
C=CI-1;
CM=CMI;

