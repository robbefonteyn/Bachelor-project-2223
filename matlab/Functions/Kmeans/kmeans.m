function [C,CM,DISTTOCM] = kmeans(A,N,I,DIM)

% [C,CM] = kmeans(A,N,I,DIM)
% This function performs k-means analysis on the column vectors (samples) of A (No preprocessing is done)
% It is advised that PCA is performed (see function pcasam) before using this function (only the first (=DIM) PC's are used!)
% N is the number of clusters that has to be found
% I is the number of iterations of the algorithm
% DIM is the number of components/dimensions that is used
% The outputs are:
% C is the class (or in this case cluster) distinction (row-)vector
% CM contains the cluster means (column vectors)

samples = size(A,2);
A=A(1:DIM,:);

CI=(1:N)';
CI=CI*ones(1,ceil(samples/N));
CI=CI(1:samples);

for loop=1:I
   CMI=[];
   for cluster=1:N
      E =(CI==cluster);
      CMI(:,cluster)=(E*A')'/sum(E);
   end
   CI=[];
   for s=1:samples
      sam=A(:,s);
      dis=sqrt(sum((CMI-sam*ones(1,N)).^2));
      [m,p]=min(dis);
      CI=[CI p];
   end
end
CMI=[];
for cluster=1:N
   E =(CI==cluster);
   CMI(:,cluster)=(E*A')'/sum(E);
end
DISTTOCM=[];
for s=1:samples
   sam=A(:,s);
   sm=CMI(:,CI(s));
   dis=sqrt(sum((sam-sm).^2));
   DISTTOCM=[DISTTOCM dis];
end
   
C=CI-1;
CM=CMI;




