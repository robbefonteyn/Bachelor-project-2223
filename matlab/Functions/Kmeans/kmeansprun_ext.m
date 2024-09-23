function [C,CM,DISTTOCM,NRCLRET] = kmeansprun_ext(A,MAXNRPOINT,EXTTRESH,QUALCRIT,RAND,PLOT,DIM,A_ORG)

% [C,CM,DISTTOCM,NRCLRET] = kmeansprun_ext(A,MAXNRPOINT,EXTTRESH,QUALCRIT,RAND,PLOT,DIM,A_ORG)

if EXTTRESH <= 0
   disp('EXTTRESH should be greater than 0 !!')
   return
end

NUMCLUS=1;
MAXITER=100;
STARTPRUN_EXT=3;
REUSESAM=0;

if nargin < 4 | isempty(QUALCRIT)
   QUALCRIT = 1;
end

if nargin < 5 | isempty(RAND)
   RAND = 1;
end

if nargin < 6 | isempty(PLOT)
   PLOT = 1;
end

if nargin < 7 | isempty(DIM)
   DIM = size(A,1);
end

samples = size(A,2);

A=A(1:DIM,:);
CONV=0;

if NUMCLUS==1
   WAIT=STARTPRUN_EXT;
   CI=ones(1,samples);
else
   WAIT=0;
   CI=(1:NUMCLUS)';
	CI=CI*ones(1,ceil(samples/NUMCLUS));
	NP=randperm(samples);
	CI=CI(NP);
end
Q=CI;

for loop=1:MAXITER
   disp(['Starting iteration ',num2str(loop),' with ',num2str(NUMCLUS),' cluster(s)'])
   WAIT=WAIT+1;
   CMI=[];
   for cluster=1:NUMCLUS
      E = (CI==cluster);
      if sum(E)==0
         CMI(:,cluster)=QM(:,cluster);
      else
         CMI(:,cluster)=(E*A')'/sum(E);
      end
   end
   
   if WAIT >= STARTPRUN_EXT
      cluster=1;
      while cluster <= NUMCLUS
         E = (CI==cluster);
         if sum(E)<=MAXNRPOINT
            NUMCLUS=NUMCLUS-1;
            CMI(:,cluster)=[];
            if REUSESAM==0
               CI=(CI-E.*(CI+1));
            else
               CI=(CI-E.*CI);
               WAIT=0;
            end
            CI=CI-(CI>cluster);
         else
            cluster=cluster+1;
         end
      end
      
      cluster=1;
      while cluster <= NUMCLUS
         ACL=A(:,find(CI==cluster));
         CLSAM = size(ACL,2);
         ACLMU=ACL-CMI(:,cluster)*ones(1,CLSAM);
         if QUALCRIT==1
            mdis=max(max(abs(ACLMU)));
         elseif QUALCRIT==2
            mdis=max(sqrt(sum(ACLMU.^2)));
         else
            mdis=max(sum(abs(ACLMU)));
         end
                     
         if mdis > EXTTRESH
            WAIT=0;
            if RAND == 1
               SEP=(1:2)';
					SEP=SEP*ones(1,ceil(CLSAM/2));
					NSEP=randperm(CLSAM);
					SEP=SEP(NSEP);
               C1=mean(ACL(:,find(SEP==1)),2);
               C2=mean(ACL(:,find(SEP==2)),2);
            else
               [U1,S1,V] = svd(ACLMU,0);
					S1=diag(S1);
   				S1=sqrt(1/(CLSAM-1))*S1(1);
               ST=S1*U1(:,1);
               C1=CMI(:,cluster)-ST;
               C2=CMI(:,cluster)+ST;
            end
            CMI=[CMI(:,1:(cluster-1)) C1 C2 CMI(:,(cluster+1):NUMCLUS)];
            NUMCLUS=NUMCLUS+1;
            CI=CI+(CI>cluster);
            cluster=cluster+2;
         else
            cluster=cluster+1;
         end
      end
   end
   
   CIC=[];
   for s=1:samples
      if CI(s)==-1
         CIC=[CIC -1];
      else
         sam=A(:,s);
      	dis=sqrt(sum((CMI-sam*ones(1,NUMCLUS)).^2));
      	[m,p]=min(dis);
         CIC=[CIC p];
      end
   end
   CI=CIC;
   
   QM=CMI;
   
   if sum((Q==CI)==0)==0 & WAIT >= STARTPRUN_EXT
      disp(['Convergence after ',num2str(loop),' iterations'])
      disp([num2str(NUMCLUS),' cluster(s) left'])
      disp(' ')
      CONV=1;
      break
   end
   
   Q=CI;
   
end

if CONV==0
   CMI=[];
	for cluster=1:NUMCLUS
   	E =(CI==cluster);
   	if sum(E)==0
      	CMI(:,cluster)=QM(:,cluster);
   	else
      	CMI(:,cluster)=(E*A')'/sum(E);
      end
   end
   disp('No convergence after maximum iterations')
   disp([num2str(NUMCLUS),' cluster(s) left'])
   disp(' ')
end

DISTTOCM=[];
for s=1:samples
   if CI(s)==-1
      DISTTOCM=[DISTTOCM -1];
   else
      sam=A(:,s);
      sm=CMI(:,CI(s));
      
      if QUALCRIT==1
         dis=max(abs(sam-sm));
      elseif QUALCRIT==2
         dis=sqrt(sum((sam-sm).^2));
      else
         dis=sum(abs(sam-sm));
      end
      DISTTOCM=[DISTTOCM dis];
   end
end

if nargin==8
   AT=A_ORG';
else
   AT=A';
end

SP=ceil(NUMCLUS/5);
f1=figure;
f2=figure;
if PLOT==1
   for loop=1:NUMCLUS
      MV=AT(find(CI==loop),:);
      if ~isempty(MV)
         MN=CMI(:,loop)';
         ST=std(MV,0,1);
         figure(f1)
         subplot(SP,5,loop)
         for iter=1:sum(CI==loop)
            plot(MV(iter,:))
            hold on
         end
         hold off
         figure(f2)
      	subplot(SP,5,loop)
         errorbar([1:size(AT,2)],MN,ST)
         title(['Nr=',num2str(sum(CI==loop))])
      end
   end
end

C=CI-(CI>-1);
CM=CMI;
NRCLRET=NUMCLUS;




