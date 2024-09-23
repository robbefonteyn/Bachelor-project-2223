function [C,CM,DISTTOCM,NRCLRET] = kmeansprun_ext_merge3(A,MAXNRPOINT,EXTTRESH,QUALCRIT,RAND,PLOT,DIM,A_ORG)

% [C,CM,DISTTOCM,NRCLRET] = kmeansprun_ext_merge3(A,MAXNRPOINT,EXTTRESH,QUALCRIT,RAND,PLOT,DIM,A_ORG)

if EXTTRESH <= 0
   disp('EXTTRESH should be greater than 0 !!')
   return
end

NUMCLUS=1;
MAXITER=400;
STARTPRUN_EXT=3;
REUSESAM=0;

if nargin < 4 | isempty(QUALCRIT)
   QUALCRIT = 2;
end

if nargin < 5 | isempty(RAND)
   RAND = 0;
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
      if sum(CI==cluster)==0
         CMI(:,cluster)=QM(:,cluster);
      else
         CMI(:,cluster)=mean_misval(A(:,find(CI==cluster)),QUALCRIT);
      end
   end
   
   if WAIT >= STARTPRUN_EXT
      DM=ones(NUMCLUS,NUMCLUS)*inf;
      for l1=1:NUMCLUS
         for l2=(l1+1):NUMCLUS
            ACL=A(:,find(CI==l1|CI==l2));
            CLSAM = size(ACL,2);
            if CLSAM >0
               DM(l1,l2)=max(dist_misval(ACL,mean_misval(ACL,QUALCRIT),QUALCRIT));
            else
               DM(l1,l2)=0;
            end            
         end
      end
      while min(min(DM)) <= EXTTRESH
         WAIT=0;
         [mdis row]=min(DM,[],1);
         [mn col]=min(mdis);
         row=row(col);
         C1=min(row,col);
         C2=max(row,col);
         disp(['         Merging cluster ',num2str(C1),' and cluster ',num2str(C2)])
         ACL=A(:,find(CI==C1|CI==C2));
         CLSAM = size(ACL,2);
			if CLSAM > 0
            CMI(:,C1)=mean_misval(ACL,QUALCRIT);
         else
            CMI(:,C1)=(CMI(:,C1)+CMI(:,C2))/2;
         end
         CMI(:,C2)=[];
         CI=CI + ((CI==C2).*(ones(1,samples)*C1-CI));
         CI=CI-(CI>C2);
         NUMCLUS=NUMCLUS-1;
         DM(C2,:)=[];
         DM(:,C2)=[];
         for l1=(C1+1):NUMCLUS
            ACL=A(:,find(CI==l1|CI==C1));
            CLSAM = size(ACL,2);
            if CLSAM > 0
               DM(C1,l1)=max(dist_misval(ACL,mean_misval(ACL,QUALCRIT),QUALCRIT));
            else
               DM(C1,l1)=0;
            end            
         end
         for l1=1:(C1-1)
            ACL=A(:,find(CI==l1|CI==C1));
            CLSAM = size(ACL,2);
            if CLSAM > 0
               DM(l1,C1)=max(dist_misval(ACL,mean_misval(ACL,QUALCRIT),QUALCRIT));
            else
               DM(l1,C1)=0;
            end            
         end
      end
      
      cluster=1;
      while cluster <= NUMCLUS
         ACL=A(:,find(CI==cluster));
         CLSAM = size(ACL,2);
         mdis=max(dist_misval(ACL,CMI(:,cluster),QUALCRIT));
         if mdis > EXTTRESH
            disp(['         Extending cluster ',num2str(cluster)]);
            WAIT=0;
            if RAND == 1
               SEP=(1:2)';
					SEP=SEP*ones(1,ceil(CLSAM/2));
					NSEP=randperm(CLSAM);
					SEP=SEP(NSEP);
               C1=mean_misval(ACL(:,find(SEP==1)),QUALCRIT);
               C2=mean_misval(ACL(:,find(SEP==2)),QUALCRIT);
            else
               ACLMU=ACL-CMI(:,cluster)*ones(1,CLSAM);
               ACLMU=replacemisval(ACLMU);
               [U1,S1,V] = svd(ACLMU',0);
					S1=diag(S1);
   				S1=sqrt(1/(CLSAM-1))*S1(1);
               ST=S1*V(:,1);
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
   
   if WAIT >= STARTPRUN_EXT
      cluster=1;
      while cluster <= NUMCLUS
         E = (CI==cluster);
         if sum(E)<=MAXNRPOINT
            disp(['        Pruning cluster ',num2str(cluster)]);
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
   end

   
   CIC=[];
   for s=1:samples
      if CI(s)==-1
         CIC=[CIC -1];
      else
         sam=A(:,s);
         dis=dist_misval(CMI,sam,QUALCRIT);
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
      if sum(CI==cluster)==0
         CMI(:,cluster)=QM(:,cluster);
      else
         CMI(:,cluster)=mean_misval(A(:,find(CI==cluster)),QUALCRIT);
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
      dis=dist_misval(sam,sm,QUALCRIT);
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
         MN=mean_misval(MV',QUALCRIT)';
         figure(f1)
         subplot(SP,5,loop)
         for iter=1:sum(CI==loop)
            plot(MV(iter,:))
            hold on
         end
         hold off
         figure(f2)
      	subplot(SP,5,loop)
         plot(MN)
         title(['Nr=',num2str(sum(CI==loop))])
      end
   end
end

C=CI-(CI>-1);
CM=CMI;
NRCLRET=NUMCLUS;




