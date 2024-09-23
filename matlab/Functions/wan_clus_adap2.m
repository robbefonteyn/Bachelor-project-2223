function [C,CM,NRCLRET] = wan_clus_adap2(A,MAXNRPOINT,EXTTRESH,CL,QUALCRIT,PLOT,DIM,A_ORG)

# Input:
# A = An array or dataset
# MAXNRPOINT = minimal number of genes in the sphere
# EXTTRESH = Rk_prelim in the paper (?), Used in the calculation of DELTARAD. It is the preliminary estimate
#            of the radius
# CL = array of cluster assignment (with collumns equal to sample, zeroes are further used to put genes into)
# QUALCRIT = an integer (0,1 or 2) used for quality criterion
# PLOT = 0 or 1, generates a plot
# DIM = the number of components/dimensions that is used
# A_ORG = variable that has no use in this algorithm, (might just be copied from the kmeans functions?)

# Outputs:
# C = the class (or in this case cluster) distinction (row-)vector
# CM = contains the cluster means (column vectors)
# NRCLRET = number of clusters found

# samples = Number of columns in the dataset
samples = size(A,2);

if EXTTRESH <= 0
   disp('EXTTRESH should be greater than 0 !!')
   return
end

# MAXITER Maximum times the algorithm will iterate,
# NRWAN is variable used for calculating the DELTARAD variable (Called DIV in paper, is also a fraction)
MAXITER=500;
NRWAN=30;

# Setting variables to standard values if they are not specified
if nargin < 4 | isempty(CL)
   C = zeros(1,samples);
else
   C = CL;
end

if nargin < 5 | isempty(QUALCRIT)
   QUALCRIT = 2;
end

if nargin < 6 | isempty(PLOT)
   PLOT = 1;
end

if nargin < 7 | isempty(DIM)
   DIM = size(A,1);
end

# Reshapes the array (or matrix?) to a matrix with rows equal to DIM and
# columns equal to the original number of columns in the array/matrix (Done by the ":")
# Only necessary if DIM is given as an input
A=A(1:DIM,:);

# Cluster should start at maximum value in C, to not interfere with existing cluster that are passed with "CL"
cluster=max(C);
CM=[];
CONV=1;
it=0;
TOOFEWPOINTS=0;

while CONV == 1 & TOOFEWPOINTS < 5
   it = it + 1
# Turns the arrays into something
   PT = find(C == 0);
   SP = A(:,PT);
   if isempty(SP)
      break
   end

# Calculates mean of the newly made matrix SP with the MODE set to the Quality criterion
   ME1 = mean_misval(SP, QUALCRIT);
# Calculates the distance using the mean, the array and the quality criterion)
   RAD = max(dist_misval(SP, ME1, QUALCRIT))
   DELTARAD = (RAD-EXTTRESH) / NRWAN;
   RAD = RAD - DELTARAD;

# Q = GENES_IN_SPHERE in the paper, variable calculated to determine profiles within sphere
   Q = find(dist_misval(SP, ME1, QUALCRIT) < RAD);
   CONV = 0;
   TEL = 0;
   for iter = 1:MAXITER
      # Recalculate mean using the new genes in sphere
      ME2 = mean_misval(SP(:, Q), QUALCRIT);

      # Check that sets the radius to the RK_PRELIM (preliminary estimate of the radius)
      if iter < NRWAN + TEL
         RAD = RAD - DELTARAD;
      else
         RAD = EXTTRESH;
      end
      # Find the new GENES_IN_SPHERE (or Q)
      Q = find(dist_misval(SP, ME2, QUALCRIT) < RAD);

      # Check if GENES_IN_SPHERE is empty and we have not reached the last iteration(?)
      if isempty(Q) & iter < NRWAN + TEL
         # Recalculate the GENES_IN_SPHERE, increase TEL and Recalculate the Radius by adding DELTARAD
         TEL = TEL + 1;
         RAD = RAD + DELTARAD;
         Q = find(dist_misval(SP, ME2, QUALCRIT) < RAD);
      end


      if iter >= NRWAN + TEL
         # Check in convergence is reached and set CONV to 1 if so
         if ((ME1==ME2)+((~isfinite(ME1))&(~isfinite(ME2))))
            CONV=1;
            break
         end
      end
      # If no convergence is reached and the check before this is not TRUE Use the recalculated
      # mean as the first mean and repeat
      ME1=ME2;
   end


   if CONV == 1
      if length(Q) > MAXNRPOINT
         cluster = cluster + 1;
         CM = [CM ME2];
         C(PT(Q)) = cluster;
         disp(['		Found cluster with ',num2str(length(Q)),' genes / Convergence after ',num2str(iter),' loops'])
         TOOFEWPOINTS=0;
      else
         C(PT(Q))=-1;
         disp(['Discarded cluster with ',num2str(length(Q)),' genes'])
         TOOFEWPOINTS = TOOFEWPOINTS + 1;
      end
   else
      disp(['Warning: no convergence: Not all clusters may have been found !'])
   end         
end

# NC is used to generate a plot
NC = 7;

if nargin==8
   AT = A_ORG';
else
   # AT is equal to the transposed array/matrix
   AT = A';
end

# Code to generate a plot of the dataset
if PLOT==1
   SP=ceil(cluster/NC);
   f=figure;
   fl=0;
   for loop=1:cluster
   	MV=AT(find(C==loop),:);
    	if ~isempty(MV)
      	MN=mean_misval(MV',QUALCRIT)';
      	hs = subplot(SP,NC,loop);
      	set(hs,'FontSize',7);
      	for iter=1:sum(C==loop)
				plot(MV(iter,:),'k:')
				hold on
      	end
        	p = plot(MN,'k');
      	set(p,'LineWidth',2.5);
      % 	t = title(['NOG=',num2str(sum(C==loop))]);
       	t = title(['Cluster ',num2str(loop),'  NOG=',num2str(sum(C==loop))]);
      	set(t,'FontSize',7);
  			set(t,'Color',[0 0 1]);
        if NC*ceil(loop/NC)>cluster & fl==0
         %  if loop==26
      		xl = xlabel('Experiments');
      		set(xl,'FontSize',7);
      		yl = ylabel('Normalized expr. level');
            set(yl,'FontSize',7);
            fl=1;
   		end
        	hold off
      end
  end
end
   
C=C-(C>-1);
NRCLRET=cluster;
