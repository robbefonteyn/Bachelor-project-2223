function [RAD,CONV,PCRET,PBRET,VARRET] = exp_max(AS,CK,QUAL)

E=size(AS,1);
D=E-2;
R=sqrt(E-1);

RD=dist_misval(AS,CK,2);
RD=RD(find(isfinite(RD)));
samples = size(RD,2);

MAXITER=1000;

CDIF=0.001;

S=0.95;

PC=sum(RD<QUAL)/samples;
PB=1-PC;
VAR=(1/D)*sum(RD(find(RD<QUAL)).^2)/sum(RD<QUAL);

CONV=0;
for em=1:MAXITER
   em
   prc=clusterdistrib(RD,VAR,D,R);
   prb=background(RD,D,R);
   prcpc=prc*PC;
   prbpb=prb*PB;
   pr=sum([prcpc;prbpb],1);
   pcr=prcpc./pr;
   pbr=prbpb./pr;
   VAR_new=(1/D)*sum((RD.^2).*pcr)/sum(pcr);
   PC_new=sum(pcr)/samples;
   PB_new=sum(pbr)/samples;
   if abs(VAR_new-VAR)<CDIF & abs(PC_new-PC)<CDIF & abs(PB_new-PB)<CDIF
      CONV=1;
      break
   end
   PC=PC_new;
   PB=PB_new;
   VAR=VAR_new;     
end

if CONV==1
   PC=PC_new;
   PB=PB_new;
   VAR=VAR_new;
   
   SD=(2*(pi^(D/2)))/GAMMA(D/2);
   SD1=(2*(pi^((D+1)/2)))/GAMMA((D+1)/2);
   CC=SD*(1/((2*pi*VAR)^(D/2)));
   CB=(SD/(SD1*(R^D)));
   LO=(S/(1-S))*((PB*CB)/(PC*CC));
   DIS=-2*VAR*log(LO);

   if DIS < 0
      CONV=0;
      disp(['Not possible to calculate radius'])
      RAD=0;
      VARRET=0;
      PCRET=0;
      PBRET=0;
      return
   end
   
   RAD=sqrt(DIS);
   VARRET=VAR;
   PCRET=PC;
   PBRET=PB;
   
else
   disp(['No Convergence'])
   RAD=0;
   VARRET=0;
   PCRET=0;
   PBRET=0;
end


   