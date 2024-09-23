function P=clusterdistrib(r,VAR,D,R)

SD=(2*(pi^(D/2)))/GAMMA(D/2);

P=SD*(1/((2*pi*VAR)^(D/2)))*(r.^(D-1)).*exp(-((r.^2)/(2*VAR)));

P(find(r>R))=0;