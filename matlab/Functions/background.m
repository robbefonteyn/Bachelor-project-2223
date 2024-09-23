function P=background(r,D,R)

SD=(2*(pi^(D/2)))/GAMMA(D/2);
SD1=(2*(pi^((D+1)/2)))/GAMMA((D+1)/2);

P=(SD/(SD1*(R^D)))*(r.^(D-1));

P(find(r>R))=1;