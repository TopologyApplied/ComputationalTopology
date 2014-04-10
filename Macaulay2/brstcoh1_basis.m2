ofilename="result_brstcoh10d1spindiffbas.txt"
ofile=ofilename<<""<<endl

--needsPackage "DGAlgebras"

D=10;Nspin=16;
R=QQ[t_1..t_16,Degrees=>toList(16:{1,0})];
 rel=matrix{{t_16*t_1-t_15*t_4-t_8*t_12+t_9*t_11,
  t_16*t_2-t_15*t_5+t_7*t_12-t_9*t_10,
  t_16*t_3-t_15*t_6-t_7*t_11+t_8*t_10,
  t_14*t_7-t_13*t_10-t_2*t_6+t_3*t_5,
  t_14*t_8-t_13*t_11+t_1*t_6-t_3*t_4,
  t_14*t_9-t_13*t_12-t_1*t_5+t_2*t_4,
  t_10*t_1+t_11*t_2+t_12*t_3+t_14*t_15,
  t_7*t_1+t_8*t_2+t_9*t_3+t_13*t_15,
  t_7*t_4+t_8*t_5+t_9*t_6+t_13*t_16,
  t_10*t_4+t_11*t_5+t_12*t_6+t_14*t_16
  }};
R=R/ideal(rel);

C=QQ[c_1..c_16,SkewCommutative =>apply(toList(1..16),i->c_i),Degrees=>apply(toList(1..(16)),i->{0,-1})]
use C;
B=basis C;

diff(c_i,basis(j,C))

use R
for i from 1 to 16 do
(
	(cbase,cdiffi)=coefficients(diff(c_i,B),Monomials=>B);
	if i==1 then cdiff=matrix{{sum(1..16,j->(t_j)*sub(cdiffi,R))}} else	cdiff=cdiff+matrix{{sum(1..16,j->(t_j)*sub(cdiffi,R))}};
)



d=map(R^(2^Nspin),R^(2^Nspin),cdiff);
M=prune homology(d,d);
NumgenM=numgens M
