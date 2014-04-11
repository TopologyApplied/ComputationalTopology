--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

ofilename="result_brstcoh10d2spindiffbas.txt"
ofile=ofilename<<""<<endl

--needsPackage "DGAlgebras"

D=10;Nspin=16;
R=QQ[t_1..t_16,s_1..s_16,Degrees=>toList(32:{1,0})];
 rel=matrix{{t_16*t_1-t_15*t_4-t_8*t_12+t_9*t_11,
  t_16*t_2-t_15*t_5+t_7*t_12-t_9*t_10,
  t_16*t_3-t_15*t_6-t_7*t_11+t_8*t_10,
  t_14*t_7-t_13*t_10-t_2*t_6+t_3*t_5,
  t_14*t_8-t_13*t_11+t_1*t_6-t_3*t_4,
  t_14*t_9-t_13*t_12-t_1*t_5+t_2*t_4,
  t_10*t_1+t_11*t_2+t_12*t_3+t_14*t_15,
  t_7*t_1+t_8*t_2+t_9*t_3+t_13*t_15,
  t_7*t_4+t_8*t_5+t_9*t_6+t_13*t_16,
  t_10*t_4+t_11*t_5+t_12*t_6+t_14*t_16,
  s_16*s_1-s_15*s_4-s_8*s_12+s_9*s_11,
  s_16*s_2-s_15*s_5+s_7*s_12-s_9*s_10,
  s_16*s_3-s_15*s_6-s_7*s_11+s_8*s_10,
  s_14*s_7-s_13*s_10-s_2*s_6+s_3*s_5,
  s_14*s_8-s_13*s_11+s_1*s_6-s_3*s_4,
  s_14*s_9-s_13*s_12-s_1*s_5+s_2*s_4,
  s_10*s_1+s_11*s_2+s_12*s_3+s_14*s_15,
  s_7*s_1+s_8*s_2+s_9*s_3+s_13*s_15,
  s_7*s_4+s_8*s_5+s_9*s_6+s_13*s_16,
  s_10*s_4+s_11*s_5+s_12*s_6+s_14*s_16
  }};
R=R/ideal(rel);

C=QQ[c_1..c_16,SkewCommutative =>apply(toList(1..16),i->c_i),Degrees=>apply(toList(1..(Nspin)),i->{0,-1})]
use C;
B=basis C;

use R
for i from 1 to Nspin do
(
	(cbase,cdiffi)=coefficients(diff(c_i,B),Monomials=>B);
	if i==1 then cdiff=matrix{{sum(1..16,j->(t_j+s_j)*sub(cdiffi,R))}} else	cdiff=cdiff+matrix{{sum(1..16,j->(t_j+s_j)*sub(cdiffi,R))}};
)

d=map(R^(2^Nspin),R^(2^Nspin),cdiff);
M=prune homology(d,d);
NumgenM=numgens M

d1=map(R^(2^(Nspin-1)),R^(2^(Nspin-1)),cdiff1);
M1=prune homology(d1,d1);
NumgenM1=numgens M1

--Q=R ** B
--use Q;
--d=map(Q^1,Q^1,matrix{{sum(1..16,i->(t_i+s_i)*b_i)}});
M=prune homology(d,d);
NumgenM=numgens M;
GenM=generators image M.cache.pruningMap;
s=reduceHilbert(hilbertSeries(M));
<<"Number of cohomology generators="<<NumgenM<<endl;
<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--	<<"The cohomology generators are:"<<endl<<GenM<<endl;
ofile<<"Number of cohomology generators="<<NumgenM<<endl;
ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
ofile<<"The cohomology generators are:"<<endl<<GenM<<endl; 


ofile<<close
<< get ofilename



evenbasis=basis({0,0},C)|basis({0,-2},C)|basis({0,-4},C)|basis({0,-6},C)|basis({0,-8},C)|basis({0,-10},C)|basis({0,-12},C)|basis({0,-14},C)|basis({0,-16},C)
oddbasis=basis({0,-1},C)|basis({0,-3},C)|basis({0,-5},C)|basis({0,-7},C)|basis({0,-9},C)|basis({0,-11},C)|basis({0,-13},C)|basis({0,-15},C)

for i from 1 to Nspin do
(
	(cbase,cdiffi)=coefficients(diff(c_i,oddbasis),Monomials=>evenbasis);
	if i==1 then cdiff1=sub(cdiffi,R) else	cdiff1=cdiff1+sub(cdiffi,R);
)
