--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

t1 = cpuTime()
needsPackage "DGAlgebras"
needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
D=10;
spindim=16;
Nsusy=1;
Nspin=spindim*Nsusy;
msusy=id_(ZZ^Nsusy);
kk= ZZ/32003;
R=kk[t_1..t_Nspin];
ofilename="result_cohgen"|D|"dN"|Nsusy|"dKoszul.txt"
ofile=ofilename<<""<<endl


--gamma matrices computed based on Clifford Algebra--
D2=sub(D/2,ZZ);
G=QQ[e_1..e_D2,SkewCommutative => true]
use G;
evenbasis=basis(0,G)|basis(2,G)|basis(4,G);
oddbasis=basis(1,G)|basis(3,G)|basis(5,G);
gamma={};
	for i from 1 to D2 do
	(
		(TxGamma, gammai)=coefficients(e_i**evenbasis,Monomials=>oddbasis); 
		gamma=join(gamma,{sub(gammai,R)});
		(TxGamma, gammai)=coefficients(diff(e_i, evenbasis),Monomials=>oddbasis); 
		gamma=join(gamma,{sub(gammai,R)});
	)
use R;

Tleft = matrix{{t_1..t_Nspin}};	Tright= matrix{{t_1}..{t_Nspin}};
Ilist=apply(toList(0..(D-1)),i->(Tleft*(gamma#i)*Tright)_(0,0))  --to get a list, not matrix--
numgens trim ideal Ilist


