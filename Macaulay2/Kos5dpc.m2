--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
needsPackage "BoijSoederberg" 
needsPackage "ChainComplexExtras"
D=5
Nspin=8
R=ZZ/32003[t_1..t_Nspin]
ofilename="result_Kos"|D|"dpc.txt"
ofile=ofilename<<""<<endl


-- DiracC1.DiracGamma;DiracB1.DiracGamma --
--gamma#i=gamma4D#i,gamma#5=gamma4d#1 gamma4d#2 gamma4d#3 gamma4d#4
----NEED to check all relations: Tleft*gamma*Tright!=0----
gamma={matrix{{0, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 
  0, 0, 0}, {0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0}, {0, 
  0, 0, -1, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0,
   0, 0}},
			matrix{{0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, -1, 
  0, 0, 0}, {0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, -1, 0, 0, 0, 0, 0}, {0, 
  0, 0, -1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 
  0, 0}},
			matrix{{0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 
  0, 0, -1}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {-1, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0, 
  0, 0}},
			matrix{{0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 
  0, 0, -1}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, -1, 0, 0, 0, 0, 0, 0}, {1, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0, 
  0, 0}},
			matrix{{0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 
  0, -1, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 
  0}, {0, 0, -1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 
  0, 0, 0, 0, 0}}}

redux = (Mg)->(
i=0;while i<=D and Mg_i!=0 do
			(
			s=reduceHilbert(hilbertSeries(Mg_i));
			<<"D="<<D<<"_s"<<i<<"="<<s<<endl;
			ofile<<"D="<<D<<"_s"<<i<<"="<<s<<endl;
			i=i+1;
			)
)

Tleft = matrix{{t_1..t_Nspin}}
Tright= matrix{{t_1}..{t_Nspin}}
I1=apply(toList(0..(D-1)),i->Tleft*(gamma#i)*Tright)
I=ideal(I1)
A = koszulComplexDGA(I) 
K = toComplex A
M = HH K
ofile<<"The cohomology of F+ (D="<<D<<"):"<<endl;
ofile<<M<<endl;
redux(M)


Mc=M
Ra=R
i=1;while i<=1 do
(
	Rc=Ra[u_i]/(t_i*u_i-1);
	f=map(Rc,Ra);
	Mc=apply(toList(0..D),j->prune(tensor(Rc,f,Mc_j)));
	if Mc==apply(toList(0..D),j->0) then break;
	Ra=Rc;
	i=i+1;
)
Mc
ofile<<"The cohomology of F by the localization of F+ (D="<<D<<"):"<<endl;
ofile<<Mc<<endl;
-- apply(toList(0..D),i->isFreeModule Mc_i)
redux(Mc)


ofile<<close
<< get ofilename
