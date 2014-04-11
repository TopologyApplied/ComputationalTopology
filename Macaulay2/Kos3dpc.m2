--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
needsPackage "BoijSoederberg" 
needsPackage "ChainComplexExtras"
D=3
Nspin=2
R=ZZ/32003[t_1..t_Nspin]
ofilename="result"|D|"dKoszul.txt"
ofile=ofilename<<""<<endl


-- DiracC1.DiracGamma;DiracB1.DiracGamma --
--gamma#i=gamma4D#i,gamma#5=gamma4d#1 gamma4d#2 gamma4d#3 gamma4d#4
----NEED to check all relations: Tleft*gamma*Tright!=0----
gamma={matrix{{1, 0}, {0, -1}},
			matrix{{1, 0}, {0, 1}},
			matrix{{0, -1}, {-1, 0}}}

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
M = prune M
ofile<<"The cohomology of F+ (D="<<D<<"):"<<endl;
ofile<<M<<endl;
redux(M)


M0=M
--Mc=M
--Ra=R
i=1;while i<=Nspin do
(
	Rc=R[u_1..u_i]/(apply(toList(1..i),i->t_i*u_i-1));
--DEBUG:
	<<"i="<<i<<endl;
	f=map(Rc,R);
	Mc=apply(toList(0..D),j->prune(tensor(Rc,f,M_j)));
--	lozcokerff=Mc;
	cokerFf=apply(toList(0..D),j->prune(M_j/ker(map(M_j,M_j,t_i^10))));
	lozcokerff=apply(toList(0..D),j->prune(tensor(Rc,f,cokerFf_j)));
	ofile<<"F"|horizontalJoin(unstack "+"^(0,Nspin-i+1))|horizontalJoin(unstack "-"^(0,i-1))<<":"<<endl<<M<<endl;
	ofile<<"F"|horizontalJoin(unstack "+"^(0,Nspin-i))|horizontalJoin(unstack "-"^(0,i-1))<<" from the localization of the cokernel of F+=>F+:"<<endl<<lozcokerff<<endl;
	ofile<<"F"|horizontalJoin(unstack "+"^(0,Nspin-i))|horizontalJoin(unstack "-"^(0,i))<<" from the localization of F+:"<<endl<<Mc<<endl;
--verify lozcokerff being equivalent to Mc--
ofile<<"--verify lozcokerff being equivalent to Mc--"<<endl;
fm=apply(toList(0..D),j->inducedMap(lozcokerff_j,Mc_j));
ofile<<"apply(toList(0..D),j->prune ker fm_j)"<<endl<<apply(toList(0..D),k->prune ker fm_k)<<endl;
ofile<<"apply(toList(0..D),j->prune coker fm_j)"<<endl<<apply(toList(0..D),k->prune coker fm_k)<<endl;

frr=map(R,Rc);
t1N=apply(toList(0..D),j->apply(toList(1..(numgens lozcokerff_j)),k->apply(toList(1..(numgens M_j)),l->1)));
Alhpa=apply(toList(0..D),j->map(frr**lozcokerff_j,M_j,t_i^10*map(frr**lozcokerff_j,M_j,t1N_j)));
cokerAlpha0=apply(toList(0..D),j->prune coker Alhpa_j);
matrixAlpha0=apply(toList(0..D),j->relations(cokerAlpha0_j));
matrixAlpha=apply(toList(0..D),j->submatrix'(matrixAlpha0_j,{(numColumns matrixAlpha0_j)-1}));
cokerAlpha=apply(toList(0..D),j->coker matrixAlpha_j);
--cokerAlpha=apply(toList(0..D),j->prune coker Alhpa_j);
kerAlpha=apply(toList(0..D),j->prune ker Alhpa_j);
apply(toList(1..D),j->prune Ext^1(cokerAlpha_j,kerAlpha_(j-1)));
ofile<<"cokerAlpha:"<<endl<<cokerAlpha<<endl;
ofile<<"kerAlpha:"<<endl<<kerAlpha<<endl;

--	if Mc==apply(toList(0..D),j->0) then break;
--	Ra=Rc;
--M is over the ring R--
	M={prune cokerAlpha_0}|apply(toList(1..D),j->prune (cokerAlpha_j ++ kerAlpha_(j-1)));
	i=i+1;
)
Mc
ofile<<"The cohomology of F"|horizontalJoin(unstack "-"^(0,Nspin))|" by the localization of F+ (D="<<D<<"):"<<endl;
ofile<<M<<endl;
-- apply(toList(0..D),i->isFreeModule Mc_i)
redux(Mc)


ofile<<close
<< get ofilename

M
