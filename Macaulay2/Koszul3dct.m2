--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
D=3
spindim=2
Nsusy=1
Nspin=spindim*Nsusy
msusy=id_(ZZ^Nsusy)
kk= ZZ/32003
--R=kk[t_1..t_Nspin,p_1..p_D,SkewCommutative =>{p_1..p_D}]
R=kk[t_1..t_Nspin,p_1..p_D,SkewCommutative =>{p_1,p_2,p_3}]
ofilename="result_cohgen"|D|"dN"|Nsusy|"dKoszul.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl


-- DiracC1.DiracGamma;DiracB1.DiracGamma --
--gamma#i=gamma4D#i,gamma#5=gamma4d#1 gamma4d#2 gamma4d#3 gamma4d#4
----NEED to check all relations: Tleft*gamma*Tright!=0----
gamma={matrix{{1, 0}, {0, -1}}**msusy,
matrix{{1, 0}, {0, 1}}**msusy,
matrix{{0, -1}, {-1, 0}}**msusy}

redux = (n0,j,i)->(
n:=n0;M:=1;while n<=j and M!=0 do
			(
			M = HH_n K;
--		M=prune M;--also improve the ambient free module --
			M=trim M;--faster--
			s = reduceHilbert(hilbertSeries(M));
--			NumgenM=numgens trim M;
--			GenM=mingens M;
			NumgenM=numgens M;
			GenM=generators M;
--			resM=resolution(M);
			<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
--			<<"The resolution of M_"<<n<<": "<<resM<<endl;
			ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
--			ofile<<"The resolution of M_"<<n<<": "<<resM<<endl;
			n=n+1;
			)
)



Tleft = matrix{{t_1..t_Nspin}}
Tright= matrix{{t_1}..{t_Nspin}}
plist= matrix{{p_1}..{p_D}}
d=sum(0..(D-1),i->Tleft*(gamma#i)*Tright*plist^{i})
I1={}
i=D-1
--I1={Tleft*(gamma#5)*Tright,Tleft*(gamma#6)*Tright,Tleft*(gamma#7)*Tright}
--for i from i0+1 to D-1 do
--for i from i0 to D-1 do --if i quit at i0+1--
--(
--	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
--	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<",NSusy="<<Nsusy<<":"<<endl;
--	I1=I1|{Tleft*(gamma#i)*Tright};
	I1=apply(toList(0..(D-1)),i->Tleft*(gamma#i)*Tright);
	I1=I1|apply(toList(0..(Nspin-1)),j->Tleft*matrix(apply(toList(0..(D-1)),k->gamma#k_j))*plist);
	I=ideal(I1);
--	rR=M0/I;
--  A = koszulComplexDGA(rR);
	A = koszulComplexDGA(I);--over the ideal--
--	apply(maxDegree A + 1, j -> numgens prune homology(j,A));
--	HA = homologyAlgebra(A);
--	HA = homologyAlgebra(A,GenDegreeLimit=>4,RelDegreeLimit=>4);
--	HA = homologyAlgebra(A,GenDegreeLimit=>i+3);
--	HAnumgen=numgens HA;
--      HAcycles=HA.cache.cycles;
--	<<"Number of total non-trivial cohomology generators="<<HAnumgen<<endl;
--	<<"The cycles are:"<<HAcycles<<endl;
--	ofile<<"Number of total non-trivial cohomology generators="<<HAnumgen<<endl;
--	ofile<<"The cycles are:"<<HAcycles<<endl;
	K = toComplex A;
-- K.dd;
	redux(0,Nspin+i,i);
--)

ofile<<close
<< get ofilename
