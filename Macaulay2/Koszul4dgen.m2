--ULIMit -v unlimited; GC_FREE_SPACE_DIVISOR=12 GC_MAXIMUM_HEAP_SIZE=7000000000 GC_NPROCS=4 GC_MARKERS=4 M2 Koszul4d.m2
needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
D=4
Nspin=4
kk= ZZ/32003
R=kk[t_1..t_Nspin];
ofilename="result_cohgen"|D|"dKoszul.txt"
ofile=ofilename<<""<<endl

-- DiracC2.DiracGamma --
----NEED to check all relations: Tleft*gamma*Tright!=0----  
gamma={matrix{{0, 0, -1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, 0, 1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, -1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
			matrix{{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}



redux = (n0,j,i)->(
n:=n0;M:=1;while n<=j and M!=0 do
			(
			M = HH_n K;
			prune M;
			s = reduceHilbert(hilbertSeries(M));
			<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<numgens M<<endl;
			n=n+1;
			)
)

--I0list={}
--M0list={}
M0=R
M1=0
i=0;while i<=D-1 and M1==0 do 
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	Tleft = matrix{{t_1..t_Nspin}};
	Tright= matrix{{t_1}..{t_Nspin}};
	I0 = ideal(Tleft*(gamma#i)*Tright);
--	I0list = I0list|{I0};
	if i>=0 then
	 (
--	rR=M0/I0;
--	A = koszulComplexDGA(rR);--to compute valid numgens--
		A = koszulComplexDGA(I0);--over the ideal--
  K = toComplex A;
  M1 = HH_1 K;
  prune M1;
     apply(maxDegree A + 1, j -> numgens prune homology(j,A));
     HA = homologyAlgebra(A,GenDegreeLimit=>i+2,RelDegreeLimit=>i+2);
--   HA = homologyAlgebra(A);
     HAnumgen=numgens HA;
     HAcycles=HA.cache.cycles;
	<<"Number of non-trivial 0-th cohomology generators="<<HAnumgen<<endl;
	ofile<<"Number of non-trivial 0-th cohomology generators="<<HAnumgen<<endl;
  <<"The cycles are:"<<HAcycles<<endl;
	ofile<<"The cycles are:"<<HAcycles<<endl;
	);
  if M1==0 then
  	(
		M0=M0/I0;--New QuotientRing applied--
		prune M0;
--	M0list=M0list | {M0};
		s0 =reduceHilbert(hilbertSeries(M0));
		<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s0="<< s0 <<endl;
		ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s0="<< s0 <<endl;
--		ofile<<"Number of 0-th cohomology generators="<<numgens M0<<endl;
		i=i+1;
		)
	else
		(
		redux(0,0,i);
		s = reduceHilbert(hilbertSeries(M1));
		<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s1="<< s <<endl;
		ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s1="<< s <<endl;
--		ofile<<"Number of non-trivial 1st cohomology generators="<<numgens M1<<endl;
		redux(2,i,i);
		);
)
--i quit at 5, i.e. D=6+4--

i0=i
Tleft = matrix{{t_1..t_Nspin}}
Tright= matrix{{t_1}..{t_Nspin}}
--Mlist={}
I1={}
--I1={Tleft*(gamma#5)*Tright,Tleft*(gamma#6)*Tright,Tleft*(gamma#7)*Tright}
for i from i0 to D-1 do
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	I1=I1|{Tleft*(gamma#i)*Tright};
	I=ideal(I1);
--	rR=M0/I;
--  A = koszulComplexDGA(rR);--to compute valid numgens--
	A = koszulComplexDGA(I);--over the ideal--
	apply(maxDegree A + 1, j -> numgens prune homology(j,A));
	HA = homologyAlgebra(A);
--	HA = homologyAlgebra(A,GenDegreeLimit=>4,RelDegreeLimit=>4);
--	HA = homologyAlgebra(A,GenDegreeLimit=>i+3);
	HAnumgen=numgens HA;
  HAcycles=HA.cache.cycles;
	<<"Number of total non-trivial cohomology generators="<<HAnumgen<<endl;
	<<"The cycles are:"<<HAcycles<<endl;
	ofile<<"Number of total non-trivial cohomology generators="<<HAnumgen<<endl;
	ofile<<"The cycles are:"<<HAcycles<<endl;
--	K = toComplex A;
-- K.dd;
--	redux(0,i,i);
)

ofile<<close
<< get ofilename
