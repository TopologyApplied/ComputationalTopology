--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

--ulimit -v unlimited; GC_FREE_SPACE_DIVISOR=12 GC_MAXIMUM_HEAP_SIZE=7000000000 GC_NPROCS=4 GC_MARKERS=4 M2 Koszul4d.m2
needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
D=4
Nspin=4
kk= ZZ/32003
R=kk[t_1..t_Nspin];
ofilename="result_"|D|"dKoszul.txt"
ofile=ofilename<<""<<endl

-- DiracC2.DiracGamma --
----NEED to check all relations: Tleft*gamma*Tright!=0----  
gamma={matrix{{0, 0, -1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, 0, 1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, -1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
			matrix{{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}


redux = (n0,j,i)->(
n:=n0;M=1;while n<=j and M!=0 do
			(
			M = HH_n K;
			M=prune M;--also improve the ambient free module --
--			M=trim M;--faster--
			s = reduceHilbert(hilbertSeries(M));
--			NumgenM=numgens trim M;
--			GenM=mingens M;
			NumgenM=numgens M;
			GenM=generators M;
			<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s"<<n<<"="<<s<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			n=n+1;
			)
)

--I0list={}
--M0list={}
M0=R
M1=0
i=10;while i<=D and M1==0 do 
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	Tleft = matrix{{t_1..t_Nspin}};
	Tright= matrix{{t_1}..{t_Nspin}};
	I0 = ideal(Tleft*(gamma#i)*Tright);
--	I0list = I0list|{I0};
	if i>=0 then
	 (
		A = koszulComplexDGA(I0);--over the ideal--
		K = toComplex A;
	  M1 = HH_1 K;
		M1=trim M1;
	);
  if M1==0 then
  	(
		M0=M0/I0;--New QuotientRing applied--
--		prune M0;
		M0=trim M0;--faster--
--	M0list=M0list | {M0};
		s0 =reduceHilbert(hilbertSeries(M0));
--		NumgenM=numgens trim module(M0);
--		GenM=mingens module(M0);
		NumgenM=numgens module(M0);
		GenM=generators module(M0);
		<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s0="<< s0 <<endl;
		<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		<<"The 0-th cohomology generators are:"<<GenM<<endl;
		ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s0="<< s0 <<endl;
		ofile<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		ofile<<"The 0-th cohomology generators are:"<<GenM<<endl;
		i=i+1;
		)
	else
		(
			redux(0,0,i);
			s = reduceHilbert(hilbertSeries(M1));
			NumgenM=numgens M1;
			GenM=generators M1;
			<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s1="<<s<<endl;
			<<"Number of 1st cohomology generators="<<NumgenM<<endl;
			<<"The 1st cohomology generators are:"<<GenM<<endl;
			ofile<<"(D="<<(i+1)<<"+"<<(D-i-1)<<")_s1="<<s<<endl;
			ofile<<"Number of 1-st cohomology generators="<<NumgenM<<endl;
			ofile<<"The 1-st cohomology generators are:"<<GenM<<endl;
			redux(2,i,i);
		);
)

i0=0
Tleft = matrix{{t_1..t_Nspin}}
Tright= matrix{{t_1}..{t_Nspin}}
--Mlist={}
I1={}
--I1={Tleft*(gamma#9)*Tright}
--for i from i0+1 to D-1 do
for i from i0 to D-1 do --if i quit at i0+1--
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(i+1)<<"+"<<(D-i-1)<<":"<<endl;
	I1=I1|{Tleft*(gamma#i)*Tright};
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
	redux(0,i,i);
)

ofile<<close
<< get ofilename
