--Author: Renjun Xu--
--E-mail: rxu@ucdavis.edu--
--Advisor: Albert Schwarz--
--Please feel free to contact me--

ofilename="result_brstcoh10d2spin.txt"
ofile=ofilename<<""<<endl
t1 = cpuTime()

needsPackage "DGAlgebras"

R=QQ[t_1..t_16,s_1..s_16];
 --rel=Simplify[Table[Table[Subscript[t,i],{i,1,16}].gamma10d[n].Table[Subscript[t,i],{i,1,16}],{n,1,10}]/2]--
 --Flatten[Table[{Simplify[(rel[[2*i-1]]+rel[[2*i]])/2],Simplify[(rel[[2*i-1]]-rel[[2*i]])/2]},{i,1,5}]]--
rel=matrix{{
   t_3*t_9-t_1*t_11-t_7*t_13+t_5*t_15,
  -t_4*t_10+t_2*t_12+t_8*t_14-t_6*t_16,
  -t_4*t_9+t_1*t_12+t_8*t_13-t_5*t_16,
  -t_3*t_10+t_2*t_11+t_7*t_14-t_6*t_15,
  t_2*t_9-t_1*t_10-t_8*t_15+t_7*t_16,
  t_4*t_11-t_3*t_12-t_6*t_13+t_5*t_14,
  -t_2*t_13+t_1*t_14+t_4*t_15-t_3*t_16,
  -t_6*t_9+t_5*t_10+t_8*t_11-t_7*t_12,
  t_2*t_5-t_1*t_6-t_4*t_7+t_3*t_8,
  t_10*t_13-t_9*t_14-t_12*t_15+t_11*t_16,
  s_3*s_9-s_1*s_11-s_7*s_13+s_5*s_15,
  -s_4*s_10+s_2*s_12+s_8*s_14-s_6*s_16,
  -s_4*s_9+s_1*s_12+s_8*s_13-s_5*s_16,
  -s_3*s_10+s_2*s_11+s_7*s_14-s_6*s_15,
  s_2*s_9-s_1*s_10-s_8*s_15+s_7*s_16,
  s_4*s_11-s_3*s_12-s_6*s_13+s_5*s_14,
  -s_2*s_13+s_1*s_14+s_4*s_15-s_3*s_16,
  -s_6*s_9+s_5*s_10+s_8*s_11-s_7*s_12,
  s_2*s_5-s_1*s_6-s_4*s_7+s_3*s_8,
  s_10*s_13-s_9*s_14-s_12*s_15+s_11*s_16
  }};
-- rel=matrix{{t_16*t_1-t_15*t_4-t_8*t_12+t_9*t_11,
--  t_16*t_2-t_15*t_5+t_7*t_12-t_9*t_10,
--  t_16*t_3-t_15*t_6-t_7*t_11+t_8*t_10,
--  t_14*t_7-t_13*t_10-t_2*t_6+t_3*t_5,
--  t_14*t_8-t_13*t_11+t_1*t_6-t_3*t_4,
--  t_14*t_9-t_13*t_12-t_1*t_5+t_2*t_4,
--  t_10*t_1+t_11*t_2+t_12*t_3+t_14*t_15,
--  t_7*t_1+t_8*t_2+t_9*t_3+t_13*t_15,
--  t_7*t_4+t_8*t_5+t_9*t_6+t_13*t_16,
--  t_10*t_4+t_11*t_5+t_12*t_6+t_14*t_16,
--  s_16*s_1-s_15*s_4-s_8*s_12+s_9*s_11,
--  s_16*s_2-s_15*s_5+s_7*s_12-s_9*s_10,
--  s_16*s_3-s_15*s_6-s_7*s_11+s_8*s_10,
--  s_14*s_7-s_13*s_10-s_2*s_6+s_3*s_5,
--  s_14*s_8-s_13*s_11+s_1*s_6-s_3*s_4,
--  s_14*s_9-s_13*s_12-s_1*s_5+s_2*s_4,
--  s_10*s_1+s_11*s_2+s_12*s_3+s_14*s_15,
--  s_7*s_1+s_8*s_2+s_9*s_3+s_13*s_15,
--  s_7*s_4+s_8*s_5+s_9*s_6+s_13*s_16,
--  s_10*s_4+s_11*s_5+s_12*s_6+s_14*s_16
--  }};
R=R/ideal(rel);

--Q=QQ[q_1..q_16,SkewCommutative =>apply(toList(1..16),i->q_i)]

--R=R ** Q
use R;
I=ideal apply(toList(1..16),i->(t_i+s_i))
A=koszulComplexDGA(I)
--A=freeDGAlgebra(R,toList(16:{1,1}))		--NOT 16:{0,1}
--A.natural
--setDiff(A,apply(toList(1..16),i->(t_i+s_i)))
--isHomogeneous(A)
--
--M=HH_0 A
--M=prune M
--s=reduceHilbert(hilbertSeries(M))
--NumgenM=numgens M
--
--HA = homologyAlgebra(A)

K = toComplex A;
n=0;M=1;while n<=5 and M!=0 do
	(
			M = HH_n K;
			M=prune M;--also improve the ambient free module --
---		M=trim M;--faster--
--			GenM=generators image M.cache.pruningMap;
--			s = reduceHilbert(hilbertSeries(M));
--			NumgenM=numgens M;
			--resM=res M;
			resM=res(M,LengthLimit=>12,SyzygyLimit=>12,Strategy=>0);
			B=betti(resM);
--			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl;
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
--			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl; 
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			n=n+1;
	)



t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename