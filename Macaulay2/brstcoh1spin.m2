--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

ofilename="result_brstcoh10d1spin.txt"
ofile=ofilename<<""<<endl
t1 = cpuTime()

needsPackage "DGAlgebras"

R=QQ[t_1..t_16];
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

--Q=QQ[q_1..q_16,SkewCommutative =>apply(toList(1..16),i->q_i)]

--R=R ** Q
use R;
I=ideal apply(toList(1..16),i->(t_i))
K = toComplex koszulComplexDGA(I)
n=0;M=1;while n<=5 and M!=0 do
	(
			M = HH_n K;
			M=prune M;--also improve the ambient free module --
--		M=trim M;--faster--
			GenM=generators image M.cache.pruningMap;
			s = reduceHilbert(hilbertSeries(M));
			NumgenM=numgens M;
			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl; 
			n=n+1;
	)



t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename