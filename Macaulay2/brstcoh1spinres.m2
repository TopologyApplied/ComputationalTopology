ofilename="result_brstcoh10d1spin_res.txt"
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

use R;
I=ideal apply(toList(1..16),i->(t_i))
K = toComplex koszulComplexDGA(I)
			M = HH_1 K
			M=prune M
			relations M
			relations coimage oo
			
			
			resM=res M
			B=betti(resM)
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl

t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename