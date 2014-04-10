--Author: Renjun Xu--
--E-mail: rxu@ucdavis.edu--
--Advisor: Albert Schwarz--
--Please feel free to contact me--

ofilename="result_brstcoh10d2spinsym.txt"
ofile=ofilename<<""<<endl
t1 = cpuTime()

needsPackage "DGAlgebras"

R=QQ[u_1..u_16,v_1..v_16,SkewCommutative =>apply(toList(1..16),i->v_i)];
 rel=matrix{{(u_16+v_16)*(u_1+v_1)-(u_15+v_15)*(u_4+v_4)-(u_8+v_8)*(u_12+v_12)+(u_9+v_9)*(u_11+v_11),
  (u_16+v_16)*(u_2+v_2)-(u_15+v_15)*(u_5+v_5)+(u_7+v_7)*(u_12+v_12)-(u_9+v_9)*(u_10+v_10),
  (u_16+v_16)*(u_3+v_3)-(u_15+v_15)*(u_6+v_6)-(u_7+v_7)*(u_11+v_11)+(u_8+v_8)*(u_10+v_10),
  (u_14+v_14)*(u_7+v_7)-(u_13+v_13)*(u_10+v_10)-(u_2+v_2)*(u_6+v_6)+(u_3+v_3)*(u_5+v_5),
  (u_14+v_14)*(u_8+v_8)-(u_13+v_13)*(u_11+v_11)+(u_1+v_1)*(u_6+v_6)-(u_3+v_3)*(u_4+v_4),
  (u_14+v_14)*(u_9+v_9)-(u_13+v_13)*(u_12+v_12)-(u_1+v_1)*(u_5+v_5)+(u_2+v_2)*(u_4+v_4),
  (u_10+v_10)*(u_1+v_1)+(u_11+v_11)*(u_2+v_2)+(u_12+v_12)*(u_3+v_3)+(u_14+v_14)*(u_15+v_15),
  (u_7+v_7)*(u_1+v_1)+(u_8+v_8)*(u_2+v_2)+(u_9+v_9)*(u_3+v_3)+(u_13+v_13)*(u_15+v_15),
  (u_7+v_7)*(u_4+v_4)+(u_8+v_8)*(u_5+v_5)+(u_9+v_9)*(u_6+v_6)+(u_13+v_13)*(u_16+v_16),
  (u_10+v_10)*(u_4+v_4)+(u_11+v_11)*(u_5+v_5)+(u_12+v_12)*(u_6+v_6)+(u_14+v_14)*(u_16+v_16),
(u_16-v_16)*(u_1-v_1)-(u_15-v_15)*(u_4-v_4)-(u_8-v_8)*(u_12-v_12)+(u_9-v_9)*(u_11-v_11),
  (u_16-v_16)*(u_2-v_2)-(u_15-v_15)*(u_5-v_5)+(u_7-v_7)*(u_12-v_12)-(u_9-v_9)*(u_10-v_10),
  (u_16-v_16)*(u_3-v_3)-(u_15-v_15)*(u_6-v_6)-(u_7-v_7)*(u_11-v_11)+(u_8-v_8)*(u_10-v_10),
  (u_14-v_14)*(u_7-v_7)-(u_13-v_13)*(u_10-v_10)-(u_2-v_2)*(u_6-v_6)+(u_3-v_3)*(u_5-v_5),
  (u_14-v_14)*(u_8-v_8)-(u_13-v_13)*(u_11-v_11)+(u_1-v_1)*(u_6-v_6)-(u_3-v_3)*(u_4-v_4),
  (u_14-v_14)*(u_9-v_9)-(u_13-v_13)*(u_12-v_12)-(u_1-v_1)*(u_5-v_5)+(u_2-v_2)*(u_4-v_4),
  (u_10-v_10)*(u_1-v_1)+(u_11-v_11)*(u_2-v_2)+(u_12-v_12)*(u_3-v_3)+(u_14-v_14)*(u_15-v_15),
  (u_7-v_7)*(u_1-v_1)+(u_8-v_8)*(u_2-v_2)+(u_9-v_9)*(u_3-v_3)+(u_13-v_13)*(u_15-v_15),
  (u_7-v_7)*(u_4-v_4)+(u_8-v_8)*(u_5-v_5)+(u_9-v_9)*(u_6-v_6)+(u_13-v_13)*(u_16-v_16),
  (u_10-v_10)*(u_4-v_4)+(u_11-v_11)*(u_5-v_5)+(u_12-v_12)*(u_6-v_6)+(u_14-v_14)*(u_16-v_16)
  }};
R=R/ideal(rel);

--Q=QQ[q_1..q_16,SkewCommutative =>apply(toList(1..16),i->q_i)]

--R=R ** Q
use R;
I=ideal apply(toList(1..16),i->(u_i))
K = toComplex koszulComplexDGA(I)
n=0;M=1;while n<=5 and M!=0 do
	(
			M = HH_n K;
			M=prune M;--also improve the ambient free module --
---		M=trim M;--faster--
			GenM=generators image M.cache.pruningMap;
			s = reduceHilbert(hilbertSeries(M));
			NumgenM=numgens M;
			resM=res M;
			B=betti(resM);
			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl;
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl; 
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			n=n+1;
	)



t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename