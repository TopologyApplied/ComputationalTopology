--Author: Renjun Xu--
--E-mail: rxu@ucdavis.edu--
--Please feel free to contact me--

ofilename="result_brstcoh10d2spin.txt"
ofile=ofilename<<""<<endl
t1 = cpuTime()

needsPackage "DGAlgebras"

R=QQ[t_1..t_16,s_1..s_16,MonomialSize=>8];
 rel=matrix{{t_16*t_1-t_15*t_4-t_8*t_12+t_9*t_11,
  t_16*t_2-t_15*t_5+t_7*t_12-t_9*t_10,
  t_16*t_3-t_15*t_6-t_7*t_11+t_8*t_10,
  t_14*t_7-t_13*t_10-t_2*t_6+t_3*t_5,
  t_14*t_8-t_13*t_11+t_1*t_6-t_3*t_4,
  t_14*t_9-t_13*t_12-t_1*t_5+t_2*t_4,
  t_10*t_1+t_11*t_2+t_12*t_3+t_14*t_15,
  t_7*t_1+t_8*t_2+t_9*t_3+t_13*t_15,
  t_7*t_4+t_8*t_5+t_9*t_6+t_13*t_16,
  t_10*t_4+t_11*t_5+t_12*t_6+t_14*t_16,
  s_16*s_1-s_15*s_4-s_8*s_12+s_9*s_11,
  s_16*s_2-s_15*s_5+s_7*s_12-s_9*s_10,
  s_16*s_3-s_15*s_6-s_7*s_11+s_8*s_10,
  s_14*s_7-s_13*s_10-s_2*s_6+s_3*s_5,
  s_14*s_8-s_13*s_11+s_1*s_6-s_3*s_4,
  s_14*s_9-s_13*s_12-s_1*s_5+s_2*s_4,
  s_10*s_1+s_11*s_2+s_12*s_3+s_14*s_15,
  s_7*s_1+s_8*s_2+s_9*s_3+s_13*s_15,
  s_7*s_4+s_8*s_5+s_9*s_6+s_13*s_16,
  s_10*s_4+s_11*s_5+s_12*s_6+s_14*s_16
  }};
R=R/ideal(rel);

use R;
I=ideal apply(toList(1..16),i->(t_i+s_i))
-- first differential
d1 = koszul(2,gens I);
-- second differential
d2 = koszul(3,gens I);
-- compute homology of a pair of maps
M = prune HH(d1,d2);
--This avoids having to generate the entire Koszul complex to simply compute the first differentia
--			GenM=generators image M.cache.pruningMap;
--			s = reduceHilbert(hilbertSeries(M));
--			NumgenM=numgens M;
--relations M;
--betti oo
--relations coimage o ;

			resM=resolution(M,LengthLimit=>4)		--,Strategy=>0
			B=betti(resM)
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
--			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl; 
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;


--K = toComplex koszulComplexDGA(I)
			M = HH_1 K;
			M=prune M;
			relations M
			relations coimage oo
			
			resM=res M
			B=betti(resM)


t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename