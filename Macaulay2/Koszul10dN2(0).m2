--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

﻿needsPackage "DGAlgebras"

--R=ZZ/32003[t_1..t_16,s_1..s_16,MonomialSize=>8];
--R=ZZ/32003[t_1..t_16,s_1..s_16,Degrees=>{16:{1},16:{-1}},MonomialSize=>8];
R=ZZ/32003[t_1..t_16,s_1..s_16,Degrees=>{16:{1,1},16:{-1,1}},MonomialSize=>8];
 --rel=Simplify[Table[Table[Subscript[t,i],{i,1,16}].gamma10d[n].Table[Subscript[t,i],{i,1,16}],{n,1,10}]/2]--
 --Flatten[Table[{Simplify[(rel[[2*i-1]]+rel[[2*i]])/2],Simplify[(rel[[2*i-1]]-rel[[2*i]])/2]},{i,1,5}]]--
rel=matrix{{
		s_12*t_2-s_10*t_4-s_16*t_6+s_14*t_8-s_4*t_10+s_2*t_12+s_8*t_14-s_6*t_16,
		-s_11*t_1+s_9*t_3+s_15*t_5-s_13*t_7+s_3*t_9-s_1*t_11-s_7*t_13+s_5*t_15,
		s_11*t_2-s_10*t_3-s_15*t_6+s_14*t_7-s_3*t_10+s_2*t_11+s_7*t_14-s_6*t_15,
		s_12*t_1-s_9*t_4-s_16*t_5+s_13*t_8-s_4*t_9+s_1*t_12+s_8*t_13-s_5*t_16,
		-s_12*t_3+s_11*t_4+s_14*t_5-s_13*t_6+s_4*t_11-s_3*t_12-s_6*t_13+s_5*t_14,
		-s_10*t_1+s_9*t_2+s_16*t_7-s_15*t_8+s_2*t_9-s_1*t_10-s_8*t_15+s_7*t_16,
		s_10*t_5-s_9*t_6-s_12*t_7+s_11*t_8-s_6*t_9+s_5*t_10+s_8*t_11-s_7*t_12,
		s_14*t_1-s_13*t_2-s_16*t_3+s_15*t_4-s_2*t_13+s_1*t_14+s_4*t_15-s_3*t_16,
		-s_14*t_9+s_13*t_10+s_16*t_11-s_15*t_12+s_10*t_13-s_9*t_14-s_12*t_15+s_11*t_16,
		-s_6*t_1+s_5*t_2+s_8*t_3-s_7*t_4+s_2*t_5-s_1*t_6-s_4*t_7+s_3*t_8
  }};
I=ideal(rel);

n=1;
-- first differential
d1 = koszul(n,gens I);
-- second differential
d2 = koszul(n+1,gens I);
-- compute homology of a pair of maps
M = prune HH(d1,d2);
--This avoids having to generate the entire Koszul complex to simply compute the first differentia
			s = reduceHilbert(hilbertSeries(M))
--			GenM=generators image M.cache.pruningMap
--			NumgenM=numgens M
			ofilename="result_Koszul10d2N"|n|".txt";
			ofile=ofilename<<""<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			ofile<<close;

resM=res M
B=betti(resM)
apply(toList(0..9),i->tally sort degrees resM#i)


cohM = (n1,n2)->(
n:=n1;M=1;while n<=n2 and M!=0 do
	(
	t1 = cpuTime();
		-- first differential
	d1 = koszul(n,gens I);
	-- second differential
	d2 = koszul(n+1,gens I);
	-- compute homology of a pair of maps
	M = prune HH(d1,d2);
	--This avoids having to generate the entire Koszul complex to simply compute the first differentia
	s = reduceHilbert(hilbertSeries(M));
--	GenM=generators image M.cache.pruningMap;
--	NumgenM=numgens M;
			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<s<<endl;
--			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			ofilename="result_brstcoh10d2N"|n|".txt";
			ofile=ofilename<<""<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			t2 = cpuTime();
			<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<close;
	n=n+1;
))

cohM(2,10)
