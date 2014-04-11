--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

﻿needsPackage "DGAlgebras"

R=QQ[t_1..t_16,MonomialSize=>8];
--R=ZZ/32003[t_1..t_16,s_1..s_16,Degrees=>{16:{1,1},16:{-1,1}},MonomialSize=>8];
 --rel=Simplify[Table[Table[Subscript[t,i],{i,1,16}].gamma10d[n].Table[Subscript[t,i],{i,1,16}],{n,1,10}]/2]--
 --Flatten[Table[{Simplify[(rel[[2*i-1]]+rel[[2*i]])/2],Simplify[(rel[[2*i-1]]-rel[[2*i]])/2]},{i,1,5}]]--
rel={
   t_3*t_9-t_1*t_11-t_7*t_13+t_5*t_15,
  -t_4*t_10+t_2*t_12+t_8*t_14-t_6*t_16,
  -t_4*t_9+t_1*t_12+t_8*t_13-t_5*t_16,
  -t_3*t_10+t_2*t_11+t_7*t_14-t_6*t_15,
  t_2*t_9-t_1*t_10-t_8*t_15+t_7*t_16,
  t_4*t_11-t_3*t_12-t_6*t_13+t_5*t_14,
  -t_2*t_13+t_1*t_14+t_4*t_15-t_3*t_16,
  -t_6*t_9+t_5*t_10+t_8*t_11-t_7*t_12,
  t_2*t_5-t_1*t_6-t_4*t_7+t_3*t_8,
  t_10*t_13-t_9*t_14-t_12*t_15+t_11*t_16
  };
--rel=matrix{rel}
--	R=R/ideal(rel);

--Q=QQ[q_1..q_16,SkewCommutative =>apply(toList(1..16),i->q_i)]

--R=R ** Q
--use R;
--I=ideal(apply(toList(1..16),i->(t_i))|rel)

A=freeDGAlgebra(R,toList(26:{1,-1})) 		--NOT 16:{0,1}, the degree of {theta,t} in the differential
A.natural
setDiff(A,toList(t_1..t_16)|rel)
isHomogeneous(A)
apply(maxDegree A+1, i->numgens prune homology(i,A))
--M=prune HH_5 A 		--if 5 is the maximum non-zero degree above
---s=reduceHilbert(hilbertSeries(M))
---NumgenM=numgens M
--degree generators M 		--the GenDegreeLimit in homologyAlgebra
--relM=relations M 		--the RelDegreeLimit in homologyAlgebra
--degrees relM
--betti relM
--relations coimage relM;
--
HA = homologyAlgebra(A,GenDegreeLimit=>8,RelDegreeLimit=>14) 	--If DGAlgebra has generators in even degrees, then one must specify the options GenDegreeLimit and RelDegreeLimit


-- first differential
d1 = koszul(3,gens I);
-- second differential
d2 = koszul(4,gens I);
-- compute homology of a pair of maps
M = prune HH(d1,d2);
--This avoids having to generate the entire Koszul complex to simply compute the first differentia
			GenM=generators image M.cache.pruningMap
			betti GenM
			s = reduceHilbert(hilbertSeries(M))
			NumgenM=numgens M
relations M
betti oo
relations coimage oo


			resM=resolution(M,LengthLimit=>4)		--,Strategy=>0
			B=betti(resM)
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
--			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl; 
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;


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
	GenM=generators image M.cache.pruningMap;
	s = reduceHilbert(hilbertSeries(M));
	NumgenM=numgens M;
			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<s<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl;
			ofilename="result_brstcoh10dcheck"|n|".txt";
			ofile=ofilename<<""<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<endl<<GenM<<endl;
			ofile<<"The degrees of "<<n<<"-th cohomology generators are:"<<(betti GenM)<<endl;
			t2 = cpuTime();
			<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<close;
	n=n+1;
))

cohM(0,10)