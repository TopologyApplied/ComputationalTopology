needsPackage "DGAlgebras"

--R=ZZ/32003[t_1..t_16,s_1..s_16,Degrees=>{16:{1,0},16:{0,1}},MonomialSize=>8];
R=ZZ/32003[t_1..t_16,s_1..s_16,MonomialSize=>8];
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
R=R/ideal(rel);

--Q=QQ[q_1..q_16,SkewCommutative =>apply(toList(1..16),i->q_i)]

--R=R ** Q
--use R;
--I=ideal apply(toList(1..16),i->(t_i,s_i))

A=freeDGAlgebra(R,toList(32:{-1,1})) 		--NOT 16:{0,1}, the degree of {theta,t} in the differential
A.natural
setDiff(A,toList(t_1..t_16)|toList(s_1..s_16))
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
relationHA=ideal HA;
betti oo
degrees relationHA
ofile<<"relations of the generators of cohomology ring:"<<endl<<(html relationHA)<<endl<<endl
ofile<<"betti ideal HA:"<<endl<<(betti relationHA)<<endl
ofile<<"degrees ideal HA:"<<endl<<(degrees relationHA)<<endl

numgens HA
apply(toList(0..((numgens HA)-1)),i->degree ((generators HA)_i))
reduceHilbert(hilbertSeries(HA))
