--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

﻿needsPackage "DGAlgebras"

R=QQ[t_1..t_16,MonomialSize=>8];
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
  t_10*t_13-t_9*t_14-t_12*t_15+t_11*t_16
  }};
R=R/ideal(rel);

use R;
A=freeDGAlgebra(R,toList(16:{1,1}))		--NOT 16:{0,1}
A.natural
setDiff(A,apply(toList(1..16),i->(t_i)))
isHomogeneous(A)
apply(maxDegree A+1, i->numgens prune homology(i,A))

HA = homologyAlgebra(A,GenDegreeLimit=>12,RelDegreeLimit=>14)
numgens HA
apply(toList(0..(numgens HA)),i->degree ((generators HA)_i))
reduceHilbert(hilbertSeries(HA))
