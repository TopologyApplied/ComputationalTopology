--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
R = ZZ/10007[t_1..t_2]
-- R = CC[t_1..t_2]
T3left=matrix{{t_1..t_2}}
T3right=matrix{{t_1}..{t_2}}
-- DiracC1.DiracGamma --
gamma3d1=matrix{{1, 0}, {0, -1}}
gamma3d2=matrix{{1, 0}, {0, 1}}
gamma3d3=matrix{{0, -1}, {-1, 0}}
-- RealGamma --
-- gamma3d1=matrix{{0, 1}, {-1, 0}}
-- gamma3d2=matrix{{0, 1}, {1, 0}}
-- gamma3d3=matrix{{1, 0}, {0, -1}}
epsilon3d = matrix{{0, 1}, {-1, 0}}

I=ideal(T3left*gamma3d1*T3right, T3left*gamma3d2*T3right, T3left*gamma3d3*T3right)
-- I=ideal(T3left*(epsilon3d*gamma3d1)*T3right, T3left*(epsilon3d*gamma3d2)*T3right, T3left*(epsilon3d*gamma3d3)*T3right)
rR = R/I
A = koszulComplexDGA(I)
K = toComplex A
K.dd
M0 = HH_0 K
M1 = HH_1 K
M2 = HH_2 K
prune M0
prune M1
prune M2
s0 = hilbertSeries M0
s1 = hilbertSeries M1
s2 = hilbertSeries M2
reduceHilbert s0
reduceHilbert s1
reduceHilbert s2