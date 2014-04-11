--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
R = ZZ/10007[t_1..t_4]
T4left=matrix{{t_1..t_4}}
T4right=matrix{{t_1}..{t_4}}
-- gamma4d1=matrix{{0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
-- gamma4d2=matrix{{0, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
-- gamma4d3=matrix{{0, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}
-- gamma4d4=matrix{{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}}
-- DiracC2.DiracGamma --
gamma4d1=matrix{{0, -1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, -1}, {0, 0, -1, 0}}
gamma4d2=matrix{{0, -1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, -1}, {0, 0, -1, 0}}
gamma4d3=matrix{{0, 1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, -1, 0}}
gamma4d4=matrix{{0, 0, -1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, -1, 0, 0}}

I=ideal(T4left*gamma4d1*T4right, T4left*gamma4d2*T4right, T4left*gamma4d3*T4right, T4left*gamma4d4*T4right)
rR = R/I
A = koszulComplexDGA(I)
K = toComplex A
K.dd
M0 = HH_0 K
M1 = HH_1 K
M2 = HH_2 K
M3 = HH_3 K
prune M0
prune M1
prune M2
prune M3
s0 = hilbertSeries M0
s1 = hilbertSeries M1
s2 = hilbertSeries M2
s3 = hilbertSeries M3
reduceHilbert s0
reduceHilbert s1
reduceHilbert s2
reduceHilbert s3