needsPackage "DGAlgebras"
needsPackage "BoijSoederberg"
R=ZZ/32003[t,u]/(t*u-1)
I=ideal(t^2)
A = koszulComplexDGA(I)
K = toComplex A
M0=trim HH_0 K
reduceHilbert(hilbertSeries(M0))
M0=trim HH_0 K
reduceHilbert(hilbertSeries(M0))



A=ZZ[t]
B=ZZ[C,SkewCommutative=>{C}]
R=tensor(A,B)
use R
d=t^2*C
