--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
needsPackage "BoijSoederberg" 
needsPackage "ChainComplexExtras"
R=ZZ/32003[t]
I=ideal(t^2)
A = koszulComplexDGA(I) 
K = toComplex A
M0 = HH_0 K
M0=prune M0
s0 = reduceHilbert(hilbertSeries(M0))
M1 = HH_1 K
M1=prune M1
s1 = reduceHilbert(hilbertSeries(M1))

Mh0 = HH^0 K
Mh0=prune Mh0
sh0 = reduceHilbert(hilbertSeries(Mh0))
Mh1 = HH^1 K
Mh1=prune Mh1
sh1 = reduceHilbert(hilbertSeries(Mh1))

Rc=R[u]/(t*u-1)
f=map(Rc,R)
Mc0=tensor(Rc,f,M0)
Mc0=prune Mc0
Mc1=tensor(Rc,f,M1)
Mc1=prune Mc1


Rc=R[u]/(t*u-1)
Kc=substitute(K,Rc)
Mc0 = HH_0 Kc
Mc0=prune Mc0
sc0 = reduceHilbert(hilbertSeries(Mc0))
Mc1 = HH_1 Kc
Mc1=prune Mc1
sc1 = reduceHilbert(hilbertSeries(Mc1))

gens Rc
f=map(R,Rc,{t})
i->map(R^i,Rc^i,f,{{t}})
map(K,Kc,i->map(R^i,Rc^i,f,t*id_(R^i)))




Rm=Rc/(t)
Km=substitute(K,Rm)
Mm0 = HH_0 Km
Mm0=prune Mm0
sm0 = reduceHilbert(hilbertSeries(Mm0))
Mm1 = HH_1 Km
Mm1=prune Mm1
sm1 = reduceHilbert(hilbertSeries(Mm1))
