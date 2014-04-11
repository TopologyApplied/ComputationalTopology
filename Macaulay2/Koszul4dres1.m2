--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"
needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
D=4; Nspin=4;
R=ZZ/32003[t_1..t_Nspin];

gamma={matrix{{0, 0, -1, 0}, {0, 0, 0, -1}, {-1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, 0, 1, 0}, {0, 0, 0, -1}, {1, 0, 0, 0}, {0, -1, 0, 0}},
			matrix{{0, -1, 0, 0}, {-1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}},
			matrix{{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}}


	Tleft = matrix{{t_1..t_Nspin}}; Tright= transpose Tleft;
	I=ideal(Tleft*(gamma#0)*Tright,Tleft*(gamma#1)*Tright);
	A=koszulComplexDGA(I)
	K = toComplex A	 --over the ideal--
	M2 = HH_0 K;
	M2=prune M2
	s2 = reduceHilbert(hilbertSeries(M2))
	NumgenM2=numgens M2
	GenM2=generators M2
	resM2=res M2
	B2=betti(resM2)

--third method--	
	Q0=R/ideal(Tleft*(gamma#0)*Tright)	--New QuotientRing applied--
	Tleft = matrix{{t_1..t_Nspin}}; Tright= transpose Tleft;
	I=ideal(Tleft*(gamma#1)*Tright)
	A=koszulComplexDGA(I)
	K = toComplex A	 --over the ideal--
	M3 = HH_0 K;
	M3=prune M3
	s3 = reduceHilbert(hilbertSeries(M3))
	NumgenM3=numgens M3
	GenM3=generators M3
--	use R;
--	resM3=resolution(map(R,ring M3)**M3)
	resM3=res pushForward(map(ring M3,R), M3)
	B3=betti(resM3)


	Q=Q0/ideal(Tleft*(gamma#1)*Tright)	--New QuotientRing applied--
--  Q=flattenRing(Q,Result=>Thing)
--	resM1 = res (R^1 / ideal (flattenRing(Q,Result=>Thing)))
	resM1=res (R^1 / flattenRing(Q,Result=>Ideal))
	B1=betti(resM1)
	Q=prune Q
	s1 =reduceHilbert(hilbertSeries(Q))
	NumgenM1=numgens Q^1
	GenM1=generators Q^1

