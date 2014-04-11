--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=2;
--kk= ZZ/32749;--249999797479451
----ringmap R->S
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}];
S =QQ[psi,Degrees=>{1}];
ofilename="result_gapevs"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

--ev{1,2}:
--{la_1=>3*ps,la_2=>2*ps^2}
--ev{1,3}:
--{la_1=>4*ps,la_2=>3*ps^2}

ev={{la_1=>3*psi,la_2=>2*psi^2,ps=>psi},{la_1=>4*psi,la_2=>3*psi^2,ps=>psi}}
gapseq={{1, 2}, {1, 3}}

Idetpsi={{2*ps^2-3*ps*la_1+la_2,-11*ps^2*la_1+6*ps*la_2,35*ps^2*la_2,3*ps*la_1^2-3*ps*la_2-la_1*la_2,ps^2*la_1^2-ps^2*la_2-ps*la_1*la_2+la_2^2,11*ps^2*la_1^2-11*ps^2*la_2-6*ps*la_1*la_2,-la_1^3+2*la_1*la_2,-3*ps*la_1^3+6*ps*la_1*la_2+la_1^2*la_2-la_2^2,-ps*la_1^3+la_1^4+2*ps*la_1*la_2-3*la_1^2*la_2+la_2^2},{-11*ps^2*la_1+6*ps*la_2,35*ps^2*la_2,ps^2*la_1^2-ps^2*la_2-ps*la_1*la_2+la_2^2,11*ps^2*la_1^2-11*ps^2*la_2-6*ps*la_1*la_2,-la_1^3+2*la_1*la_2,-3*ps*la_1^3+6*ps*la_1*la_2+la_1^2*la_2-la_2^2,-ps*la_1^3+la_1^4+2*ps*la_1*la_2-3*la_1^2*la_2+la_2^2}};

DetM={-11*ps^2*la_1+6*ps*la_2,35*ps^2*la_2,11*ps^2*la_1^2-11*ps^2*la_2-6*ps*la_1*la_2,ps^2*la_1^2-ps^2*la_2-ps*la_1*la_2+la_2^2,-la_1^3+2*la_1*la_2,-3*ps*la_1^3+6*ps*la_1*la_2+la_1^2*la_2-la_2^2,-ps*la_1^3+la_1^4+2*ps*la_1*la_2-3*la_1^2*la_2+la_2^2};
--relationlambda={-la_3^2,la_2^2-2*la_1*la_3,-la_1^2+2*la_2}

--compare with faber's result--
relationsfaber=QQ[kappa]/ideal(kappa^(g-1))
sfaber=reduceHilbert(hilbertSeries relationsfaber)


 ofile<<"g="<<g<<":"<<endl;
--for i from 0 to #ev-1 do
--(
-- f = map(S,R,ev#i);
-- I=kernel f;
-- genI=gens prune I;
-- numgenI=numgens prune I;
-- ofile<<"Considering gap sequence "<<gapseq#i<<":"<<endl;
-- ofile<<"Generators of the kernel of the evaluation map ev_H:"<<endl<<genI<<endl;
-- ofile<<"Number of Generators:"<<numgenI<<endl<<endl<<"**************"<<endl<<endl;
--)

I=trim intersect(apply(#ev-1,i->ideal(kernel map(S,R,ev#i),ideal(Idetpsi#i))))
s=reduceHilbert(hilbertSeries(R/I))
genI=gens prune I
genQ=gens prune (R/I)

I'=trim ideal(DetM)
IsubI'=isSubset(I',I)
s'=reduceHilbert(hilbertSeries(R/I'))
R'=R/I'
genI'=gens prune I'
genQ'=gens prune (R/I')
--intersect(I',I)==I'

ofile<<"The intersection of the kernels of evaluation maps corresponding to all gap sequences -- I=ker ev:"<<endl<<(toExternalString I)<<endl<<endl;
ofile<<"The ideal of k*Omega_mu -- I'=Det(Psi_mu):"<<endl<<(toExternalString I')<<endl<<endl;
ofile<<"Is I' a subset of I: "<<IsubI'<<"!"<<endl<<endl;
ofile<<"Generators of the intersection of all ideals I:"<<endl<<(html genI)<<endl<<endl;
ofile<<"Generators of ideals I'=Det(Psi_mu):"<<endl<<(html genI')<<endl<<endl;
ofile<<"The Hilbert Series of R/I:"<<endl<<(html s)<<endl<<endl;
ofile<<"The Hilbert Series of R/I':"<<endl<<(html s')<<endl<<endl;
ofile<<"The Hilbert Series of Faber's relations - "<<(toExternalString relationsfaber)<<":"<<endl<<(html sfaber)<<endl<<endl;
ofile<<"Generators of R/I:"<<endl<<genQ<<endl<<endl;
ofile<<"Generators of R/I':"<<endl<<genQ'<<endl<<endl;

t2 = cpuTime()
<<"The time takes in this computation process: "<<(t2-t1)<<endl;
ofile<<"The time takes in this computation process: "<<(t2-t1)<<endl;

ofile<<close
<< get ofilename

