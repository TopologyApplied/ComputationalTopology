--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=2;
kk=QQ;-- ZZ/32003;--77362086801983
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>apply(toList(1..D),i->c_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10}]
ofilename="result_gapseq"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

DetM={2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2}
--DetM={2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2,-11*ps^2*la_1+6*ps*la_2} --add{3},same resolution
--{
--g=2--{2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2}s
relationlambda={la_2^2,-la_1^2+2*la_2}
--{
--{la_2^2,-la_1^2+2*la_2},
--{-la_3^2,la_2^2-2*la_1*la_3,-la_1^2+2*la_2},
--{la_4^2,-la_3^2+2*la_2*la_4,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--{-la_5^2,la_4^2-2*la_3*la_5,-la_3^2+2*la_2*la_4-2*la_1*la_5,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--{la_6^2,-la_5^2+2*la_4*la_6,la_4^2-2*la_3*la_5+2*la_2*la_6,-la_3^2+2*la_2*la_4-2*la_1*la_5+2*la_6,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2}
--}
I=ideal(DetM|relationlambda)
resI=resolution(I)
--resI=resolution(I,LengthLimit=>2)
B=betti(resI)
I=prune I
genI=mingens I
ofile<<"g="<<g<<":"<<endl
ofile<<"Simplified relations:"<<endl<<genI<<endl
ofile<<resI<<endl
ofile<<resI.dd<<endl
ofile<<B<<endl


ofile<<close
<< get ofilename