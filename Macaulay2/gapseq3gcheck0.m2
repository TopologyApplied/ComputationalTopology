--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=3;
maxkappa=4;
--kk=QQ;-- 77362086801983
--R=QQ[la_1..la_g,ps]
R=QQ[kappa_1..kappa_(maxkappa),Degrees=>apply(toList(1..maxkappa),i->{i})]
relationsfaber={kappa_1^2,kappa_2,kappa_3,kappa_4}
--relationsoursimp={-(5/9)*kappa_1^2+10*kappa_2,1/12*kappa_1*kappa_2,kappa_1^4/20736}
DetM={-(25/9)*kappa_1^2+50*kappa_2+1/3*(-(kappa_1^3/3456)+kappa_3/120),(85*kappa_1^3)/288-75/4*kappa_1*kappa_2-20*(kappa_1^3/3456-kappa_3/120),245/96*kappa_1^2*kappa_2-175/3*kappa_1*(kappa_1^3/3456-kappa_3/120),kappa_1^3/36+kappa_1^4/82944-1/2*kappa_1*kappa_2+4*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+4*kappa_3,(125*kappa_1^3)/864-35/12*kappa_1*kappa_2+12*(kappa_1^3/3456-kappa_3/120)+1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+50*kappa_3,-((85*kappa_1^4)/10368)+25/96*kappa_1^2*kappa_2+77/9*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-35/6*kappa_1*kappa_3+100*kappa_4,-((85*kappa_1^4)/3456)+155/144*kappa_1^2*kappa_2+25*kappa_1*(kappa_1^3/3456-kappa_3/120)-75/4*kappa_1*kappa_3,1/3*(-(kappa_1^3/3456)+kappa_3/120),kappa_1^4/82944-8*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120),kappa_1^4/6912-kappa_1^5/995328-8/3*kappa_1*(kappa_1^3/3456-kappa_3/120)+1/288*kappa_1^2*(kappa_1^3/3456-kappa_3/120),(5*kappa_1^4)/10368-115/9*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120),-(kappa_1^4/82944)+1/18*kappa_1*(kappa_1^3/3456-kappa_3/120),-(kappa_1^4/3456)+4/3*kappa_1*(kappa_1^3/3456-kappa_3/120),-(kappa_1^4/20736)+kappa_1^5/995328+2/9*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/216*kappa_1^2*(kappa_1^3/3456-kappa_3/120)}
relationlambda={kappa_1^4/82944-1/18*kappa_1*(kappa_1^3/3456-kappa_3/120),-(1/9)*(kappa_1^3/3456-kappa_3/120)^2}


sfaber=reduceHilbert(hilbertSeries(ideal relationsfaber))
soursimp=reduceHilbert(hilbertSeries(ideal (DetM|relationlambda)))
scombined=reduceHilbert(hilbertSeries(ideal (relationsfaber|DetM|relationlambda)))

ofilename="result_gapseq"|g|"gcheck.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ofile<<"g="<<g<<":"<<endl
ofile<<"hilbertSeries with relations from Faber:"<<endl<<sfaber<<endl
ofile<<"hilbertSeries with combined relations from Faber and us:"<<endl<<scombined<<endl
ofile<<"hilbertSeries with relations from us:"<<endl<<soursimp<<endl

ofile<<close
<< get ofilename