--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=3;
maxkappa=3*g-1;--g-1
--kk=QQ;-- 77362086801983
--R=QQ[la_1..la_g,ps]
R=QQ[kappa_1..kappa_(maxkappa),Degrees=>apply(toList(1..maxkappa),i->{i})]
relationtails=apply(toList((g-1)..maxkappa),i->kappa_i)
relationsfaber={kappa_1^2,kappa_2}|relationtails
relationsoursimp={-(5/9)*kappa_1^2+10*kappa_2,1/12*kappa_1*kappa_2,kappa_1^4/20736}
relationkappalambda={-(5/72)*kappa_1^2+kappa_2,kappa_1^3/288-1/12*kappa_1*kappa_2-4/3*(kappa_1^3/3456-kappa_3/120)+kappa_3,1/288*kappa_1^2*kappa_2-11/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(121*kappa_1*kappa_3)/1440+kappa_4,-(1/3)*kappa_2*(kappa_1^3/3456-kappa_3/120)+1/288*kappa_1^2*kappa_3+1/4*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/12*kappa_1*kappa_4-4/5*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+kappa_5,-(1/3)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/4*kappa_2*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+1/288*kappa_1^2*kappa_4-1/5*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)-1/12*kappa_1*kappa_5+2/3*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)+kappa_6,1/4*kappa_3*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/3*(kappa_1^3/3456-kappa_3/120)*kappa_4-1/5*kappa_2*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+1/288*kappa_1^2*kappa_5+1/6*kappa_1*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)-1/12*kappa_1*kappa_6-4/7*(-(1/480)*kappa_3*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+(kappa_1^2*kappa_5)/72576+1/72*kappa_1*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)-kappa_7/240)+kappa_7,1/4*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)*kappa_4-1/5*kappa_3*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)-1/3*(kappa_1^3/3456-kappa_3/120)*kappa_5+1/6*kappa_2*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)+1/288*kappa_1^2*kappa_6-1/7*kappa_1*(-(1/480)*kappa_3*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+(kappa_1^2*kappa_5)/72576+1/72*kappa_1*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)-kappa_7/240)-1/12*kappa_1*kappa_7+1/2*(-(1/600)*kappa_3*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+1/756*(kappa_1^3/3456-kappa_3/120)*kappa_5+1/84*kappa_1*(-(1/480)*kappa_3*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+(kappa_1^2*kappa_5)/72576+1/72*kappa_1*(-(1/360)*(kappa_1^3/3456-kappa_3/120)*kappa_3+1/60*kappa_1*(-((kappa_1^2*kappa_3)/34560)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+kappa_5/252)+(kappa_1*kappa_5)/3024)-kappa_7/240)-(kappa_1*kappa_7)/2880)+kappa_8}

sfaber=reduceHilbert(hilbertSeries(ideal relationsfaber))
scombined=reduceHilbert(hilbertSeries(ideal (relationsfaber|relationsoursimp|relationkappalambda)))
soursimp=reduceHilbert(hilbertSeries(ideal (relationsoursimp|relationtails)))
soursimp2=reduceHilbert(hilbertSeries(ideal (relationsoursimp|relationkappalambda|relationtails)))

ofilename="result_gapseq"|g|"gcheck.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ofile<<"g="<<g<<":"<<endl
ofile<<"hilbertSeries with relations from Faber:"<<endl<<sfaber<<endl
ofile<<"hilbertSeries with combined relations from Faber and us:"<<endl<<scombined<<endl
ofile<<"hilbertSeries with relations from us:"<<endl<<soursimp<<endl
ofile<<"hilbertSeries with relations from us+:"<<endl<<soursimp2<<endl

ofile<<close
<< get ofilename