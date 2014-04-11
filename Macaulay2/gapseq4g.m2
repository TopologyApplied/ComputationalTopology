--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=4;
kk=QQ;-- ZZ/32749;--77362086801983
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>apply(toList(1..D),i->c_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10}]
ofilename="result_gapseq"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

DetM={1624*ps^4-735*ps^3*la_1+175*ps^2*la_2-21*ps*la_3+la_4,-6769*ps^4*la_1+1960*ps^3*la_2-322*ps^2*la_3+28*ps*la_4,85*ps^4-60*ps^3*la_1+25*ps^2*la_1^2-14*ps^2*la_2-6*ps*la_1*la_2+la_2^2+6*ps*la_3-la_1*la_3,401*ps^4-255*ps^3*la_1+85*ps^2*la_1^2-40*ps^2*la_2-15*ps*la_1*la_2+12*ps*la_3+la_1*la_3-la_4,831*ps^5-661*ps^4*la_1+285*ps^3*la_1^2-120*ps^3*la_2-90*ps^2*la_1*la_2+15*ps*la_2^2+79*ps^2*la_3-9*ps*la_1*la_3-la_2*la_3-6*ps*la_4+la_1*la_4,4872*ps^5-2205*ps^4*la_1+735*ps^3*la_1^2-210*ps^3*la_2-175*ps^2*la_1*la_2+112*ps^2*la_3+21*ps*la_1*la_3-18*ps*la_4-la_1*la_4,9*ps^4-10*ps^3*la_1+10*ps^2*la_1^2-10*ps*la_1^3-9*ps^2*la_2+19*ps*la_1*la_2+la_1^2*la_2-la_2^2-9*ps*la_3-la_1*la_3+la_4,19*ps^5-24*ps^4*la_1+25*ps^3*la_1^2-25*ps^2*la_1^3-20*ps^3*la_2+44*ps^2*la_1*la_2+6*ps*la_1^2*la_2-5*ps*la_2^2-la_1*la_2^2-19*ps^2*la_3-7*ps*la_1*la_3+la_1^2*la_3+la_2*la_3+6*ps*la_4-la_1*la_4,-49*ps^5-85*ps^4*la_1+85*ps^3*la_1^2-85*ps^2*la_1^3-70*ps^3*la_2+155*ps^2*la_1*la_2+15*ps*la_1^2*la_2-15*ps*la_2^2-71*ps^2*la_3-14*ps*la_1*la_3-la_1^2*la_3+la_2*la_3+14*ps*la_4+la_1*la_4,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3-la_4,10*ps*la_1^4-30*ps*la_1^2*la_2-la_1^3*la_2+10*ps*la_2^2+2*la_1*la_2^2+20*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3-10*ps*la_4-la_1*la_4,-la_1^5+4*la_1^3*la_2-3*la_1*la_2^2-3*la_1^2*la_3+2*la_2*la_3+2*la_1*la_4}
--DetM={2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2,-11*ps^2*la_1+6*ps*la_2} --add{3},same resolution

relationlambda={la_4^2,-la_3^2+2*la_2*la_4,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2}
--{
--g=2--{la_2^2,-la_1^2+2*la_2},
--g=3--{-la_3^2,la_2^2-2*la_1*la_3,-la_1^2+2*la_2},
--g=4--{la_4^2,-la_3^2+2*la_2*la_4,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--g=5--{-la_5^2,la_4^2-2*la_3*la_5,-la_3^2+2*la_2*la_4-2*la_1*la_5,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--g=6--{la_6^2,-la_5^2+2*la_4*la_6,la_4^2-2*la_3*la_5+2*la_2*la_6,-la_3^2+2*la_2*la_4-2*la_1*la_5+2*la_6,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2}
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