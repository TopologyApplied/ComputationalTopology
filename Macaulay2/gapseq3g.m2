--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=3;
--kk=QQ;-- 77362086801983
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>apply(toList(1..D),i->c_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10}]
ofilename="result_gapseq"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

--relations from degreelimit Not add more simplified relations--
degreelimit={la_1*ps^2*la_2^2,la_1*ps*la_2*la_3,la_1^2*ps^3*la_2,la_1^2*ps^2*la_3,la_1^2*ps*la_2^2,la_1^3*ps^4,la_1^3*ps^2*la_2,la_1^3*ps*la_3,la_1*la_2^3,la_1*la_2^2*la_3,la_1*la_2*la_3^2,la_1^2*la_2*la_3,la_1^3*la_2^2,la_1*la_3^2}
DetM={50*ps^3-35*ps^2*la_1+10*ps*la_2-la_3,-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,4*ps^4-6*ps^3*la_1+7*ps^2*la_1^2-5*ps^2*la_2-3*ps*la_1*la_2+la_2^2+3*ps*la_3-la_1*la_3,50*ps^4-35*ps^3*la_1+35*ps^2*la_1^2-25*ps^2*la_2-10*ps*la_1*la_2+9*ps*la_3+la_1*la_3,-la_1^3+2*la_1*la_2-la_3,-6*ps*la_1^3+12*ps*la_1*la_2+la_1^2*la_2-la_2^2-6*ps*la_3-la_1*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3}
--sumLamda[g,g+2]--{50*ps^3-35*ps^2*la_1+10*ps*la_2-la_3,-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,4*ps^4-6*ps^3*la_1+7*ps^2*la_1^2-5*ps^2*la_2-3*ps*la_1*la_2+la_2^2+3*ps*la_3-la_1*la_3,50*ps^4-35*ps^3*la_1+35*ps^2*la_1^2-25*ps^2*la_2-10*ps*la_1*la_2+9*ps*la_3+la_1*la_3,100*ps^5-70*ps^4*la_1+55*ps^3*la_1^2-35*ps^3*la_2-30*ps^2*la_1*la_2+10*ps*la_2^2+28*ps^2*la_3-7*ps*la_1*la_3-la_2*la_3,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-la_1^3+2*la_1*la_2-la_3,-6*ps*la_1^3+12*ps*la_1*la_2+la_1^2*la_2-la_2^2-6*ps*la_3-la_1*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3}
--DetM={2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2,-11*ps^2*la_1+6*ps*la_2} --add{3},same resolution

relationlambda={-la_3^2,la_2^2-2*la_1*la_3,-la_1^2+2*la_2}
--{
--{la_2^2,-la_1^2+2*la_2},
--{-la_3^2,la_2^2-2*la_1*la_3,-la_1^2+2*la_2},
--{la_4^2,-la_3^2+2*la_2*la_4,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--{-la_5^2,la_4^2-2*la_3*la_5,-la_3^2+2*la_2*la_4-2*la_1*la_5,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2},
--{la_6^2,-la_5^2+2*la_4*la_6,la_4^2-2*la_3*la_5+2*la_2*la_6,-la_3^2+2*la_2*la_4-2*la_1*la_5+2*la_6,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2}
--}
I=ideal(DetM|relationlambda|degreelimit)
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