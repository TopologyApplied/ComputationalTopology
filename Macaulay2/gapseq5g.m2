--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=5;
--kk= ZZ/32749;--249999797479451
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>apply(toList(1..D),i->c_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10}]
ofilename="result_gapseq"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

DetM={67284*ps^5-22449*ps^4*la_1+4536*ps^3*la_2-546*ps^2*la_3+36*ps*la_4-la_5,-269325*ps^5*la_1+63273*ps^4*la_2-9450*ps^3*la_3+870*ps^2*la_4-45*ps*la_5,9485*ps^5-4501*ps^4*la_1+1015*ps^3*la_1^2-280*ps^3*la_2-210*ps^2*la_1*la_2+21*ps*la_2^2+175*ps^2*la_3-11*ps*la_1*la_3-la_2*la_3-10*ps*la_4+la_1*la_4,27335*ps^6-15015*ps^5*la_1+3850*ps^4*la_1^2-749*ps^4*la_2-1050*ps^3*la_1*la_2+140*ps^2*la_2^2+825*ps^3*la_3-55*ps^2*la_1*la_3-15*ps*la_2*la_3+la_3^2-85*ps^2*la_4+15*ps*la_1*la_4-la_2*la_4,27482*ps^5-11760*ps^4*la_1+1960*ps^3*la_1^2-28*ps^3*la_2-322*ps^2*la_1*la_2+154*ps^2*la_3+28*ps*la_1*la_3-22*ps*la_4-la_1*la_4+la_5,105595*ps^6-55468*ps^5*la_1+12831*ps^4*la_1^2-1561*ps^4*la_2-3220*ps^3*la_1*la_2+322*ps^2*la_2^2+2240*ps^3*la_3-42*ps^2*la_1*la_3-28*ps*la_2*la_3-245*ps^2*la_4+18*ps*la_1*la_4+la_2*la_4+10*ps*la_5-la_1*la_5,403704*ps^6-134694*ps^5*la_1+22449*ps^4*la_1^2+4767*ps^4*la_2-4536*ps^3*la_1*la_2+1260*ps^3*la_3+546*ps^2*la_1*la_3-330*ps^2*la_4-36*ps*la_1*la_4+30*ps*la_5+la_1*la_5,665*ps^5-424*ps^4*la_1+195*ps^3*la_1^2-65*ps^2*la_1^3-140*ps^3*la_2+100*ps^2*la_1*la_2+10*ps*la_1^2*la_2-7*ps*la_2^2-la_1*la_2^2-35*ps^2*la_3-13*ps*la_1*la_3+la_1^2*la_3+la_2*la_3+10*ps*la_4-la_1*la_4,575*ps^6-444*ps^5*la_1+239*ps^4*la_1^2-90*ps^3*la_1^3-154*ps^4*la_2+120*ps^3*la_1*la_2+25*ps^2*la_1^2*la_2-14*ps^2*la_2^2-6*ps*la_1*la_2^2+la_2^3-30*ps^3*la_3-36*ps^2*la_1*la_3+6*ps*la_1^2*la_3+6*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2+25*ps^2*la_4-6*ps*la_1*la_4+la_1^2*la_4-la_2*la_4,2037*ps^5-1225*ps^4*la_1+525*ps^3*la_1^2-175*ps^2*la_1^3-378*ps^3*la_2+287*ps^2*la_1*la_2+21*ps*la_1^2*la_2-21*ps*la_2^2-119*ps^2*la_3-18*ps*la_1*la_3-la_1^2*la_3+la_2*la_3+18*ps*la_4+la_1*la_4-la_5,9345*ps^6-6517*ps^5*la_1+3045*ps^4*la_1^2-1015*ps^3*la_1^3-1890*ps^4*la_2+1400*ps^3*la_1*la_2+210*ps^2*la_1^2*la_2-147*ps^2*la_2^2-21*ps*la_1*la_2^2-440*ps^3*la_3-243*ps^2*la_1*la_3+11*ps*la_1^2*la_3+28*ps*la_2*la_3+la_1*la_2*la_3-la_3^2+180*ps^2*la_4-8*ps*la_1*la_4-la_1^2*la_4-10*ps*la_5+la_1*la_5,7987*ps^6-13720*ps^5*la_1+5880*ps^4*la_1^2-1960*ps^3*la_1^3-3626*ps^4*la_2+2954*ps^3*la_1*la_2+322*ps^2*la_1^2*la_2-322*ps^2*la_2^2-1190*ps^3*la_3-238*ps^2*la_1*la_3-28*ps*la_1^2*la_3+28*ps*la_2*la_3+245*ps^2*la_4+25*ps*la_1*la_4+la_1^2*la_4-la_2*la_4-25*ps*la_5-la_1*la_5,14*ps^5-15*ps^4*la_1+15*ps^3*la_1^2-15*ps^2*la_1^3+15*ps*la_1^4-14*ps^3*la_2+29*ps^2*la_1*la_2-44*ps*la_1^2*la_2-la_1^3*la_2+14*ps*la_2^2+2*la_1*la_2^2-14*ps^2*la_3+29*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3-14*ps*la_4-la_1*la_4+la_5,55*ps^6-64*ps^5*la_1+65*ps^4*la_1^2-65*ps^3*la_1^3+65*ps^2*la_1^4-56*ps^4*la_2+120*ps^3*la_1*la_2-185*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+56*ps^2*la_2^2+19*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3-55*ps^3*la_3+119*ps^2*la_1*la_3+11*ps*la_1^2*la_3-la_1^3*la_3-19*ps*la_2*la_3-55*ps^2*la_4-11*ps*la_1*la_4+la_1^2*la_4+la_2*la_4+10*ps*la_5-la_1*la_5,875*ps^6-175*ps^5*la_1+175*ps^4*la_1^2-175*ps^3*la_1^3+175*ps^2*la_1^4-154*ps^4*la_2+329*ps^3*la_1*la_2-504*ps^2*la_1^2*la_2-21*ps*la_1^3*la_2+154*ps^2*la_2^2+42*ps*la_1*la_2^2-155*ps^3*la_3+330*ps^2*la_1*la_3+20*ps*la_1^2*la_3+la_1^3*la_3-41*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2-155*ps^2*la_4-20*ps*la_1*la_4-la_1^2*la_4+la_2*la_4+20*ps*la_5+la_1*la_5,-la_1^5+4*la_1^3*la_2-3*la_1*la_2^2-3*la_1^2*la_3+2*la_2*la_3+2*la_1*la_4-la_5,-15*ps*la_1^5+60*ps*la_1^3*la_2+la_1^4*la_2-45*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3-45*ps*la_1^2*la_3-la_1^3*la_3+30*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2+30*ps*la_1*la_4+la_1^2*la_4-2*la_2*la_4-15*ps*la_5-la_1*la_5,la_1^6-5*la_1^4*la_2+6*la_1^2*la_2^2-la_2^3+4*la_1^3*la_3-6*la_1*la_2*la_3+la_3^2-3*la_1^2*la_4+2*la_2*la_4+2*la_1*la_5}
--DetM={2*ps^2-3*ps*la_1+la_2,la_1^2-la_2,-la_1^3+2*la_1*la_2,-11*ps^2*la_1+6*ps*la_2} --add{3},same resolution

relationlambda={-la_5^2,la_4^2-2*la_3*la_5,-la_3^2+2*la_2*la_4-2*la_1*la_5,la_2^2-2*la_1*la_3+2*la_4,-la_1^2+2*la_2}
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