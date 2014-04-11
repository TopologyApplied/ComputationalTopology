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
maxkappa=4;
--kk=QQ;-- 77362086801983
--R=QQ[la_1..la_g,ps]
R=QQ[kappa_1..kappa_(maxkappa),Degrees=>apply(toList(1..maxkappa),i->{i})]
relationsfaber={kappa_1^2-32/3*kappa_2,kappa_1^(g-1)}|apply(toList((g-1)..maxkappa),i->kappa_i)
--relationsoursimp={-(5/9)*kappa_1^2+10*kappa_2,1/12*kappa_1*kappa_2,kappa_1^4/20736}
DetM={(175*kappa_1^3)/288-245/4*kappa_1*kappa_2-42*(kappa_1^3/3456-kappa_3/120)+1624*kappa_3+1/4*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),245/36*kappa_1^2*kappa_2-322/3*kappa_1*(kappa_1^3/3456-kappa_3/120)-6769/12*kappa_1*kappa_3+42*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),(11*kappa_1^3)/96+kappa_1^4/82944-5*kappa_1*kappa_2+12*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+85*kappa_3,(245*kappa_1^3)/576-85/4*kappa_1*kappa_2+24*(kappa_1^3/3456-kappa_3/120)+1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+401*kappa_3+1/4*(-(1/36)*kappa_1*(kappa_1^3/3456-kappa_3/120)+(kappa_1*kappa_3)/1440),-((115*kappa_1^4)/4608)+25/16*kappa_1^2*kappa_2+149/6*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-661/12*kappa_1*kappa_3-9*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+831*kappa_4,-((175*kappa_1^4)/3456)+35/8*kappa_1^2*kappa_2+245/6*kappa_1*(kappa_1^3/3456-kappa_3/120)-735/4*kappa_1*kappa_3-27*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+4872*kappa_4,(7*kappa_1^3)/192+kappa_1^4/82944-5/6*kappa_1*kappa_2-18*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+9*kappa_3+1/4*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),-((17*kappa_1^4)/13824)-kappa_1^5/995328+5/48*kappa_1^2*kappa_2-15/2*kappa_1*(kappa_1^3/3456-kappa_3/120)+1/288*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-2*kappa_1*kappa_3+9*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+19*kappa_4,-((5*kappa_1^4)/1536)+25/72*kappa_1^2*kappa_2-26*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-85/12*kappa_1*kappa_3+21*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-49*kappa_4,-(kappa_1^4/82944)+1/18*kappa_1*(kappa_1^3/3456-kappa_3/120)+1/4*(-(1/36)*kappa_1*(kappa_1^3/3456-kappa_3/120)+(kappa_1*kappa_3)/1440),-((5*kappa_1^4)/6912)+10/3*kappa_1*(kappa_1^3/3456-kappa_3/120)-15*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),kappa_1^5/995328-1/216*kappa_1^2*(kappa_1^3/3456-kappa_3/120)+1/24*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)}
relationlambda={kappa_1^4/82944-1/18*kappa_1*(kappa_1^3/3456-kappa_3/120)+1/2*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),-(1/9)*(kappa_1^3/3456-kappa_3/120)^2+1/576*kappa_1^2*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),1/16*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)^2}
--relationtest={(161*kappa_1^3)/144-490/3*kappa_1*kappa_2-56*(kappa_1^3/3456-kappa_3/120)+6769*kappa_3+1/4*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),63/4*kappa_1^2*kappa_2-182*kappa_1*(kappa_1^3/3456-kappa_3/120)-7483/4*kappa_1*kappa_3+54*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),(95*kappa_1^3)/288+kappa_1^4/82944-25*kappa_1*kappa_2+20*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+725*kappa_3,(581*kappa_1^3)/576-175/2*kappa_1*kappa_2+30*(kappa_1^3/3456-kappa_3/120)+1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+2786*kappa_3+1/4*(-(1/36)*kappa_1*(kappa_1^3/3456-kappa_3/120)+(kappa_1*kappa_3)/1440),-((91*kappa_1^4)/1536)+875/144*kappa_1^2*kappa_2+113/2*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-4501/12*kappa_1*kappa_3-15*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+9485*kappa_4,-((161*kappa_1^4)/1728)+973/72*kappa_1^2*kappa_2+56*kappa_1*(kappa_1^3/3456-kappa_3/120)-980*kappa_1*kappa_3-33*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+40614*kappa_4,(101*kappa_1^3)/576+kappa_1^4/82944-35/4*kappa_1*kappa_2-24*(kappa_1^3/3456-kappa_3/120)-1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)+194*kappa_3+1/4*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440),-((107*kappa_1^4)/13824)-kappa_1^5/995328+125/144*kappa_1^2*kappa_2-83/6*kappa_1*(kappa_1^3/3456-kappa_3/120)+1/288*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-106/3*kappa_1*kappa_3+15*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+665*kappa_4,-((77*kappa_1^4)/4608)+7/3*kappa_1^2*kappa_2-128/3*kappa_1*(kappa_1^3/3456-kappa_3/120)-1/864*kappa_1^2*(kappa_1^3/3456-kappa_3/120)-1225/12*kappa_1*kappa_3+27*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+273*kappa_4,kappa_1^3/288-kappa_1^4/82944-1/12*kappa_1*kappa_2-2*(kappa_1^3/3456-kappa_3/120)+1/18*kappa_1*(kappa_1^3/3456-kappa_3/120)+kappa_3+1/4*(-(1/36)*kappa_1*(kappa_1^3/3456-kappa_3/120)+(kappa_1*kappa_3)/1440),-(1/768)*kappa_1^4+1/18*kappa_1^2*kappa_2+1/6*kappa_1*(kappa_1^3/3456-kappa_3/120)-5/4*kappa_1*kappa_3-21*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)-1/48*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)+134*kappa_4,kappa_1^5/995328-1/216*kappa_1^2*(kappa_1^3/3456-kappa_3/120)+1/24*kappa_1*(1/36*kappa_1*(kappa_1^3/3456-kappa_3/120)-(kappa_1*kappa_3)/1440)}

sfaber=reduceHilbert(hilbertSeries(ideal relationsfaber))
soursimp=reduceHilbert(hilbertSeries(ideal (DetM|relationlambda)))
scombined=reduceHilbert(hilbertSeries(ideal (relationsfaber|DetM|relationlambda)))
--scombined2=reduceHilbert(hilbertSeries(ideal (relationsfaber|DetM)))
--scombined3=reduceHilbert(hilbertSeries(ideal (relationsfaber|relationlambda)))
--stest=reduceHilbert(hilbertSeries(ideal (relationsfaber|relationtest)))

ofilename="result_gapseq"|g|"gcheck.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ofile<<"g="<<g<<":"<<endl
ofile<<"hilbertSeries with relations from Faber:"<<endl<<sfaber<<endl
ofile<<"hilbertSeries with relations from us:"<<endl<<soursimp<<endl
ofile<<"hilbertSeries with combined relations from Faber and us:"<<endl<<scombined<<endl
--ofile<<"hilbertSeries with combined relations from Faber and test:"<<endl<<stest<<endl

ofile<<close
<< get ofilename