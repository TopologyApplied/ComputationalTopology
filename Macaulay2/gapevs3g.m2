--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=3;
--kk= ZZ/32749;--249999797479451
----ringmap R->S
R=QQ[la_1..la_g,ps,Degrees=>apply(toList(1..g),i->{i})|{{1}}];
S =QQ[psi,Degrees=>{1}];
ofilename="result_gapevs"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ev={{la_1=>6*psi,la_2=>11*psi^2,la_3=>6*psi^3,ps=>psi},{la_1=>7*psi,la_2=>14*psi^2,la_3=>8*psi^3,ps=>psi},{la_1=>8*psi,la_2=>17*psi^2,la_3=>10*psi^3,ps=>psi},{la_1=>9*psi,la_2=>23*psi^2,la_3=>15*psi^3,ps=>psi}}
gapseq={{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 3, 5}}
--Idetpsi={{-la_1^3+2*la_1*la_2-la_3,-6*ps*la_1^3+12*ps*la_1*la_2+la_1^2*la_2-la_2^2-6*ps*la_3-la_1*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2},{-6*ps*la_1^3+12*ps*la_1*la_2+la_1^2*la_2-la_2^2-6*ps*la_3-la_1*la_3,-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2},{-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2},{-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2}}
Idetpsi={{11*ps^2-6*ps*la_1+la_2,50*ps^3-35*ps^2*la_1+10*ps*la_2-la_3,-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,-1960*ps^3*la_3,5*ps^3-6*ps^2*la_1+6*ps*la_1^2-5*ps*la_2-la_1*la_2+la_3,4*ps^4-6*ps^3*la_1+7*ps^2*la_1^2-5*ps^2*la_2-3*ps*la_1*la_2+la_2^2+3*ps*la_3-la_1*la_3,50*ps^4-35*ps^3*la_1+35*ps^2*la_1^2-25*ps^2*la_2-10*ps*la_1*la_2+9*ps*la_3+la_1*la_3,100*ps^5-70*ps^4*la_1+55*ps^3*la_1^2-35*ps^3*la_2-30*ps^2*la_1*la_2+10*ps*la_2^2+28*ps^2*la_3-7*ps*la_1*la_3-la_2*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-450*ps^5*la_1+675*ps^4*la_1^2-505*ps^4*la_2-255*ps^3*la_1*la_2+85*ps^2*la_2^2+225*ps^3*la_3-40*ps^2*la_1*la_3-15*ps*la_2*la_3,735*ps^4*la_2-735*ps^3*la_1*la_2+560*ps^3*la_3+175*ps^2*la_1*la_3,-6*ps*la_1^3+12*ps*la_1*la_2+la_1^2*la_2-la_2^2-6*ps*la_3-la_1*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2},{50*ps^3-35*ps^2*la_1+10*ps*la_2-la_3,-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,-1960*ps^3*la_3,4*ps^4-6*ps^3*la_1+7*ps^2*la_1^2-5*ps^2*la_2-3*ps*la_1*la_2+la_2^2+3*ps*la_3-la_1*la_3,50*ps^4-35*ps^3*la_1+35*ps^2*la_1^2-25*ps^2*la_2-10*ps*la_1*la_2+9*ps*la_3+la_1*la_3,100*ps^5-70*ps^4*la_1+55*ps^3*la_1^2-35*ps^3*la_2-30*ps^2*la_1*la_2+10*ps*la_2^2+28*ps^2*la_3-7*ps*la_1*la_3-la_2*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-450*ps^5*la_1+675*ps^4*la_1^2-505*ps^4*la_2-255*ps^3*la_1*la_2+85*ps^2*la_2^2+225*ps^3*la_3-40*ps^2*la_1*la_3-15*ps*la_2*la_3,735*ps^4*la_2-735*ps^3*la_1*la_2+560*ps^3*la_3+175*ps^2*la_1*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-35*ps^2*la_1^3+70*ps^2*la_1*la_2+10*ps*la_1^2*la_2-10*ps*la_2^2-35*ps^2*la_3-10*ps*la_1*la_3-la_1^2*la_3+la_2*la_3,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2},{-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,-1960*ps^3*la_3,4*ps^4-6*ps^3*la_1+7*ps^2*la_1^2-5*ps^2*la_2-3*ps*la_1*la_2+la_2^2+3*ps*la_3-la_1*la_3,100*ps^5-70*ps^4*la_1+55*ps^3*la_1^2-35*ps^3*la_2-30*ps^2*la_1*la_2+10*ps*la_2^2+28*ps^2*la_3-7*ps*la_1*la_3-la_2*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-450*ps^5*la_1+675*ps^4*la_1^2-505*ps^4*la_2-255*ps^3*la_1*la_2+85*ps^2*la_2^2+225*ps^3*la_3-40*ps^2*la_1*la_3-15*ps*la_2*la_3,735*ps^4*la_2-735*ps^3*la_1*la_2+560*ps^3*la_3+175*ps^2*la_1*la_3,-7*ps^2*la_1^3+14*ps^2*la_1*la_2+3*ps*la_1^2*la_2-3*ps*la_2^2-la_1*la_2^2-7*ps^2*la_3-3*ps*la_1*la_3+la_1^2*la_3+la_2*la_3,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-55*ps^3*la_1^3+110*ps^3*la_1*la_2+30*ps^2*la_1^2*la_2-30*ps^2*la_2^2-10*ps*la_1*la_2^2-55*ps^3*la_3-30*ps^2*la_1*la_3+7*ps*la_1^2*la_3+13*ps*la_2*la_3+la_1*la_2*la_3-la_3^2,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2},{-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,-1960*ps^3*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-450*ps^5*la_1+675*ps^4*la_1^2-505*ps^4*la_2-255*ps^3*la_1*la_2+85*ps^2*la_2^2+225*ps^3*la_3-40*ps^2*la_1*la_3-15*ps*la_2*la_3,735*ps^4*la_2-735*ps^3*la_1*la_2+560*ps^3*la_3+175*ps^2*la_1*la_3,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2}};

--DetM={-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,735*ps^3*la_2-175*ps^2*la_3,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2}
DetM={-225*ps^3*la_1+85*ps^2*la_2-15*ps*la_3,735*ps^3*la_2-175*ps^2*la_3,-225*ps^4*la_1+225*ps^3*la_1^2-140*ps^3*la_2-85*ps^2*la_1*la_2+70*ps^2*la_3+15*ps*la_1*la_3,-1960*ps^3*la_3,735*ps^4*la_2-735*ps^3*la_1*la_2+560*ps^3*la_3+175*ps^2*la_1*la_3,-450*ps^5*la_1+675*ps^4*la_1^2-505*ps^4*la_2-255*ps^3*la_1*la_2+85*ps^2*la_2^2+225*ps^3*la_3-40*ps^2*la_1*la_3-15*ps*la_2*la_3,-225*ps^3*la_1^3+450*ps^3*la_1*la_2+85*ps^2*la_1^2*la_2-85*ps^2*la_2^2-225*ps^3*la_3-85*ps^2*la_1*la_3-15*ps*la_1^2*la_3+15*ps*la_2*la_3,300*ps^6-210*ps^5*la_1+85*ps^4*la_1^2-25*ps^4*la_2-60*ps^3*la_1*la_2+25*ps^2*la_2^2+54*ps^3*la_3-14*ps^2*la_1*la_3-6*ps*la_2*la_3+la_3^2,-ps^3*la_1^3+2*ps^3*la_1*la_2+ps^2*la_1^2*la_2-ps^2*la_2^2-ps*la_1*la_2^2+la_2^3-ps^3*la_3-ps^2*la_1*la_3+ps*la_1^2*la_3+ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,la_1^4-3*la_1^2*la_2+la_2^2+2*la_1*la_3,6*ps*la_1^4-18*ps*la_1^2*la_2-la_1^3*la_2+6*ps*la_2^2+2*la_1*la_2^2+12*ps*la_1*la_3+la_1^2*la_3-2*la_2*la_3,ps*la_1^4-la_1^5-3*ps*la_1^2*la_2+4*la_1^3*la_2+ps*la_2^2-3*la_1*la_2^2+2*ps*la_1*la_3-3*la_1^2*la_3+2*la_2*la_3,35*ps^2*la_1^4-105*ps^2*la_1^2*la_2-10*ps*la_1^3*la_2+35*ps^2*la_2^2+20*ps*la_1*la_2^2+70*ps^2*la_1*la_3+10*ps*la_1^2*la_3+la_1^3*la_3-20*ps*la_2*la_3-2*la_1*la_2*la_3+la_3^2,7*ps^2*la_1^4-21*ps^2*la_1^2*la_2-3*ps*la_1^3*la_2+7*ps^2*la_2^2+6*ps*la_1*la_2^2+la_1^2*la_2^2-la_2^3+14*ps^2*la_1*la_3+3*ps*la_1^2*la_3-la_1^3*la_3-6*ps*la_2*la_3,6*ps^2*la_1^4-6*ps*la_1^5-18*ps^2*la_1^2*la_2+23*ps*la_1^3*la_2+la_1^4*la_2+6*ps^2*la_2^2-16*ps*la_1*la_2^2-3*la_1^2*la_2^2+la_2^3+12*ps^2*la_1*la_3-17*ps*la_1^2*la_3-la_1^3*la_3+10*ps*la_2*la_3+4*la_1*la_2*la_3-la_3^2,3*ps^2*la_1^4-3*ps*la_1^5+la_1^6-9*ps^2*la_1^2*la_2+12*ps*la_1^3*la_2-5*la_1^4*la_2+3*ps^2*la_2^2-9*ps*la_1*la_2^2+6*la_1^2*la_2^2-la_2^3+6*ps^2*la_1*la_3-9*ps*la_1^2*la_3+4*la_1^3*la_3+6*ps*la_2*la_3-6*la_1*la_2*la_3+la_3^2};

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