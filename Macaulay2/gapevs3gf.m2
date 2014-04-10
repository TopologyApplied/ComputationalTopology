--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=3;
--kk= ZZ/32749;--249999797479451
----ringmap R->S
R=QQ[a,b,c,x,Degrees=>apply(toList(1..g),i->{i})|{{1}}];
S =QQ[y,Degrees=>{1}];
ofilename="result_gapevs"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ev={{a=>6*y,b=>11*y^2,c=>6*y^3,x=>y},{a=>7*y,b=>14*y^2,c=>8*y^3,x=>y},{a=>8*y,b=>17*y^2,c=>10*y^3,x=>y},{a=>9*y,b=>23*y^2,c=>15*y^3,x=>y}}
gapSeq={{1, 2, 3}, {1, 2, 4}, {1, 2, 5}, {1, 3, 5}}

--Idety=evsfromDety[g]--
IdetPsi={{11*x^2-6*x*a+b,50*x^3-35*x^2*a+10*x*b-c,-225*x^3*a+85*x^2*b-15*x*c,735*x^3*b-175*x^2*c,-1960*x^3*c,5*x^3-6*x^2*a+6*x*a^2-5*x*b-a*b+c,4*x^4-6*x^3*a+7*x^2*a^2-5*x^2*b-3*x*a*b+b^2+3*x*c-a*c,50*x^4-35*x^3*a+35*x^2*a^2-25*x^2*b-10*x*a*b+9*x*c+a*c,100*x^5-70*x^4*a+55*x^3*a^2-35*x^3*b-30*x^2*a*b+10*x*b^2+28*x^2*c-7*x*a*c-b*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,-450*x^5*a+675*x^4*a^2-505*x^4*b-255*x^3*a*b+85*x^2*b^2+225*x^3*c-40*x^2*a*c-15*x*b*c,735*x^4*b-735*x^3*a*b+560*x^3*c+175*x^2*a*c,-6*x*a^3+12*x*a*b+a^2*b-b^2-6*x*c-a*c,-7*x^2*a^3+14*x^2*a*b+3*x*a^2*b-3*x*b^2-a*b^2-7*x^2*c-3*x*a*c+a^2*c+b*c,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,-35*x^2*a^3+70*x^2*a*b+10*x*a^2*b-10*x*b^2-35*x^2*c-10*x*a*c-a^2*c+b*c,-55*x^3*a^3+110*x^3*a*b+30*x^2*a^2*b-30*x^2*b^2-10*x*a*b^2-55*x^3*c-30*x^2*a*c+7*x*a^2*c+13*x*b*c+a*b*c-c^2,-225*x^3*a^3+450*x^3*a*b+85*x^2*a^2*b-85*x^2*b^2-225*x^3*c-85*x^2*a*c-15*x*a^2*c+15*x*b*c,a^4-3*a^2*b+b^2+2*a*c,6*x*a^4-18*x*a^2*b-a^3*b+6*x*b^2+2*a*b^2+12*x*a*c+a^2*c-2*b*c,7*x^2*a^4-21*x^2*a^2*b-3*x*a^3*b+7*x^2*b^2+6*x*a*b^2+a^2*b^2-b^3+14*x^2*a*c+3*x*a^2*c-a^3*c-6*x*b*c,35*x^2*a^4-105*x^2*a^2*b-10*x*a^3*b+35*x^2*b^2+20*x*a*b^2+70*x^2*a*c+10*x*a^2*c+a^3*c-20*x*b*c-2*a*b*c+c^2,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,6*x^2*a^4-6*x*a^5-18*x^2*a^2*b+23*x*a^3*b+a^4*b+6*x^2*b^2-16*x*a*b^2-3*a^2*b^2+b^3+12*x^2*a*c-17*x*a^2*c-a^3*c+10*x*b*c+4*a*b*c-c^2,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2},{50*x^3-35*x^2*a+10*x*b-c,-225*x^3*a+85*x^2*b-15*x*c,735*x^3*b-175*x^2*c,-1960*x^3*c,4*x^4-6*x^3*a+7*x^2*a^2-5*x^2*b-3*x*a*b+b^2+3*x*c-a*c,50*x^4-35*x^3*a+35*x^2*a^2-25*x^2*b-10*x*a*b+9*x*c+a*c,100*x^5-70*x^4*a+55*x^3*a^2-35*x^3*b-30*x^2*a*b+10*x*b^2+28*x^2*c-7*x*a*c-b*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,-450*x^5*a+675*x^4*a^2-505*x^4*b-255*x^3*a*b+85*x^2*b^2+225*x^3*c-40*x^2*a*c-15*x*b*c,735*x^4*b-735*x^3*a*b+560*x^3*c+175*x^2*a*c,-7*x^2*a^3+14*x^2*a*b+3*x*a^2*b-3*x*b^2-a*b^2-7*x^2*c-3*x*a*c+a^2*c+b*c,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,-35*x^2*a^3+70*x^2*a*b+10*x*a^2*b-10*x*b^2-35*x^2*c-10*x*a*c-a^2*c+b*c,-55*x^3*a^3+110*x^3*a*b+30*x^2*a^2*b-30*x^2*b^2-10*x*a*b^2-55*x^3*c-30*x^2*a*c+7*x*a^2*c+13*x*b*c+a*b*c-c^2,-225*x^3*a^3+450*x^3*a*b+85*x^2*a^2*b-85*x^2*b^2-225*x^3*c-85*x^2*a*c-15*x*a^2*c+15*x*b*c,a^4-3*a^2*b+b^2+2*a*c,6*x*a^4-18*x*a^2*b-a^3*b+6*x*b^2+2*a*b^2+12*x*a*c+a^2*c-2*b*c,7*x^2*a^4-21*x^2*a^2*b-3*x*a^3*b+7*x^2*b^2+6*x*a*b^2+a^2*b^2-b^3+14*x^2*a*c+3*x*a^2*c-a^3*c-6*x*b*c,35*x^2*a^4-105*x^2*a^2*b-10*x*a^3*b+35*x^2*b^2+20*x*a*b^2+70*x^2*a*c+10*x*a^2*c+a^3*c-20*x*b*c-2*a*b*c+c^2,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,6*x^2*a^4-6*x*a^5-18*x^2*a^2*b+23*x*a^3*b+a^4*b+6*x^2*b^2-16*x*a*b^2-3*a^2*b^2+b^3+12*x^2*a*c-17*x*a^2*c-a^3*c+10*x*b*c+4*a*b*c-c^2,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2},{-225*x^3*a+85*x^2*b-15*x*c,735*x^3*b-175*x^2*c,-1960*x^3*c,4*x^4-6*x^3*a+7*x^2*a^2-5*x^2*b-3*x*a*b+b^2+3*x*c-a*c,100*x^5-70*x^4*a+55*x^3*a^2-35*x^3*b-30*x^2*a*b+10*x*b^2+28*x^2*c-7*x*a*c-b*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,-450*x^5*a+675*x^4*a^2-505*x^4*b-255*x^3*a*b+85*x^2*b^2+225*x^3*c-40*x^2*a*c-15*x*b*c,735*x^4*b-735*x^3*a*b+560*x^3*c+175*x^2*a*c,-7*x^2*a^3+14*x^2*a*b+3*x*a^2*b-3*x*b^2-a*b^2-7*x^2*c-3*x*a*c+a^2*c+b*c,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,-55*x^3*a^3+110*x^3*a*b+30*x^2*a^2*b-30*x^2*b^2-10*x*a*b^2-55*x^3*c-30*x^2*a*c+7*x*a^2*c+13*x*b*c+a*b*c-c^2,-225*x^3*a^3+450*x^3*a*b+85*x^2*a^2*b-85*x^2*b^2-225*x^3*c-85*x^2*a*c-15*x*a^2*c+15*x*b*c,a^4-3*a^2*b+b^2+2*a*c,6*x*a^4-18*x*a^2*b-a^3*b+6*x*b^2+2*a*b^2+12*x*a*c+a^2*c-2*b*c,7*x^2*a^4-21*x^2*a^2*b-3*x*a^3*b+7*x^2*b^2+6*x*a*b^2+a^2*b^2-b^3+14*x^2*a*c+3*x*a^2*c-a^3*c-6*x*b*c,35*x^2*a^4-105*x^2*a^2*b-10*x*a^3*b+35*x^2*b^2+20*x*a*b^2+70*x^2*a*c+10*x*a^2*c+a^3*c-20*x*b*c-2*a*b*c+c^2,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,6*x^2*a^4-6*x*a^5-18*x^2*a^2*b+23*x*a^3*b+a^4*b+6*x^2*b^2-16*x*a*b^2-3*a^2*b^2+b^3+12*x^2*a*c-17*x*a^2*c-a^3*c+10*x*b*c+4*a*b*c-c^2,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2},{-225*x^3*a+85*x^2*b-15*x*c,735*x^3*b-175*x^2*c,-1960*x^3*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,-450*x^5*a+675*x^4*a^2-505*x^4*b-255*x^3*a*b+85*x^2*b^2+225*x^3*c-40*x^2*a*c-15*x*b*c,735*x^4*b-735*x^3*a*b+560*x^3*c+175*x^2*a*c,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,-225*x^3*a^3+450*x^3*a*b+85*x^2*a^2*b-85*x^2*b^2-225*x^3*c-85*x^2*a*c-15*x*a^2*c+15*x*b*c,a^4-3*a^2*b+b^2+2*a*c,6*x*a^4-18*x*a^2*b-a^3*b+6*x*b^2+2*a*b^2+12*x*a*c+a^2*c-2*b*c,7*x^2*a^4-21*x^2*a^2*b-3*x*a^3*b+7*x^2*b^2+6*x*a*b^2+a^2*b^2-b^3+14*x^2*a*c+3*x*a^2*c-a^3*c-6*x*b*c,35*x^2*a^4-105*x^2*a^2*b-10*x*a^3*b+35*x^2*b^2+20*x*a*b^2+70*x^2*a*c+10*x*a^2*c+a^3*c-20*x*b*c-2*a*b*c+c^2,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,6*x^2*a^4-6*x*a^5-18*x^2*a^2*b+23*x*a^3*b+a^4*b+6*x^2*b^2-16*x*a*b^2-3*a^2*b^2+b^3+12*x^2*a*c-17*x*a^2*c-a^3*c+10*x*b*c+4*a*b*c-c^2,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2}};
--Idety=value get "6g.txt";
--for i in `find . -iname "*.dat"`; do sed -e '{s/\n//g;s/\\\[y\]/x/g;s/\\\[Lambda\]\[/la_/g;s/ //g;s/\]//g;}' $i >${i%.*}.txt;done;
--for i in `find . -iname "*.txt"`; do tr -d '\n' <$i >${i%.*}.lst;done;

--DetM={-225*x^3*a+85*x^2*b-15*x*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,a^4-3*a^2*b+b^2+2*a*c,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,735*x^3*b-175*x^2*c,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2}
DetM={-225*x^3*a+85*x^2*b-15*x*c,735*x^3*b-175*x^2*c,-225*x^4*a+225*x^3*a^2-140*x^3*b-85*x^2*a*b+70*x^2*c+15*x*a*c,-1960*x^3*c,735*x^4*b-735*x^3*a*b+560*x^3*c+175*x^2*a*c,-450*x^5*a+675*x^4*a^2-505*x^4*b-255*x^3*a*b+85*x^2*b^2+225*x^3*c-40*x^2*a*c-15*x*b*c,-225*x^3*a^3+450*x^3*a*b+85*x^2*a^2*b-85*x^2*b^2-225*x^3*c-85*x^2*a*c-15*x*a^2*c+15*x*b*c,300*x^6-210*x^5*a+85*x^4*a^2-25*x^4*b-60*x^3*a*b+25*x^2*b^2+54*x^3*c-14*x^2*a*c-6*x*b*c+c^2,-x^3*a^3+2*x^3*a*b+x^2*a^2*b-x^2*b^2-x*a*b^2+b^3-x^3*c-x^2*a*c+x*a^2*c+x*b*c-2*a*b*c+c^2,a^4-3*a^2*b+b^2+2*a*c,6*x*a^4-18*x*a^2*b-a^3*b+6*x*b^2+2*a*b^2+12*x*a*c+a^2*c-2*b*c,x*a^4-a^5-3*x*a^2*b+4*a^3*b+x*b^2-3*a*b^2+2*x*a*c-3*a^2*c+2*b*c,35*x^2*a^4-105*x^2*a^2*b-10*x*a^3*b+35*x^2*b^2+20*x*a*b^2+70*x^2*a*c+10*x*a^2*c+a^3*c-20*x*b*c-2*a*b*c+c^2,7*x^2*a^4-21*x^2*a^2*b-3*x*a^3*b+7*x^2*b^2+6*x*a*b^2+a^2*b^2-b^3+14*x^2*a*c+3*x*a^2*c-a^3*c-6*x*b*c,6*x^2*a^4-6*x*a^5-18*x^2*a^2*b+23*x*a^3*b+a^4*b+6*x^2*b^2-16*x*a*b^2-3*a^2*b^2+b^3+12*x^2*a*c-17*x*a^2*c-a^3*c+10*x*b*c+4*a*b*c-c^2,3*x^2*a^4-3*x*a^5+a^6-9*x^2*a^2*b+12*x*a^3*b-5*a^4*b+3*x^2*b^2-9*x*a*b^2+6*a^2*b^2-b^3+6*x^2*a*c-9*x*a^2*c+4*a^3*c+6*x*b*c-6*a*b*c+c^2};

--relationlambda={-c^2,b^2-2*a*c,-a^2+2*b}

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
-- ofile<<"Considering gap sequence "<<gapSeq#i<<":"<<endl;
-- ofile<<"Generators of the kernel of the evaluation map ev_H:"<<endl<<genI<<endl;
-- ofile<<"Number of Generators:"<<numgenI<<endl<<endl<<"**************"<<endl<<endl;
--)

I=trim intersect(apply(#ev-1,i->ideal(kernel map(S,R,ev#i),ideal(IdetPsi#i))))
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
ofile<<"The ideal of k*Omega_mu -- I'=Det(y_mu):"<<endl<<(toExternalString I')<<endl<<endl;
ofile<<"Is I' a subset of I: "<<IsubI'<<"!"<<endl<<endl;
ofile<<"Generators of the intersection of all ideals I:"<<endl<<(html genI)<<endl<<endl;
ofile<<"Generators of ideals I'=Det(y_mu):"<<endl<<(html genI')<<endl<<endl;
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
