--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
(*read data from external files*)
t1 = cpuTime()
--needsPackage "SymmetricPolynomials"
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
g=6;
--kk= ZZ/32749;--249999797479451
----ringmap R->S
R=QQ[a,b,c,d,e,f,x,Degrees=>apply(toList(1..g),i->{i})|{{1}}];
S =QQ[y,Degrees=>{1}];
ofilename="result_gapevs"|g|"g.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl

ev={{a=>21*y,b=>175*y^2,c=>735*y^3,d=>1624*y^4,e=>1764*y^5,f=>720*y^6,x=>y},{a=>22*y,b=>190*y^2,c=>820*y^3,d=>1849*y^4,e=>2038*y^5,f=>840*y^6,x=>y},{a=>23*y,b=>205*y^2,c=>905*y^3,d=>2074*y^4,e=>2312*y^5,f=>960*y^6,x=>y},{a=>24*y,b=>220*y^2,c=>990*y^3,d=>2299*y^4,e=>2586*y^5,f=>1080*y^6,x=>y},{a=>25*y,b=>235*y^2,c=>1075*y^3,d=>2524*y^4,e=>2860*y^5,f=>1200*y^6,x=>y},{a=>26*y,b=>250*y^2,c=>1160*y^3,d=>2749*y^4,e=>3134*y^5,f=>1320*y^6,x=>y},{a=>23*y,b=>207*y^2,c=>925*y^3,d=>2144*y^4,e=>2412*y^5,f=>1008*y^6,x=>y},{a=>24*y,b=>223*y^2,c=>1020*y^3,d=>2404*y^4,e=>2736*y^5,f=>1152*y^6,x=>y},{a=>25*y,b=>239*y^2,c=>1115*y^3,d=>2664*y^4,e=>3060*y^5,f=>1296*y^6,x=>y},{a=>27*y,b=>271*y^2,c=>1305*y^3,d=>3184*y^4,e=>3708*y^5,f=>1584*y^6,x=>y},{a=>25*y,b=>241*y^2,c=>1135*y^3,d=>2734*y^4,e=>3160*y^5,f=>1344*y^6,x=>y},{a=>26*y,b=>258*y^2,c=>1240*y^3,d=>3029*y^4,e=>3534*y^5,f=>1512*y^6,x=>y},{a=>27*y,b=>277*y^2,c=>1365*y^3,d=>3394*y^4,e=>4008*y^5,f=>1728*y^6,x=>y},{a=>24*y,b=>226*y^2,c=>1056*y^3,d=>2545*y^4,e=>2952*y^5,f=>1260*y^6,x=>y},{a=>26*y,b=>260*y^2,c=>1270*y^3,d=>3159*y^4,e=>3744*y^5,f=>1620*y^6,x=>y},{a=>27*y,b=>277*y^2,c=>1377*y^3,d=>3466*y^4,e=>4140*y^5,f=>1800*y^6,x=>y},{a=>27*y,b=>280*y^2,c=>1410*y^3,d=>3589*y^4,e=>4323*y^5,f=>1890*y^6,x=>y},{a=>29*y,b=>316*y^2,c=>1646*y^3,d=>4285*y^4,e=>5237*y^5,f=>2310*y^6,x=>y},{a=>30*y,b=>340*y^2,c=>1842*y^3,d=>4951*y^4,e=>6192*y^5,f=>2772*y^6,x=>y},{a=>27*y,b=>285*y^2,c=>1485*y^3,d=>3954*y^4,e=>4968*y^5,f=>2240*y^6,x=>y},{a=>29*y,b=>323*y^2,c=>1751*y^3,d=>4796*y^4,e=>6140*y^5,f=>2800*y^6,x=>y},{a=>31*y,b=>365*y^2,c=>2065*y^3,d=>5834*y^4,e=>7624*y^5,f=>3520*y^6,x=>y},{a=>36*y,b=>505*y^2,c=>3480*y^3,d=>12139*y^4,e=>19524*y^5,f=>10395*y^6,x=>y}};
gapSeq={{1,2,3,4,5,6},{1,2,3,4,5,7},{1,2,3,4,5,8},{1,2,3,4,5,9},{1,2,3,4,5,10},{1,2,3,4,5,11},{1,2,3,4,6,7},{1,2,3,4,6,8},{1,2,3,4,6,9},{1,2,3,4,6,11},{1,2,3,4,7,8},{1,2,3,4,7,9},{1,2,3,4,8,9},{1,2,3,5,6,7},{1,2,3,5,6,9},{1,2,3,5,6,10},{1,2,3,5,7,9},{1,2,3,5,7,11},{1,2,3,6,7,11},{1,2,4,5,7,8},{1,2,4,5,7,10},{1,2,4,5,8,11},{1,3,5,7,9,11}}

--IdetPsi=evsfromDety[g]--
--IdetPsi=value get "IdetPsi6g.dat";
--$ for i in `find . -iname "*.lst"`; do sed -e '{s/\\\[Psi\]/x/g;s/\\\[Lambda\]\[1\]/a/g;s/\\\[Lambda\]\[2\]/b/g;s/\\\[Lambda\]\[3\]/c/g;s/\\\[Lambda\]\[4\]/d/g;s/\\\[Lambda\]\[5\]/e/g;s/\\\[Lambda\]\[6\]/f/g; s/ //g;}' $i >${i%.*}.dat;done;
--$ for i in `find . -iname "*.dat"`; do tr -d '\n' <$i >${i%.*}.txt;done;
--$ for i in `find . -iname "*.txt"`; do sed -e '{s/\\\[Psi\]/x/g;s/\\\[Lambda\]\[1\]/a/g;s/\\\[Lambda\]\[2\]/b/g;s/\\\[Lambda\]\[3\]/c/g;s/\\\[Lambda\]\[4\]/d/g;s/\\\[Lambda\]\[5\]/e/g;s/\\\[Lambda\]\[6\]/f/g; s/ //g;}' $i >${i%.*}.dat;done;


DetM=value get "DetM6g.dat";
--relationlambda={-c^2,b^2-2*a*c,-a^2+2*b}

--compare with faber's result--
maxkappa=2;
relationsfaber=QQ[kappa_1..kappa_(maxkappa),Degrees=>apply(toList(1..maxkappa),i->{i})]/ideal({127*kappa_1^3-2304*kappa_1*kappa_2,113*kappa_1^4-36864*kappa_2^2})
sfaber=reduceHilbert(hilbertSeries relationsfaber)


 ofile<<"g="<<g<<":"<<endl;
--for i from 0 to #ev-1 do
--(
-- mapf = map(S,R,ev#i);
-- I=kernel mapf;
-- genI=gens prune I;
-- numgenI=numgens prune I;
-- ofile<<"Considering gap sequence "<<gapSeq#i<<":"<<endl;
-- ofile<<"Generators of the kernel of the evaluation map ev_H:"<<endl<<genI<<endl;
-- ofile<<"Number of Generators:"<<numgenI<<endl<<endl<<"**************"<<endl<<endl;
--)

I=trim intersect(apply(#ev-1,i->ideal(kernel map(S,R,ev#i),ideal(value get ("./6g/IdetPsi6g_"|(i+1)|".dat")))))
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
