--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
-- needsPackage "DGAlgebras"
-- needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
D=3;
spindim=2;
Nsusy=1;
Nspin=spindim*Nsusy;
msusy=id_(ZZ^Nsusy);
kk= ZZ/32003;
----c:=partial/partial C, e:=partial/partial theta, p:=partial/partial x----
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1..c_D}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>{c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10}]
--R=kk[t_1..t_Nspin,c_1..c_D,SkewCommutative =>apply(toList(1..D),i->c_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})]
R=kk[t_1..t_Nspin,c_1..c_D,e_1..e_Nspin,SkewCommutative =>apply(toList(1..D),i->c_i)|apply(toList(1..Nspin),i->e_i),Degrees=>apply(toList(1..Nspin),i->{1,0})|apply(toList(1..D),i->{-1,1})|apply(toList(1..Nspin),i->{0,1})];
ofilename="result_cohgen"|D|"dN"|Nsusy|"dt0d0.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl


-- DiracC1.DiracGamma;DiracB1.DiracGamma --
---SO(6)xSO(3) with complexification---
----NEED to check all relations: Tleft*gamma*Tright!=0----  
gamma={matrix{{1, 0}, {0, -1}}**msusy,
matrix{{1, 0}, {0, 1}}**msusy,
matrix{{0, -1}, {-1, 0}}**msusy}


Tleft = matrix{{t_1..t_Nspin}};
--Tright= matrix{{t_1}..{t_Nspin}};
--Tright=transpose Tleft;
clist= matrix{{c_1}..{c_D}};
thetalist= matrix{{e_1}..{e_Nspin}};
d1=sum(0..(D-1),i->Tleft*(gamma#i)*(transpose Tleft)*clist^{i})+Tleft*thetalist
--d1=sum(0..(D-1),i->Tleft*(gamma#i)*(transpose Tleft)*clist^{i})
--isHomogeneous d1
d1=matrix{{d1_(0,0)}}	--to be homogeneous
--isHomogeneous d1
--f = map(R^1,R^1,d1, Degree => {1,1})
--M=prune homology(f,f)
d2=map(source d1,,d1)
C=new ChainComplex; C.ring = R;
C#0 = target d1; C#1 = source d1; C#2 = source d2;
C.dd#1 = d1; C.dd#2 = d2;
--C
--isHomogeneous C
--M=(ker d1)/(image d2);
gamma=symbol gamma; spindim=symbol spindim; Nsusy=symbol Nsusy; Nspin=symbol Nspin; msusy=symbol msusy; Tleft=symbol Tleft; 
clist=symbol clist; thetalist=symbol thetalist;
--d1=symbol d1; d2=symbol d2; 
M = HH_1 C; 
M=trim M;  --faster--
GenM=generators M
M=prune M;  --also improve the ambient free module --
s = reduceHilbert(hilbertSeries(M))
NumgenM=numgens M
--resM=resolution(M)
resM=resolution(M,LengthLimit=>8)
--	relations(coimage relations(M))
B=betti(resM)
t2 = cpuTime()
<<"(D="<<D<<")_s"<<n<<"="<<s<<endl
<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl
<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl
 <<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
 <<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl
<<"The time takes in this computation process: "<<t2-t1<<endl
ofile<<"(D="<<D<<")_s"<<n<<"="<<s<<endl
ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl
ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl
 ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
 ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl

ofile<<"The time takes in this computation process: "<<t2-t1<<endl

ofile<<close
<< get ofilename
