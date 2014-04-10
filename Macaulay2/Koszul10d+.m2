--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
--crontab -e */5 * * * * renice -0 -p $(pidof M2)

t1 = cpuTime()
needsPackage "DGAlgebras"
needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
D=10;
spindim=16;
Nsusy=1;
Nspin=spindim*Nsusy;
msusy=id_(ZZ^Nsusy);
kk= ZZ/32003;
--R=kk[t_1..t_Nspin];
R=QQ[t_1..t_Nspin];
ofilename="result_cohgen"|D|"dN"|Nsusy|"dKoszul.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl


-- gamma=WeylC1[10].WeylGamma[10, i];gamma[[1 ;; 16, 1 ;; 16]]--
----NEED to check all relations: Tleft*gamma*Tright!=0, and #(trim ideal(Ilist))=D----  
gamma0={matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 
  0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 
  0}, {0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 
  0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
  0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 
  0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 
  0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
  0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 
  0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 
  0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 
  0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, -1,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, -1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
  0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 
  0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 
  0, 0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0,
   0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, {0, 
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
  0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
  0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 
  0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 
  0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 
  0}}**msusy,
--
matrix{{0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0,
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, -1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0,
   0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 
  0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}}**msusy};

DiracC2=matrix{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0}, {0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, -1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
   0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 
  0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 
  0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0,
   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0}, {-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}


--calculate the cohomologies from n0 to nmax at reduced dimension iDim--
redux = (n0,nmax,iDim)->(
n:=n0;M=1;while n<=nmax and M!=0 do
	(
			M=symbol M;
			M = HH_n K;
			M=prune M;--also improve the ambient free module --
--		M=trim M;--faster--
			s = reduceHilbert(hilbertSeries(M));
			GenM=generators image M.cache.pruningMap;
			NumgenM=numgens M;
--!Resolution can't be derived after RingMap since we'll over QuotientRing, NOT PolynomialRing!--
--		use R;
--		M=map(R,ring M)**M;	--map M from the QuotientRing to PolynomialRing--
--		resM=resolution(map(R,ring M)**M);
			resM=res pushForward(map(ring M,R), M);
			B=betti(resM);
--		use ring M;
			<<"(D="<<(iDim+1)<<"+"<<(D-iDim-1)<<")_s"<<n<<"="<<s<<endl;
			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			ofile<<"(D="<<(iDim+1)<<"+"<<(D-iDim-1)<<")_s"<<n<<"="<<s<<endl;
			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
			ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			s=symbol s; NumgenM=symbol NumgenM; GenM=symbol GenM; resM=symbol resM; B=symbol B; 
			n=n+1;
	)
)

--gamma matrices computed based on Clifford Algebra--
D2=sub(D/2,ZZ);
G=QQ[e_1..e_D2,SkewCommutative => true]
use G;
--evenbasis=basis(0,G,Degree=>{0})|basis(2,G,Degree=>{0})|basis(4,G,Degree=>{0});
--oddbasis=basis(1,G,Degree=>{0})|basis(3,G,Degree=>{0})|basis(5,G,Degree=>{0});
evenbasis=basis(0,G)|basis(2,G)|basis(4,G);
oddbasis=basis(1,G)|basis(3,G)|basis(5,G);
gammaD={};gamma={};
--monomialeven={1,e_1*e_2,}
	for i from 1 to D2 do
	(--need to be careful about the sign below, sign(diff(e_i, evenbasis))*abs(quotientRemainder)--
--		(TxGamma, gammai)=coefficients(matrix{(flatten entries (e_i**evenbasis))/(j->(itmp=(quotientRemainder(product(gens G),j))_0;itmp*(((coefficients(itmp))_1++id_(G^1))_(0,0))))},Monomials=>evenbasis); 
--		gamma=join(gamma,{sub(gammai,R)});
--		(TxGamma, gammai)=coefficients(matrix{(flatten entries diff(e_i, evenbasis))/(j->(itmp=(quotientRemainder(product(gens G),j))_0;itmp*(((coefficients(itmp))_1++id_(G^1))_(0,0))))},Monomials=>evenbasis); 
--		gammaD=join(gammaD,{sub(gammai,R)});
		(TxGamma, gammai)=coefficients(matrix{(flatten entries (e_i**evenbasis))/(j->(itmp=(quotientRemainder(product(gens G),j))_0;itmp))},Monomials=>evenbasis); 
		gamma=join(gamma,{sub(gammai,R)});
		(TxGamma, gammai)=coefficients(matrix{(flatten entries diff(e_i,evenbasis))/(j->(itmp=(quotientRemainder(product(gens G),j))_0;itmp*(-1)^((degree (j*id_(G^1)))_0)))},Monomials=>evenbasis); 
		gamma=join(gamma,{sub(gammai,R)});
--		(TxGamma, gammai)=coefficients(matrix{(flatten entries (e_i**evenbasis))/(j->(itmp=(quotientRemainder(product(gens G),j))_0;itmp*(-1)^((degree (itmp*id_(G^1)))_0)))},Monomials=>evenbasis); 
--		gamma=join(gamma,{sub(gammai,R)*(0*id_(R^8)||id_(R^8))|(-1*id_(R^8)||0*id_(R^8))});
--		(TxGamma, gammai)=coefficients(e_i**oddbasis,Monomials=>evenbasis); 
--		gamma=join(gamma,{sub(gammai,R)});
--		(TxGamma, gammai)=coefficients(diff(e_i, oddbasis),Monomials=>evenbasis); 
--		gamma=join(gamma,{sub(gammai,R)});
	)
	gammaD={};
	for i from 1 to D2 do
	(
		(TxGamma, gammai)=coefficients(e_i**(basis G),Monomials=>(basis G)); 
		--the differential in M2: diff(e_1,e_1)=1,inconsistent with diff(e_1,e_1*e_2)=-e_2--
		(TxGamma, gammaj)=coefficients(matrix{(flatten entries diff(e_i,(basis G)))/(j->j*(-1)^((degree (j*id_(G^1)))_0))},Monomials=>(basis G));
		gammaD=join(gammaD,{sub(gammai+gammaj,R),sub(gammai-gammaj,R)});
		--gammaD=join(gammaD,{sub(gammaj,R)});
	)
	--charge conjugation matrix, B1=product(gamma#3..gamma#(D-1))--
	--verify by B*gamma*inverse(B)=-transpose(gamma)--
	--apply(0..9,i->sub(DiracB1*gammaD#i*inverse(DiracB1),ZZ)==-sub(gammaD#i,ZZ))
	DiracB1=sub(sub(product(apply(toList(0..(D-1)),i->gammaD#i)),ZZ),R)
--DiracB1=product(flatten {apply(toList(1..4),i->gammaD#(2*i)+gammaD#(2*i+1)),apply(toList(2..4),i->gammaD#(2*i)-gammaD#(2*i+1))})
--use the Dirac conjugation matrix to raise one index--
gamma={};
	  --use even basis--
gamma=apply(toList(0..(D-1)),i->submatrix(gammaD#i*DiracB1,positions(flatten (degrees basis G)_1,even),positions(flatten (degrees basis G)_1,even)));
--gamma=apply(toList(0..4),i->submatrix((gammaD#(2*i)+gammaD#(2*i+1))*DiracB1,positions(flatten (degrees basis G)_1,even),positions(flatten (degrees basis G)_1,even)));
--gamma=gamma|apply(toList(0..4),i->submatrix((gammaD#(2*i)-gammaD#(2*i+1))*DiracB1,positions(flatten (degrees basis G)_1,even),positions(flatten (degrees basis G)_1,even)));
--gamma=flatten apply(toList(0..4),i->{sub(gamma#i+gamma#(i+1),ZZ),sub(gamma#i-gamma#(i+1),ZZ)});
-----TxGamma=symbol TxGamma; gammai=symbol gammai; G=symbol G;
use R;
--gamma=gamma0;

Tright= matrix{{t_1}..{t_Nspin}};Tleft =transpose Tright;
--Tleft = matrix{{t_1..t_Nspin}};	--Wrong! Not keep homogeneity in matrix multiplication!
--Ilist=apply(toList(0..(D-1)),i->Tleft*(gamma#i)*Tright)
--Ilist=apply(toList(0..(D-1)),i->(Tleft*sub(gamma#i,ZZ)*Tright)_(0,0))  --to get a list, not matrix--
Ilist={t_16*t_1-t_15*t_4-t_8*t_12+t_9*t_11,
  t_16*t_2-t_15*t_5+t_7*t_12-t_9*t_10,
  t_16*t_3-t_15*t_6-t_7*t_11+t_8*t_10,
  t_14*t_7-t_13*t_10-t_2*t_6+t_3*t_5,
  t_14*t_8-t_13*t_11+t_1*t_6-t_3*t_4,
  t_14*t_9-t_13*t_12-t_1*t_5+t_2*t_4,
  t_10*t_1+t_11*t_2+t_12*t_3+t_14*t_15,
  t_7*t_1+t_8*t_2+t_9*t_3+t_13*t_15,
  t_7*t_4+t_8*t_5+t_9*t_6+t_13*t_16,
  t_10*t_4+t_11*t_5+t_12*t_6+t_14*t_16}
--numgens trim ideal Ilist
----gamma=symbol gamma; spindim=symbol spindim; Nspin=symbol Nspin; msusy=symbol msusy; Tleft=symbol Tleft; Tright=symbol Tright; 
--M0list={}
M0=R;
--resolution: set isubDim=10--
isubDim=0;deltas0=0;while isubDim<=D-1 and deltas0==0 do 
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<",NSusy="<<Nsusy<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<",NSusy="<<Nsusy<<":"<<endl;
--	Tleft = matrix{{t_1..t_Nspin}};	Tright= matrix{{t_1}..{t_Nspin}};
--	I0 = ideal(Tleft*(gamma#isubDim)*Tright);
	I0 = ideal((map(M0,R)) (Ilist_isubDim));
--	A = koszulComplexDGA(I0);--over the ideal--
--  K = toComplex koszulComplexDGA(I0);--over the ideal--
--
		M00=M0;
		M0=M0/I0;--New QuotientRing applied--
--	M0=prune M0; --NO prune, otherwise would change the ring--
		M0=trim M0;--faster--
--	M0list=M0list | {M0};
		s0 =reduceHilbert(hilbertSeries(M0));
		use ring (numerator s0);
		deltas0=(numerator s0)-(1+T)^(isubDim+1);	 --isubDim+1 since count from 0--
		T=symbol T;
  if deltas0==0 then
  	(
--		use M0;	--current step still need to work over M00--
--	M0^1=module(M0)--
		NumgenM=numgens M0^1;
		GenM=generators M0^1;
		resM=res (R^1 / flattenRing(M0,Result=>Ideal));
--	resM=res (R^1 / ideal (flattenRing(M0,Result=>Thing)));
		B=betti(resM);
		use M0;
		<<"(D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<")_s0="<<s0<<endl;
		<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		<<"The 0-th cohomology generators are:"<<GenM<<endl;
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
		ofile<<"(D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<")_s0="<<s0<<endl;
		ofile<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		ofile<<"The 0-th cohomology generators are:"<<GenM<<endl;
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
		isubDim=isubDim+1;
		)
)
--isubDim quit at 5, i.e. D=6+4--

use M00;
idimS1=isubDim	 --record where isubDim exit--
--Tleft = matrix{{t_1..t_Nspin}};Tright= matrix{{t_1}..{t_Nspin}};
--Mlist={}
I1={}
--I1={Tleft*(gamma#5)*Tright,Tleft*(gamma#6)*Tright,Tleft*(gamma#7)*Tright}
--!Resolution can't be easy derived with induction method since we'll over QuotientRing, NOT PolynomialRing!--
--!resolution: for i from 0!--
--for isubDim from idimS1+1 to D-1 do
for isubDim from idimS1 to D-1 do --if isubDim quit at idimS1--
(
	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<":"<<endl;
	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<",NSusy="<<Nsusy<<":"<<endl;
--	I1=I1|{Tleft*(gamma#isubDim)*Tright};
	I1=I1|{(map(M00,R)) (Ilist_isubDim)};
	I=ideal(I1);
--	rR=M0/I;
--  A = koszulComplexDGA(rR);
--	A = koszulComplexDGA(I);--over the ideal--
--	apply(maxDegree A + 1, j -> numgens prune homology(j,A));
	K = toComplex koszulComplexDGA(I);--over the ideal--
-- K.dd;
	redux(1,1,isubDim);
)
---- first differential
--d1 = koszul(1,gens I)
---- second differential
--d2 = koszul(2,gens I)
---- compute homology of a pair of maps
--H1 = prune HH(d1,d2)
----This avoids having to generate the entire Koszul complex to simply compute the first differentia


t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename
