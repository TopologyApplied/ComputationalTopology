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
R=QQ[t_1..t_Nspin];
--R=QQ[t_1..t_Nspin,Degrees=>apply(toList(1..Nspin),i->{1,0})]
ofilename="result_Koszul"|D|"t+d.txt"
ofile=ofilename<<""<<endl

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
        
			I=ideal(Ilist);
			K = toComplex koszulComplexDGA(I);
			
	reachzero={};reachzeroatdelta0={};reachzeroearly={};nonezeroatdelta0={}; --numpt0={};numptearly={};numptnot={};
for n from 0 to sub(D/2,ZZ) do	--could change for loop to parallel apply or scan on a function--
(
			M = HH_n K;
			resM=res pushForward(map(ring M,R), M);
			B=betti(resM);
			<<endl<<endl<<"***********************************"<<endl;
			<<"The resolution of the module of the "<<n<<"-th cohomology:"<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
			ofile<<endl<<endl<<"***********************************"<<endl;
			ofile<<"The resolution of the module of the "<<n<<"-th cohomology:"<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
--NO empty line--
			Theta=QQ[p_1..p_Nspin,SkewCommutative =>apply(toList(1..Nspin),i->p_i)];
			S=tensor(R,Theta);	--same as R**Theta except options, where p:=theta,default set degree p=(1,0),and degree t=(0,1)--
			use S;
	for j from 1 to length resM do (	--could change for loop to parallel apply or scan--
	--scan(1..length(resM),j->(		--use too much memory
			<<endl<<"reachzero: {n-th cohomology, j-th component of the resolution, degree of mu, {number of t and theta}}="<<endl;
			<<endl<<"d(mu_"<<j<<")=0:"<<endl;
			ofile<<endl<<"reachzero: {n-th cohomology, j-th component of the resolution, degree of mu, {number of t and theta}}="<<endl;
			ofile<<endl<<"d(mu_"<<j<<")=0:"<<endl;
			--for igen from 0 to (numgens resM_j)-1 do
			for igen in {0,(numgens resM_j)-1} do
			(
				dideltaa=sub((resM_j)^{igen},S);	--consider one generator, a_i=d^-1 delta a_i+1--
					--DEBUG--
					--<<"igen="<<igen<<", dideltaa50="<<dideltaa<<endl;
				idelta=j; while idelta>=0 and dideltaa!=0 do
				(
					deltaa=sub(resM.dd_idelta,S)*(transpose dideltaa);	--apply delta to a--
						--DEBUG--
						--<<"igen="<<igen<<", idelta="<<idelta<<", deltaa55="<<deltaa<<endl;
					lastdideltaa=dideltaa;
						----dideltaa=d^-1 not change the total degree of t and theta----
					dideltaa=sum(1..Nspin,i'->p_i'*diff(t_i',matrix{(flatten entries transpose deltaa)/(k'->(if sum(degree(k'+1)) ==0 then 0 else k'/sum(degree(k'+1)) ))}));
						--DEBUG--
						--<<"igen="<<igen<<", idelta="<<idelta<<", dideltaa60="<<dideltaa<<endl;
					idelta=idelta-1;
				);
				if dideltaa==0 then
				(
					reachzero=sort unique(reachzero|{{n,j,sum ((degrees (gens resM_j))_0)_igen,select(sort unique ((flatten entries lastdideltaa)/(j'->degree(j'))),k'->#k'>=2)}});
					if idelta==-1 then 
							(
							--reachzeroatdelta0=unique(reachzeroatdelta0|{{n,j,degree(sum flatten entries lastdideltaa)}}); --numpt0=unique(numpt0|{degree(lastdideltaa_(0,0))});
							--<<"reachzeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<reachzeroatdelta0<<endl;
							<<"Reach zero at delta0!"<<endl<<reachzero<<endl;
							ofile<<"Reach zero at delta0!"<<endl<<reachzero<<endl;
							)
					else 
							(
							--reachzeroearly=reachzeroearly|{{n,j,idelta,igen,degree(sum flatten entries lastdideltaa)}}; --numptearly=unique(numptearly|{degree(lastdideltaa_(0,0))});
							--<<"reachzeroearly: {n-th cohomology, j-th component of the resolution, reach zero at delta_i, act on k-th generator, {number of t and theta}}="<<endl<<reachzeroearly<<endl;
							<<"Reach zero at early at delta_"<<idelta<<endl<<reachzero<<endl;
							ofile<<"Reach zero at delta0!"<<endl<<reachzero<<endl;
							)
				)
				else
				(
						--nonezeroatdelta0=unique(nonezeroatdelta0|{{n,j,degree(sum flatten entries lastdideltaa)}}); --numptnot=unique(numptnot|{degree(lastdideltaa_(0,0))});
						--<<"nonezeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<nonezeroatdelta0<<endl;
						<<"Nonzero when reach at delta0!"<<endl<<reachzero<<endl;
						ofile<<"Nonzero when reach at delta0!"<<endl<<reachzero<<endl;
				);
				--reachzero=unique(reachzero|{sub(resM.dd_0,S)*(transpose dideltaa)});
				--numpt=unique(numpt|{degree (dideltaa_(0,0))});
			);
		);
		<<endl<<"**reachzero: {n-th cohomology, j-th component of the resolution, degree of mu, {number of t and theta}}="<<endl<<reachzero<<endl;
		ofile<<endl<<"**reachzero: {n-th cohomology, j-th component of the resolution, degree of mu, {number of t and theta}}="<<endl<<reachzero<<endl;
		--<<"reachzeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<reachzeroatdelta0<<endl;
		--<<"reachzeroearly: {n-th cohomology, j-th component of the resolution, reach zero at delta_i, act on k-th generator, {number of t and theta}}="<<endl<<reachzeroearly<<endl;
		--<<"nonezeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<nonezeroatdelta0<<endl;
		--ofile<<"reachzeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<reachzeroatdelta0<<endl;
		--ofile<<"reachzeroearly: {n-th cohomology, j-th component of the resolution, reach zero at delta_i, act on k-th generator, {number of t and theta}}="<<endl<<reachzeroearly<<endl;
		--ofile<<"nonezeroatdelta0: {n-th cohomology, j-th component of the resolution, {number of t and theta}}="<<endl<<nonezeroatdelta0<<endl;
)

t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename
	
