--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

needsPackage "DGAlgebras"

R=ZZ/32003[t_1..t_32,MonomialSize=>8];
--R=ZZ/32003[t_1..t_16,s_1..s_16,Degrees=>{16:{1},16:{-1}},MonomialSize=>8];
 --rel=Simplify[Table[Table[Subscript[t,i],{i,1,16}].gamma10d[n].Table[Subscript[t,i],{i,1,16}],{n,1,10}]/2]--
 --Flatten[Table[{Simplify[(rel[[2*i-1]]+rel[[2*i]])/2],Simplify[(rel[[2*i-1]]-rel[[2*i]])/2]},{i,1,5}]]--
rel=matrix{{
		t_15*t_17-t_13*t_19+t_11*t_21-t_9*t_23-t_7*t_25+t_5*t_27-t_3*t_29+t_1*t_31,
		-t_16*t_18+t_14*t_20-t_12*t_22+t_10*t_24+t_8*t_26-t_6*t_28+t_4*t_30-t_2*t_32,
		t_14*t_17-t_13*t_18+t_10*t_21-t_9*t_22-t_6*t_25+t_5*t_26-t_2*t_29+t_1*t_30,
		t_16*t_19-t_15*t_20+t_12*t_23-t_11*t_24-t_8*t_27+t_7*t_28-t_4*t_31+t_3*t_32,
		t_12*t_17-t_11*t_18+t_10*t_19-t_9*t_20-t_4*t_25+t_3*t_26-t_2*t_27+t_1*t_28,
		-t_16*t_21+t_15*t_22-t_14*t_23+t_13*t_24+t_8*t_29-t_7*t_30+t_6*t_31-t_5*t_32,
		t_8*t_17-t_7*t_18+t_6*t_19-t_5*t_20-t_4*t_21+t_3*t_22-t_2*t_23+t_1*t_24,
		t_16*t_25-t_15*t_26+t_14*t_27-t_13*t_28-t_12*t_29+t_11*t_30-t_10*t_31+t_9*t_32,
		t_8*t_9-t_7*t_10+t_6*t_11-t_5*t_12-t_4*t_13+t_3*t_14-t_2*t_15+t_1*t_16,
		-t_24*t_25+t_23*t_26-t_22*t_27+t_21*t_28+t_20*t_29-t_19*t_30+t_18*t_31-t_17*t_32,
		-t_16*t_17-t_15*t_18+t_14*t_19+t_13*t_20-t_12*t_21-t_11*t_22+t_10*t_23+t_9*t_24+t_8*t_25+t_7*t_26-t_6*t_27-t_5*t_28+t_4*t_29+t_3*t_30-t_2*t_31-t_1*t_32
  }};
I=ideal(rel);

n=1;
-- first differential
d1 = koszul(n,gens I);
-- second differential
d2 = koszul(n+1,gens I);
-- compute homology of a pair of maps
M = prune HH(d1,d2);
--This avoids having to generate the entire Koszul complex to simply compute the first differentia
--			s = reduceHilbert(hilbertSeries(M))
--			GenM=generators image M.cache.pruningMap
--			NumgenM=numgens M
			ofilename="result_Koszul11dres"|n|".txt";
			ofile=ofilename<<""<<endl;
--			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;

resM=res M
B=betti(resM)
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl
			ofile<<close;
apply(toList(0..9),i->tally sort degrees resM#i)


cohM = (n1,n2)->(
n:=n1;M=1;while n<=n2 and M!=0 do
	(
	t1 = cpuTime();
		-- first differential
	d1 = koszul(n,gens I);
	-- second differential
	d2 = koszul(n+1,gens I);
	-- compute homology of a pair of maps
	M = prune HH(d1,d2);
	--This avoids having to generate the entire Koszul complex to simply compute the first differentia
	s = reduceHilbert(hilbertSeries(M));
--	GenM=generators image M.cache.pruningMap;
--	NumgenM=numgens M;
			<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<s<<endl;
--			<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			ofilename="result_brstcoh10d2N"|n|".txt";
			ofile=ofilename<<""<<endl;
			ofile<<endl<<"Hilbert Series s_"<<n<<"="<<endl<<(html s)<<endl;
--			ofile<<"Number of "<<n<<"-th cohomology generators="<<NumgenM<<endl;
--			ofile<<"The "<<n<<"-th cohomology generators are:"<<GenM<<endl;
			t2 = cpuTime();
			<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<"The time takes in this computation process: "<<t2-t1<<endl;
			ofile<<close;
	n=n+1;
))

cohM(2,10)
