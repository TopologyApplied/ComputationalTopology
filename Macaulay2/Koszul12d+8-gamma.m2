--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Koszul10d.m2
t1 = cpuTime()
needsPackage "DGAlgebras"
needsPackage "BoijSoederberg"
-- needsPackage "ChainComplexExtras"
D=12;
spindim=64;
Nsusy=1;
Nspin=spindim*Nsusy;
msusy=id_(ZZ^Nsusy);
kk= ZZ/32003;
R=kk[t_1..t_Nspin];
ofilename=D|"d8.txt"
ofile=ofilename<<""<<endl
--ofile<<run "hostname"<<run "date"<<endl


---Tleft = matrix{{t_1..t_Nspin}};	Tright= matrix{{t_1}..{t_Nspin}};
--Ilist=apply(toList(0..(D-1)),i->Tleft*(gamma#i)*Tright)
---Ilist=apply(toList(0..(D-1)),i->(Tleft*(gamma#i)*Tright)_(0,0))	--to get a list, not matrix--
Ilist={-2*t_31*t_33-2*t_32*t_34-2*t_29*t_35-2*t_30*t_36+2*t_27*t_37+2*t_28*t_38+2*t_25*t_39+2*t_26*t_40-2*t_23*t_41-2*t_24*t_42-2*t_21*t_43-2*t_22*t_44+2*t_19*t_45+2*t_20*t_46+2*t_17*t_47+2*t_18*t_48+2*t_15*t_49+2*t_16*t_50+2*t_13*t_51+2*t_14*t_52-2*t_11*t_53-2*t_12*t_54-2*t_9*t_55-2*t_10*t_56+2*t_7*t_57+2*t_8*t_58+2*t_5*t_59+2*t_6*t_60-2*t_3*t_61-2*t_4*t_62-2*t_1*t_63-2*t_2*t_64,-2*t_31*t_33+2*t_32*t_34-2*t_29*t_35+2*t_30*t_36+2*t_27*t_37-2*t_28*t_38+2*t_25*t_39-2*t_26*t_40-2*t_23*t_41+2*t_24*t_42-2*t_21*t_43+2*t_22*t_44+2*t_19*t_45-2*t_20*t_46+2*t_17*t_47-2*t_18*t_48+2*t_15*t_49-2*t_16*t_50+2*t_13*t_51-2*t_14*t_52-2*t_11*t_53+2*t_12*t_54-2*t_9*t_55+2*t_10*t_56+2*t_7*t_57-2*t_8*t_58+2*t_5*t_59-2*t_6*t_60-2*t_3*t_61+2*t_4*t_62-2*t_1*t_63+2*t_2*t_64,-2*t_30*t_33-2*t_29*t_34+2*t_32*t_35+2*t_31*t_36+2*t_26*t_37+2*t_25*t_38-2*t_28*t_39-2*t_27*t_40-2*t_22*t_41-2*t_21*t_42+2*t_24*t_43+2*t_23*t_44+2*t_18*t_45+2*t_17*t_46-2*t_20*t_47-2*t_19*t_48+2*t_14*t_49+2*t_13*t_50-2*t_16*t_51-2*t_15*t_52-2*t_10*t_53-2*t_9*t_54+2*t_12*t_55+2*t_11*t_56+2*t_6*t_57+2*t_5*t_58-2*t_8*t_59-2*t_7*t_60-2*t_2*t_61-2*t_1*t_62+2*t_4*t_63+2*t_3*t_64,-2*t_30*t_33-2*t_29*t_34-2*t_32*t_35-2*t_31*t_36+2*t_26*t_37+2*t_25*t_38+2*t_28*t_39+2*t_27*t_40-2*t_22*t_41-2*t_21*t_42-2*t_24*t_43-2*t_23*t_44+2*t_18*t_45+2*t_17*t_46+2*t_20*t_47+2*t_19*t_48+2*t_14*t_49+2*t_13*t_50+2*t_16*t_51+2*t_15*t_52-2*t_10*t_53-2*t_9*t_54-2*t_12*t_55-2*t_11*t_56+2*t_6*t_57+2*t_5*t_58+2*t_8*t_59+2*t_7*t_60-2*t_2*t_61-2*t_1*t_62-2*t_4*t_63-2*t_3*t_64,-2*t_28*t_33-2*t_27*t_34+2*t_26*t_35+2*t_25*t_36-2*t_32*t_37-2*t_31*t_38+2*t_30*t_39+2*t_29*t_40-2*t_20*t_41-2*t_19*t_42+2*t_18*t_43+2*t_17*t_44-2*t_24*t_45-2*t_23*t_46+2*t_22*t_47+2*t_21*t_48+2*t_12*t_49+2*t_11*t_50-2*t_10*t_51-2*t_9*t_52+2*t_16*t_53+2*t_15*t_54-2*t_14*t_55-2*t_13*t_56+2*t_4*t_57+2*t_3*t_58-2*t_2*t_59-2*t_1*t_60+2*t_8*t_61+2*t_7*t_62-2*t_6*t_63-2*t_5*t_64,-2*t_28*t_33-2*t_27*t_34+2*t_26*t_35+2*t_25*t_36+2*t_32*t_37+2*t_31*t_38-2*t_30*t_39-2*t_29*t_40-2*t_20*t_41-2*t_19*t_42+2*t_18*t_43+2*t_17*t_44+2*t_24*t_45+2*t_23*t_46-2*t_22*t_47-2*t_21*t_48+2*t_12*t_49+2*t_11*t_50-2*t_10*t_51-2*t_9*t_52-2*t_16*t_53-2*t_15*t_54+2*t_14*t_55+2*t_13*t_56+2*t_4*t_57+2*t_3*t_58-2*t_2*t_59-2*t_1*t_60-2*t_8*t_61-2*t_7*t_62+2*t_6*t_63+2*t_5*t_64,-2*t_24*t_33-2*t_23*t_34+2*t_22*t_35+2*t_21*t_36-2*t_20*t_37-2*t_19*t_38+2*t_18*t_39+2*t_17*t_40+2*t_32*t_41+2*t_31*t_42-2*t_30*t_43-2*t_29*t_44+2*t_28*t_45+2*t_27*t_46-2*t_26*t_47-2*t_25*t_48+2*t_8*t_49+2*t_7*t_50-2*t_6*t_51-2*t_5*t_52+2*t_4*t_53+2*t_3*t_54-2*t_2*t_55-2*t_1*t_56-2*t_16*t_57-2*t_15*t_58+2*t_14*t_59+2*t_13*t_60-2*t_12*t_61-2*t_11*t_62+2*t_10*t_63+2*t_9*t_64,-2*t_24*t_33-2*t_23*t_34+2*t_22*t_35+2*t_21*t_36-2*t_20*t_37-2*t_19*t_38+2*t_18*t_39+2*t_17*t_40-2*t_32*t_41-2*t_31*t_42+2*t_30*t_43+2*t_29*t_44-2*t_28*t_45-2*t_27*t_46+2*t_26*t_47+2*t_25*t_48+2*t_8*t_49+2*t_7*t_50-2*t_6*t_51-2*t_5*t_52+2*t_4*t_53+2*t_3*t_54-2*t_2*t_55-2*t_1*t_56+2*t_16*t_57+2*t_15*t_58-2*t_14*t_59-2*t_13*t_60+2*t_12*t_61+2*t_11*t_62-2*t_10*t_63-2*t_9*t_64,-2*t_16*t_33-2*t_15*t_34+2*t_14*t_35+2*t_13*t_36-2*t_12*t_37-2*t_11*t_38+2*t_10*t_39+2*t_9*t_40+2*t_8*t_41+2*t_7*t_42-2*t_6*t_43-2*t_5*t_44+2*t_4*t_45+2*t_3*t_46-2*t_2*t_47-2*t_1*t_48-2*t_32*t_49-2*t_31*t_50+2*t_30*t_51+2*t_29*t_52-2*t_28*t_53-2*t_27*t_54+2*t_26*t_55+2*t_25*t_56+2*t_24*t_57+2*t_23*t_58-2*t_22*t_59-2*t_21*t_60+2*t_20*t_61+2*t_19*t_62-2*t_18*t_63-2*t_17*t_64,-2*t_16*t_33-2*t_15*t_34+2*t_14*t_35+2*t_13*t_36-2*t_12*t_37-2*t_11*t_38+2*t_10*t_39+2*t_9*t_40+2*t_8*t_41+2*t_7*t_42-2*t_6*t_43-2*t_5*t_44+2*t_4*t_45+2*t_3*t_46-2*t_2*t_47-2*t_1*t_48+2*t_32*t_49+2*t_31*t_50-2*t_30*t_51-2*t_29*t_52+2*t_28*t_53+2*t_27*t_54-2*t_26*t_55-2*t_25*t_56-2*t_24*t_57-2*t_23*t_58+2*t_22*t_59+2*t_21*t_60-2*t_20*t_61-2*t_19*t_62+2*t_18*t_63+2*t_17*t_64,-2*t_16*t_17-2*t_15*t_18+2*t_14*t_19+2*t_13*t_20-2*t_12*t_21-2*t_11*t_22+2*t_10*t_23+2*t_9*t_24+2*t_8*t_25+2*t_7*t_26-2*t_6*t_27-2*t_5*t_28+2*t_4*t_29+2*t_3*t_30-2*t_2*t_31-2*t_1*t_32+2*t_48*t_49+2*t_47*t_50-2*t_46*t_51-2*t_45*t_52+2*t_44*t_53+2*t_43*t_54-2*t_42*t_55-2*t_41*t_56-2*t_40*t_57-2*t_39*t_58+2*t_38*t_59+2*t_37*t_60-2*t_36*t_61-2*t_35*t_62+2*t_34*t_63+2*t_33*t_64,-2*t_16*t_17-2*t_15*t_18+2*t_14*t_19+2*t_13*t_20-2*t_12*t_21-2*t_11*t_22+2*t_10*t_23+2*t_9*t_24+2*t_8*t_25+2*t_7*t_26-2*t_6*t_27-2*t_5*t_28+2*t_4*t_29+2*t_3*t_30-2*t_2*t_31-2*t_1*t_32-2*t_48*t_49-2*t_47*t_50+2*t_46*t_51+2*t_45*t_52-2*t_44*t_53-2*t_43*t_54+2*t_42*t_55+2*t_41*t_56+2*t_40*t_57+2*t_39*t_58-2*t_38*t_59-2*t_37*t_60+2*t_36*t_61+2*t_35*t_62-2*t_34*t_63-2*t_33*t_64}
gamma=symbol gamma; spindim=symbol spindim; Nspin=symbol Nspin; msusy=symbol msusy; Tleft=symbol Tleft; Tright=symbol Tright; 
--M0list={}
M0=R;
--resolution: set isubDim=10--
isubDim=0;deltas0=0;
I1={};
while isubDim<=8-1 and deltas0==0 do --isubDim<=D-1
(
---	ofile<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<",NSusy="<<Nsusy<<":"<<endl;
---	<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl<<"D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<",NSusy="<<Nsusy<<":"<<endl;
--	Tleft = matrix{{t_1..t_Nspin}};	Tright= matrix{{t_1}..{t_Nspin}};
--I1=I1|{(map(M0,R)) (Ilist_isubDim)};
I1=I1|{(Ilist_isubDim)};
isubDim=isubDim+1;
)
(
--	I0 = ideal(Tleft*(gamma#isubDim)*Tright);
--	I0 = I0|ideal((map(M0,R)) (Ilist_isubDim));
	I0 = ideal(I1);
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
		<<"The Hilber series:(D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<")_s0="<<s0<<endl;
		<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		<<"The 0-th cohomology generators are:"<<GenM<<endl;
			<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
		ofile<<"The Hilber series:(D="<<(isubDim+1)<<"+"<<(D-isubDim-1)<<")_s0="<<s0<<endl;
		ofile<<"Number of 0-th cohomology generators="<<NumgenM<<endl;
		ofile<<"The 0-th cohomology generators are:"<<GenM<<endl;
			ofile<<"The resolution of M_"<<n<<": "<<endl<<resM<<endl;
			ofile<<"The Betti diagram of the free resolution of M_"<<n<<": "<<endl<<B<<endl;
---		isubDim=isubDim+1;
		)
)

t2 = cpuTime()
<<"The time takes in this computation process: "<<t2-t1<<endl;
ofile<<"The time takes in this computation process: "<<t2-t1<<endl;

ofile<<close
<< get ofilename
