--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Wei0.m2
needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
g=3
R=ZZ/32003[psi,lambda_1..lambda_g]
ofilename="result_Weierstrass"|g|"g0.txt"
ofile=ofilename<<""<<endl

Ri=R/ideal(9*psi^4-10*psi^3*lambda_1+10*psi^2*lambda_1^2-10*psi *lambda_1^3-9*psi^2*lambda_2+19*psi *lambda_1*lambda_2+lambda_1^2*lambda_2-lambda_2^2-9*psi *lambda_3-lambda_1*lambda_3,85*psi^4-60*psi^3*lambda_1+25*psi^2*lambda_1^2-14*psi^2*lambda_2-6*psi *lambda_1*lambda_2+lambda_2^2+6*psi *lambda_3-lambda_1*lambda_3,401*psi^4-255*psi^3*lambda_1+85*psi^2*lambda_1^2-40*psi^2*lambda_2-15*psi *lambda_1*lambda_2+12*psi *lambda_3+lambda_1*lambda_3,201*psi^6-271*psi^5*lambda_1+285*psi^4*lambda_1^2-285*psi^3*lambda_1^3-210*psi^4*lambda_2+480*psi^3*lambda_1*lambda_2+90*psi^2*lambda_1^2*lambda_2-75*psi^2*lambda_2^2-15*psi *lambda_1*lambda_2^2-200*psi^3*lambda_3-99*psi^2*lambda_1*lambda_3+9*psi *lambda_1^2*lambda_3+20*psi *lambda_2*lambda_3+lambda_1*lambda_2*lambda_3-lambda_3^2)
--cohW_{2,1,1}= 9*psi^4-10*psi^3*lambda_1+10*psi^2*lambda_1^2-10*psi *lambda_1^3-9*psi^2*lambda_2+19*psi *lambda_1*lambda_2+lambda_1^2*lambda_2-lambda_2^2-9*psi *lambda_3-lambda_1*lambda_3,
--cohW_{2,2}= 85*psi^4-60*psi^3*lambda_1+25*psi^2*lambda_1^2-14*psi^2*lambda_2-6*psi *lambda_1*lambda_2+lambda_2^2+6*psi *lambda_3-lambda_1*lambda_3,
--cohW_{3,1}= 401*psi^4-255*psi^3*lambda_1+85*psi^2*lambda_1^2-40*psi^2*lambda_2-15*psi *lambda_1*lambda_2+12*psi *lambda_3+lambda_1*lambda_3

s = reduceHilbert(hilbertSeries(Ri))
ofile<<"(g="<<g<<")_s="<<s<<endl

ofile<<close
<< get ofilename