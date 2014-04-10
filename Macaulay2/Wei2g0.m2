--ulimit -v unlimited;ulimit -d unlimited; GC_NPROCS=4 GC_MARKERS=4 M2 Wei0.m2
needsPackage "DGAlgebras"
-- needsPackage "ChainComplexExtras"
g=2
R=ZZ/32003[psi,lambda_1..lambda_g]
ofilename="result_Weierstrass"|g|"g0.txt"
ofile=ofilename<<""<<endl

Ri=R/ideal(lambda_1^2-2*lambda_2,11*psi^2-6*psi *lambda_1+lambda_2,psi^2-psi *lambda_1+lambda_1^2-lambda_2,50*psi^3-35*psi^2*lambda_1+10*psi *lambda_2,5*psi^3-6*psi^2*lambda_1+6*psi *lambda_1^2-5*psi *lambda_2-lambda_1*lambda_2,26*psi^4-35*psi^3*lambda_1+35*psi^2*lambda_1^2-25*psi^2*lambda_2-10*psi *lambda_1*lambda_2,4*psi^4-6*psi^3*lambda_1+7*psi^2*lambda_1^2-5*psi^2*lambda_2-3*psi *lambda_1*lambda_2+lambda_2^2,28*psi^5-46*psi^4*lambda_1+55*psi^3*lambda_1^2-35*psi^3*lambda_2-30*psi^2*lambda_1*lambda_2+10*psi *lambda_2^2)
--cohW_{1,1}= psi^2-psi*lambda_1+lambda_1^2-lambda_2,
--cohW_{2}= 11*psi^2-6*psi*lambda_1+lambda_2
--cohW_{2,1}=5*psi^3-6*psi^2*lambda_1+6*psi *lambda_1^2-5*psi *lambda_2-lambda_1*lambda_2
--cohW_{2}=11*psi^2-6*psi *lambda_1+lambda_2
--cohW_{3}=50*psi^3-35*psi^2*lambda_1+10*psi *lambda_2
--cohW_{1,1}=psi^2-psi *lambda_1+lambda_1^2-lambda_2
--cohW_{2,1}=5*psi^3-6*psi^2*lambda_1+6*psi *lambda_1^2-5*psi *lambda_2-lambda_1*lambda_2
--cohW_{3,1}=26*psi^4-35*psi^3*lambda_1+35*psi^2*lambda_1^2-25*psi^2*lambda_2-10*psi *lambda_1*lambda_2
--cohW_{2,2}=4*psi^4-6*psi^3*lambda_1+7*psi^2*lambda_1^2-5*psi^2*lambda_2-3*psi *lambda_1*lambda_2+lambda_2^2
--cohW_{3,2}=28*psi^5-46*psi^4*lambda_1+55*psi^3*lambda_1^2-35*psi^3*lambda_2-30*psi^2*lambda_1*lambda_2+10*psi *lambda_2^2

s = reduceHilbert(hilbertSeries(Ri))
ofile<<"(g="<<g<<")_s="<<s<<endl

ofile<<close
<< get ofilename