(* ::Package:: *)

(* ::Input:: *)
(*(*Weyl Dimension Formula for the highest weight vector[Subscript[x, 1],Subscript[x, 2],...,Subscript[x, l]]-group:D_l=so (2l)*)*)
(*(*http://arxiv.org/abs/1012.5256*)*)


(* ::Input:: *)
(*(*For the the tori like T0 and T1 all irreducible representations are one-dimensional.*)*)


(* ::Input:: *)
(*(*Weyl Dimension Formula for the highest weight vector[Subscript[x, 1],Subscript[x, 2],...,Subscript[x, l]]-group:B_l=so (2l+1)*)*)


(* ::Input:: *)
(*euler[k_]:=Sum[(-1)^n*Binomial[s-1+k-2n,s-1]*Binomial[v,n],{n,0,v}]*)


(* ::Input:: *)
(*euler10p1[k_]:=Factor[FunctionExpand[euler[k]/.s->32/.v->10]]*)


(* ::Input:: *)
(*dimB[l_]:=Product[1+Sum[Subscript[x,k],{k,i,j-1}]/(j-i),{i,1,l-1},{j,i+1,l}]*Product[1+(2*Sum[Subscript[x,k],{k,i,l-1}]+Subscript[x,l])/(2*l+1-2*i),{i,1,l}]*Product[1+(Sum[Subscript[x,k],{k,i,j-1}]+2*Sum[Subscript[x,k],{k,j,l-1}]+Subscript[x,l])/(2*l+1-i-j),{i,1,l-1},{j,i+1,l}]*)


(* ::Input:: *)
(*dimB[5]/.Subscript[x, 1]->a/.Subscript[x, 2]->b/.Subscript[x, 3]->c/.Subscript[x, 4]->d/.Subscript[x, 5]->e*)


(* ::Input:: *)
(*dim11[a_,b_,c_,d_,e_]:=dimB[5]/.Subscript[x, 1]->a/.Subscript[x, 2]->b/.Subscript[x, 3]->c/.Subscript[x, 4]->d/.Subscript[x, 5]->e*)


(* ::Input:: *)
(*dimD[l_]:=Product[1+Sum[Subscript[x,k],{k,i,j-1}]/(j-i),{i,1,l-1},{j,i+1,l}]*Product[1+(Sum[Subscript[x,k],{k,i,l-2}]+Subscript[x,l])/(l-i),{i,1,l-1}]*Product[1+(Sum[Subscript[x,k],{k,i,j-1}]+2*Sum[Subscript[x,k],{k,j,l-2}]+Subscript[x,l-1]+Subscript[x,l])/(2*l-i-j),{i,1,l-2},{j,i+1,l-1}]*)


(* ::Input:: *)
(*dim10[a_,b_,c_,d_,e_]:=dimD[5]/.Subscript[x, 1]->a/.Subscript[x, 2]->b/.Subscript[x, 3]->c/.Subscript[x, 4]->d/.Subscript[x, 5]->e*)


(* ::Input:: *)
(*(*D=11+0:Subscript[d, 0]*)*)


(* ::Input:: *)
(*Series[(1+9T+34T^2+66T^3+66T^4+34T^5+9T^6+T^7)/(1-T)^23,{T,0,16}]*)
(**)


(* ::Input:: *)
(*SeriesCoefficient[(1+9T+34T^2+66T^3+66T^4+34T^5+9T^6+T^7)/(1-T)^23,{T,0,k}]*)


(* ::Input:: *)
(*Factor[1/25545471085854720000 (504+191 k+24 k^2+k^3)^2 (100566385920000+329922471168000 k+463748770281600 k^2+377683285086720 k^3+202232210685072 k^4+76310043229184 k^5+21171752556352 k^6+4435537544000 k^7+712827303181 k^8+88504603712 k^9+8484496756 k^10+622603840 k^11+34331942 k^12+1376704 k^13+37892 k^14+640 k^15+5 k^16)]*)


(* ::Input:: *)
(*(*the dimension of conjectured cohomologies for D=11+0: H^(k,0)*)*)


(* ::Input:: *)
(*Factor[Sum[dim11[0,i,0,0,k-2i],{i,0,k/2}]]*)


(* ::Input:: *)
(*h11dim0[k_]:=((1+k) (2+k) (3+k) (4+k) (5+k) (6+k) (7+k)^2 (8+k)^2 (9+k)^2 (10+k) (11+k) (12+k) (13+k) (14+k) (15+k) (38760+13232 k+2107 k^2+160 k^3+5 k^4))/25545471085854720000*)


(* ::Input:: *)
(*(*D=11+0:Subscript[d, 1]*)*)


(* ::Input:: *)
(*Expand[(1+T)^11 (1-T)^2]*)


(* ::Input:: *)
(*Series[(11T^4+67T^5+142T^6+142T^7+67T^8+11T^9)/(1-T)^23,{T,0,24}]*)


(* ::Input:: *)
(*SeriesCoefficient[(11T^4+67T^5+142T^6+142T^7+67T^8+11T^9)/(1-T)^23,{T,0,k}]*)


(* ::Input:: *)
(*Factor[1/2554547108585472000 k (105+71 k+15 k^2+k^3)^2 (-20638126080-17398430208 k+21102900480 k^2+22011381024 k^3+1215706720 k^4-4521453376 k^5-1779634880 k^6-127707622 k^7+92502550 k^8+35440751 k^9+6592820 k^10+766597 k^11+58310 k^12+2833 k^13+80 k^14+k^15)]*)


(* ::Input:: *)
(*(*the dimension of conjectured cohomologies for D=11+0: H^(k,1)*)*)


(* ::Input:: *)
(*Factor[Sum[dim11[1,i,0,0,k-4-2i],{i,0,(k-4)/2}]]*)


(* ::Input:: *)
(*h11dim1[k_]:=1/2554547108585472000 (-3+k) (-2+k) (-1+k) k (1+k) (2+k) (3+k)^2 (4+k) (5+k)^2 (6+k) (7+k)^2 (8+k) (9+k) (10+k) (11+k) (12+k) (13+k) (58+10 k+k^2)*)


(* ::Input:: *)
(*(*D=11+0:Subscript[d, 2]*)*)


(* ::Input:: *)
(*Series[T^6(1+9T+34T^2+66T^3+66T^4+34T^5+9T^6+T^7)/(1-T)^23,{T,0,24}]*)


(* ::Input:: *)
(*SeriesCoefficient[T^6(1+9T+34T^2+66T^3+66T^4+34T^5+9T^6+T^7)/(1-T)^23,{T,0,k}]*)


(* ::Input:: *)
(*Factor[1/25545471085854720000 k (6+11 k+6 k^2+k^3)^2 (-51819264000+60137683200 k+5845772160 k^2-13940746128 k^3-1420032064 k^4+1009175152 k^5+194467280 k^6-1128899 k^7-4008112 k^8-1532444 k^9-374000 k^10-22378 k^11+8176 k^12+1892 k^13+160 k^14+5 k^15)]*)


(* ::Input:: *)
(*(**This is exactly T^4 times those of D=10+1:Subscript[d, 1] **)*)


(* ::Input:: *)
(*(*the dimension of conjectured cohomologies for D=11+0: H^(k,2)*)*)


(* ::Input:: *)
(*Factor[Sum[dim11[0,i,0,0,k-6-2i],{i,0,(k-6)/2}]]*)


(* ::Input:: *)
(*h11dim2[k_]:=((-5+k) (-4+k) (-3+k) (-2+k) (-1+k) k (1+k)^2 (2+k)^2 (3+k)^2 (4+k) (5+k) (6+k) (7+k) (8+k) (9+k) (7140+908 k+307 k^2+40 k^3+5 k^4))/25545471085854720000*)


(* ::Input:: *)
(*(*D=10+1:Subscript[d, 0]*)*)


(* ::Input:: *)
(*Series[(1+9T+35T^2+75T^3+91T^4+51T^5-8T^6-24T^7-9T^8-T^9)/(1-T)^23,{T,0,16}]*)


(* ::Input:: *)
(*Expand[(1+T)^10(1-T)]*)


(* ::Input:: *)
(*SeriesCoefficient[(1+9T+35T^2+75T^3+91T^4+51T^5-8T^6-24T^7-9T^8-T^9)/(1-T)^23,{T,0,k}]*)


(* ::Input:: *)
(*Factor[1/25545471085854720000 (25545471085854720000+103143893445992448000 k+187347710958353971200 k^2+206181520017351782400 k^3+155654661280400509824 k^4+86321378863373173632 k^5+36710405639733048720 k^6+12323920722308733440 k^7+3332512980642777904 k^8+736252237424534152 k^9+134195169088712285 k^10+20302428009907680 k^11+2556868819153194 k^12+267987827355752 k^13+23293588422475 k^14+1666871422720 k^15+97031056604 k^16+4511767512 k^17+163065315 k^18+4393760 k^19+82474 k^20+952 k^21+5 k^22)]*)


(* ::Input:: *)
(*h10p1dim0[k_]:=((1+k) (2+k) (3+k) (4+k) (5+k) (6+k) (7+k) (8+k) (9+k) (10+k) (11+k) (12+k) (13+k) (4102358400+3517874640 k+1377034344 k^2+323465858 k^3+50304233 k^4+5370505 k^5+392126 k^6+18592 k^7+497 k^8+5 k^9))/25545471085854720000*)


(* ::Input:: *)
(*(*the dimension of conjectured cohomologies for D=10+1: H^(k,0)*)*)


h10p1dim0[k_]:=(Sum[A[0,i,0,j,k-2i-j-2l],{i,0,k/2},{j,0,k-2i},{l,0,(k-2i-j)/2}]-Sum[A[0,i,0,j,k-2i-j-4],{i,0,k/2},{j,0,k-2i-4}]+Sum[A[i,j,0,l,k-2i-2j-l],{i,1,k/2},{j,0,(k-2i)/2},{l,0,k-2i-2j}]+Sum[A[0,j,0,l,k-4-2j-l],{j,0,(k-4)/2},{l,0,k-4-2j}])


h10p1dim0[5]


(* ::InheritFromParent:: *)
(*A[0,0,0,0,1]+A[0,0,0,0,3]+A[0,0,0,0,5]+A[0,0,0,1,0]+A[0,0,0,1,2]+A[0,0,0,1,4]+A[0,0,0,2,1]+A[0,0,0,2,3]+A[0,0,0,3,0]+A[0,0,0,3,2]+A[0,0,0,4,1]+A[0,0,0,5,0]+A[0,1,0,0,1]+A[0,1,0,0,3]+A[0,1,0,1,0]+A[0,1,0,1,2]+A[0,1,0,2,1]+A[0,1,0,3,0]+A[0,2,0,0,1]+A[0,2,0,1,0]+A[1,0,0,0,3]+A[1,0,0,1,2]+A[1,0,0,2,1]+A[1,0,0,3,0]+A[1,1,0,0,1]+A[1,1,0,1,0]+A[2,0,0,0,1]+A[2,0,0,1,0]/.A->dim10*)


h10p1dim0[4]


(* ::InheritFromParent:: *)
(*A[0,0,0,0,0]+A[0,0,0,0,2]+A[0,0,0,0,4]+A[0,0,0,1,1]+A[0,0,0,1,3]+A[0,0,0,2,0]+A[0,0,0,2,2]+A[0,0,0,3,1]+A[0,0,0,4,0]+A[0,1,0,0,0]+A[0,1,0,0,2]+A[0,1,0,1,1]+A[0,1,0,2,0]+A[0,2,0,0,0]+A[1,0,0,0,2]+A[1,0,0,1,1]+A[1,0,0,2,0]+A[1,1,0,0,0]+A[2,0,0,0,0]/.A->dim10*)


(* ::Input:: *)
(*Factor[Sum[dim10[0,i,0,j,k-2i-j-2l],{i,0,k/2},{j,0,k-2i},{l,0,1}]+Sum[dim10[0,i,0,j,k-2i-j-2l],{i,0,k/2},{j,0,k-2i},{l,3,(k-2i-j)/2}]+Sum[dim10[i,j,0,l,k-2i-2j-l],{i,1,k/2},{j,0,(k-2i)/2},{l,0,k-2i-2j}]]*)


(* ::Input:: *)
(*(*D=10+1:Subscript[d, 1]*)*)


(* ::Input:: *)
(*Series[(T^4+9T^5+34T^6+66T^7+66T^8+34T^9+9T^10+T^11)/(1-T)^23,{T,0,16}]*)


(* ::Input:: *)
(*(SeriesData[T, 0, {1, 32, 518, 5664, 47126, 318624, 1825341, 9121408, 40615576, 163802848, 606087053, 2078937984, 6667246210, 20134523840, 57602006955, 156915041760, 408829922805}, 0, 17, 1])-(SeriesData[T, 0, {1, 32, 517, 5632, 46618, 313248, 1782417, 8844352, 39104692, 156625216, 575710355, 1962444640, 6256996590}, 4, 17, 1])*)


(* ::Input:: *)
(*SeriesCoefficient[(T^4+9T^5+34T^6+66T^7+66T^8+34T^9+9T^10+T^11)/(1-T)^23,{T,0,k}]*)


(* ::Input:: *)
(*Factor[1/25545471085854720000 k (60+47 k+12 k^2+k^3)^2 (-42247941120-28115724672 k+46586882304 k^2+35622631056 k^3-2075583488 k^4-7388382464 k^5-2391844448 k^6-171418547 k^7+116993056 k^8+51278644 k^9+11342624 k^10+1606886 k^11+150752 k^12+9092 k^13+320 k^14+5 k^15)]*)


(* ::Input:: *)
(*h10p1dim1[k_]:=((-3+k) (-2+k) (-1+k) k (1+k) (2+k) (3+k)^2 (4+k)^2 (5+k)^2 (6+k) (7+k) (8+k) (9+k) (10+k) (11+k) (10584+2776 k+667 k^2+80 k^3+5 k^4))/25545471085854720000*)


(* ::Input:: *)
(*(*the dimension of conjectured cohomologies for D=10+1: H^(k,1)*)*)


Factor[Sum[dim10[i,j,0,l,k-4-2i-2j-l],{i,0,(k-4)/2},{j,0,(k-4)/2-i},{l,0,k-4-2i-2j}]]


(* ::Input:: *)
(*Factor[Sum[dim10[j,i-j,0,l,k-4-2i-l],{i,0,(k-4)/2},{j,0,i},{l,0,k-4-2i}]]*)


dim10[1,0,0,0,0]


(* ::Input:: *)
(*(*Euler character number*)*)


(* ::Input:: *)
(*Factor[h10p1dim0[k]-h10p1dim1[k]]*)


(* ::Input:: *)
(*euler10p1[k]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=9+2:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^9/(1-T)^23*)


(* ::Input:: *)
(*Series[(1+T)^9/(1-T)^23,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[(1+T)^9/(1-T)^23,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=8+3:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^8/(1-T)^24*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=7+4:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^7/(1-T)^25*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=6+5:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^6/(1-T)^26*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=5+6:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^5/(1-T)^27*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=4+7:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^4/(1-T)^28*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=3+8:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^3/(1-T)^29*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=2+9:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)^2/(1-T)^30*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*D=1+10:Subscript[d, 0]*)*)


(* ::Input:: *)
(*hilbertseries=(1+T)/(1-T)^31*)


(* ::Input:: *)
(*Series[hilbertseries,{T,0,16}]*)


(* ::Input:: *)
(*Factor[SeriesCoefficient[hilbertseries,{T,0,k},Assumptions-> (k>= 0)]]*)
