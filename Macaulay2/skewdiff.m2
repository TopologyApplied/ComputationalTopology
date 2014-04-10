D=10;Nspin=16;
R=QQ[t_1..t_16,s_1..s_16,Degrees=>toList(32:{1,0})];
 --rel=Simplify[Table[Table[Subscript[t,i],{i,1,16}].gamma10d[n].Table[Subscript[t,i],{i,1,16}],{n,1,10}]/2]--
 --Flatten[Table[{Simplify[(rel[[2*i-1]]+rel[[2*i]])/2],Simplify[(rel[[2*i-1]]-rel[[2*i]])/2]},{i,1,5}]]--
rel=matrix{{
   t_3*t_9-t_1*t_11-t_7*t_13+t_5*t_15,
  -t_4*t_10+t_2*t_12+t_8*t_14-t_6*t_16,
  -t_4*t_9+t_1*t_12+t_8*t_13-t_5*t_16,
  -t_3*t_10+t_2*t_11+t_7*t_14-t_6*t_15,
  t_2*t_9-t_1*t_10-t_8*t_15+t_7*t_16,
  t_4*t_11-t_3*t_12-t_6*t_13+t_5*t_14,
  -t_2*t_13+t_1*t_14+t_4*t_15-t_3*t_16,
  -t_6*t_9+t_5*t_10+t_8*t_11-t_7*t_12,
  t_2*t_5-t_1*t_6-t_4*t_7+t_3*t_8,
  t_10*t_13-t_9*t_14-t_12*t_15+t_11*t_16,
  s_3*s_9-s_1*s_11-s_7*s_13+s_5*s_15,
  -s_4*s_10+s_2*s_12+s_8*s_14-s_6*s_16,
  -s_4*s_9+s_1*s_12+s_8*s_13-s_5*s_16,
  -s_3*s_10+s_2*s_11+s_7*s_14-s_6*s_15,
  s_2*s_9-s_1*s_10-s_8*s_15+s_7*s_16,
  s_4*s_11-s_3*s_12-s_6*s_13+s_5*s_14,
  -s_2*s_13+s_1*s_14+s_4*s_15-s_3*s_16,
  -s_6*s_9+s_5*s_10+s_8*s_11-s_7*s_12,
  s_2*s_5-s_1*s_6-s_4*s_7+s_3*s_8,
  s_10*s_13-s_9*s_14-s_12*s_15+s_11*s_16
  }};
R=R/ideal(rel);

C=QQ[e_1..e_16,SkewCommutative =>apply(toList(1..16),i->e_i)]

Q=R**C
use Q

ex1=24*(-4*(s_12*t_11+s_11*t_12)*e_1*e_2+4*(s_11*t_10+s_10*t_11)*e_1*e_4+(s_14*t_9-s_13*t_10+3*s_16*t_11+s_15*t_12-s_10*t_13+s_9*t_14+s_12*t_15+3*s_11*t_16)*e_1*e_6-4*(s_14*t_11+s_11*t_14)*e_1*e_8+(s_12*t_3+3*s_11*t_4-s_14*t_5+s_13*t_6+3*s_4*t_11+s_3*t_12+s_6*t_13-s_5*t_14)*e_1*e_10+(-3*s_11*t_2-s_10*t_3-s_15*t_6+s_14*t_7-s_3*t_10-3*s_2*t_11+s_7*t_14-s_6*t_15)*e_1*e_12+(s_10*t_5-s_9*t_6-s_12*t_7-3*s_11*t_8-s_6*t_9+s_5*t_10-3*s_8*t_11-s_7*t_12)*e_1*e_14+4*(s_11*t_6+s_6*t_11)*e_1*e_16-4*(s_12*t_9+s_9*t_12)*e_2*e_3+(s_14*t_9-s_13*t_10-s_16*t_11-3*s_15*t_12-s_10*t_13+s_9*t_14-3*s_12*t_15-s_11*t_16)*e_2*e_5+4*(s_13*t_12+s_12*t_13)*e_2*e_7+(-3*s_12*t_3-s_11*t_4-s_14*t_5+s_13*t_6-s_4*t_11-3*s_3*t_12+s_6*t_13-s_5*t_14)*e_2*e_9+(3*s_12*t_1+s_9*t_4+s_16*t_5-s_13*t_8+s_4*t_9+3*s_1*t_12-s_8*t_13+s_5*t_16)*e_2*e_11+(s_10*t_5-s_9*t_6+3*s_12*t_7+s_11*t_8-s_6*t_9+s_5*t_10+s_8*t_11+3*s_7*t_12)*e_2*e_13-4*(s_12*t_5+s_5*t_12)*e_2*e_15-4*(s_10*t_9+s_9*t_10)*e_3*e_4-4*(s_16*t_9+s_9*t_16)*e_3*e_6+(3*s_14*t_9+s_13*t_10+s_16*t_11-s_15*t_12+s_10*t_13+3*s_9*t_14-s_12*t_15+s_11*t_16)*e_3*e_8+(-s_12*t_1-3*s_9*t_4+s_16*t_5-s_13*t_8-3*s_4*t_9-s_1*t_12-s_8*t_13+s_5*t_16)*e_3*e_10+(s_10*t_1+3*s_9*t_2-s_16*t_7+s_15*t_8+3*s_2*t_9+s_1*t_10+s_8*t_15-s_7*t_16)*e_3*e_12+4*(s_9*t_8+s_8*t_9)*e_3*e_14+(-s_10*t_5-3*s_9*t_6+s_12*t_7-s_11*t_8-3*s_6*t_9-s_5*t_10-s_8*t_11+s_7*t_12)*e_3*e_16+4*(s_15*t_10+s_10*t_15)*e_4*e_5+(-s_14*t_9-3*s_13*t_10+s_16*t_11-s_15*t_12-3*s_10*t_13-s_9*t_14-s_12*t_15+s_11*t_16)*e_4*e_7+(s_11*t_2+3*s_10*t_3-s_15*t_6+s_14*t_7+3*s_3*t_10+s_2*t_11+s_7*t_14-s_6*t_15)*e_4*e_9+(-3*s_10*t_1-s_9*t_2-s_16*t_7+s_15*t_8-s_2*t_9-3*s_1*t_10+s_8*t_15-s_7*t_16)*e_4*e_11-4*(s_10*t_7+s_7*t_10)*e_4*e_13+(3*s_10*t_5+s_9*t_6+s_12*t_7-s_11*t_8+s_6*t_9+3*s_5*t_10-s_8*t_11+s_7*t_12)*e_4*e_15-4*(s_16*t_15+s_15*t_16)*e_5*e_6+4*(s_15*t_14+s_14*t_15)*e_5*e_8+(s_14*t_1-s_13*t_2-s_16*t_3-3*s_15*t_4-s_2*t_13+s_1*t_14-3*s_4*t_15-s_3*t_16)*e_5*e_10+4*(s_15*t_2+s_2*t_15)*e_5*e_12+(-s_10*t_1+s_9*t_2+s_16*t_7+3*s_15*t_8+s_2*t_9-s_1*t_10+3*s_8*t_15+s_7*t_16)*e_5*e_14+(-s_11*t_2+s_10*t_3-3*s_15*t_6-s_14*t_7+s_3*t_10-s_2*t_11-s_7*t_14-3*s_6*t_15)*e_5*e_16-4*(s_16*t_13+s_13*t_16)*e_6*e_7+(s_14*t_1-s_13*t_2+3*s_16*t_3+s_15*t_4-s_2*t_13+s_1*t_14+s_4*t_15+3*s_3*t_16)*e_6*e_9-4*(s_16*t_1+s_1*t_16)*e_6*e_11+(-s_10*t_1+s_9*t_2-3*s_16*t_7-s_15*t_8+s_2*t_9-s_1*t_10-s_8*t_15-3*s_7*t_16)*e_6*e_13+(s_12*t_1-s_9*t_4+3*s_16*t_5+s_13*t_8-s_4*t_9+s_1*t_12+s_8*t_13+3*s_5*t_16)*e_6*e_15-4*(s_14*t_13+s_13*t_14)*e_7*e_8+4*(s_13*t_4+s_4*t_13)*e_7*e_10+(-s_14*t_1-3*s_13*t_2+s_16*t_3-s_15*t_4-3*s_2*t_13-s_1*t_14-s_4*t_15+s_3*t_16)*e_7*e_12+(s_12*t_1-s_9*t_4-s_16*t_5-3*s_13*t_8-s_4*t_9+s_1*t_12-3*s_8*t_13-s_5*t_16)*e_7*e_14+(-s_12*t_3+s_11*t_4+s_14*t_5+3*s_13*t_6+s_4*t_11-s_3*t_12+3*s_6*t_13+s_5*t_14)*e_7*e_16-4*(s_14*t_3+s_3*t_14)*e_8*e_9+(3*s_14*t_1+s_13*t_2+s_16*t_3-s_15*t_4+s_2*t_13+3*s_1*t_14-s_4*t_15+s_3*t_16)*e_8*e_11+(-s_11*t_2+s_10*t_3+s_15*t_6+3*s_14*t_7+s_3*t_10-s_2*t_11+3*s_7*t_14+s_6*t_15)*e_8*e_13+(-s_12*t_3+s_11*t_4-3*s_14*t_5-s_13*t_6+s_4*t_11-s_3*t_12-s_6*t_13-3*s_5*t_14)*e_8*e_15-4*(s_4*t_3+s_3*t_4)*e_9*e_10+4*(s_3*t_2+s_2*t_3)*e_9*e_12+(s_6*t_1-s_5*t_2+3*s_8*t_3+s_7*t_4-s_2*t_5+s_1*t_6+s_4*t_7+3*s_3*t_8)*e_9*e_14-4*(s_6*t_3+s_3*t_6)*e_9*e_16-4*(s_4*t_1+s_1*t_4)*e_10*e_11+(s_6*t_1-s_5*t_2-s_8*t_3-3*s_7*t_4-s_2*t_5+s_1*t_6-3*s_4*t_7-s_3*t_8)*e_10*e_13+4*(s_5*t_4+s_4*t_5)*e_10*e_15-4*(s_2*t_1+s_1*t_2)*e_11*e_12-4*(s_8*t_1+s_1*t_8)*e_11*e_14+(3*s_6*t_1+s_5*t_2+s_8*t_3-s_7*t_4+s_2*t_5+3*s_1*t_6-s_4*t_7+s_3*t_8)*e_11*e_16+4*(s_7*t_2+s_2*t_7)*e_12*e_13+(-s_6*t_1-3*s_5*t_2+s_8*t_3-s_7*t_4-3*s_2*t_5-s_1*t_6-s_4*t_7+s_3*t_8)*e_12*e_15-4*(s_8*t_7+s_7*t_8)*e_13*e_14+4*(s_7*t_6+s_6*t_7)*e_13*e_16-4*(s_8*t_5+s_5*t_8)*e_14*e_15-4*(s_6*t_5+s_5*t_6)*e_15*e_16);

ex2=(-4*s_14*t_9+4*s_13*t_10+4*s_16*t_11-4*s_15*t_12+4*s_10*t_13-4*s_9*t_14-4*s_12*t_15+4*s_11*t_16)*e_1*e_6+(-4*s_12*t_3+4*s_11*t_4+4*s_14*t_5-4*s_13*t_6+4*s_4*t_11-4*s_3*t_12-4*s_6*t_13+4*s_5*t_14)*e_1*e_10+(-4*s_11*t_2+4*s_10*t_3+4*s_15*t_6-4*s_14*t_7+4*s_3*t_10-4*s_2*t_11-4*s_7*t_14+4*s_6*t_15)*e_1*e_12+(-4*s_10*t_5+4*s_9*t_6+4*s_12*t_7-4*s_11*t_8+4*s_6*t_9-4*s_5*t_10-4*s_8*t_11+4*s_7*t_12)*e_1*e_14+(-4*s_14*t_9+4*s_13*t_10+4*s_16*t_11-4*s_15*t_12+4*s_10*t_13-4*s_9*t_14-4*s_12*t_15+4*s_11*t_16)*e_2*e_5+(-4*s_12*t_3+4*s_11*t_4+4*s_14*t_5-4*s_13*t_6+4*s_4*t_11-4*s_3*t_12-4*s_6*t_13+4*s_5*t_14)*e_2*e_9+(4*s_12*t_1-4*s_9*t_4-4*s_16*t_5+4*s_13*t_8-4*s_4*t_9+4*s_1*t_12+4*s_8*t_13-4*s_5*t_16)*e_2*e_11+(-4*s_10*t_5+4*s_9*t_6+4*s_12*t_7-4*s_11*t_8+4*s_6*t_9-4*s_5*t_10-4*s_8*t_11+4*s_7*t_12)*e_2*e_13+(4*s_14*t_9-4*s_13*t_10-4*s_16*t_11+4*s_15*t_12-4*s_10*t_13+4*s_9*t_14+4*s_12*t_15-4*s_11*t_16)*e_3*e_8+(4*s_12*t_1-4*s_9*t_4-4*s_16*t_5+4*s_13*t_8-4*s_4*t_9+4*s_1*t_12+4*s_8*t_13-4*s_5*t_16)*e_3*e_10+(-4*s_10*t_1+4*s_9*t_2+4*s_16*t_7-4*s_15*t_8+4*s_2*t_9-4*s_1*t_10-4*s_8*t_15+4*s_7*t_16)*e_3*e_12+(4*s_10*t_5-4*s_9*t_6-4*s_12*t_7+4*s_11*t_8-4*s_6*t_9+4*s_5*t_10+4*s_8*t_11-4*s_7*t_12)*e_3*e_16+(4*s_14*t_9-4*s_13*t_10-4*s_16*t_11+4*s_15*t_12-4*s_10*t_13+4*s_9*t_14+4*s_12*t_15-4*s_11*t_16)*e_4*e_7+(-4*s_11*t_2+4*s_10*t_3+4*s_15*t_6-4*s_14*t_7+4*s_3*t_10-4*s_2*t_11-4*s_7*t_14+4*s_6*t_15)*e_4*e_9+(-4*s_10*t_1+4*s_9*t_2+4*s_16*t_7-4*s_15*t_8+4*s_2*t_9-4*s_1*t_10-4*s_8*t_15+4*s_7*t_16)*e_4*e_11+(4*s_10*t_5-4*s_9*t_6-4*s_12*t_7+4*s_11*t_8-4*s_6*t_9+4*s_5*t_10+4*s_8*t_11-4*s_7*t_12)*e_4*e_15+(-4*s_14*t_1+4*s_13*t_2+4*s_16*t_3-4*s_15*t_4+4*s_2*t_13-4*s_1*t_14-4*s_4*t_15+4*s_3*t_16)*e_5*e_10+(4*s_10*t_1-4*s_9*t_2-4*s_16*t_7+4*s_15*t_8-4*s_2*t_9+4*s_1*t_10+4*s_8*t_15-4*s_7*t_16)*e_5*e_14+(4*s_11*t_2-4*s_10*t_3-4*s_15*t_6+4*s_14*t_7-4*s_3*t_10+4*s_2*t_11+4*s_7*t_14-4*s_6*t_15)*e_5*e_16+(-4*s_14*t_1+4*s_13*t_2+4*s_16*t_3-4*s_15*t_4+4*s_2*t_13-4*s_1*t_14-4*s_4*t_15+4*s_3*t_16)*e_6*e_9+(4*s_10*t_1-4*s_9*t_2-4*s_16*t_7+4*s_15*t_8-4*s_2*t_9+4*s_1*t_10+4*s_8*t_15-4*s_7*t_16)*e_6*e_13+(-4*s_12*t_1+4*s_9*t_4+4*s_16*t_5-4*s_13*t_8+4*s_4*t_9-4*s_1*t_12-4*s_8*t_13+4*s_5*t_16)*e_6*e_15+(4*s_14*t_1-4*s_13*t_2-4*s_16*t_3+4*s_15*t_4-4*s_2*t_13+4*s_1*t_14+4*s_4*t_15-4*s_3*t_16)*e_7*e_12+(-4*s_12*t_1+4*s_9*t_4+4*s_16*t_5-4*s_13*t_8+4*s_4*t_9-4*s_1*t_12-4*s_8*t_13+4*s_5*t_16)*e_7*e_14+(4*s_12*t_3-4*s_11*t_4-4*s_14*t_5+4*s_13*t_6-4*s_4*t_11+4*s_3*t_12+4*s_6*t_13-4*s_5*t_14)*e_7*e_16+(4*s_14*t_1-4*s_13*t_2-4*s_16*t_3+4*s_15*t_4-4*s_2*t_13+4*s_1*t_14+4*s_4*t_15-4*s_3*t_16)*e_8*e_11+(4*s_11*t_2-4*s_10*t_3-4*s_15*t_6+4*s_14*t_7-4*s_3*t_10+4*s_2*t_11+4*s_7*t_14-4*s_6*t_15)*e_8*e_13+(4*s_12*t_3-4*s_11*t_4-4*s_14*t_5+4*s_13*t_6-4*s_4*t_11+4*s_3*t_12+4*s_6*t_13-4*s_5*t_14)*e_8*e_15+(-4*s_6*t_1+4*s_5*t_2+4*s_8*t_3-4*s_7*t_4+4*s_2*t_5-4*s_1*t_6-4*s_4*t_7+4*s_3*t_8)*e_9*e_14+(-4*s_6*t_1+4*s_5*t_2+4*s_8*t_3-4*s_7*t_4+4*s_2*t_5-4*s_1*t_6-4*s_4*t_7+4*s_3*t_8)*e_10*e_13+(4*s_6*t_1-4*s_5*t_2-4*s_8*t_3+4*s_7*t_4-4*s_2*t_5+4*s_1*t_6+4*s_4*t_7-4*s_3*t_8)*e_11*e_16+(4*s_6*t_1-4*s_5*t_2-4*s_8*t_3+4*s_7*t_4-4*s_2*t_5+4*s_1*t_6+4*s_4*t_7-4*s_3*t_8)*e_12*e_15;

ox1=sum(1..16,i->((t_i+s_i)*diff(e_i,ex1)))
ox2=sum(1..16,i->((t_i+s_i)*diff(e_i,ex2)))
