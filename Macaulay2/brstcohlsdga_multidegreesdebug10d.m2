--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

﻿needsPackage "DGAlgebras"

R = ZZ/101[x,y]
A = freeDGAlgebra(R,{{1,1},{1,1},{1,1},{3,3}})
A.natural
setDiff(A,{x,y,z,x*T_2*T_3-y*T_1*T_3+z*T_1*T_2})
isHomogeneous(A)


R = ZZ/101[x,y,z]
A = freeDGAlgebra(R,{{1,1},{1,1},{1,1},{3,3}})
A.natural
setDiff(A,{x,y,z,x*T_2*T_3-y*T_1*T_3+z*T_1*T_2})
isHomogeneous(A)


R = ZZ/101[x,y,z]
A = freeDGAlgebra(R,{{1},{1},{1},{3}})
A.natural
setDiff(A,{x,y,z,x*T_2*T_3-y*T_1*T_3+z*T_1*T_2})
isHomogeneous(A)


A.diff
gensList=gens A.natural
diffList=apply(gensList, f->A.diff(f))
degree first gensList
#oo
homDegreeShift := {1} | (toList ((#(degree first gensList)-1):0))
degree diffList#1
degree diffList#0
degree diffList#2
degree diffList#3
apply(#diffList, i->degree gensList#i - homDegreeShift == degree diffList#i)


i16 : degree diffList#1

o16 = {0, 1}

o16 : List

i17 : degree diffList#0

o17 = {0, 1}

o17 : List

i18 : degree diffList#2

o18 = {0, 1}

o18 : List

i19 : degree diffList#3

o19 = {2, 3}

o19 : List

i20 : apply(#diffList, i->degree gensList#i - homDegreeShift == degree diffList#i)

o20 = {true, true, true, true}