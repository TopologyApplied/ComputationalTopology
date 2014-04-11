--Author: Renjun Xu
--E-mail: rxu@ucdavis.edu
--This projected is licensed under the terms of the MIT license.

i1 : needsPackage "ChainComplexExtras"

o1 = ChainComplexExtras

o1 : Package

i2 : R = ZZ/10007[t_1..t_2]/(t_1^2-t_2^2, 2*t_1*t_2)

o2 = R

o2 : QuotientRing

i3 : K = koszulComplex(ideal vars R)

      1      2      1
o3 = R  <-- R  <-- R

     0      1      2

o3 : ChainComplex

i4 : M = HH K

o4 = 0 : cokernel | t_1 t_2 |

     1 : subquotient ({1} | t_2 t_1  0   |, {1} | -t_2 |)
                      {1} | 0   -t_2 t_1 |  {1} | t_1  |

     2 : image {2} | t_2^2 |

o4 : GradedModule

i5 : s = hilbertSeries M
stdio:5:5:(3): error: no method found for applying hilbertSeries to:
     argument   :  0 : cokernel | t_1 t_2 |           . (of class GradedModule)
                                                      .
                   1 : subquotient ({1} | t_2 t_1  0  .
                                    {1} | 0   -t_2 t_1.
                   ....................................

i6 : prune M

o6 = 0 : cokernel | t_2 t_1 |

     1 : cokernel {2} | t_2 t_1 0   0   |
                  {2} | 0   0   t_2 t_1 |

     2 : cokernel {4} | t_2 t_1 |

o6 : GradedModule

i7 : M1 = HH_1 K

o7 = subquotient ({1} | t_2 t_1  0   |, {1} | -t_2 |)
                  {1} | 0   -t_2 t_1 |  {1} | t_1  |

                               2
o7 : R-module, subquotient of R

i8 : s1 = hilbertSeries M1

       2     3     4
     2T  - 4T  + 2T
o8 = ---------------
                2
         (1 - T)

o8 : Expression of class Divide

i9 :  M2 = HH_2 K

o9 = image {2} | t_2^2 |

                             1
o9 : R-module, submodule of R

i10 : s2 = hilbertSeries M2

       4     5    6
      T  - 2T  + T
o10 = -------------
                2
         (1 - T)

o10 : Expression of class Divide

i11 : M0 = HH_0 K

o11 = cokernel | t_1 t_2 |

                             1
o11 : R-module, quotient of R

i12 : s0 = hilbertSeries M0

                2
      1 - 2T + T
o12 = -----------
               2
        (1 - T)

o12 : Expression of class Divide

i13 :  R = ZZ/10007[t_1..t_2]/(t_1^2-t_2^2, t_1^2+t_2^2, 2*t_1*t_2)

o13 = R

o13 : QuotientRing

i14 : K = koszulComplex(ideal vars R)

       1      2      1
o14 = R  <-- R  <-- R

      0      1      2

o14 : ChainComplex

i15 : M0 = HH_0 K

o15 = cokernel | t_1 t_2 |

                             1
o15 : R-module, quotient of R

i16 : s0 = hilbertSeries M0

                2
      1 - 2T + T
o16 = -----------
               2
        (1 - T)

o16 : Expression of class Divide

i17 : M1 = HH_1 K

o17 = subquotient ({1} | t_2 t_1 0   0   |, {1} | -t_2 |)
                   {1} | 0   0   t_2 t_1 |  {1} | t_1  |

                                2
o17 : R-module, subquotient of R

i18 :  s1 = hilbertSeries M1

     	   2     3     4
      3T  - 6T  + 3T
o18 = ---------------
                	  2
          (1 - T)

o18 : Expression of class Divide

i19 : M2 = HH_2 K

o19 = image {2} | t_2 t_1 |

                             							 1
o19 : R-module, submodule of R

i20 : s2 = hilbertSeries M2

        	3     4     	5
      2T  - 4T  + 2T
o20 = ---------------
                 2
          (1 - T)

o20 : Expression of class Divide
