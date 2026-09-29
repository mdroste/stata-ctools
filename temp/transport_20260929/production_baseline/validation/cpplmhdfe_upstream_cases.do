* Numerical examples adapted from the upstream ppmlhdfe test suite.
* Source: https://github.com/sergiocorreia/ppmlhdfe/tree/master/test

clear
input byte(y id1 id2)
0 1 1
1 1 1
0 2 1
0 2 2
1 2 2
end
ppml_native_compare y, absorb(id1 id2) name("upstream hard1")

clear
input byte FTA int(ikt_fe jkt_fe ijk_fe)
1 1 3 1
1 2 4 1
1 1 5 2
1 2 6 2
1 3 1 3
1 4 2 3
1 3 7 4
1 4 8 4
1 5 1 5
1 6 2 5
0 5 3 6
1 6 4 6
0 5 5 7
1 6 6 7
1 5 7 8
1 6 8 8
end
ppml_native_compare FTA, absorb(ijk_fe jkt_fe ikt_fe) extra("keepsingletons") name("upstream hard2")

clear
input byte(y x1 x2 x3 x4)
0  5 11  2  2
0  5  2 11  2
0  5  2  2 11
0  0 -1 -1 -1
1  5  5  5  5
2  4  4  4  4
3  3  3  3  3
4  2  2  2  2
5  1  1  1  1
end
ppml_native_compare y x1 x2 x3 x4, name("upstream hard4")

clear
input byte(y x1 x2 x3 x4 x5)
0  0  0  0  0 1
0  0  0  0  0 2
0  5 11  2  2 0
0  5  2 11  2 0 
0  5  2  2 11 0
0  0 -1 -1 -1 0
1  5  5  5  5 0
2  4  4  4  4 0
3  3  3  3  3 0
4  2  2  2  2 0
5  1  1  1  1 0
end
ppml_native_compare y x1 x2 x3 x4 x5, name("upstream hard5")

clear
input double trade byte(fta exporter trend c)
35412.88183544995 0  8 15 1
619124.7429669127 0  8 16 1
 645947.355868591 0 35  0 1
2113023.937040806 1 35 12 1
end
ppml_native_compare trade fta, absorb(exporter##c.trend c) referenceopts("guess(variable trade)") name("upstream slopes1")

clear
input double(id1 y) byte(x z c)
1   67793.6092504498 0 -1 1
1              68774 0  0 1
2 29460066.818006903 1 -1 1
2  40370377.65078497 1 10 1
3  4215.552981376648 0  7 1
3  4288.904005765915 0 15 1
4  483882.9337809086 0  7 1
4 306164.57000368834 0  8 1
5 29557.571687340736 0  6 1
5  170794.3896008134 0 12 1
6 1963894.9483253956 0  6 1
6 3364687.2930059433 0 12 1
end
ppml_native_compare y x, absorb(id1##c.z c) name("upstream slopes2")

clear
input byte(index id1 id2 x z) double y
	1 1 1 0 -1  3369071263
	2 1 1 0  2  2535746605
	3 1 1 1 15 15241548895
	4 1 2 0 11        1040
	5 1 2 0 13         166
	6 1 2 0 15    21682500
	7 1 3 0 -1           0
	8 1 3 0  7      764800
end
* At tight tolerance, Stata rejects the upstream covariance as nonsymmetric
* or highly singular and posts zeros. Use an analytic sandwich oracle. The rows
* in id2=1 are exactly fitted (so the slope variance is zero). In id2=2 the
* fitted means form (a, a*r, a*r^2); solve the count and trend score equations.
* The normalized constant's variance is 6/5 * sum((y-mu)^2) / sum(y)^2.
scalar slopes3_total = y[4]+y[5]+y[6]
scalar slopes3_m = (y[5]+2*y[6])/slopes3_total
scalar slopes3_r = (slopes3_m-1 + sqrt((1-slopes3_m)^2 + ///
    4*(2-slopes3_m)*slopes3_m))/(2*(2-slopes3_m))
scalar slopes3_a = slopes3_total/(1+slopes3_r+slopes3_r^2)
scalar slopes3_v = (6/5)*((y[4]-slopes3_a)^2 + ///
    (y[5]-slopes3_a*slopes3_r)^2 + (y[6]-slopes3_a*slopes3_r^2)^2) / ///
    (y[1]+y[2]+y[3]+slopes3_total)^2
matrix slopes3_V = (0,0 \ 0,slopes3_v)
ppml_native_compare y x, absorb(id1 id2##c.z) name("upstream slopes3") vceref(slopes3_V)
matrix drop slopes3_V
scalar drop slopes3_total slopes3_m slopes3_r slopes3_a slopes3_v

clear
input long y byte(x z i j)
	      0 0 -1  1  1
	   5211 0 -1  1  4
	  10498 0  2  2  4
	   1310 0  2  2  5
	    189 0  3  3  2
	      0 0  3  3  6
	   1468 0  6  4  5
	     34 0  6  4  7
	  12000 0  7  5  3
	   1036 0  7  5  5
	      0 0  8  6  3
	   1247 0  8  6  7
	    460 0  9  7  1
	     70 0  9  7  7
	     28 0 14  8  2
	      4 0 14  8  3
	     17 0 14  8  6
	      0 0 15  9  1
	      8 0 15  9  6
	 203886 0  3 10  8
	3840369 1  3 10  9
	3961747 1  8 11  9
	 710771 0  8 11 10
	 810565 0 13 12  8
	 249014 0 13 12 10
	 819845 0 16 13  8
	5589636 1 16 13  9
end
* Match the reference's default stopping rule and retained rows, including
* the ill-conditioned separation direction that it treats as numerical zero.
ppml_native_compare y x, absorb(i j##c.z) name("upstream slopes4 default sample")
ppml_native_compare y x, absorb(i j##c.z) extra("relu_accelerate(1)") name("upstream slopes4 accelerated ReLU")
capture noisily cpplmhdfe y x, absorb(i j##c.z) separation(relu) relu_maxiter(1) relu_strict(1)
if _rc==9010 test_pass "native strict ReLU iteration limit"
else test_fail "native strict ReLU iteration limit" "expected r(9010), got r(`=_rc')"


clear
set obs 6
gen y=max(_n-3,0)
gen double x1=100*(_n<3)+(_n==6)
ppml_native_compare y x1, name("upstream ill_conditioned")
clear
drawnorm x1, n(1000) seed(101010) double
gen double u=rpoisson(1)
gen double y=exp(1+10*x1)*u
gen double x2=(y==0)
ppml_native_compare y x1 x2, name("upstream collinear2")
clear
drawnorm u x1 x2, n(1000) seed(101010) double
gen double y=exp(40+x1+x2+u)
ppml_native_compare y x1 x2, name("upstream large_y")
