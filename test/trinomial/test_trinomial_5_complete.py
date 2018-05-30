import sys
from sympy.external import import_module
matchpy = import_module("matchpy")

if not matchpy:
    #bin/test will not execute any tests now
    disabled = True

if sys.version_info[:2] < (3, 6):
    disabled = True

from sympy.integrals.rubi.rubi import rubi_integrate
from sympy.functions import log, sqrt, exp, cos, sin, tan, sec, csc, cot
from sympy.functions.elementary.hyperbolic import atanh as arctanh
from sympy.functions.elementary.hyperbolic import asinh as arcsinh
from sympy.functions.elementary.hyperbolic import acosh as arccosh
from sympy.functions.elementary.trigonometric import atan as arctan
from sympy.functions.elementary.trigonometric import asin as arcsin
from sympy.functions.elementary.trigonometric import acos as arccos
from sympy.integrals.rubi.utility_function import (EllipticE, EllipticF,
    hypergeom, rubi_test, AppellF1, EllipticPi, Log, Sqrt, ArcTan, ArcTanh, ArcSin, Hypergeometric2F1)
from sympy import pi as Pi
from sympy import S, hyper, I, simplify, exp_polar, symbols


a, b, c, d, e, f, m, n, x, u , k, p= symbols('a b c d e f m n x u k p')
A, B, a, b, c, d, e, f, g, h, y, z, m, n, p, q, u, v, w, F = symbols('A B a b c d e f g h y z m n p q u v w F', real=True, imaginary=False)

def test_1():
	test = [
    
    ## ::Package:: *)

## ::Title:: *)
##Integration problems of the form x**m (a x**q+b x**n+c x**(S(2) n-q))**p*)


## ::Section::Closed:: *)
##Integrands of the form x**m (a x**S(2)+b x**S(3)+c x**S(4))**p*)


## ::Subsection::Closed:: *)
##x**m (a x**S(2)+b x**S(3)+c x**S(4))**p*)


## ::Subsubsection::Closed:: *)
##p>S(0)*)


[x**S(2)*(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(2), (a*x**S(5))/S(5) + (b*x**S(6))/S(6) + (c*x**S(7))/S(7)],
[x*(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(2), (a*x**S(4))/S(4) + (b*x**S(5))/S(5) + (c*x**S(6))/S(6)],
[a*x**S(2) + b*x**S(3) + c*x**S(4), x, S(1), (a*x**S(3))/S(3) + (b*x**S(4))/S(4) + (c*x**S(5))/S(5)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))/x, x, S(2), (a*x**S(2))/S(2) + (b*x**S(3))/S(3) + (c*x**S(4))/S(4)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(2), x, S(2), a*x + (b*x**S(2))/S(2) + (c*x**S(3))/S(3)],


[x**S(2)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(3), (a**S(2)*x**S(7))/S(7) + (a*b*x**S(8))/S(4) + ((b**S(2) + S(2)*a*c)*x**S(9))/S(9) + (b*c*x**S(10))/S(5) + (c**S(2)*x**S(11))/S(11)],
[x*(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(3), (a**S(2)*x**S(6))/S(6) + (S(2)*a*b*x**S(7))/S(7) + ((b**S(2) + S(2)*a*c)*x**S(8))/S(8) + (S(2)*b*c*x**S(9))/S(9) + (c**S(2)*x**S(10))/S(10)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(3), (a**S(2)*x**S(5))/S(5) + (a*b*x**S(6))/S(3) + ((b**S(2) + S(2)*a*c)*x**S(7))/S(7) + (b*c*x**S(8))/S(4) + (c**S(2)*x**S(9))/S(9)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2)/x, x, S(3), (a**S(2)*x**S(4))/S(4) + (S(2)*a*b*x**S(5))/S(5) + ((b**S(2) + S(2)*a*c)*x**S(6))/S(6) + (S(2)*b*c*x**S(7))/S(7) + (c**S(2)*x**S(8))/S(8)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2)/x**S(2), x, S(3), (a**S(2)*x**S(3))/S(3) + (a*b*x**S(4))/S(2) + ((b**S(2) + S(2)*a*c)*x**S(5))/S(5) + (b*c*x**S(6))/S(3) + (c**S(2)*x**S(7))/S(7)],


## ::Subsubsection::Closed:: *)
##p<S(0)*)


[x**S(5)/(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(8), -((b*x)/c**S(2)) + x**S(2)/(S(2)*c) + (b*(b**S(2) - S(3)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(c**S(3)*Sqrt(b**S(2) - S(4)*a*c)) + ((b**S(2) - a*c)*Log(a + b*x + c*x**S(2)))/(S(2)*c**S(3))],
[x**S(4)/(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(7), x/c - ((b**S(2) - S(2)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(c**S(2)*Sqrt(b**S(2) - S(4)*a*c)) - (b*Log(a + b*x + c*x**S(2)))/(S(2)*c**S(2))],
[x**S(3)/(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(6), (b*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(c*Sqrt(b**S(2) - S(4)*a*c)) + Log(a + b*x + c*x**S(2))/(S(2)*c)],
[x**S(2)/(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(3), (-S(2)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/Sqrt(b**S(2) - S(4)*a*c)],
[x/(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(8), (b*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a*Sqrt(b**S(2) - S(4)*a*c)) + Log(x)/a - Log(a + b*x + c*x**S(2))/(S(2)*a)],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(-S(1)), x, S(9), -(S(1)/(a*x)) - ((b**S(2) - S(2)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(2)*Sqrt(b**S(2) - S(4)*a*c)) - (b*Log(x))/a**S(2) + (b*Log(a + b*x + c*x**S(2)))/(S(2)*a**S(2))],
[S(1)/(x*(a*x**S(2) + b*x**S(3) + c*x**S(4))), x, S(9), -S(1)/(S(2)*a*x**S(2)) + b/(a**S(2)*x) + (b*(b**S(2) - S(3)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(3)*Sqrt(b**S(2) - S(4)*a*c)) + ((b**S(2) - a*c)*Log(x))/a**S(3) - ((b**S(2) - a*c)*Log(a + b*x + c*x**S(2)))/(S(2)*a**S(3))],
[S(1)/(x**S(2)*(a*x**S(2) + b*x**S(3) + c*x**S(4))), x, S(9), -S(1)/(S(3)*a*x**S(3)) + b/(S(2)*a**S(2)*x**S(2)) - (b**S(2) - a*c)/(a**S(3)*x) - ((b**S(4) - S(4)*a*b**S(2)*c + S(2)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(4)*Sqrt(b**S(2) - S(4)*a*c)) - (b*(b**S(2) - S(2)*a*c)*Log(x))/a**S(4) + (b*(b**S(2) - S(2)*a*c)*Log(a + b*x + c*x**S(2)))/(S(2)*a**S(4))],


[x**S(8)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), (S(2)*(b**S(2) - S(3)*a*c)*x)/(c**S(2)*(b**S(2) - S(4)*a*c)) - (b*x**S(2))/(c*(b**S(2) - S(4)*a*c)) + (x**S(3)*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2))) - (S(2)*(b**S(4) - S(6)*a*b**S(2)*c + S(6)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(c**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) - (b*Log(a + b*x + c*x**S(2)))/c**S(3)],
[x**S(7)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), -((b*x)/(c*(b**S(2) - S(4)*a*c))) + (x**S(2)*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2))) + (b*(b**S(2) - S(6)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(c**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + Log(a + b*x + c*x**S(2))/(S(2)*c**S(2))],
[x**S(6)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(4), (x*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2))) + (S(4)*a*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(5)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(4), (S(2)*a + b*x)/((b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2))) - (S(2)*b*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(4)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(4), -((b + S(2)*c*x)/((b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2)))) + (S(4)*c*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(3)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), (b**S(2) - S(2)*a*c + b*c*x)/(a*(b**S(2) - S(4)*a*c)*(a + b*x + c*x**S(2))) + (b*(b**S(2) - S(6)*a*c)*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + Log(x)/a**S(2) - Log(a + b*x + c*x**S(2))/(S(2)*a**S(2))],
[x**S(2)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), (-S(2)*(b**S(2) - S(3)*a*c))/(a**S(2)*(b**S(2) - S(4)*a*c)*x) + (b**S(2) - S(2)*a*c + b*c*x)/(a*(b**S(2) - S(4)*a*c)*x*(a + b*x + c*x**S(2))) - (S(2)*(b**S(4) - S(6)*a*b**S(2)*c + S(6)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) - (S(2)*b*Log(x))/a**S(3) + (b*Log(a + b*x + c*x**S(2)))/a**S(3)],
[x**S(1)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), -((S(3)*b**S(2) - S(8)*a*c)/(S(2)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(2))) + (b*(S(3)*b**S(2) - S(11)*a*c))/(a**S(3)*(b**S(2) - S(4)*a*c)*x) + (b**S(2) - S(2)*a*c + b*c*x)/(a*(b**S(2) - S(4)*a*c)*x**S(2)*(a + b*x + c*x**S(2))) + (b*(S(3)*b**S(4) - S(20)*a*b**S(2)*c + S(30)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(4)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + ((S(3)*b**S(2) - S(2)*a*c)*Log(x))/a**S(4) - ((S(3)*b**S(2) - S(2)*a*c)*Log(a + b*x + c*x**S(2)))/(S(2)*a**S(4))],
[x**S(0)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2), x, S(9), -((S(2)*(S(2)*b**S(2) - S(5)*a*c))/(S(3)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(3))) + (b*(S(2)*b**S(2) - S(7)*a*c))/(a**S(3)*(b**S(2) - S(4)*a*c)*x**S(2)) - (S(2)*(S(2)*b**S(4) - S(9)*a*b**S(2)*c + S(5)*a**S(2)*c**S(2)))/(a**S(4)*(b**S(2) - S(4)*a*c)*x) + (b**S(2) - S(2)*a*c + b*c*x)/(a*(b**S(2) - S(4)*a*c)*x**S(3)*(a + b*x + c*x**S(2))) - (S(2)*(S(2)*b**S(6) - S(15)*a*b**S(4)*c + S(30)*a**S(2)*b**S(2)*c**S(2) - S(10)*a**S(3)*c**S(3))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(5)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) - (S(2)*b*(S(2)*b**S(2) - S(3)*a*c)*Log(x))/a**S(5) + (b*(S(2)*b**S(2) - S(3)*a*c)*Log(a + b*x + c*x**S(2)))/a**S(5)],
[S(1)/(x*(a*x**S(2) + b*x**S(3) + c*x**S(4))**S(2)), x, S(9), -((S(5)*b**S(2) - S(12)*a*c)/(S(4)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(4))) + (b*(S(5)*b**S(2) - S(17)*a*c))/(S(3)*a**S(3)*(b**S(2) - S(4)*a*c)*x**S(3)) - (S(5)*b**S(4) - S(22)*a*b**S(2)*c + S(12)*a**S(2)*c**S(2))/(S(2)*a**S(4)*(b**S(2) - S(4)*a*c)*x**S(2)) + (b*(S(5)*b**S(4) - S(27)*a*b**S(2)*c + S(29)*a**S(2)*c**S(2)))/(a**S(5)*(b**S(2) - S(4)*a*c)*x) + (b**S(2) - S(2)*a*c + b*c*x)/(a*(b**S(2) - S(4)*a*c)*x**S(4)*(a + b*x + c*x**S(2))) + (b*(S(5)*b**S(6) - S(42)*a*b**S(4)*c + S(105)*a**S(2)*b**S(2)*c**S(2) - S(70)*a**S(3)*c**S(3))*ArcTanh((b + S(2)*c*x)/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(6)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + ((S(5)*b**S(4) - S(12)*a*b**S(2)*c + S(3)*a**S(2)*c**S(2))*Log(x))/a**S(6) - ((S(5)*b**S(4) - S(12)*a*b**S(2)*c + S(3)*a**S(2)*c**S(2))*Log(a + b*x + c*x**S(2)))/(S(2)*a**S(6))],


## ::Subsection::Closed:: *)
##x**m (a x**S(2)+b x**S(3)+c x**S(4))**(p/S(2))*)


## ::Subsubsection::Closed:: *)
##p>S(0)*)


[x**S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(8), (b*(S(35)*b**S(2) - S(116)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(960)*c**S(3)) - ((S(105)*b**S(4) - S(460)*a*b**S(2)*c + S(256)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(1920)*c**S(4)*x) - ((S(7)*b**S(2) - S(16)*a*c)*x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(240)*c**S(2)) + (x**S(2)*(b + S(8)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(40)*c) + (b*(S(7)*b**S(2) - S(12)*a*c)*(b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(256)*c**(S(9)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(7), -(((S(5)*b**S(2) - S(12)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(96)*c**S(2))) + (b*(S(15)*b**S(2) - S(52)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(192)*c**S(3)*x) + (x*(b + S(6)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(24)*c) - ((b**S(2) - S(4)*a*c)*(S(5)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(128)*c**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(5), -((b*(b + S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(8)*c**S(2)*x)) + ((a + b*x + c*x**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(3)*c*x) + (b*(b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(16)*c**(S(5)/S(2))*x*Sqrt(a + b*x + c*x**S(2)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x, x, S(4), ((b + S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*c*x) - ((b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*c**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(2), x, S(7), Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x - (Sqrt(a)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)) + (b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*Sqrt(c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(3), x, S(7), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(2)) - (b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*Sqrt(a)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) + (Sqrt(c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(4), x, S(6), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(S(2)*x**S(3))) - (b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*a*x**S(2)) + ((b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*a**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(5), x, S(7), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(S(3)*x**S(4))) - (b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(12)*a*x**S(3)) + ((S(3)*b**S(2) - S(8)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(24)*a**S(2)*x**S(2)) - (b*(b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(16)*a**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/x**S(6), x, S(8), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(S(4)*x**S(5))) - (b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(24)*a*x**S(4)) + ((S(5)*b**S(2) - S(12)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(96)*a**S(2)*x**S(3)) - (b*(S(15)*b**S(2) - S(52)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(192)*a**S(3)*x**S(2)) + ((b**S(2) - S(4)*a*c)*(S(5)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(128)*a**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],


[x*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(10), ((S(1155)*b**S(6) - S(8988)*a*b**S(4)*c + S(18896)*a**S(2)*b**S(2)*c**S(2) - S(6720)*a**S(3)*c**S(3))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(286720)*c**S(5)) - (b*(S(3465)*b**S(6) - S(30660)*a*b**S(4)*c + S(81648)*a**S(2)*b**S(2)*c**S(2) - S(58816)*a**S(3)*c**S(3))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(573440)*c**S(6)*x) - (b*(S(231)*b**S(4) - S(1560)*a*b**S(2)*c + S(2416)*a**S(2)*c**S(2))*x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(71680)*c**S(4)) + ((S(99)*b**S(4) - S(568)*a*b**S(2)*c + S(560)*a**S(2)*c**S(2))*x**S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(35840)*c**S(3)) - (x**S(3)*(b*(S(11)*b**S(2) + S(68)*a*c) + S(10)*c*(S(11)*b**S(2) - S(28)*a*c)*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4480)*c**S(2)) + (x*(S(3)*b + S(14)*c*x)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)))/(S(112)*c) + (S(3)*(b**S(2) - S(4)*a*c)**S(2)*(S(33)*b**S(4) - S(72)*a*b**S(2)*c + S(16)*a**S(2)*c**S(2))*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(32768)*c**(S(13)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(10), -((b*(S(105)*b**S(4) - S(728)*a*b**S(2)*c + S(1168)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(17920)*c**S(4))) + ((S(315)*b**S(6) - S(2520)*a*b**S(4)*c + S(5488)*a**S(2)*b**S(2)*c**S(2) - S(2048)*a**S(3)*c**S(3))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(35840)*c**S(5)*x) + ((S(7)*b**S(2) - S(32)*a*c)*(S(3)*b**S(2) - S(4)*a*c)*x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4480)*c**S(3)) - (b*(S(9)*b**S(2) - S(44)*a*c)*x**S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(2240)*c**S(2)) + (x**S(3)*(b**S(2) + S(24)*a*c + S(10)*b*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(280)*c) + (S(1)/S(7))*x*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)) - (S(3)*b*(b**S(2) - S(4)*a*c)**S(2)*(S(3)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(2048)*c**(S(11)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x, x, S(8), ((S(35)*b**S(4) - S(216)*a*b**S(2)*c + S(240)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(3840)*c**S(3)) - (b*(S(105)*b**S(4) - S(760)*a*b**S(2)*c + S(1296)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(7680)*c**S(4)*x) - (x*(b*(S(7)*b**S(2) + S(12)*a*c) + S(6)*c*(S(7)*b**S(2) - S(20)*a*c)*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(960)*c**S(2)) + ((S(3)*b + S(10)*c*x)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)))/(S(60)*c*x) + ((b**S(2) - S(4)*a*c)**S(2)*(S(7)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(1024)*c**(S(9)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(2), x, S(6), (S(3)*b*(b**S(2) - S(4)*a*c)*(b + S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(128)*c**S(3)*x) - (b*(b + S(2)*c*x)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)))/(S(16)*c**S(2)*x**S(3)) + (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(5)/S(2))/(S(5)*c*x**S(5)) - (S(3)*b*(b**S(2) - S(4)*a*c)**S(2)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(256)*c**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(3), x, S(5), -((S(3)*(b**S(2) - S(4)*a*c)*(b + S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(64)*c**S(2)*x)) + ((b + S(2)*c*x)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)))/(S(8)*c*x**S(3)) + (S(3)*(b**S(2) - S(4)*a*c)**S(2)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(128)*c**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(4), x, S(8), ((b**S(2) + S(8)*a*c + S(2)*b*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(8)*c*x) + (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/(S(3)*x**S(3)) - (a**(S(3)/S(2))*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)) - (b*(b**S(2) - S(12)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(16)*c**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(5), x, S(8), (S(3)*(S(3)*b + S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*x) - (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(4) - (S(3)*Sqrt(a)*b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) + (S(3)*(b**S(2) + S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*Sqrt(c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(6), x, S(8), -((S(3)*(b - S(2)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*x**S(2))) - (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/(S(2)*x**S(5)) - (S(3)*(b**S(2) + S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*Sqrt(a)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) + (S(3)*b*Sqrt(c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(7), x, S(9), ((b**S(2) - S(8)*a*c + S(2)*b*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(8)*a*x**S(2)) - (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/(S(3)*x**S(6)) - (b*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)))/(S(4)*a*x**S(5)) + (b*(b**S(2) - S(12)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(16)*a**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) + (c**(S(3)/S(2))*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(8), x, S(8), -(((b**S(2) - S(12)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(32)*a*x**S(3))) + (b*(S(3)*b**S(2) - S(20)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(64)*a**S(2)*x**S(2)) - ((b + S(6)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(8)*x**S(4)) - (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/(S(4)*x**S(7)) - (S(3)*(b**S(2) - S(4)*a*c)**S(2)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(128)*a**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/x**S(9), x, S(9), -(((b**S(2) - S(8)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(80)*a*x**S(4))) + (b*(S(5)*b**S(2) - S(28)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(320)*a**S(2)*x**S(3)) - ((S(15)*b**S(4) - S(100)*a*b**S(2)*c + S(128)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(640)*a**S(3)*x**S(2)) - (S(3)*(b + S(4)*c*x)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(40)*x**S(5)) - (a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))/(S(5)*x**S(8)) + (S(3)*b*(b**S(2) - S(4)*a*c)**S(2)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(256)*a**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],


## ::Subsubsection::Closed:: *)
##p<S(0)*)


[x**S(3)/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(6), -((S(3)*b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*c**S(2)*x)) + Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(S(2)*c) + ((S(3)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*c**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x**S(2)/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(4), Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(c*x) - (b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*c**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(3), (x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(Sqrt(c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[S(1)/Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)), x, S(3), -((x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(Sqrt(a)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))))],
[S(1)/(x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))), x, S(4), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(a*x**S(2))) + (b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*a**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[S(1)/(x**S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))), x, S(6), -(Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))/(S(2)*a*x**S(3))) + (S(3)*b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*a**S(2)*x**S(2)) - ((S(3)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*a**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],


[x**S(7)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(8), (S(2)*x**S(4)*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) + ((S(5)*b**S(2) - S(12)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(2)*c**S(2)*(b**S(2) - S(4)*a*c)) - (b*(S(15)*b**S(2) - S(52)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*c**S(3)*(b**S(2) - S(4)*a*c)*x) - (S(2)*b*x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(c*(b**S(2) - S(4)*a*c)) + (S(3)*(S(5)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*c**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x**S(6)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(7), (S(2)*x**S(3)*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - (S(2)*b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(c*(b**S(2) - S(4)*a*c)) + ((S(3)*b**S(2) - S(8)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(c**S(2)*(b**S(2) - S(4)*a*c)*x) - (S(3)*b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*c**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x**S(5)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(6), (S(2)*x**S(2)*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - (S(2)*b*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(c*(b**S(2) - S(4)*a*c)*x) + (x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((b + S(2)*c*x)/(S(2)*Sqrt(c)*Sqrt(a + b*x + c*x**S(2)))))/(c**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x**S(4)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(1), (S(2)*x*(S(2)*a + b*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x**S(3)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(1), -((S(2)*x*(b + S(2)*c*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))))],
[x**S(2)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(4), (S(2)*x*(b**S(2) - S(2)*a*c + b*c*x))/(a*(b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - (x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(a**(S(3)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[x/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(6), (S(2)*(b**S(2) - S(2)*a*c + b*c*x))/(a*(b**S(2) - S(4)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - ((S(3)*b**S(2) - S(8)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(a**S(2)*(b**S(2) - S(4)*a*c)*x**S(2)) + (S(3)*b*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(2)*a**(S(5)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[S(1)/(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2)), x, S(7), (S(2)*(b**S(2) - S(2)*a*c + b*c*x))/(a*(b**S(2) - S(4)*a*c)*x*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - ((S(5)*b**S(2) - S(12)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(2)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(3)) + (b*(S(15)*b**S(2) - S(52)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*a**S(3)*(b**S(2) - S(4)*a*c)*x**S(2)) - (S(3)*(S(5)*b**S(2) - S(4)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(8)*a**(S(7)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[S(1)/(x*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))), x, S(8), (S(2)*(b**S(2) - S(2)*a*c + b*c*x))/(a*(b**S(2) - S(4)*a*c)*x**S(2)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - ((S(7)*b**S(2) - S(16)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(3)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(4)) + (b*(S(35)*b**S(2) - S(116)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(12)*a**S(3)*(b**S(2) - S(4)*a*c)*x**S(3)) - ((S(105)*b**S(4) - S(460)*a*b**S(2)*c + S(256)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(24)*a**S(4)*(b**S(2) - S(4)*a*c)*x**S(2)) + (S(5)*b*(S(7)*b**S(2) - S(12)*a*c)*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(16)*a**(S(9)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],
[S(1)/(x**S(2)*(a*x**S(2) + b*x**S(3) + c*x**S(4))**(S(3)/S(2))), x, S(9), (S(2)*(b**S(2) - S(2)*a*c + b*c*x))/(a*(b**S(2) - S(4)*a*c)*x**S(3)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4))) - ((S(9)*b**S(2) - S(20)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(4)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(5)) + (b*(S(21)*b**S(2) - S(68)*a*c)*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(8)*a**S(3)*(b**S(2) - S(4)*a*c)*x**S(4)) - ((S(105)*b**S(4) - S(448)*a*b**S(2)*c + S(240)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(32)*a**S(4)*(b**S(2) - S(4)*a*c)*x**S(3)) + (b*(S(315)*b**S(4) - S(1680)*a*b**S(2)*c + S(1808)*a**S(2)*c**S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))/(S(64)*a**S(5)*(b**S(2) - S(4)*a*c)*x**S(2)) - (S(15)*(S(21)*b**S(4) - S(56)*a*b**S(2)*c + S(16)*a**S(2)*c**S(2))*x*Sqrt(a + b*x + c*x**S(2))*ArcTanh((S(2)*a + b*x)/(S(2)*Sqrt(a)*Sqrt(a + b*x + c*x**S(2)))))/(S(128)*a**(S(11)/S(2))*Sqrt(a*x**S(2) + b*x**S(3) + c*x**S(4)))],


## ::Section::Closed:: *)
##Integrands of the form x**m (a x+b x**S(3)+c x**S(5))**p*)


## ::Subsection::Closed:: *)
##x**m (a x+b x**S(3)+c x**S(5))**p*)


## ::Subsubsection::Closed:: *)
##p>S(0)*)


[x**m*(a*x + b*x**S(3) + c*x**S(5)), x, S(2), (a*x**(S(2) + m))/(S(2) + m) + (b*x**(S(4) + m))/(S(4) + m) + (c*x**(S(6) + m))/(S(6) + m)],

[x**S(2)*(a*x + b*x**S(3) + c*x**S(5)), x, S(2), (a*x**S(4))/S(4) + (b*x**S(6))/S(6) + (c*x**S(8))/S(8)],
[x**S(1)*(a*x + b*x**S(3) + c*x**S(5)), x, S(2), (a*x**S(3))/S(3) + (b*x**S(5))/S(5) + (c*x**S(7))/S(7)],
[x**S(0)*(a*x + b*x**S(3) + c*x**S(5)), x, S(1), (a*x**S(2))/S(2) + (b*x**S(4))/S(4) + (c*x**S(6))/S(6)],
[(a*x + b*x**S(3) + c*x**S(5))/x**S(1), x, S(2), a*x + (b*x**S(3))/S(3) + (c*x**S(5))/S(5)],
[(a*x + b*x**S(3) + c*x**S(5))/x**S(2), x, S(2), (b*x**S(2))/S(2) + (c*x**S(4))/S(4) + a*Log(x)],
[(a*x + b*x**S(3) + c*x**S(5))/x**S(3), x, S(2), -(a/x) + b*x + (c*x**S(3))/S(3)],


[x**m*(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(3), (a**S(2)*x**(S(3) + m))/(S(3) + m) + (S(2)*a*b*x**(S(5) + m))/(S(5) + m) + ((b**S(2) + S(2)*a*c)*x**(S(7) + m))/(S(7) + m) + (S(2)*b*c*x**(S(9) + m))/(S(9) + m) + (c**S(2)*x**(S(11) + m))/(S(11) + m)],

[x**S(2)*(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(3), (a**S(2)*x**S(5))/S(5) + (S(2)*a*b*x**S(7))/S(7) + ((b**S(2) + S(2)*a*c)*x**S(9))/S(9) + (S(2)*b*c*x**S(11))/S(11) + (c**S(2)*x**S(13))/S(13)],
[x**S(1)*(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(4), (a**S(2)*x**S(4))/S(4) + (a*b*x**S(6))/S(3) + ((b**S(2) + S(2)*a*c)*x**S(8))/S(8) + (b*c*x**S(10))/S(5) + (c**S(2)*x**S(12))/S(12)],
[x**S(0)*(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(3), (a**S(2)*x**S(3))/S(3) + (S(2)*a*b*x**S(5))/S(5) + ((b**S(2) + S(2)*a*c)*x**S(7))/S(7) + (S(2)*b*c*x**S(9))/S(9) + (c**S(2)*x**S(11))/S(11)],
[(a*x + b*x**S(3) + c*x**S(5))**S(2)/x**S(1), x, S(4), (a**S(2)*x**S(2))/S(2) + (a*b*x**S(4))/S(2) + ((b**S(2) + S(2)*a*c)*x**S(6))/S(6) + (b*c*x**S(8))/S(4) + (c**S(2)*x**S(10))/S(10)],
[(a*x + b*x**S(3) + c*x**S(5))**S(2)/x**S(2), x, S(3), a**S(2)*x + (S(2)*a*b*x**S(3))/S(3) + ((b**S(2) + S(2)*a*c)*x**S(5))/S(5) + (S(2)*b*c*x**S(7))/S(7) + (c**S(2)*x**S(9))/S(9)],


## ::Subsubsection::Closed:: *)
##p<S(0)*)


[x**S(8)/(a*x + b*x**S(3) + c*x**S(5)), x, S(9), -((b*x**S(2))/(S(2)*c**S(2))) + x**S(4)/(S(4)*c) + (b*(b**S(2) - S(3)*a*c)*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*c**S(3)*Sqrt(b**S(2) - S(4)*a*c)) + ((b**S(2) - a*c)*Log(a + b*x**S(2) + c*x**S(4)))/(S(4)*c**S(3))],
[x**S(7)/(a*x + b*x**S(3) + c*x**S(5)), x, S(6), -((b*x)/c**S(2)) + x**S(3)/(S(3)*c) + ((b**S(2) - a*c - (b*(b**S(2) - S(3)*a*c))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*c**(S(5)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) + ((b**S(2) - a*c + (b*(b**S(2) - S(3)*a*c))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*c**(S(5)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(6)/(a*x + b*x**S(3) + c*x**S(5)), x, S(8), x**S(2)/(S(2)*c) - ((b**S(2) - S(2)*a*c)*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*c**S(2)*Sqrt(b**S(2) - S(4)*a*c)) - (b*Log(a + b*x**S(2) + c*x**S(4)))/(S(4)*c**S(2))],
[x**S(5)/(a*x + b*x**S(3) + c*x**S(5)), x, S(5), x/c - ((b - (b**S(2) - S(2)*a*c)/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*c**(S(3)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - ((b + (b**S(2) - S(2)*a*c)/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*c**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(4)/(a*x + b*x**S(3) + c*x**S(5)), x, S(7), (b*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*c*Sqrt(b**S(2) - S(4)*a*c)) + Log(a + b*x**S(2) + c*x**S(4))/(S(4)*c)],
[x**S(3)/(a*x + b*x**S(3) + c*x**S(5)), x, S(4), -((Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*Sqrt(c)*Sqrt(b**S(2) - S(4)*a*c))) + (Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*Sqrt(c)*Sqrt(b**S(2) - S(4)*a*c))],
[x**S(2)/(a*x + b*x**S(3) + c*x**S(5)), x, S(4), -(ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c))/Sqrt(b**S(2) - S(4)*a*c))],
[x/(a*x + b*x**S(3) + c*x**S(5)), x, S(4), (Sqrt(S(2))*Sqrt(c)*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(b**S(2) - S(4)*a*c)*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - (Sqrt(S(2))*Sqrt(c)*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(b**S(2) - S(4)*a*c)*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[S(1)/(a*x + b*x**S(3) + c*x**S(5)), x, S(9), (b*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*a*Sqrt(b**S(2) - S(4)*a*c)) + Log(x)/a - Log(a + b*x**S(2) + c*x**S(4))/(S(4)*a)],
[S(1)/(x*(a*x + b*x**S(3) + c*x**S(5))), x, S(5), -(S(1)/(a*x)) - (Sqrt(c)*(S(1) + b/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*a*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - (Sqrt(c)*(S(1) - b/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*a*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[S(1)/(x**S(2)*(a*x + b*x**S(3) + c*x**S(5))), x, S(10), -S(1)/(S(2)*a*x**S(2)) - ((b**S(2) - S(2)*a*c)*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*a**S(2)*Sqrt(b**S(2) - S(4)*a*c)) - (b*Log(x))/a**S(2) + (b*Log(a + b*x**S(2) + c*x**S(4)))/(S(4)*a**S(2))],


[x**S(11)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(10), ((b**S(2) - S(3)*a*c)*x**S(2))/(c**S(2)*(b**S(2) - S(4)*a*c)) - (b*x**S(4))/(S(2)*c*(b**S(2) - S(4)*a*c)) + (x**S(6)*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) - ((b**S(4) - S(6)*a*b**S(2)*c + S(6)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(c**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) - (b*Log(a + b*x**S(2) + c*x**S(4)))/(S(2)*c**S(3))],
[x**S(10)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(7), ((S(3)*b**S(2) - S(10)*a*c)*x)/(S(2)*c**S(2)*(b**S(2) - S(4)*a*c)) - (b*x**S(3))/(S(2)*c*(b**S(2) - S(4)*a*c)) + (x**S(5)*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) - ((S(3)*b**S(3) - S(13)*a*b*c - (S(3)*b**S(4) - S(19)*a*b**S(2)*c + S(20)*a**S(2)*c**S(2))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*c**(S(5)/S(2))*(b**S(2) - S(4)*a*c)*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - ((S(3)*b**S(3) - S(13)*a*b*c + (S(3)*b**S(4) - S(19)*a*b**S(2)*c + S(20)*a**S(2)*c**S(2))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*c**(S(5)/S(2))*(b**S(2) - S(4)*a*c)*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(9)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(10), -((b*x**S(2))/(S(2)*c*(b**S(2) - S(4)*a*c))) + (x**S(4)*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + (b*(b**S(2) - S(6)*a*c)*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*c**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + Log(a + b*x**S(2) + c*x**S(4))/(S(4)*c**S(2))],
[x**S(8)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(6), -((b*x)/(S(2)*c*(b**S(2) - S(4)*a*c))) + (x**S(3)*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + ((b**S(2) - S(6)*a*c - (b*(b**S(2) - S(8)*a*c))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*c**(S(3)/S(2))*(b**S(2) - S(4)*a*c)*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) + ((b**S(2) - S(6)*a*c + (b*(b**S(2) - S(8)*a*c))/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*c**(S(3)/S(2))*(b**S(2) - S(4)*a*c)*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(7)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), (x**S(2)*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + (S(2)*a*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(6)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), (x*(S(2)*a + b*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + ((b - (b**S(2) + S(4)*a*c)/Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*Sqrt(c)*(b**S(2) - S(4)*a*c)*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) + ((b**S(2) + S(4)*a*c + b*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*Sqrt(c)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(5)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), (S(2)*a + b*x**S(2))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) - (b*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(4)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), -((x*(b + S(2)*c*x**S(2)))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4)))) + (Sqrt(c)*(S(2)*b - Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - (Sqrt(c)*(S(2)*b + Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(Sqrt(S(2))*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(3)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), -((b + S(2)*c*x**S(2))/(S(2)*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4)))) + (S(2)*c*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(b**S(2) - S(4)*a*c)**(S(3)/S(2))],
[x**S(2)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(5), (x*(b**S(2) - S(2)*a*c + b*c*x**S(2)))/(S(2)*a*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + (Sqrt(c)*(b**S(2) - S(12)*a*c + b*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - (Sqrt(c)*(b**S(2) - S(12)*a*c - b*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[x**S(1)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(10), (b**S(2) - S(2)*a*c + b*c*x**S(2))/(S(2)*a*(b**S(2) - S(4)*a*c)*(a + b*x**S(2) + c*x**S(4))) + (b*(b**S(2) - S(6)*a*c)*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*a**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + Log(x)/a**S(2) - Log(a + b*x**S(2) + c*x**S(4))/(S(4)*a**S(2))],
[x**S(0)/(a*x + b*x**S(3) + c*x**S(5))**S(2), x, S(6), -((S(3)*b**S(2) - S(10)*a*c)/(S(2)*a**S(2)*(b**S(2) - S(4)*a*c)*x)) + (b**S(2) - S(2)*a*c + b*c*x**S(2))/(S(2)*a*(b**S(2) - S(4)*a*c)*x*(a + b*x**S(2) + c*x**S(4))) - (Sqrt(c)*(S(3)*b**S(3) - S(16)*a*b*c + (S(3)*b**S(2) - S(10)*a*c)*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) + (Sqrt(c)*(S(3)*b**S(3) - S(16)*a*b*c - (S(3)*b**S(2) - S(10)*a*c)*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a**S(2)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[S(1)/(x**S(1)*(a*x + b*x**S(3) + c*x**S(5))**S(2)), x, S(10), -((b**S(2) - S(3)*a*c)/(a**S(2)*(b**S(2) - S(4)*a*c)*x**S(2))) + (b**S(2) - S(2)*a*c + b*c*x**S(2))/(S(2)*a*(b**S(2) - S(4)*a*c)*x**S(2)*(a + b*x**S(2) + c*x**S(4))) - ((b**S(4) - S(6)*a*b**S(2)*c + S(6)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(a**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) - (S(2)*b*Log(x))/a**S(3) + (b*Log(a + b*x**S(2) + c*x**S(4)))/(S(2)*a**S(3))],
[S(1)/(x**S(2)*(a*x + b*x**S(3) + c*x**S(5))**S(2)), x, S(7), -((S(5)*b**S(2) - S(14)*a*c)/(S(6)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(3))) + (b*(S(5)*b**S(2) - S(19)*a*c))/(S(2)*a**S(3)*(b**S(2) - S(4)*a*c)*x) + (b**S(2) - S(2)*a*c + b*c*x**S(2))/(S(2)*a*(b**S(2) - S(4)*a*c)*x**S(3)*(a + b*x**S(2) + c*x**S(4))) + (Sqrt(c)*(S(5)*b**S(4) - S(29)*a*b**S(2)*c + S(28)*a**S(2)*c**S(2) + b*(S(5)*b**S(2) - S(19)*a*c)*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b - Sqrt(b**S(2) - S(4)*a*c))) - (Sqrt(c)*(S(5)*b**S(4) - S(29)*a*b**S(2)*c + S(28)*a**S(2)*c**S(2) - b*(S(5)*b**S(2) - S(19)*a*c)*Sqrt(b**S(2) - S(4)*a*c))*ArcTan((Sqrt(S(2))*Sqrt(c)*x)/Sqrt(b + Sqrt(b**S(2) - S(4)*a*c))))/(S(2)*Sqrt(S(2))*a**S(3)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))*Sqrt(b + Sqrt(b**S(2) - S(4)*a*c)))],
[S(1)/(x**S(3)*(a*x + b*x**S(3) + c*x**S(5))**S(2)), x, S(10), -((S(3)*b**S(2) - S(8)*a*c)/(S(4)*a**S(2)*(b**S(2) - S(4)*a*c)*x**S(4))) + (b*(S(3)*b**S(2) - S(11)*a*c))/(S(2)*a**S(3)*(b**S(2) - S(4)*a*c)*x**S(2)) + (b**S(2) - S(2)*a*c + b*c*x**S(2))/(S(2)*a*(b**S(2) - S(4)*a*c)*x**S(4)*(a + b*x**S(2) + c*x**S(4))) + (b*(S(3)*b**S(4) - S(20)*a*b**S(2)*c + S(30)*a**S(2)*c**S(2))*ArcTanh((b + S(2)*c*x**S(2))/Sqrt(b**S(2) - S(4)*a*c)))/(S(2)*a**S(4)*(b**S(2) - S(4)*a*c)**(S(3)/S(2))) + ((S(3)*b**S(2) - S(2)*a*c)*Log(x))/a**S(4) - ((S(3)*b**S(2) - S(2)*a*c)*Log(a + b*x**S(2) + c*x**S(4)))/(S(4)*a**S(4))],


## ::Subsection::Closed:: *)
##x**m (a x+b x**S(3)+c x**S(5))**(p/S(2))*)


[x/Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(3), (S(2)*x**S(2)*Sqrt(S(1) + (S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c)))*Sqrt(S(1) + (S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))*AppellF1(S(3)/S(4), S(1)/S(2), S(1)/S(2), S(7)/S(4), -((S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c))), -((S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))))/(S(3)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


## ::Subsection::Closed:: *)
##x**(m/S(2)) (a x+b x**S(3)+c x**S(5))**(p/S(2))*)


## ::Subsubsection::Closed:: *)
##p>S(0)*)


[x**(S(3)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(5), -((S(2)*(b**S(2) - S(3)*a*c)*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(S(15)*c**(S(3)/S(2))*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))) + (Sqrt(x)*(b + S(3)*c*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(15)*c) + (S(2)*a**(S(1)/S(4))*(b**S(2) - S(3)*a*c)*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(15)*c**(S(7)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (a**(S(1)/S(4))*(S(2)*b**S(2) + Sqrt(a)*b*Sqrt(c) - S(6)*a*c)*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(30)*c**(S(7)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(x)*Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(5), ((b + S(2)*c*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(8)*c*Sqrt(x)) - ((b**S(2) - S(4)*a*c)*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((b + S(2)*c*x**S(2))/(S(2)*Sqrt(c)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(16)*c**(S(3)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(a*x + b*x**S(3) + c*x**S(5))/Sqrt(x), x, S(5), (b*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(S(3)*Sqrt(c)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (S(1)/S(3))*Sqrt(x)*Sqrt(a*x + b*x**S(3) + c*x**S(5)) - (a**(S(1)/S(4))*b*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(3)*c**(S(3)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (a**(S(1)/S(4))*(b + S(2)*Sqrt(a)*Sqrt(c))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(6)*c**(S(3)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(a*x + b*x**S(3) + c*x**S(5))/x**(S(3)/S(2)), x, S(8), Sqrt(a*x + b*x**S(3) + c*x**S(5))/(S(2)*Sqrt(x)) - (Sqrt(a)*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((S(2)*a + b*x**S(2))/(S(2)*Sqrt(a)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(2)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (b*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((b + S(2)*c*x**S(2))/(S(2)*Sqrt(c)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(4)*Sqrt(c)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


[x**(S(3)/S(2))*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)), x, S(8), ((S(15)*b**S(4) - S(100)*a*b**S(2)*c + S(128)*a**S(2)*c**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(1280)*c**S(3)*Sqrt(x)) - (x**(S(3)/S(2))*(b*(S(5)*b**S(2) - S(4)*a*c) + S(4)*c*(S(5)*b**S(2) - S(16)*a*c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(640)*c**S(2)) + (Sqrt(x)*(S(3)*b + S(8)*c*x**S(2))*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)))/(S(80)*c) - (S(3)*b*(b**S(2) - S(4)*a*c)**S(2)*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((b + S(2)*c*x**S(2))/(S(2)*Sqrt(c)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(512)*c**(S(7)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(x)*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)), x, S(6), ((S(8)*b**S(4) - S(57)*a*b**S(2)*c + S(84)*a**S(2)*c**S(2))*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(S(315)*c**(S(5)/S(2))*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (Sqrt(x)*(b*(S(4)*b**S(2) - S(9)*a*c) + S(6)*c*(S(2)*b**S(2) - S(7)*a*c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(315)*c**S(2)) + ((S(3)*b + S(7)*c*x**S(2))*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)))/(S(63)*c*Sqrt(x)) - (a**(S(1)/S(4))*(S(8)*b**S(4) - S(57)*a*b**S(2)*c + S(84)*a**S(2)*c**S(2))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(315)*c**(S(11)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (a**(S(1)/S(4))*(S(8)*b**S(4) - S(57)*a*b**S(2)*c + S(84)*a**S(2)*c**S(2) + S(4)*Sqrt(a)*b*Sqrt(c)*(b**S(2) - S(6)*a*c))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(630)*c**(S(11)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2))/Sqrt(x), x, S(6), -((S(3)*(b**S(2) - S(4)*a*c)*(b + S(2)*c*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(128)*c**S(2)*Sqrt(x))) + ((b + S(2)*c*x**S(2))*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)))/(S(16)*c*x**(S(3)/S(2))) + (S(3)*(b**S(2) - S(4)*a*c)**S(2)*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((b + S(2)*c*x**S(2))/(S(2)*Sqrt(c)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(256)*c**(S(5)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2))/x**(S(3)/S(2)), x, S(6), -((S(2)*b*(b**S(2) - S(8)*a*c)*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(S(35)*c**(S(3)/S(2))*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))) + (Sqrt(x)*(b**S(2) + S(10)*a*c + S(3)*b*c*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(35)*c) + (a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2))/(S(7)*Sqrt(x)) + (S(2)*a**(S(1)/S(4))*b*(b**S(2) - S(8)*a*c)*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(35)*c**(S(7)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (a**(S(1)/S(4))*(Sqrt(a)*Sqrt(c)*(b**S(2) - S(20)*a*c) + S(2)*b*(b**S(2) - S(8)*a*c))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(70)*c**(S(7)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


## ::Subsubsection::Closed:: *)
##p<S(0)*)


[x**(S(3)/S(2))/Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(4), (Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((b + S(2)*c*x**S(2))/(S(2)*Sqrt(c)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(2)*Sqrt(c)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(x)/Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(2), (Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(2)*a**(S(1)/S(4))*c**(S(1)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[S(1)/(Sqrt(x)*Sqrt(a*x + b*x**S(3) + c*x**S(5))), x, S(4), -((Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((S(2)*a + b*x**S(2))/(S(2)*Sqrt(a)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(2)*Sqrt(a)*Sqrt(a*x + b*x**S(3) + c*x**S(5))))],
[S(1)/(x**(S(3)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))), x, S(6), (Sqrt(c)*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(a*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - Sqrt(a*x + b*x**S(3) + c*x**S(5))/(a*x**(S(3)/S(2))) - (c**(S(1)/S(4))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(a**(S(3)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (c**(S(1)/S(4))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(2)*a**(S(3)/S(4))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


[x**(S(3)/S(2))/(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)), x, S(5), (x**(S(3)/S(2))*(b**S(2) - S(2)*a*c + b*c*x**S(2)))/(a*(b**S(2) - S(4)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (b*Sqrt(c)*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(a*(b**S(2) - S(4)*a*c)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (b*c**(S(1)/S(4))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(a**(S(3)/S(4))*(b**S(2) - S(4)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (c**(S(1)/S(4))*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(2)*a**(S(3)/S(4))*(b - S(2)*Sqrt(a)*Sqrt(c))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[Sqrt(x)/(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2)), x, S(5), (Sqrt(x)*(b**S(2) - S(2)*a*c + b*c*x**S(2)))/(a*(b**S(2) - S(4)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((S(2)*a + b*x**S(2))/(S(2)*Sqrt(a)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(2)*a**(S(3)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[S(1)/(Sqrt(x)*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2))), x, S(6), (b**S(2) - S(2)*a*c + b*c*x**S(2))/(a*(b**S(2) - S(4)*a*c)*Sqrt(x)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (S(2)*Sqrt(c)*(b**S(2) - S(3)*a*c)*x**(S(3)/S(2))*(a + b*x**S(2) + c*x**S(4)))/(a**S(2)*(b**S(2) - S(4)*a*c)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - (S(2)*(b**S(2) - S(3)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(a**S(2)*(b**S(2) - S(4)*a*c)*x**(S(3)/S(2))) - (S(2)*c**(S(1)/S(4))*(b**S(2) - S(3)*a*c)*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticE(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(a**(S(7)/S(4))*(b**S(2) - S(4)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (c**(S(1)/S(4))*(S(2)*b**S(2) + Sqrt(a)*b*Sqrt(c) - S(6)*a*c)*Sqrt(x)*(Sqrt(a) + Sqrt(c)*x**S(2))*Sqrt((a + b*x**S(2) + c*x**S(4))/(Sqrt(a) + Sqrt(c)*x**S(2))**S(2))*EllipticF(S(2)*ArcTan((c**(S(1)/S(4))*x)/a**(S(1)/S(4))), (S(1)/S(4))*(S(2) - b/(Sqrt(a)*Sqrt(c)))))/(S(2)*a**(S(7)/S(4))*(b**S(2) - S(4)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],
[S(1)/(x**(S(3)/S(2))*(a*x + b*x**S(3) + c*x**S(5))**(S(3)/S(2))), x, S(7), (b**S(2) - S(2)*a*c + b*c*x**S(2))/(a*(b**S(2) - S(4)*a*c)*x**(S(3)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5))) - ((S(3)*b**S(2) - S(8)*a*c)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))/(S(2)*a**S(2)*(b**S(2) - S(4)*a*c)*x**(S(5)/S(2))) + (S(3)*b*Sqrt(x)*Sqrt(a + b*x**S(2) + c*x**S(4))*ArcTanh((S(2)*a + b*x**S(2))/(S(2)*Sqrt(a)*Sqrt(a + b*x**S(2) + c*x**S(4)))))/(S(4)*a**(S(5)/S(2))*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


[x**((S(3)*(n - S(1)))/S(2))/(a*x**(n - S(1)) + b*x**n + c*x**(n + S(1)))**(S(3)/S(2)), x, S(1), -((S(2)*x**((S(1)/S(2))*(-S(1) + n))*(b + S(2)*c*x))/((b**S(2) - S(4)*a*c)*Sqrt(a*x**(-S(1) + n) + b*x**n + c*x**(S(1) + n))))],


## ::Section::Closed:: *)
##Integrands of the form x**m (d+e x**S(2))**n (a x+b x**S(3)+c x**S(5))**p*)


## ::Subsection:: *)
##x**m (d+e x**S(2))**n (a x+b x**S(3)+c x**S(5))**p*)


## ::Subsection::Closed:: *)
##x**m (d+e x**S(2))**n (a x+b x**S(3)+c x**S(5))**(p/S(2))*)


[x*(d + e*x**S(2))/Sqrt(a*x + b*x**S(3) + c*x**S(5)), x, S(7), (S(2)*d*x**S(2)*Sqrt(S(1) + (S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c)))*Sqrt(S(1) + (S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))*AppellF1(S(3)/S(4), S(1)/S(2), S(1)/S(2), S(7)/S(4), -((S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c))), -((S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))))/(S(3)*Sqrt(a*x + b*x**S(3) + c*x**S(5))) + (S(2)*e*x**S(4)*Sqrt(S(1) + (S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c)))*Sqrt(S(1) + (S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))*AppellF1(S(7)/S(4), S(1)/S(2), S(1)/S(2), S(11)/S(4), -((S(2)*c*x**S(2))/(b - Sqrt(b**S(2) - S(4)*a*c))), -((S(2)*c*x**S(2))/(b + Sqrt(b**S(2) - S(4)*a*c)))))/(S(7)*Sqrt(a*x + b*x**S(3) + c*x**S(5)))],


## ::Subsection:: *)
##x**(m/S(2)) (d+e x**S(2))**n (a x+b x**S(3)+c x**S(5))**(p/S(2))*)


    ]

    for i in test:
	    try:
	        r = rubi_integrate(nsimplify(i[0]), i[1])

	        if len(i) == 5:
	            assert rubi_test(r, i[1], i[3], expand=True, _diff=True, _numerical=True) or rubi_test(r, i[1], i[4], expand=True, _diff=True, _numerical=True)
	        else:
	            assert rubi_test(r, i[1], i[3], expand=True, _diff=True, _numerical=True)
	    except:   
	        print(i)
	        print("Fail************\n")
