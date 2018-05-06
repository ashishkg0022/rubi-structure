from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Sum, Set, With, Module, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec, ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma, FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False
    from sympy import And, Or

def binomial_products(rubi):

    def cons_f1(b, n, x, p):
        return FreeQ(List(b, n, p), x)

    cons1 = CustomConstraint(cons_f1)

    def cons_f2(b, x):
        return FreeQ(b, x)

    cons2 = CustomConstraint(cons_f2)

    def cons_f3(n, x):
        return FreeQ(n, x)

    cons3 = CustomConstraint(cons_f3)

    def cons_f4(p, x):
        return FreeQ(p, x)

    cons4 = CustomConstraint(cons_f4)
    pattern1 = Pattern(Integral((x_**n_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons1)
    rule1 = ReplacementRule(pattern1, lambda b, n, x, p : b**IntPart(p)*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(n*p), x))
    rubi.add(rule1)


    def cons_f5(n, p):
        return ZeroQ(p + S(1) + S(1)/n)

    cons5 = CustomConstraint(cons_f5)

    def cons_f6(a, x):
        return FreeQ(a, x)

    cons6 = CustomConstraint(cons_f6)
    pattern2 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons3, cons4, cons5)
    rule2 = ReplacementRule(pattern2, lambda n, a, x, p, b : x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule2)


    def cons_f7(n, p):
        return NegativeIntegerQ(p + S(1) + S(1)/n)

    cons7 = CustomConstraint(cons_f7)

    def cons_f8(p):
        return NonzeroQ(p + S(1))

    cons8 = CustomConstraint(cons_f8)
    pattern3 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons3, cons4, cons7, cons8)
    rule3 = ReplacementRule(pattern3, lambda n, a, x, p, b : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule3)


    def cons_f9(n):
        return NonzeroQ(S(3)*n + S(1))

    cons9 = CustomConstraint(cons_f9)
    pattern4 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**S(2), x_), cons6, cons2, cons3, cons9)
    rule4 = ReplacementRule(pattern4, lambda b, n, a, x : Int(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n), x))
    rubi.add(rule4)


    def cons_f10(n):
        return RationalQ(n)

    cons10 = CustomConstraint(cons_f10)

    def cons_f11(n):
        return Less(n, S(0))

    cons11 = CustomConstraint(cons_f11)

    def cons_f12(p):
        return IntegerQ(p)

    cons12 = CustomConstraint(cons_f12)
    pattern5 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons10, cons11, cons12)
    rule5 = ReplacementRule(pattern5, lambda n, a, x, p, b : Int(x**(n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule5)


    def cons_f13(n, p):
        return PositiveIntegerQ(n, p)

    cons13 = CustomConstraint(cons_f13)
    pattern6 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons13)
    rule6 = ReplacementRule(pattern6, lambda n, a, x, p, b : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule6)


    def cons_f14(n):
        return PositiveIntegerQ(n)

    cons14 = CustomConstraint(cons_f14)

    def cons_f15(p):
        return RationalQ(p)

    cons15 = CustomConstraint(cons_f15)

    def cons_f16(p):
        return Greater(p, S(0))

    cons16 = CustomConstraint(cons_f16)

    def cons_f17(n, p):
        return Or(IntegerQ(S(2)*p), And(Equal(n, S(2)), IntegerQ(S(4)*p)), And(Equal(n, S(2)), IntegerQ(S(3)*p)), Less(Denominator(p + S(1)/n), Denominator(p)))

    cons17 = CustomConstraint(cons_f17)
    pattern7 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons16, cons17)
    rule7 = ReplacementRule(pattern7, lambda n, a, x, p, b : a*n*p*Int((a + b*x**n)**(p + S(-1)), x)/(n*p + S(1)) + x*(a + b*x**n)**p/(n*p + S(1)))
    rubi.add(rule7)


    def cons_f18(b, a):
        return PosQ(b/a)

    cons18 = CustomConstraint(cons_f18)

    def cons_f19(a):
        return PositiveQ(a)

    cons19 = CustomConstraint(cons_f19)
    pattern8 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons6, cons2, cons18, cons19)
    rule8 = ReplacementRule(pattern8, lambda b, a, x : S(2)*EllipticE(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(5)/4)*Rt(b/a, S(2))))
    rubi.add(rule8)


    def cons_f20(a):
        return Not(PositiveQ(a))

    cons20 = CustomConstraint(cons_f20)
    pattern9 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_), cons6, cons2, cons18, cons20)
    rule9 = ReplacementRule(pattern9, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-5)/4), x)/(a*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule9)


    def cons_f21(a, x, b):
        return FreeQ(List(a, b), x)

    cons21 = CustomConstraint(cons_f21)
    pattern10 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-7)/6), x_), cons6, cons2, cons21)
    rule10 = ReplacementRule(pattern10, lambda b, a, x : Subst(Int((-b*x**S(2) + S(1))**(S(-1)/3), x), x, x/sqrt(a + b*x**S(2)))/((a/(a + b*x**S(2)))**(S(2)/3)*(a + b*x**S(2))**(S(2)/3)))
    rubi.add(rule10)


    def cons_f22(p):
        return Less(p, S(-1))

    cons22 = CustomConstraint(cons_f22)
    pattern11 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons22, cons17)
    rule11 = ReplacementRule(pattern11, lambda n, a, x, p, b : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(S(1)/(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons21)
    rule12 = ReplacementRule(pattern12, lambda b, a, x : Int((-x*Rt(b, S(3)) + S(2)*Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))**S(2)) + Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))**S(2)))
    rubi.add(rule12)


    def cons_f23(n):
        return PositiveIntegerQ(n/S(2) + S(-3)/2)

    cons23 = CustomConstraint(cons_f23)

    def cons_f24(a, b):
        return PosQ(a/b)

    cons24 = CustomConstraint(cons_f24)
    pattern13 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons23, cons24, )
    def With13(b, n, a, x):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r + s*x), x)/(a*n)
    rule13 = ReplacementRule(pattern13, lambda b, n, a, x : With13(b, n, a, x))
    rubi.add(rule13)


    def cons_f25(a, b):
        return NegQ(a/b)

    cons25 = CustomConstraint(cons_f25)
    pattern14 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons23, cons25, )
    def With14(b, n, a, x):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r - s*x), x)/(a*n)
    rule14 = ReplacementRule(pattern14, lambda b, n, a, x : With14(b, n, a, x))
    rubi.add(rule14)


    def cons_f26(a, b):
        return Or(PositiveQ(a), PositiveQ(b))

    cons26 = CustomConstraint(cons_f26)
    pattern15 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons24, cons26)
    rule15 = ReplacementRule(pattern15, lambda b, a, x : ArcTan(x*Rt(b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(b, S(2))))
    rubi.add(rule15)


    def cons_f27(a, b):
        return Or(NegativeQ(a), NegativeQ(b))

    cons27 = CustomConstraint(cons_f27)
    pattern16 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons24, cons27)
    rule16 = ReplacementRule(pattern16, lambda b, a, x : -ArcTan(x*Rt(-b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(-b, S(2))))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons24)
    rule17 = ReplacementRule(pattern17, lambda b, a, x : ArcTan(x/Rt(a/b, S(2)))*Rt(a/b, S(2))/a)
    rubi.add(rule17)


    def cons_f28(a, b):
        return Or(PositiveQ(a), NegativeQ(b))

    cons28 = CustomConstraint(cons_f28)
    pattern18 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons25, cons28)
    rule18 = ReplacementRule(pattern18, lambda b, a, x : atanh(x*Rt(-b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(-b, S(2))))
    rubi.add(rule18)


    def cons_f29(a, b):
        return Or(NegativeQ(a), PositiveQ(b))

    cons29 = CustomConstraint(cons_f29)
    pattern19 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons25, cons29)
    rule19 = ReplacementRule(pattern19, lambda b, a, x : -atanh(x*Rt(b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(b, S(2))))
    rubi.add(rule19)

    pattern20 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons25)
    rule20 = ReplacementRule(pattern20, lambda b, a, x : Rt(-a/b, S(2))*atanh(x/Rt(-a/b, S(2)))/a)
    rubi.add(rule20)


    def cons_f30(n):
        return PositiveIntegerQ(n/S(4) + S(-1)/2)

    cons30 = CustomConstraint(cons_f30)
    pattern21 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons30, cons24, )
    def With21(b, n, a, x):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        v = Symbol('v')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n)
    rule21 = ReplacementRule(pattern21, lambda b, n, a, x : With21(b, n, a, x))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons30, cons25, )
    def With22(b, n, a, x):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(2*Pi*k/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r + s*x*cos(2*Pi*k/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    rule22 = ReplacementRule(pattern22, lambda b, n, a, x : With22(b, n, a, x))
    rubi.add(rule22)


    def cons_f31(a, b):
        return Or(PositiveQ(a/b), And(PosQ(a/b), AtomQ(SplitProduct(SumBaseQ, a)), AtomQ(SplitProduct(SumBaseQ, b))))

    cons31 = CustomConstraint(cons_f31)
    pattern23 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons31, )
    def With23(b, a, x):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r)
    rule23 = ReplacementRule(pattern23, lambda b, a, x : With23(b, a, x))
    rubi.add(rule23)


    def cons_f32(a, b):
        return Not(PositiveQ(a/b))

    cons32 = CustomConstraint(cons_f32)
    pattern24 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons32, )
    def With24(b, a, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*a)
    rule24 = ReplacementRule(pattern24, lambda b, a, x : With24(b, a, x))
    rubi.add(rule24)


    def cons_f33(n):
        return PositiveIntegerQ(n/S(4) + S(-1))

    cons33 = CustomConstraint(cons_f33)

    def cons_f34(a, b):
        return PositiveQ(a/b)

    cons34 = CustomConstraint(cons_f34)
    pattern25 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons33, cons34, )
    def With25(b, n, a, x):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*r*Int((sqrt(S(2))*r - s*x**(n/S(4)))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a) + sqrt(S(2))*r*Int((sqrt(S(2))*r + s*x**(n/S(4)))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a)
    rule25 = ReplacementRule(pattern25, lambda b, n, a, x : With25(b, n, a, x))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons33, cons32, )
    def With26(b, n, a, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**(n/S(2))), x)/(S(2)*a)
    rule26 = ReplacementRule(pattern26, lambda b, n, a, x : With26(b, n, a, x))
    rubi.add(rule26)


    def cons_f35(b):
        return PosQ(b)

    cons35 = CustomConstraint(cons_f35)
    pattern27 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons19, cons35)
    rule27 = ReplacementRule(pattern27, lambda b, a, x : asinh(x*Rt(b, S(2))/sqrt(a))/Rt(b, S(2)))
    rubi.add(rule27)


    def cons_f36(b):
        return NegQ(b)

    cons36 = CustomConstraint(cons_f36)
    pattern28 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons19, cons36)
    rule28 = ReplacementRule(pattern28, lambda b, a, x : asin(x*Rt(-b, S(2))/sqrt(a))/Rt(-b, S(2)))
    rubi.add(rule28)

    pattern29 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons20)
    rule29 = ReplacementRule(pattern29, lambda b, a, x : Subst(Int(S(1)/(-b*x**S(2) + S(1)), x), x, x/sqrt(a + b*x**S(2))))
    rubi.add(rule29)


    def cons_f37(a):
        return PosQ(a)

    cons37 = CustomConstraint(cons_f37)
    pattern30 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons37, )
    def With30(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(-sqrt(S(3)) + S(1)))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(S(3)*r*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3)))
    rule30 = ReplacementRule(pattern30, lambda b, a, x : With30(b, a, x))
    rubi.add(rule30)


    def cons_f38(a):
        return NegQ(a)

    cons38 = CustomConstraint(cons_f38)
    pattern31 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons38, )
    def With31(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(-sqrt(S(3)) + S(1)))), S(-7) + S(4)*sqrt(S(3)))/(S(3)*r*sqrt(-s*(r*x + s)/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(a + b*x**S(3)))
    rule31 = ReplacementRule(pattern31, lambda b, a, x : With31(b, a, x))
    rubi.add(rule31)

    pattern32 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons18, )
    def With32(b, a, x):
        q = Rt(b/a, S(4))
        return sqrt((a + b*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), S(1)/2)/(S(2)*q*sqrt(a + b*x**S(4)))
    rule32 = ReplacementRule(pattern32, lambda b, a, x : With32(b, a, x))
    rubi.add(rule32)


    def cons_f39(b, a):
        return NegQ(b/a)

    cons39 = CustomConstraint(cons_f39)
    pattern33 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons39, cons19)
    rule33 = ReplacementRule(pattern33, lambda b, a, x : EllipticF(asin(x*Rt(-b, S(4))/Rt(a, S(4))), S(-1))/(Rt(a, S(4))*Rt(-b, S(4))))
    rubi.add(rule33)


    def cons_f40(a):
        return NegativeQ(a)

    cons40 = CustomConstraint(cons_f40)

    def cons_f41(b):
        return PositiveQ(b)

    cons41 = CustomConstraint(cons_f41)
    pattern34 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons40, cons41, )
    def With34(b, a, x):
        q = Rt(-a*b, S(2))
        if IntegerQ(q):
            return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(-a)*sqrt(a + b*x**S(4)))
        print("Unable to Integrate")
    rule34 = ReplacementRule(pattern34, lambda b, a, x : With34(b, a, x))
    rubi.add(rule34)

    pattern35 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons40, cons41, )
    def With35(b, a, x):
        q = Rt(-a*b, S(2))
        return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(a/(a + q*x**S(2)))*sqrt(a + b*x**S(4)))
    rule35 = ReplacementRule(pattern35, lambda b, a, x : With35(b, a, x))
    rubi.add(rule35)

    pattern36 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons39, cons20)
    rule36 = ReplacementRule(pattern36, lambda b, a, x : sqrt(S(1) + b*x**S(4)/a)*Int(S(1)/sqrt(S(1) + b*x**S(4)/a), x)/sqrt(a + b*x**S(4)))
    rubi.add(rule36)

    pattern37 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons6, cons2, cons21, )
    def With37(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(3)**(S(3)/4)*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticF(acos((r*x**S(2)*(-sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(6)*s*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6)))
    rule37 = ReplacementRule(pattern37, lambda b, a, x : With37(b, a, x))
    rubi.add(rule37)

    pattern38 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons6, cons2, cons21)
    rule38 = ReplacementRule(pattern38, lambda b, a, x : Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons6, cons2, cons18)
    rule39 = ReplacementRule(pattern39, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-5)/4), x) + S(2)*x/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons6, cons2, cons39, cons19)
    rule40 = ReplacementRule(pattern40, lambda b, a, x : S(2)*EllipticE(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(1)/4)*Rt(-b/a, S(2))))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_), cons6, cons2, cons39, cons20)
    rule41 = ReplacementRule(pattern41, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-1)/4), x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons6, cons2, cons19, cons18)
    rule42 = ReplacementRule(pattern42, lambda b, a, x : S(2)*EllipticF(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(b/a, S(2))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons6, cons2, cons19, cons39)
    rule43 = ReplacementRule(pattern43, lambda b, a, x : S(2)*EllipticF(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(-b/a, S(2))))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_), cons6, cons2, cons20)
    rule44 = ReplacementRule(pattern44, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(3)/4)*Int((S(1) + b*x**S(2)/a)**(S(-3)/4), x)/(a + b*x**S(2))**(S(3)/4))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/3), x_), cons6, cons2, cons21)
    rule45 = ReplacementRule(pattern45, lambda b, a, x : S(3)*sqrt(b*x**S(2))*Subst(Int(x/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-2)/3), x_), cons6, cons2, cons21)
    rule46 = ReplacementRule(pattern46, lambda b, a, x : S(3)*sqrt(b*x**S(2))*Subst(Int(S(1)/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(-3)/4), x_), cons6, cons2, cons21)
    rule47 = ReplacementRule(pattern47, lambda b, a, x : x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)), x)/(a + b*x**S(4))**(S(3)/4))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/6), x_), cons6, cons2, cons21)
    rule48 = ReplacementRule(pattern48, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-7)/6), x)/S(2) + S(3)*x/(S(2)*(a + b*x**S(2))**(S(1)/6)))
    rubi.add(rule48)


    def cons_f42(p):
        return Less(S(-1), p, S(0))

    cons42 = CustomConstraint(cons_f42)

    def cons_f43(p):
        return Unequal(p, S(-1)/2)

    cons43 = CustomConstraint(cons_f43)

    def cons_f44(n, p):
        return IntegerQ(p + S(1)/n)

    cons44 = CustomConstraint(cons_f44)
    pattern49 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons42, cons43, cons44)
    rule49 = ReplacementRule(pattern49, lambda n, a, x, p, b : a**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule49)


    def cons_f45(n, p):
        return Less(Denominator(p + S(1)/n), Denominator(p))

    cons45 = CustomConstraint(cons_f45)
    pattern50 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons42, cons43, cons45)
    rule50 = ReplacementRule(pattern50, lambda n, a, x, p, b : (a/(a + b*x**n))**(p + S(1)/n)*(a + b*x**n)**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule50)


    def cons_f46(n):
        return NegativeIntegerQ(n)

    cons46 = CustomConstraint(cons_f46)
    pattern51 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons4, cons46)
    rule51 = ReplacementRule(pattern51, lambda n, a, x, p, b : -Subst(Int((a + b*x**(-n))**p/x**S(2), x), x, S(1)/x))
    rubi.add(rule51)


    def cons_f47(n):
        return FractionQ(n)

    cons47 = CustomConstraint(cons_f47)
    pattern52 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons4, cons47, )
    def With52(n, a, x, p, b):
        k = Denominator(n)
        return k*Subst(Int(x**(k + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    rule52 = ReplacementRule(pattern52, lambda n, a, x, p, b : With52(n, a, x, p, b))
    rubi.add(rule52)


    def cons_f48(p):
        return PositiveIntegerQ(p)

    cons48 = CustomConstraint(cons_f48)
    pattern53 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons3, cons48)
    rule53 = ReplacementRule(pattern53, lambda n, a, x, p, b : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule53)


    def cons_f49(p):
        return Not(PositiveIntegerQ(p))

    cons49 = CustomConstraint(cons_f49)

    def cons_f50(n):
        return Not(IntegerQ(S(1)/n))

    cons50 = CustomConstraint(cons_f50)

    def cons_f51(n, p):
        return Not(NegativeIntegerQ(p + S(1)/n))

    cons51 = CustomConstraint(cons_f51)

    def cons_f52(a, p):
        return Or(IntegerQ(p), PositiveQ(a))

    cons52 = CustomConstraint(cons_f52)
    pattern54 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons3, cons4, cons49, cons50, cons51, cons52)
    rule54 = ReplacementRule(pattern54, lambda n, a, x, p, b : a**p*x*Hypergeometric2F1(-p, S(1)/n, S(1) + S(1)/n, -b*x**n/a))
    rubi.add(rule54)


    def cons_f53(a, p):
        return Not(Or(IntegerQ(p), PositiveQ(a)))

    cons53 = CustomConstraint(cons_f53)
    pattern55 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons3, cons4, cons49, cons50, cons51, cons53)
    rule55 = ReplacementRule(pattern55, lambda n, a, x, p, b : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p, x))
    rubi.add(rule55)


    def cons_f54(u, x):
        return LinearQ(u, x)

    cons54 = CustomConstraint(cons_f54)

    def cons_f55(u, x):
        return NonzeroQ(u - x)

    cons55 = CustomConstraint(cons_f55)
    pattern56 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons6, cons2, cons3, cons4, cons54, cons55)
    rule56 = ReplacementRule(pattern56, lambda u, n, a, x, p, b : Subst(Int((a + b*x**n)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule56)


    def cons_f56(a1, a2, b1, b2):
        return ZeroQ(a1*b2 + a2*b1)

    cons56 = CustomConstraint(cons_f56)

    def cons_f57(a1, a2, p):
        return Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2)))

    cons57 = CustomConstraint(cons_f57)

    def cons_f58(a1, x):
        return FreeQ(a1, x)

    cons58 = CustomConstraint(cons_f58)

    def cons_f59(b1, x):
        return FreeQ(b1, x)

    cons59 = CustomConstraint(cons_f59)

    def cons_f60(a2, x):
        return FreeQ(a2, x)

    cons60 = CustomConstraint(cons_f60)

    def cons_f61(b2, x):
        return FreeQ(b2, x)

    cons61 = CustomConstraint(cons_f61)
    pattern57 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**WC('p', S(1))*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**WC('p', S(1)), x_), cons58, cons59, cons60, cons61, cons3, cons4, cons56, cons57)
    rule57 = ReplacementRule(pattern57, lambda a1, a2, n, p, x, b2, b1 : Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule57)


    def cons_f62(n):
        return PositiveIntegerQ(S(2)*n)

    cons62 = CustomConstraint(cons_f62)

    def cons_f63(n, p):
        return Or(IntegerQ(S(2)*p), Less(Denominator(p + S(1)/n), Denominator(p)))

    cons63 = CustomConstraint(cons_f63)
    pattern58 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_), cons58, cons59, cons60, cons61, cons56, cons62, cons15, cons16, cons63)
    rule58 = ReplacementRule(pattern58, lambda a1, a2, n, p, x, b2, b1 : S(2)*a1*a2*n*p*Int((a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(S(2)*n*p + S(1)) + x*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(S(2)*n*p + S(1)))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons56, cons62, cons15, cons22, cons63)
    rule59 = ReplacementRule(pattern59, lambda a1, a2, n, p, x, b2, b1 : -x*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*n*(p + S(1))) + (S(2)*n*(p + S(1)) + S(1))*Int((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))))
    rubi.add(rule59)


    def cons_f64(n):
        return NegativeIntegerQ(S(2)*n)

    cons64 = CustomConstraint(cons_f64)
    pattern60 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons4, cons56, cons64)
    rule60 = ReplacementRule(pattern60, lambda a1, a2, n, p, x, b2, b1 : -Subst(Int((a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p/x**S(2), x), x, S(1)/x))
    rubi.add(rule60)


    def cons_f65(n):
        return FractionQ(S(2)*n)

    cons65 = CustomConstraint(cons_f65)
    pattern61 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons4, cons56, cons65, )
    def With61(a1, a2, n, p, x, b2, b1):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    rule61 = ReplacementRule(pattern61, lambda a1, a2, n, p, x, b2, b1 : With61(a1, a2, n, p, x, b2, b1))
    rubi.add(rule61)


    def cons_f66(p):
        return Not(IntegerQ(p))

    cons66 = CustomConstraint(cons_f66)
    pattern62 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**p_*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**p_, x_), cons58, cons59, cons60, cons61, cons3, cons4, cons56, cons66)
    rule62 = ReplacementRule(pattern62, lambda a1, a2, n, p, x, b2, b1 : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule62)


    def cons_f67(c, x):
        return FreeQ(c, x)

    cons67 = CustomConstraint(cons_f67)

    def cons_f68(m, x):
        return FreeQ(m, x)

    cons68 = CustomConstraint(cons_f68)
    pattern63 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons57)
    rule63 = ReplacementRule(pattern63, lambda a1, a2, n, m, p, x, b2, c, b1 : Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule63)


    def cons_f69(c, m):
        return Or(IntegerQ(m), PositiveQ(c))

    cons69 = CustomConstraint(cons_f69)

    def cons_f70(n, m):
        return IntegerQ((m + S(1))/n)

    cons70 = CustomConstraint(cons_f70)
    pattern64 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_, x_), cons2, cons67, cons68, cons3, cons4, cons69, cons70)
    rule64 = ReplacementRule(pattern64, lambda n, m, x, p, b, c : b**(S(1) - (m + S(1))/n)*c**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n), x), x, x**n)/n)
    rubi.add(rule64)


    def cons_f71(n, m):
        return Not(IntegerQ((m + S(1))/n))

    cons71 = CustomConstraint(cons_f71)
    pattern65 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons67, cons68, cons3, cons4, cons69, cons71)
    rule65 = ReplacementRule(pattern65, lambda n, m, x, p, b, c : b**IntPart(p)*c**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p), x))
    rubi.add(rule65)


    def cons_f72(m):
        return Not(IntegerQ(m))

    cons72 = CustomConstraint(cons_f72)
    pattern66 = Pattern(Integral((c_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_), cons2, cons67, cons68, cons3, cons4, cons72)
    rule66 = ReplacementRule(pattern66, lambda n, m, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(b*x**n)**p, x))
    rubi.add(rule66)


    def cons_f73(n):
        return NegQ(n)

    cons73 = CustomConstraint(cons_f73)
    pattern67 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons3, cons12, cons73)
    rule67 = ReplacementRule(pattern67, lambda n, m, a, x, p, b : Int(x**(m + n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule67)


    def cons_f74(n, m, p):
        return ZeroQ(p + S(1) + (m + S(1))/n)

    cons74 = CustomConstraint(cons_f74)

    def cons_f75(m):
        return NonzeroQ(m + S(1))

    cons75 = CustomConstraint(cons_f75)
    pattern68 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons74, cons75)
    rule68 = ReplacementRule(pattern68, lambda n, m, a, x, p, b, c : (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule68)


    def cons_f76(n, m, p):
        return ZeroQ(p + S(1) + (m + S(1))/(S(2)*n))

    cons76 = CustomConstraint(cons_f76)
    pattern69 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons76, cons75)
    rule69 = ReplacementRule(pattern69, lambda a1, a2, n, m, p, x, b2, c, b1 : (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule69)

    pattern70 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons6, cons2, cons68, cons3, cons4, cons70)
    rule70 = ReplacementRule(pattern70, lambda n, m, a, x, p, b : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p, x), x, x**n)/n)
    rubi.add(rule70)


    def cons_f77(n, m):
        return IntegerQ((m + S(1))/(S(2)*n))

    cons77 = CustomConstraint(cons_f77)
    pattern71 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons3, cons4, cons56, cons77)
    rule71 = ReplacementRule(pattern71, lambda a1, a2, n, m, p, x, b2, b1 : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a1 + b1*x)**p*(a2 + b2*x)**p, x), x, x**n)/n)
    rubi.add(rule71)

    pattern72 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons70)
    rule72 = ReplacementRule(pattern72, lambda n, m, a, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons77)
    rule73 = ReplacementRule(pattern73, lambda a1, a2, n, m, p, x, b2, c, b1 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons68, cons3, cons48)
    rule74 = ReplacementRule(pattern74, lambda n, m, a, x, p, b, c : Int(ExpandIntegrand((c*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule74)


    def cons_f78(n, m, p):
        return NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n)

    cons78 = CustomConstraint(cons_f78)
    pattern75 = Pattern(Integral(x_**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons3, cons4, cons78, cons75)
    rule75 = ReplacementRule(pattern75, lambda n, m, a, x, p, b : -b*(m + n*(p + S(1)) + S(1))*Int(x**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + x**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*(m + S(1))))
    rubi.add(rule75)


    def cons_f79(n, m, p):
        return NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n))

    cons79 = CustomConstraint(cons_f79)
    pattern76 = Pattern(Integral(x_**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons3, cons4, cons56, cons79, cons75)
    rule76 = ReplacementRule(pattern76, lambda a1, a2, n, m, p, x, b2, b1 : -b1*b2*(m + S(2)*n*(p + S(1)) + S(1))*Int(x**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + x**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*(m + S(1))))
    rubi.add(rule76)

    pattern77 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons78, cons8)
    rule77 = ReplacementRule(pattern77, lambda n, m, a, x, p, b, c : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule77)

    pattern78 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons79, cons8)
    rule78 = ReplacementRule(pattern78, lambda a1, a2, n, m, p, x, b2, c, b1 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule78)


    def cons_f80(m):
        return IntegerQ(m)

    cons80 = CustomConstraint(cons_f80)
    pattern79 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons4, cons14, cons80, )
    def With79(n, m, a, x, p, b):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p, x), x, x**k)/k
        print("Unable to Integrate")
    rule79 = ReplacementRule(pattern79, lambda n, m, a, x, p, b : With79(n, m, a, x, p, b))
    rubi.add(rule79)

    pattern80 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons4, cons56, cons62, cons80, )
    def With80(a1, a2, n, m, p, x, b2, b1):
        k = GCD(m + S(1), S(2)*n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a1 + b1*x**(n/k))**p*(a2 + b2*x**(n/k))**p, x), x, x**k)/k
        print("Unable to Integrate")
    rule80 = ReplacementRule(pattern80, lambda a1, a2, n, m, p, x, b2, b1 : With80(a1, a2, n, m, p, x, b2, b1))
    rubi.add(rule80)


    def cons_f81(m, p):
        return RationalQ(m, p)

    cons81 = CustomConstraint(cons_f81)

    def cons_f82(m):
        return Less(m, S(-1))

    cons82 = CustomConstraint(cons_f82)

    def cons_f83(n, m, p):
        return Not(NegativeIntegerQ((m + n*p + n + S(1))/n))

    cons83 = CustomConstraint(cons_f83)

    def cons_f84(n, m, a, x, p, b, c):
        return IntBinomialQ(a, b, c, n, m, p, x)

    cons84 = CustomConstraint(cons_f84)
    pattern81 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons14, cons81, cons16, cons82, cons83, cons84)
    rule81 = ReplacementRule(pattern81, lambda n, m, a, x, p, b, c : -b*c**(-n)*n*p*Int((c*x)**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + S(1))))
    rubi.add(rule81)


    def cons_f85(n, m, p):
        return NonzeroQ(m + S(2)*n*p + S(1))

    cons85 = CustomConstraint(cons_f85)

    def cons_f86(a1, a2, n, m, p, x, b2, c, b1):
        return IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)

    cons86 = CustomConstraint(cons_f86)
    pattern82 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons56, cons62, cons81, cons16, cons85, cons86)
    rule82 = ReplacementRule(pattern82, lambda a1, a2, n, m, p, x, b2, c, b1 : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule82)


    def cons_f87(n, m, p):
        return NonzeroQ(m + n*p + S(1))

    cons87 = CustomConstraint(cons_f87)
    pattern83 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons14, cons81, cons16, cons87, cons84)
    rule83 = ReplacementRule(pattern83, lambda n, m, a, x, p, b, c : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule83)

    pattern84 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons18)
    rule84 = ReplacementRule(pattern84, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule84)


    def cons_f88(m):
        return PositiveIntegerQ(m/S(4) + S(-1)/2)

    cons88 = CustomConstraint(cons_f88)
    pattern85 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons18, cons88)
    rule85 = ReplacementRule(pattern85, lambda b, m, a, x : -a*(m + S(-3))*Int(x**(m + S(-4))/(a + b*x**S(4))**(S(5)/4), x)/(b*(m + S(-4))) + x**(m + S(-3))/(b*(a + b*x**S(4))**(S(1)/4)*(m + S(-4))))
    rubi.add(rule85)


    def cons_f89(m):
        return NegativeIntegerQ(m/S(4) + S(-1)/2)

    cons89 = CustomConstraint(cons_f89)
    pattern86 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons18, cons89)
    rule86 = ReplacementRule(pattern86, lambda b, m, a, x : -b*m*Int(x**(m + S(4))/(a + b*x**S(4))**(S(5)/4), x)/(a*(m + S(1))) + x**(m + S(1))/(a*(a + b*x**S(4))**(S(1)/4)*(m + S(1))))
    rubi.add(rule86)

    pattern87 = Pattern(Integral(sqrt(x_*WC('c', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons67, cons18)
    rule87 = ReplacementRule(pattern87, lambda b, c, a, x : sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule87)


    def cons_f90(m):
        return IntegerQ(S(2)*m)

    cons90 = CustomConstraint(cons_f90)

    def cons_f91(m):
        return Greater(m, S(3)/2)

    cons91 = CustomConstraint(cons_f91)
    pattern88 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons67, cons18, cons90, cons91)
    rule88 = ReplacementRule(pattern88, lambda m, a, x, b, c : -S(2)*a*c**S(2)*(m + S(-1))*Int((c*x)**(m + S(-2))/(a + b*x**S(2))**(S(5)/4), x)/(b*(S(2)*m + S(-3))) + S(2)*c*(c*x)**(m + S(-1))/(b*(a + b*x**S(2))**(S(1)/4)*(S(2)*m + S(-3))))
    rubi.add(rule88)

    pattern89 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons67, cons18, cons90, cons82)
    rule89 = ReplacementRule(pattern89, lambda m, a, x, b, c : -b*(S(2)*m + S(1))*Int((c*x)**(m + S(2))/(a + b*x**S(2))**(S(5)/4), x)/(S(2)*a*c**S(2)*(m + S(1))) + (c*x)**(m + S(1))/(a*c*(a + b*x**S(2))**(S(1)/4)*(m + S(1))))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_), cons6, cons2, cons39)
    rule90 = ReplacementRule(pattern90, lambda b, a, x : -Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/b - S(1)/(b*x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule90)


    def cons_f92(n, m):
        return Greater(m + S(1), n)

    cons92 = CustomConstraint(cons_f92)

    def cons_f93(n, m, p):
        return Not(NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n))

    cons93 = CustomConstraint(cons_f93)
    pattern91 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons14, cons81, cons22, cons92, cons93, cons84)
    rule91 = ReplacementRule(pattern91, lambda n, m, a, x, p, b, c : -c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule91)


    def cons_f94(n, m):
        return Greater(m + S(1), S(2)*n)

    cons94 = CustomConstraint(cons_f94)

    def cons_f95(n, m, p):
        return Not(NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)))

    cons95 = CustomConstraint(cons_f95)
    pattern92 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons56, cons62, cons81, cons22, cons94, cons95, cons86)
    rule92 = ReplacementRule(pattern92, lambda a1, a2, n, m, p, x, b2, c, b1 : -c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*b1*b2*n*(p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))))
    rubi.add(rule92)

    pattern93 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons14, cons81, cons22, cons84)
    rule93 = ReplacementRule(pattern93, lambda n, m, a, x, p, b, c : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule93)

    pattern94 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons56, cons62, cons81, cons22, cons86)
    rule94 = ReplacementRule(pattern94, lambda a1, a2, n, m, p, x, b2, c, b1 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(x_/(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons21)
    rule95 = ReplacementRule(pattern95, lambda b, a, x : Int((x*Rt(b, S(3)) + Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))) - Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))))
    rubi.add(rule95)


    def cons_f96(n):
        return PositiveIntegerQ(n/S(2) + S(-1)/2)

    cons96 = CustomConstraint(cons_f96)

    def cons_f97(m):
        return PositiveIntegerQ(m)

    cons97 = CustomConstraint(cons_f97)

    def cons_f98(n, m):
        return Less(m, n + S(-1))

    cons98 = CustomConstraint(cons_f98)
    pattern96 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons96, cons97, cons98, cons24, )
    def With96(n, m, a, x, b):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) - s**(-m)*(-r)**(m + S(1))*Int(S(1)/(r + s*x), x)/(a*n)
    rule96 = ReplacementRule(pattern96, lambda n, m, a, x, b : With96(n, m, a, x, b))
    rubi.add(rule96)


    def cons_f99(n, m):
        return PositiveIntegerQ(m, n/S(2) + S(-1)/2)

    cons99 = CustomConstraint(cons_f99)
    pattern97 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons99, cons97, cons98, cons25, )
    def With97(n, m, a, x, b):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return -Dist(S(2)*s**(-m)*(-r)**(m + S(1))/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r**(m + S(1))*s**(-m)*Int(S(1)/(r - s*x), x)/(a*n)
    rule97 = ReplacementRule(pattern97, lambda n, m, a, x, b : With97(n, m, a, x, b))
    rubi.add(rule97)


    def cons_f100(n, m):
        return PositiveIntegerQ(m, n/S(4) + S(-1)/2)

    cons100 = CustomConstraint(cons_f100)
    pattern98 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons100, cons97, cons98, cons24, )
    def With98(n, m, a, x, b):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return S(2)*(S(-1))**(m/S(2))*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n) + Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x)
    rule98 = ReplacementRule(pattern98, lambda n, m, a, x, b : With98(n, m, a, x, b))
    rubi.add(rule98)

    pattern99 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons100, cons97, cons98, cons25, )
    def With99(n, m, a, x, b):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(2*Pi*k*m/n) - s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r*cos(2*Pi*k*m/n) + s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    rule99 = ReplacementRule(pattern99, lambda n, m, a, x, b : With99(n, m, a, x, b))
    rubi.add(rule99)

    pattern100 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons31, )
    def With100(b, a, x):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return -Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s)
    rule100 = ReplacementRule(pattern100, lambda b, a, x : With100(b, a, x))
    rubi.add(rule100)

    pattern101 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons32, )
    def With101(b, a, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*b) + s*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*b)
    rule101 = ReplacementRule(pattern101, lambda b, a, x : With101(b, a, x))
    rubi.add(rule101)


    def cons_f101(n, m):
        return PositiveIntegerQ(m, n/S(4))

    cons101 = CustomConstraint(cons_f101)
    pattern102 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons101, cons97, cons98, cons34, )
    def With102(n, m, a, x, b):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r) - sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r)
    rule102 = ReplacementRule(pattern102, lambda n, m, a, x, b : With102(n, m, a, x, b))
    rubi.add(rule102)


    def cons_f102(n, m):
        return Less(m, n/S(2))

    cons102 = CustomConstraint(cons_f102)
    pattern103 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons101, cons97, cons102, cons32, )
    def With103(n, m, a, x, b):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(x**m/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(x**m/(r + s*x**(n/S(2))), x)/(S(2)*a)
    rule103 = ReplacementRule(pattern103, lambda n, m, a, x, b : With103(n, m, a, x, b))
    rubi.add(rule103)


    def cons_f103(n, m):
        return Inequality(n/S(2), LessEqual, m, Less, n)

    cons103 = CustomConstraint(cons_f103)
    pattern104 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons101, cons97, cons103, cons32, )
    def With104(n, m, a, x, b):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(x**(m - n/S(2))/(r - s*x**(n/S(2))), x)/(S(2)*b) + s*Int(x**(m - n/S(2))/(r + s*x**(n/S(2))), x)/(S(2)*b)
    rule104 = ReplacementRule(pattern104, lambda n, m, a, x, b : With104(n, m, a, x, b))
    rubi.add(rule104)


    def cons_f104(n, m):
        return PositiveIntegerQ(m, n)

    cons104 = CustomConstraint(cons_f104)

    def cons_f105(n, m):
        return Greater(m, S(2)*n + S(-1))

    cons105 = CustomConstraint(cons_f105)
    pattern105 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons104, cons105)
    rule105 = ReplacementRule(pattern105, lambda n, m, a, x, b : Int(PolynomialDivide(x**m, a + b*x**n, x), x))
    rubi.add(rule105)

    pattern106 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons37, )
    def With106(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(sqrt(S(3)) + S(2))) + Int((r*x + s*(-sqrt(S(3)) + S(1)))/sqrt(a + b*x**S(3)), x)/r
    rule106 = ReplacementRule(pattern106, lambda b, a, x : With106(b, a, x))
    rubi.add(rule106)

    pattern107 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons38, )
    def With107(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return -sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(-sqrt(S(3)) + S(2))) + Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x)/r
    rule107 = ReplacementRule(pattern107, lambda b, a, x : With107(b, a, x))
    rubi.add(rule107)

    pattern108 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons18, )
    def With108(b, a, x):
        q = Rt(b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule108 = ReplacementRule(pattern108, lambda b, a, x : With108(b, a, x))
    rubi.add(rule108)

    pattern109 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons40, cons41, )
    def With109(b, a, x):
        q = Rt(-b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule109 = ReplacementRule(pattern109, lambda b, a, x : With109(b, a, x))
    rubi.add(rule109)

    pattern110 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons39, )
    def With110(b, a, x):
        q = Rt(-b/a, S(2))
        return Int((q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q - Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule110 = ReplacementRule(pattern110, lambda b, a, x : With110(b, a, x))
    rubi.add(rule110)

    pattern111 = Pattern(Integral(x_**S(4)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_), cons6, cons2, cons21, )
    def With111(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return s**S(2)*(S(-1) + sqrt(S(3)))*Int(S(1)/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2)) - Int((-S(2)*r**S(2)*x**S(4) + s**S(2)*(S(-1) + sqrt(S(3))))/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2))
    rule111 = ReplacementRule(pattern111, lambda b, a, x : With111(b, a, x))
    rubi.add(rule111)

    pattern112 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_), cons6, cons2, cons21)
    rule112 = ReplacementRule(pattern112, lambda b, a, x : -Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))))
    rubi.add(rule112)

    pattern113 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons6, cons2, cons18)
    rule113 = ReplacementRule(pattern113, lambda b, a, x : -a*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x)/S(2) + x**S(3)/(S(2)*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule113)

    pattern114 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_), cons6, cons2, cons39)
    rule114 = ReplacementRule(pattern114, lambda b, a, x : a*Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/(S(2)*b) + (a + b*x**S(4))**(S(3)/4)/(S(2)*b*x))
    rubi.add(rule114)

    pattern115 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons6, cons2, cons18)
    rule115 = ReplacementRule(pattern115, lambda b, a, x : -b*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x) - S(1)/(x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule115)

    pattern116 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_), cons6, cons2, cons39)
    rule116 = ReplacementRule(pattern116, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(1)/4)), x)/(a + b*x**S(4))**(S(1)/4))
    rubi.add(rule116)

    pattern117 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons6, cons2, cons67, cons18)
    rule117 = ReplacementRule(pattern117, lambda b, c, a, x : -a*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/S(2) + x*sqrt(c*x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule117)

    pattern118 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_), cons6, cons2, cons67, cons39)
    rule118 = ReplacementRule(pattern118, lambda b, c, a, x : a*c**S(2)*Int(S(1)/((c*x)**(S(3)/2)*(a + b*x**S(2))**(S(1)/4)), x)/(S(2)*b) + c*(a + b*x**S(2))**(S(3)/4)/(b*sqrt(c*x)))
    rubi.add(rule118)

    pattern119 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons6, cons2, cons67, cons18)
    rule119 = ReplacementRule(pattern119, lambda b, c, a, x : -b*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/c**S(2) - S(2)/(c*sqrt(c*x)*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule119)

    pattern120 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_), cons6, cons2, cons67, cons39)
    rule120 = ReplacementRule(pattern120, lambda b, c, a, x : sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(1)/4)), x)/(c**S(2)*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule120)


    def cons_f106(m):
        return RationalQ(m)

    cons106 = CustomConstraint(cons_f106)

    def cons_f107(n, m):
        return Greater(m, n + S(-1))

    cons107 = CustomConstraint(cons_f107)
    pattern121 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons4, cons14, cons106, cons107, cons87, cons84)
    rule121 = ReplacementRule(pattern121, lambda n, m, a, x, p, b, c : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule121)


    def cons_f108(n, m):
        return SumSimplerQ(m, -n)

    cons108 = CustomConstraint(cons_f108)

    def cons_f109(n, m, p):
        return NegativeIntegerQ((m + n*p + S(1))/n)

    cons109 = CustomConstraint(cons_f109)
    pattern122 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons4, cons14, cons108, cons87, cons109)
    rule122 = ReplacementRule(pattern122, lambda n, m, a, x, p, b, c : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons4, cons56, cons62, cons106, cons105, cons85, cons86)
    rule123 = ReplacementRule(pattern123, lambda a1, a2, n, m, p, x, b2, c, b1 : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule123)


    def cons_f110(n, m):
        return SumSimplerQ(m, -S(2)*n)

    cons110 = CustomConstraint(cons_f110)

    def cons_f111(n, m, p):
        return NegativeIntegerQ((m + S(2)*n*p + S(1))/(S(2)*n))

    cons111 = CustomConstraint(cons_f111)
    pattern124 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons4, cons56, cons62, cons110, cons85, cons111)
    rule124 = ReplacementRule(pattern124, lambda a1, a2, n, m, p, x, b2, c, b1 : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons4, cons14, cons106, cons82, cons84)
    rule125 = ReplacementRule(pattern125, lambda n, m, a, x, p, b, c : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule125)


    def cons_f112(n, m):
        return SumSimplerQ(m, n)

    cons112 = CustomConstraint(cons_f112)
    pattern126 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons4, cons14, cons112, cons109)
    rule126 = ReplacementRule(pattern126, lambda n, m, a, x, p, b, c : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons4, cons56, cons62, cons106, cons82, cons86)
    rule127 = ReplacementRule(pattern127, lambda a1, a2, n, m, p, x, b2, c, b1 : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule127)


    def cons_f113(n, m):
        return SumSimplerQ(m, S(2)*n)

    cons113 = CustomConstraint(cons_f113)
    pattern128 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons4, cons56, cons62, cons113, cons111)
    rule128 = ReplacementRule(pattern128, lambda a1, a2, n, m, p, x, b2, c, b1 : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule128)


    def cons_f114(m):
        return FractionQ(m)

    cons114 = CustomConstraint(cons_f114)
    pattern129 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons4, cons14, cons114, cons84, )
    def With129(n, m, a, x, p, b, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    rule129 = ReplacementRule(pattern129, lambda n, m, a, x, p, b, c : With129(n, m, a, x, p, b, c))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons4, cons56, cons62, cons114, cons86, )
    def With130(a1, a2, n, m, p, x, b2, c, b1):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(k*n))**p*(a2 + b2*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    rule130 = ReplacementRule(pattern130, lambda a1, a2, n, m, p, x, b2, c, b1 : With130(a1, a2, n, m, p, x, b2, c, b1))
    rubi.add(rule130)


    def cons_f115(n, m, p):
        return IntegersQ(m, p + (m + S(1))/n)

    cons115 = CustomConstraint(cons_f115)
    pattern131 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons42, cons43, cons115)
    rule131 = ReplacementRule(pattern131, lambda n, m, a, x, p, b : a**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule131)


    def cons_f116(n, m, p):
        return IntegersQ(m, p + (m + S(1))/(S(2)*n))

    cons116 = CustomConstraint(cons_f116)
    pattern132 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons56, cons62, cons15, cons42, cons43, cons116)
    rule132 = ReplacementRule(pattern132, lambda a1, a2, n, m, p, x, b2, b1 : (a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule132)


    def cons_f117(n, m, p):
        return Less(Denominator(p + (m + S(1))/n), Denominator(p))

    cons117 = CustomConstraint(cons_f117)
    pattern133 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons14, cons15, cons42, cons43, cons80, cons117)
    rule133 = ReplacementRule(pattern133, lambda n, m, a, x, p, b : (a/(a + b*x**n))**(p + (m + S(1))/n)*(a + b*x**n)**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule133)


    def cons_f118(n, m, p):
        return Less(Denominator(p + (m + S(1))/(S(2)*n)), Denominator(p))

    cons118 = CustomConstraint(cons_f118)
    pattern134 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons56, cons62, cons15, cons42, cons43, cons80, cons118)
    rule134 = ReplacementRule(pattern134, lambda a1, a2, n, m, p, x, b2, b1 : (a1/(a1 + b1*x**n))**(p + (m + S(1))/(S(2)*n))*(a2/(a2 + b2*x**n))**(p + (m + S(1))/(S(2)*n))*(a1 + b1*x**n)**(p + (m + S(1))/(S(2)*n))*(a2 + b2*x**n)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule134)

    pattern135 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons4, cons46, cons80)
    rule135 = ReplacementRule(pattern135, lambda n, m, a, x, p, b : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule135)

    pattern136 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons4, cons56, cons64, cons80)
    rule136 = ReplacementRule(pattern136, lambda a1, a2, n, m, p, x, b2, b1 : -Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule136)

    pattern137 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons4, cons46, cons114, )
    def With137(n, m, a, x, p, b, c):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    rule137 = ReplacementRule(pattern137, lambda n, m, a, x, p, b, c : With137(n, m, a, x, p, b, c))
    rubi.add(rule137)

    pattern138 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons4, cons56, cons64, cons114, )
    def With138(a1, a2, n, m, p, x, b2, c, b1):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(-k*n))**p*(a2 + b2*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    rule138 = ReplacementRule(pattern138, lambda a1, a2, n, m, p, x, b2, c, b1 : With138(a1, a2, n, m, p, x, b2, c, b1))
    rubi.add(rule138)


    def cons_f119(m):
        return Not(RationalQ(m))

    cons119 = CustomConstraint(cons_f119)
    pattern139 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons4, cons46, cons119)
    rule139 = ReplacementRule(pattern139, lambda n, m, a, x, p, b, c : -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons4, cons56, cons64, cons119)
    rule140 = ReplacementRule(pattern140, lambda a1, a2, n, m, p, x, b2, c, b1 : -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons4, cons47, )
    def With141(n, m, a, x, p, b):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    rule141 = ReplacementRule(pattern141, lambda n, m, a, x, p, b : With141(n, m, a, x, p, b))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons4, cons56, cons65, )
    def With142(a1, a2, n, m, p, x, b2, b1):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    rule142 = ReplacementRule(pattern142, lambda a1, a2, n, m, p, x, b2, b1 : With142(a1, a2, n, m, p, x, b2, b1))
    rubi.add(rule142)

    pattern143 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons4, cons47)
    rule143 = ReplacementRule(pattern143, lambda n, m, a, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule143)

    pattern144 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons4, cons56, cons65)
    rule144 = ReplacementRule(pattern144, lambda a1, a2, n, m, p, x, b2, c, b1 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule144)


    def cons_f120(n, m):
        return IntegerQ(n/(m + S(1)))

    cons120 = CustomConstraint(cons_f120)

    def cons_f121(n):
        return Not(IntegerQ(n))

    cons121 = CustomConstraint(cons_f121)
    pattern145 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons3, cons4, cons120, cons121)
    rule145 = ReplacementRule(pattern145, lambda n, m, a, x, p, b : Subst(Int((a + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule145)


    def cons_f122(n, m):
        return IntegerQ(S(2)*n/(m + S(1)))

    cons122 = CustomConstraint(cons_f122)

    def cons_f123(n):
        return Not(IntegerQ(S(2)*n))

    cons123 = CustomConstraint(cons_f123)
    pattern146 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons3, cons4, cons56, cons122, cons123)
    rule146 = ReplacementRule(pattern146, lambda a1, a2, n, m, p, x, b2, b1 : Subst(Int((a1 + b1*x**(n/(m + S(1))))**p*(a2 + b2*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons120, cons121)
    rule147 = ReplacementRule(pattern147, lambda n, m, a, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons122, cons123)
    rule148 = ReplacementRule(pattern148, lambda a1, a2, n, m, p, x, b2, c, b1 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule148)


    def cons_f124(n, m, p):
        return ZeroQ(p + (m + S(1))/n)

    cons124 = CustomConstraint(cons_f124)
    pattern149 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons3, cons124, cons15, cons16)
    rule149 = ReplacementRule(pattern149, lambda n, m, a, x, p, b : -b*n*p*Int(x**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*x**n)**p/(m + S(1)))
    rubi.add(rule149)


    def cons_f125(n, m, p):
        return ZeroQ(p + (m + S(1))/(S(2)*n))

    cons125 = CustomConstraint(cons_f125)
    pattern150 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons3, cons56, cons125, cons15, cons16)
    rule150 = ReplacementRule(pattern150, lambda a1, a2, n, m, p, x, b2, b1 : -S(2)*b1*b2*n*p*Int(x**(m + n)*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(m + S(1)))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons124, cons15, cons16)
    rule151 = ReplacementRule(pattern151, lambda n, m, a, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons56, cons125, cons15, cons16)
    rule152 = ReplacementRule(pattern152, lambda a1, a2, n, m, p, x, b2, c, b1 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule152)


    def cons_f126(n, m, p):
        return IntegerQ(p + (m + S(1))/n)

    cons126 = CustomConstraint(cons_f126)
    pattern153 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons126, cons15, cons16, cons87)
    rule153 = ReplacementRule(pattern153, lambda n, m, a, x, p, b, c : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule153)


    def cons_f127(n, m, p):
        return IntegerQ(p + (m + S(1))/(S(2)*n))

    cons127 = CustomConstraint(cons_f127)
    pattern154 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons56, cons127, cons15, cons16, cons85)
    rule154 = ReplacementRule(pattern154, lambda a1, a2, n, m, p, x, b2, c, b1 : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule154)

    pattern155 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons68, cons3, cons126, cons15, cons42, )
    def With155(n, m, a, x, p, b):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(-b*x**k + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    rule155 = ReplacementRule(pattern155, lambda n, m, a, x, p, b : With155(n, m, a, x, p, b))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons68, cons3, cons56, cons127, cons15, cons42, )
    def With156(a1, a2, n, m, p, x, b2, b1):
        k = Denominator(p)
        return k*(a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**(k*(m + S(1))/(S(2)*n) + S(-1))*(-b1*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x**(S(2)*n/k)*(a1 + b1*x**n)**(-S(1)/k)*(a2 + b2*x**n)**(-S(1)/k))/(S(2)*n)
    rule156 = ReplacementRule(pattern156, lambda a1, a2, n, m, p, x, b2, b1 : With156(a1, a2, n, m, p, x, b2, b1))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons126, cons15, cons42)
    rule157 = ReplacementRule(pattern157, lambda n, m, a, x, p, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons56, cons127, cons15, cons42)
    rule158 = ReplacementRule(pattern158, lambda a1, a2, n, m, p, x, b2, c, b1 : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons126, cons15, cons22)
    rule159 = ReplacementRule(pattern159, lambda n, m, a, x, p, b, c : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons56, cons126, cons15, cons22)
    rule160 = ReplacementRule(pattern160, lambda a1, a2, n, m, p, x, b2, c, b1 : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule160)


    def cons_f128(n, m):
        return FractionQ((m + S(1))/n)

    cons128 = CustomConstraint(cons_f128)
    pattern161 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons68, cons3, cons128, cons108, )
    def With161(n, m, a, x, b):
        mn = m - n
        return -a*Int(x**mn/(a + b*x**n), x)/b + x**(mn + S(1))/(b*(mn + S(1)))
    rule161 = ReplacementRule(pattern161, lambda n, m, a, x, b : With161(n, m, a, x, b))
    rubi.add(rule161)

    pattern162 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons68, cons3, cons128, cons112)
    rule162 = ReplacementRule(pattern162, lambda n, m, a, x, b : -b*Int(x**(m + n)/(a + b*x**n), x)/a + x**(m + S(1))/(a*(m + S(1))))
    rubi.add(rule162)


    def cons_f129(n, m):
        return Or(SumSimplerQ(m, n), SumSimplerQ(m, -n))

    cons129 = CustomConstraint(cons_f129)
    pattern163 = Pattern(Integral((c_*x_)**m_/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons67, cons68, cons3, cons128, cons129)
    rule163 = ReplacementRule(pattern163, lambda n, m, a, x, b, c : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m/(a + b*x**n), x))
    rubi.add(rule163)


    def cons_f130(a, p):
        return Or(NegativeIntegerQ(p), PositiveQ(a))

    cons130 = CustomConstraint(cons_f130)
    pattern164 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons49, cons130)
    rule164 = ReplacementRule(pattern164, lambda n, m, a, x, p, b, c : a**p*(c*x)**(m + S(1))*Hypergeometric2F1(-p, (m + S(1))/n, S(1) + (m + S(1))/n, -b*x**n/a)/(c*(m + S(1))))
    rubi.add(rule164)


    def cons_f131(a, p):
        return Not(Or(NegativeIntegerQ(p), PositiveQ(a)))

    cons131 = CustomConstraint(cons_f131)
    pattern165 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_), cons6, cons2, cons67, cons68, cons3, cons4, cons49, cons131)
    rule165 = ReplacementRule(pattern165, lambda n, m, a, x, p, b, c : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((c*x)**m*(S(1) + b*x**n/a)**p, x))
    rubi.add(rule165)


    def cons_f132(x, v):
        return LinearQ(v, x)

    cons132 = CustomConstraint(cons_f132)

    def cons_f133(x, v):
        return NonzeroQ(v - x)

    cons133 = CustomConstraint(cons_f133)
    pattern166 = Pattern(Integral(x_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons6, cons2, cons3, cons4, cons132, cons80, cons133)
    rule166 = ReplacementRule(pattern166, lambda n, m, a, v, p, x, b : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule166)


    def cons_f134(u, x, v):
        return LinearPairQ(u, v, x)

    cons134 = CustomConstraint(cons_f134)
    pattern167 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_), cons6, cons2, cons68, cons3, cons4, cons134)
    rule167 = ReplacementRule(pattern167, lambda u, n, m, a, v, p, x, b : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_), cons58, cons59, cons60, cons61, cons67, cons68, cons3, cons4, cons56, cons66)
    rule168 = ReplacementRule(pattern168, lambda a1, a2, n, m, p, x, b2, c, b1 : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule168)


    def cons_f135(d, b, c, a):
        return NonzeroQ(-a*d + b*c)

    cons135 = CustomConstraint(cons_f135)

    def cons_f136(q, p):
        return PositiveIntegerQ(p, q)

    cons136 = CustomConstraint(cons_f136)

    def cons_f137(d, x):
        return FreeQ(d, x)

    cons137 = CustomConstraint(cons_f137)
    pattern169 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons135, cons136)
    rule169 = ReplacementRule(pattern169, lambda q, n, a, x, p, d, b, c : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule169)


    def cons_f138(q, p):
        return IntegersQ(p, q)

    cons138 = CustomConstraint(cons_f138)
    pattern170 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons135, cons138, cons73)
    rule170 = ReplacementRule(pattern170, lambda q, n, a, x, p, d, b, c : Int(x**(n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule170)


    def cons_f139(q, x):
        return FreeQ(q, x)

    cons139 = CustomConstraint(cons_f139)
    pattern171 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons4, cons139, cons135, cons46)
    rule171 = ReplacementRule(pattern171, lambda q, n, a, x, p, d, b, c : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q/x**S(2), x), x, S(1)/x))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons4, cons139, cons135, cons47, )
    def With172(q, n, a, x, p, d, b, c):
        g = Denominator(n)
        return g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    rule172 = ReplacementRule(pattern172, lambda q, n, a, x, p, d, b, c : With172(q, n, a, x, p, d, b, c))
    rubi.add(rule172)


    def cons_f140(n, p):
        return ZeroQ(n*p + S(1))

    cons140 = CustomConstraint(cons_f140)

    def cons_f141(n):
        return IntegerQ(n)

    cons141 = CustomConstraint(cons_f141)
    pattern173 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons140, cons141)
    rule173 = ReplacementRule(pattern173, lambda n, a, x, p, d, b, c : Subst(Int(S(1)/(c - x**n*(-a*d + b*c)), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule173)


    def cons_f142(q, n, p):
        return ZeroQ(n*(p + q + S(1)) + S(1))

    cons142 = CustomConstraint(cons_f142)

    def cons_f143(q):
        return RationalQ(q)

    cons143 = CustomConstraint(cons_f143)

    def cons_f144(q):
        return Greater(q, S(0))

    cons144 = CustomConstraint(cons_f144)
    pattern174 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons135, cons142, cons143, cons144, cons8)
    rule174 = ReplacementRule(pattern174, lambda q, n, a, x, p, d, b, c : -c*q*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1)), x)/(a*(p + S(1))) - x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))))
    rubi.add(rule174)


    def cons_f145(p):
        return NegativeIntegerQ(p)

    cons145 = CustomConstraint(cons_f145)
    pattern175 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons139, cons135, cons142, cons145)
    rule175 = ReplacementRule(pattern175, lambda q, n, a, x, p, d, b, c : a**p*c**(-p + S(-1))*x*(c + d*x**n)**(-S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n))))
    rubi.add(rule175)

    pattern176 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons135, cons142)
    rule176 = ReplacementRule(pattern176, lambda q, n, a, x, p, d, b, c : x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p)*(a + b*x**n)**p*(c + d*x**n)**(-p - S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/c)
    rubi.add(rule176)


    def cons_f146(q, n, p):
        return ZeroQ(n*(p + q + S(2)) + S(1))

    cons146 = CustomConstraint(cons_f146)

    def cons_f147(q, a, p, d, b, c):
        return ZeroQ(a*d*(p + S(1)) + b*c*(q + S(1)))

    cons147 = CustomConstraint(cons_f147)
    pattern177 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons135, cons146, cons147)
    rule177 = ReplacementRule(pattern177, lambda q, n, a, x, p, d, b, c : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c))
    rubi.add(rule177)


    def cons_f148(q, p):
        return Or(And(RationalQ(p), Less(p, S(-1))), Not(And(RationalQ(q), Less(q, S(-1)))))

    cons148 = CustomConstraint(cons_f148)
    pattern178 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons139, cons135, cons146, cons148, cons8)
    rule178 = ReplacementRule(pattern178, lambda q, n, a, x, p, d, b, c : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + (b*c + n*(p + S(1))*(-a*d + b*c))*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q, x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule178)


    def cons_f149(n, a, p, d, b, c):
        return ZeroQ(a*d - b*c*(n*(p + S(1)) + S(1)))

    cons149 = CustomConstraint(cons_f149)
    pattern179 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons135, cons149)
    rule179 = ReplacementRule(pattern179, lambda n, a, x, p, d, b, c : c*x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule179)


    def cons_f150(n, p):
        return Or(And(RationalQ(p), Less(p, S(-1))), NegativeIntegerQ(p + S(1)/n))

    cons150 = CustomConstraint(cons_f150)
    pattern180 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons135, cons150)
    rule180 = ReplacementRule(pattern180, lambda n, a, x, p, d, b, c : -x*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*n*(p + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((c_ + x_**n_*WC('d', S(1)))/(a_ + x_**n_*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons3, cons135, cons10, cons11)
    rule181 = ReplacementRule(pattern181, lambda n, a, x, d, b, c : c*x/a - (-a*d + b*c)*Int(S(1)/(a*x**(-n) + b), x)/a)
    rubi.add(rule181)


    def cons_f151(n, p):
        return NonzeroQ(n*(p + S(1)) + S(1))

    cons151 = CustomConstraint(cons_f151)
    pattern182 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons3, cons135, cons151)
    rule182 = ReplacementRule(pattern182, lambda n, a, x, p, d, b, c : d*x*(a + b*x**n)**(p + S(1))/(b*(n*(p + S(1)) + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**p, x)/(b*(n*(p + S(1)) + S(1))))
    rubi.add(rule182)


    def cons_f152(q):
        return NegativeIntegerQ(q)

    cons152 = CustomConstraint(cons_f152)

    def cons_f153(q, p):
        return GreaterEqual(p, -q)

    cons153 = CustomConstraint(cons_f153)
    pattern183 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons135, cons13, cons152, cons153)
    rule183 = ReplacementRule(pattern183, lambda q, n, a, x, p, d, b, c : Int(PolynomialDivide((a + b*x**n)**p, (c + d*x**n)**(-q), x), x))
    rubi.add(rule183)

    pattern184 = Pattern(Integral(S(1)/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons3, cons135)
    rule184 = ReplacementRule(pattern184, lambda n, a, x, d, b, c : b*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule184)


    def cons_f154(d, b, c, a):
        return ZeroQ(S(3)*a*d + b*c)

    cons154 = CustomConstraint(cons_f154)
    pattern185 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons154, cons18)
    rule185 = ReplacementRule(pattern185, lambda a, x, d, b, c : sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(-x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c) + sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c))
    rubi.add(rule185)

    pattern186 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons154, cons39)
    rule186 = ReplacementRule(pattern186, lambda a, x, d, b, c : Int((-x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6) + Int((x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6))
    rubi.add(rule186)

    pattern187 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(2)/3)/(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons154)
    rule187 = ReplacementRule(pattern187, lambda a, x, d, b, c : b*Int((a + b*x**S(2))**(S(-1)/3), x)/d - (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/d)
    rubi.add(rule187)

    pattern188 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135)
    rule188 = ReplacementRule(pattern188, lambda a, x, d, b, c : sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(1)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135)
    rule189 = ReplacementRule(pattern189, lambda a, x, d, b, c : sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(3)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule189)


    def cons_f155(p):
        return Or(Equal(p, S(1)/2), Equal(Denominator(p), S(4)))

    cons155 = CustomConstraint(cons_f155)
    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**WC('p', S(1))/(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons15, cons16, cons155)
    rule190 = ReplacementRule(pattern190, lambda a, x, p, d, b, c : b*Int((a + b*x**S(2))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(2))**(p + S(-1))/(c + d*x**S(2)), x)/d)
    rubi.add(rule190)


    def cons_f156(p):
        return Equal(Denominator(p), S(4))

    cons156 = CustomConstraint(cons_f156)

    def cons_f157(p):
        return Or(Equal(p, S(-5)/4), Equal(p, S(-7)/4))

    cons157 = CustomConstraint(cons_f157)
    pattern191 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**p_/(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons15, cons22, cons156, cons157)
    rule191 = ReplacementRule(pattern191, lambda a, x, p, d, b, c : b*Int((a + b*x**S(2))**p, x)/(-a*d + b*c) - d*Int((a + b*x**S(2))**(p + S(1))/(c + d*x**S(2)), x)/(-a*d + b*c))
    rubi.add(rule191)


    def cons_f158(d, b, c, a):
        return ZeroQ(a*d + b*c)

    cons158 = CustomConstraint(cons_f158)

    def cons_f159(a, b):
        return PosQ(a*b)

    cons159 = CustomConstraint(cons_f159)
    pattern192 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons158, cons159)
    rule192 = ReplacementRule(pattern192, lambda a, x, d, b, c : a*Subst(Int(S(1)/(-S(4)*a*b*x**S(4) + S(1)), x), x, x/sqrt(a + b*x**S(4)))/c)
    rubi.add(rule192)


    def cons_f160(a, b):
        return NegQ(a*b)

    cons160 = CustomConstraint(cons_f160)
    pattern193 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons158, cons160, )
    def With193(a, x, d, b, c):
        q = Rt(-a*b, S(4))
        return a*ArcTan(q*x*(a + q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q) + a*atanh(q*x*(a - q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q)
    rule193 = ReplacementRule(pattern193, lambda a, x, d, b, c : With193(a, x, d, b, c))
    rubi.add(rule193)

    pattern194 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135)
    rule194 = ReplacementRule(pattern194, lambda a, x, d, b, c : b*Int(S(1)/sqrt(a + b*x**S(4)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(4))*(c + d*x**S(4))), x)/d)
    rubi.add(rule194)

    pattern195 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)/(c_ + x_**S(4)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135)
    rule195 = ReplacementRule(pattern195, lambda a, x, d, b, c : sqrt(a/(a + b*x**S(4)))*sqrt(a + b*x**S(4))*Subst(Int(S(1)/((c - x**S(4)*(-a*d + b*c))*sqrt(-b*x**S(4) + S(1))), x), x, x/(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule195)


    def cons_f161(p):
        return Or(Equal(p, S(3)/4), Equal(p, S(5)/4))

    cons161 = CustomConstraint(cons_f161)
    pattern196 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**p_/(c_ + x_**S(4)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons15, cons161)
    rule196 = ReplacementRule(pattern196, lambda a, x, p, d, b, c : b*Int((a + b*x**S(4))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(4))**(p + S(-1))/(c + d*x**S(4)), x)/d)
    rubi.add(rule196)

    pattern197 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('b', S(1)))*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135)
    rule197 = ReplacementRule(pattern197, lambda a, x, d, b, c : Int(S(1)/(sqrt(a + b*x**S(4))*(-x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c) + Int(S(1)/(sqrt(a + b*x**S(4))*(x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c))
    rubi.add(rule197)

    pattern198 = Pattern(Integral(S(1)/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135)
    rule198 = ReplacementRule(pattern198, lambda a, x, d, b, c : b*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - d*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule198)


    def cons_f162(d, c):
        return PosQ(d/c)

    cons162 = CustomConstraint(cons_f162)
    pattern199 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_), cons6, cons2, cons67, cons137, cons18, cons162)
    rule199 = ReplacementRule(pattern199, lambda a, x, d, b, c : sqrt(a + b*x**S(2))*EllipticE(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(c*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule199)


    def cons_f163(q, p):
        return RationalQ(p, q)

    cons163 = CustomConstraint(cons_f163)

    def cons_f164(q):
        return Less(S(0), q, S(1))

    cons164 = CustomConstraint(cons_f164)

    def cons_f165(q, n, a, x, p, d, b, c):
        return IntBinomialQ(a, b, c, d, n, p, q, x)

    cons165 = CustomConstraint(cons_f165)
    pattern200 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons135, cons163, cons22, cons164, cons165)
    rule200 = ReplacementRule(pattern200, lambda q, n, a, x, p, d, b, c : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(n*(p + S(1)) + S(1)) + d*x**n*(n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))))
    rubi.add(rule200)


    def cons_f166(q):
        return Greater(q, S(1))

    cons166 = CustomConstraint(cons_f166)
    pattern201 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons135, cons163, cons22, cons166, cons165)
    rule201 = ReplacementRule(pattern201, lambda q, n, a, x, p, d, b, c : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*n*(p + S(1))) - Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(a*d - b*c*(n*(p + S(1)) + S(1))) + d*x**n*(a*d*(n*(q + S(-1)) + S(1)) - b*c*(n*(p + q) + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule201)


    def cons_f167(q, p):
        return Not(And(Not(IntegerQ(p)), IntegerQ(q), Less(q, S(-1))))

    cons167 = CustomConstraint(cons_f167)
    pattern202 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons139, cons135, cons15, cons22, cons167, cons165)
    rule202 = ReplacementRule(pattern202, lambda q, n, a, x, p, d, b, c : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c + b*d*x**n*(n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule202)


    def cons_f168(q, p):
        return Greater(p + q, S(0))

    cons168 = CustomConstraint(cons_f168)
    pattern203 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons135, cons14, cons138, cons168)
    rule203 = ReplacementRule(pattern203, lambda q, n, a, x, p, d, b, c : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule203)


    def cons_f169(q, n, p):
        return NonzeroQ(n*(p + q) + S(1))

    cons169 = CustomConstraint(cons_f169)

    def cons_f170(p):
        return Not(And(IntegerQ(p), Greater(p, S(1))))

    cons170 = CustomConstraint(cons_f170)
    pattern204 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons135, cons143, cons166, cons169, cons170, cons165)
    rule204 = ReplacementRule(pattern204, lambda q, n, a, x, p, d, b, c : d*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*(n*(p + q) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(-a*d + b*c*(n*(p + q) + S(1))) + d*x**n*(-a*d*(n*(q + S(-1)) + S(1)) + b*c*(n*(p + S(2)*q + S(-1)) + S(1))), x), x)/(b*(n*(p + q) + S(1))))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons135, cons163, cons144, cons16, cons165)
    rule205 = ReplacementRule(pattern205, lambda q, n, a, x, p, d, b, c : n*Int((a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(n*(p + q) + S(1)) + x*(a + b*x**n)**p*(c + d*x**n)**q/(n*(p + q) + S(1)))
    rubi.add(rule205)


    def cons_f171(d, b, c, a):
        return Not(SimplerSqrtQ(b/a, d/c))

    cons171 = CustomConstraint(cons_f171)
    pattern206 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons162, cons18, cons171)
    rule206 = ReplacementRule(pattern206, lambda a, x, d, b, c : sqrt(a + b*x**S(2))*EllipticF(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(a*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule206)


    def cons_f172(d, c):
        return NegQ(d/c)

    cons172 = CustomConstraint(cons_f172)

    def cons_f173(c):
        return PositiveQ(c)

    cons173 = CustomConstraint(cons_f173)

    def cons_f174(d, b, c, a):
        return Not(And(NegQ(b/a), SimplerSqrtQ(-b/a, -d/c)))

    cons174 = CustomConstraint(cons_f174)
    pattern207 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons172, cons173, cons19, cons174)
    rule207 = ReplacementRule(pattern207, lambda a, x, d, b, c : EllipticF(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(a)*sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule207)


    def cons_f175(d, a, c, b):
        return PositiveQ(a - b*c/d)

    cons175 = CustomConstraint(cons_f175)
    pattern208 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons172, cons173, cons175)
    rule208 = ReplacementRule(pattern208, lambda a, x, d, b, c : -EllipticF(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*sqrt(a - b*c/d)*Rt(-d/c, S(2))))
    rubi.add(rule208)


    def cons_f176(c):
        return Not(PositiveQ(c))

    cons176 = CustomConstraint(cons_f176)
    pattern209 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons176)
    rule209 = ReplacementRule(pattern209, lambda a, x, d, b, c : sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*sqrt(a + b*x**S(2))), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons162, cons18)
    rule210 = ReplacementRule(pattern210, lambda a, x, d, b, c : a*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x) + b*Int(x**S(2)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x))
    rubi.add(rule210)

    pattern211 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons162, cons39)
    rule211 = ReplacementRule(pattern211, lambda a, x, d, b, c : b*Int(sqrt(c + d*x**S(2))/sqrt(a + b*x**S(2)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/d)
    rubi.add(rule211)

    pattern212 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons172, cons173, cons19)
    rule212 = ReplacementRule(pattern212, lambda a, x, d, b, c : sqrt(a)*EllipticE(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule212)

    pattern213 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons172, cons173, cons175)
    rule213 = ReplacementRule(pattern213, lambda a, x, d, b, c : -sqrt(a - b*c/d)*EllipticE(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule213)

    pattern214 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons172, cons173, cons20)
    rule214 = ReplacementRule(pattern214, lambda a, x, d, b, c : sqrt(a + b*x**S(2))*Int(sqrt(S(1) + b*x**S(2)/a)/sqrt(c + d*x**S(2)), x)/sqrt(S(1) + b*x**S(2)/a))
    rubi.add(rule214)

    pattern215 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons172, cons176)
    rule215 = ReplacementRule(pattern215, lambda a, x, d, b, c : sqrt(S(1) + d*x**S(2)/c)*Int(sqrt(a + b*x**S(2))/sqrt(S(1) + d*x**S(2)/c), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons139, cons135, cons48)
    rule216 = ReplacementRule(pattern216, lambda q, n, a, x, p, d, b, c : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule216)


    def cons_f177(n):
        return NonzeroQ(n + S(1))

    cons177 = CustomConstraint(cons_f177)
    pattern217 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons135, cons177, cons19, cons173)
    rule217 = ReplacementRule(pattern217, lambda q, n, a, x, p, d, b, c : a**p*c**q*x*AppellF1(S(1)/n, -p, -q, S(1) + S(1)/n, -b*x**n/a, -d*x**n/c))
    rubi.add(rule217)

    pattern218 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons135, cons177, cons20)
    rule218 = ReplacementRule(pattern218, lambda q, n, a, x, p, d, b, c : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule218)


    def cons_f178(mn, n):
        return EqQ(mn, -n)

    cons178 = CustomConstraint(cons_f178)

    def cons_f179(q):
        return IntegerQ(q)

    cons179 = CustomConstraint(cons_f179)

    def cons_f180(n, p):
        return Or(PosQ(n), Not(IntegerQ(p)))

    cons180 = CustomConstraint(cons_f180)
    pattern219 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons178, cons179, cons180)
    rule219 = ReplacementRule(pattern219, lambda q, mn, n, a, x, p, d, b, c : Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule219)


    def cons_f181(q):
        return Not(IntegerQ(q))

    cons181 = CustomConstraint(cons_f181)
    pattern220 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons178, cons181, cons66)
    rule220 = ReplacementRule(pattern220, lambda q, mn, n, a, x, p, d, b, c : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule220)

    pattern221 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons54, cons55)
    rule221 = ReplacementRule(pattern221, lambda q, u, n, a, x, p, d, b, c : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule221)


    def cons_f182(u, x, v):
        return PseudoBinomialPairQ(u, v, x)

    cons182 = CustomConstraint(cons_f182)
    pattern222 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_), cons4, cons139, cons182)
    rule222 = ReplacementRule(pattern222, lambda q, u, x, v, p : Int(NormalizePseudoBinomial(u, x)**p*NormalizePseudoBinomial(v, x)**q, x))
    rubi.add(rule222)


    def cons_f183(m, p):
        return IntegersQ(p, m/p)

    cons183 = CustomConstraint(cons_f183)

    def cons_f184(u, m, v, x, p):
        return PseudoBinomialPairQ(u*x**(m/p), v, x)

    cons184 = CustomConstraint(cons_f184)
    pattern223 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_), cons4, cons139, cons183, cons184)
    rule223 = ReplacementRule(pattern223, lambda q, u, m, v, p, x : Int(NormalizePseudoBinomial(v, x)**q*NormalizePseudoBinomial(u*x**(m/p), x)**p, x))
    rubi.add(rule223)


    def cons_f185(e, m):
        return Or(IntegerQ(m), PositiveQ(e))

    cons185 = CustomConstraint(cons_f185)

    def cons_f186(e, x):
        return FreeQ(e, x)

    cons186 = CustomConstraint(cons_f186)
    pattern224 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons185, cons70)
    rule224 = ReplacementRule(pattern224, lambda q, e, n, m, x, p, d, b, c : b**(S(1) - (m + S(1))/n)*e**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule224)

    pattern225 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons185, cons71)
    rule225 = ReplacementRule(pattern225, lambda q, e, n, m, x, p, d, b, c : b**IntPart(p)*e**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q, x))
    rubi.add(rule225)

    pattern226 = Pattern(Integral((e_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons72)
    rule226 = ReplacementRule(pattern226, lambda q, e, n, m, x, p, d, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule226)


    def cons_f187(n, m):
        return ZeroQ(m - n + S(1))

    cons187 = CustomConstraint(cons_f187)
    pattern227 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons139, cons135, cons187)
    rule227 = ReplacementRule(pattern227, lambda q, n, m, a, x, p, d, b, c : Subst(Int((a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule227)

    pattern228 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons135, cons138, cons73)
    rule228 = ReplacementRule(pattern228, lambda q, n, m, a, x, p, d, b, c : Int(x**(m + n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule228)

    pattern229 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons139, cons135, cons70)
    rule229 = ReplacementRule(pattern229, lambda q, n, m, a, x, p, d, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule229)

    pattern230 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons135, cons70)
    rule230 = ReplacementRule(pattern230, lambda q, e, n, m, a, x, p, d, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons136)
    rule231 = ReplacementRule(pattern231, lambda q, e, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule231)


    def cons_f188(n, m, a, p, d, b, c):
        return ZeroQ(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))

    cons188 = CustomConstraint(cons_f188)
    pattern232 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons135, cons188, cons75)
    rule232 = ReplacementRule(pattern232, lambda e, n, m, a, x, p, d, b, c : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule232)


    def cons_f189(non2, n):
        return ZeroQ(-n/S(2) + non2)

    cons189 = CustomConstraint(cons_f189)

    def cons_f190(a1, a2, n, m, p, b2, d, c, b1):
        return ZeroQ(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))

    cons190 = CustomConstraint(cons_f190)
    pattern233 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons68, cons3, cons4, cons189, cons56, cons190, cons75)
    rule233 = ReplacementRule(pattern233, lambda non2, e, b2, a1, a2, n, m, x, p, d, c, b1 : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))))
    rubi.add(rule233)


    def cons_f191(n, m, p):
        return ZeroQ(m + n*(p + S(1)) + S(1))

    cons191 = CustomConstraint(cons_f191)

    def cons_f192(e, n):
        return Or(IntegerQ(n), PositiveQ(e))

    cons192 = CustomConstraint(cons_f192)

    def cons_f193(n, m):
        return RationalQ(m, n)

    cons193 = CustomConstraint(cons_f193)

    def cons_f194(n, m):
        return Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1))))

    cons194 = CustomConstraint(cons_f194)
    pattern234 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons191, cons192, cons193, cons194)
    rule234 = ReplacementRule(pattern234, lambda e, n, m, a, x, p, d, b, c : d*e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p, x) + c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons135, cons191, cons75)
    rule235 = ReplacementRule(pattern235, lambda e, n, m, a, x, p, d, b, c : d*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/b + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*(m + S(1))))
    rubi.add(rule235)


    def cons_f195(p):
        return Not(And(IntegerQ(p), Less(p, S(-1))))

    cons195 = CustomConstraint(cons_f195)
    pattern236 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons192, cons193, cons194, cons195)
    rule236 = ReplacementRule(pattern236, lambda e, n, m, a, x, p, d, b, c : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) + e**(-n)*(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons4, cons189, cons56, cons192, cons193, cons194, cons195)
    rule237 = ReplacementRule(pattern237, lambda non2, e, b2, a1, a2, n, m, x, p, d, c, b1 : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))) + e**(-n)*(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(a1*a2*(m + S(1))))
    rubi.add(rule237)


    def cons_f196(m):
        return PositiveIntegerQ(m/S(2))

    cons196 = CustomConstraint(cons_f196)

    def cons_f197(m, p):
        return Or(IntegerQ(p), Equal(m + S(2)*p + S(1), S(0)))

    cons197 = CustomConstraint(cons_f197)
    pattern238 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons15, cons22, cons196, cons197)
    rule238 = ReplacementRule(pattern238, lambda m, a, x, p, d, b, c : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int((a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*x**S(2)*(p + S(1))*(b**(m/S(2))*x**(m + S(-2))*(c + d*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule238)


    def cons_f198(m):
        return NegativeIntegerQ(m/S(2))

    cons198 = CustomConstraint(cons_f198)
    pattern239 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons15, cons22, cons198, cons197)
    rule239 = ReplacementRule(pattern239, lambda m, a, x, p, d, b, c : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int(x**m*(a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*(p + S(1))*(b**(m/S(2))*(c + d*x**S(2)) - x**(-m + S(2))*(-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - x**(-m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule239)


    def cons_f199(n, m, p):
        return Or(IntegerQ(p), Not(RationalQ(m)), And(PositiveIntegerQ(n), NegativeIntegerQ(p + S(1)/2), LessEqual(S(-1), m, -n*(p + S(1)))))

    cons199 = CustomConstraint(cons_f199)
    pattern240 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons15, cons22, cons199)
    rule240 = ReplacementRule(pattern240, lambda e, n, m, a, x, p, d, b, c : -(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons68, cons3, cons189, cons56, cons15, cons22, cons199)
    rule241 = ReplacementRule(pattern241, lambda non2, e, b2, a1, a2, n, m, x, p, d, c, b1 : -(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1)), x)/(a1*a2*b1*b2*n*(p + S(1))) - (e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))*(-a1*a2*d + b1*b2*c)/(a1*a2*b1*b2*e*n*(p + S(1))))
    rubi.add(rule241)


    def cons_f200(n, m, p):
        return NonzeroQ(m + n*(p + S(1)) + S(1))

    cons200 = CustomConstraint(cons_f200)
    pattern242 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons135, cons200)
    rule242 = ReplacementRule(pattern242, lambda e, n, m, a, x, p, d, b, c : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(b*e*(m + n*(p + S(1)) + S(1))) - (a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**p, x)/(b*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons68, cons3, cons4, cons189, cons56, cons200)
    rule243 = ReplacementRule(pattern243, lambda non2, e, b2, a1, a2, n, m, x, p, d, c, b1 : d*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(b1*b2*e*(m + n*(p + S(1)) + S(1))) - (a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(b1*b2*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule243)


    def cons_f201(m):
        return Or(IntegerQ(m), PositiveIntegerQ(S(2)*m + S(2)), Not(RationalQ(m)))

    cons201 = CustomConstraint(cons_f201)
    pattern244 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons135, cons14, cons48, cons201)
    rule244 = ReplacementRule(pattern244, lambda e, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p/(c + d*x**n), x), x))
    rubi.add(rule244)


    def cons_f202(n):
        return Greater(n, S(0))

    cons202 = CustomConstraint(cons_f202)
    pattern245 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons14, cons193, cons82, cons202)
    rule245 = ReplacementRule(pattern245, lambda e, n, m, a, x, p, d, b, c : c**S(2)*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*Simp(-a*d**S(2)*x**n*(m + S(1)) + b*c**S(2)*n*(p + S(1)) + c*(m + S(1))*(-S(2)*a*d + b*c), x), x)/(a*(m + S(1))))
    rubi.add(rule245)

    pattern246 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons14, cons15, cons22)
    rule246 = ReplacementRule(pattern246, lambda e, n, m, a, x, p, d, b, c : Int((e*x)**m*(a + b*x**n)**(p + S(1))*Simp(a*b*d**S(2)*n*x**n*(p + S(1)) + b**S(2)*c**S(2)*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)**S(2), x), x)/(a*b**S(2)*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)**S(2)/(a*b**S(2)*e*n*(p + S(1))))
    rubi.add(rule246)


    def cons_f203(n, m, p):
        return NonzeroQ(m + n*(p + S(2)) + S(1))

    cons203 = CustomConstraint(cons_f203)
    pattern247 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons135, cons14, cons203)
    rule247 = ReplacementRule(pattern247, lambda e, n, m, a, x, p, d, b, c : d**S(2)*e**(-n + S(-1))*(e*x)**(m + n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*(p + S(2)) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*Simp(b*c**S(2)*(m + n*(p + S(2)) + S(1)) + d*x**n*(S(2)*b*c*n*(p + S(1)) + (-a*d + S(2)*b*c)*(m + n + S(1))), x), x)/(b*(m + n*(p + S(2)) + S(1))))
    rubi.add(rule247)

    pattern248 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons4, cons139, cons135, cons14, cons80, )
    def With248(q, n, m, a, x, p, d, b, c):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q, x), x, x**k)/k
        print("Unable to Integrate")
    rule248 = ReplacementRule(pattern248, lambda q, n, m, a, x, p, d, b, c : With248(q, n, m, a, x, p, d, b, c))
    rubi.add(rule248)

    pattern249 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons139, cons135, cons14, cons114, cons12, )
    def With249(q, e, n, m, a, x, p, d, b, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(k*n))**p*(c + d*e**(-n)*x**(k*n))**q, x), x, (e*x)**(S(1)/k))/e
    rule249 = ReplacementRule(pattern249, lambda q, e, n, m, a, x, p, d, b, c : With249(q, e, n, m, a, x, p, d, b, c))
    rubi.add(rule249)


    def cons_f204(q, m, p):
        return RationalQ(m, p, q)

    cons204 = CustomConstraint(cons_f204)

    def cons_f205(n, m):
        return Greater(m - n + S(1), S(0))

    cons205 = CustomConstraint(cons_f205)

    def cons_f206(q, e, n, m, a, x, p, d, b, c):
        return IntBinomialQ(a, b, c, d, e, m, n, p, q, x)

    cons206 = CustomConstraint(cons_f206)
    pattern250 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons135, cons14, cons204, cons22, cons144, cons205, cons206)
    rule250 = ReplacementRule(pattern250, lambda q, e, n, m, a, x, p, d, b, c : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(q + S(-1)) + S(1)), x), x)/(b*n*(p + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*n*(p + S(1))))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons135, cons14, cons163, cons22, cons166, cons206)
    rule251 = ReplacementRule(pattern251, lambda q, e, n, m, a, x, p, d, b, c : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule251)

    pattern252 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons135, cons14, cons163, cons22, cons164, cons206)
    rule252 = ReplacementRule(pattern252, lambda q, e, n, m, a, x, p, d, b, c : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule252)


    def cons_f207(n, m):
        return Greater(m - n + S(1), n)

    cons207 = CustomConstraint(cons_f207)
    pattern253 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons139, cons135, cons14, cons81, cons22, cons207, cons206)
    rule253 = ReplacementRule(pattern253, lambda q, e, n, m, a, x, p, d, b, c : -a*e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*n*(p + S(1))*(-a*d + b*c)) + e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*q - n + S(1)) + b*c*n*(p + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule253)


    def cons_f208(n, m):
        return Inequality(n, GreaterEqual, m - n + S(1), Greater, S(0))

    cons208 = CustomConstraint(cons_f208)
    pattern254 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons139, cons135, cons14, cons81, cons22, cons208, cons206)
    rule254 = ReplacementRule(pattern254, lambda q, e, n, m, a, x, p, d, b, c : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(n*(p + S(1))*(-a*d + b*c)) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons139, cons135, cons14, cons15, cons22, cons206)
    rule255 = ReplacementRule(pattern255, lambda q, e, n, m, a, x, p, d, b, c : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons135, cons14, cons204, cons144, cons82, cons16, cons206)
    rule256 = ReplacementRule(pattern256, lambda q, e, n, m, a, x, p, d, b, c : -e**(-n)*n*Int((e*x)**(m + n)*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*d*q + b*c*p + b*d*x**n*(p + q), x), x)/(m + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + S(1))))
    rubi.add(rule256)


    def cons_f209(q, m):
        return RationalQ(m, q)

    cons209 = CustomConstraint(cons_f209)
    pattern257 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons14, cons209, cons166, cons82, cons206)
    rule257 = ReplacementRule(pattern257, lambda q, e, n, m, a, x, p, d, b, c : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*n*(a*d*(q + S(-1)) + b*c*(p + S(1))) + c*(m + S(1))*(-a*d + b*c) + d*x**n*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)), x), x)/(a*(m + S(1))))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons14, cons209, cons164, cons82, cons206)
    rule258 = ReplacementRule(pattern258, lambda q, e, n, m, a, x, p, d, b, c : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(b*c*(m + S(1)) + d*x**n*(b*n*(p + q + S(1)) + b*(m + S(1))) + n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*(m + S(1))))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons135, cons14, cons163, cons144, cons16, cons206)
    rule259 = ReplacementRule(pattern259, lambda q, e, n, m, a, x, p, d, b, c : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule259)

    pattern260 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons4, cons135, cons14, cons143, cons166, cons206)
    rule260 = ReplacementRule(pattern260, lambda q, e, n, m, a, x, p, d, b, c : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule260)

    pattern261 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons135, cons14, cons209, cons144, cons205, cons206)
    rule261 = ReplacementRule(pattern261, lambda q, e, n, m, a, x, p, d, b, c : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(a*c*(m - n + S(1)) + x**n*(a*d*(m - n + S(1)) - n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule261)

    pattern262 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons139, cons135, cons14, cons106, cons207, cons206)
    rule262 = ReplacementRule(pattern262, lambda q, e, n, m, a, x, p, d, b, c : -e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*(q + S(-1)) + S(1)) + b*c*(m + n*(p + S(-1)) + S(1))), x), x)/(b*d*(m + n*(p + q) + S(1))) + e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q) + S(1))))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons139, cons135, cons14, cons106, cons82, cons206)
    rule263 = ReplacementRule(pattern263, lambda q, e, n, m, a, x, p, d, b, c : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(a*d*q + b*c*p) + (a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*e*(m + S(1))))
    rubi.add(rule263)


    def cons_f210(n, m):
        return LessEqual(n, m, S(2)*n + S(-1))

    cons210 = CustomConstraint(cons_f210)
    pattern264 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons135, cons14, cons106, cons210)
    rule264 = ReplacementRule(pattern264, lambda e, n, m, a, x, d, b, c : -a*e**n*Int((e*x)**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*e**n*Int((e*x)**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons135, cons14)
    rule265 = ReplacementRule(pattern265, lambda e, n, m, a, x, d, b, c : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule265)


    def cons_f211(n, m):
        return IntegersQ(m/S(2), n/S(2))

    cons211 = CustomConstraint(cons_f211)

    def cons_f212(n, m):
        return Less(S(0), m - n + S(1), n)

    cons212 = CustomConstraint(cons_f212)

    def cons_f213(n):
        return LessEqual(n, S(4))

    cons213 = CustomConstraint(cons_f213)
    pattern266 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons211, cons212, cons213)
    rule266 = ReplacementRule(pattern266, lambda n, m, a, x, d, b, c : -a*Int(x**(m - n)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(x**(m - n)/sqrt(c + d*x**n), x)/b)
    rubi.add(rule266)

    pattern267 = Pattern(Integral(x_**S(2)/((a_ + x_**S(4)*WC('b', S(1)))*sqrt(c_ + x_**S(4)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, )
    def With267(a, x, d, b, c):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(sqrt(c + d*x**S(4))*(r - s*x**S(2))), x)/(S(2)*b) + s*Int(S(1)/(sqrt(c + d*x**S(4))*(r + s*x**S(2))), x)/(S(2)*b)
    rule267 = ReplacementRule(pattern267, lambda a, x, d, b, c : With267(a, x, d, b, c))
    rubi.add(rule267)


    def cons_f214(d, b, c, a):
        return ZeroQ(-a*d + S(4)*b*c)

    cons214 = CustomConstraint(cons_f214)
    pattern268 = Pattern(Integral(x_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons214, )
    def With268(a, x, d, b, c):
        q = Rt(d/c, S(3))
        return -S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) - sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) + sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) - sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) - S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) + sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) + S(2)**(S(1)/3)*q*atanh(sqrt(c + d*x**S(3))/sqrt(c))/(S(18)*b*sqrt(c))
    rule268 = ReplacementRule(pattern268, lambda a, x, d, b, c : With268(a, x, d, b, c))
    rubi.add(rule268)


    def cons_f215(m):
        return PositiveIntegerQ(m/S(3) + S(-1)/3)

    cons215 = CustomConstraint(cons_f215)
    pattern269 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons214, cons215)
    rule269 = ReplacementRule(pattern269, lambda m, a, x, d, b, c : -a*Int(x**(m + S(-3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b + Int(x**(m + S(-3))/sqrt(c + d*x**S(3)), x)/b)
    rubi.add(rule269)


    def cons_f216(m):
        return NegativeIntegerQ(m/S(3) + S(-1)/3)

    cons216 = CustomConstraint(cons_f216)
    pattern270 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons214, cons216)
    rule270 = ReplacementRule(pattern270, lambda m, a, x, d, b, c : -b*Int(x**(m + S(3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/a + Int(x**m/sqrt(c + d*x**S(3)), x)/a)
    rubi.add(rule270)

    pattern271 = Pattern(Integral(x_**S(2)*sqrt(c_ + x_**S(4)*WC('d', S(1)))/(a_ + x_**S(4)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons135)
    rule271 = ReplacementRule(pattern271, lambda a, x, d, b, c : d*Int(x**S(2)/sqrt(c + d*x**S(4)), x)/b + (-a*d + b*c)*Int(x**S(2)/((a + b*x**S(4))*sqrt(c + d*x**S(4))), x)/b)
    rubi.add(rule271)


    def cons_f217(m):
        return IntegerQ(m/S(3) + S(-1)/3)

    cons217 = CustomConstraint(cons_f217)
    pattern272 = Pattern(Integral(x_**WC('m', S(1))*sqrt(c_ + x_**S(3)*WC('d', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons135, cons214, cons217)
    rule272 = ReplacementRule(pattern272, lambda m, a, x, d, b, c : d*Int(x**m/sqrt(c + d*x**S(3)), x)/b + (-a*d + b*c)*Int(x**m/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b)
    rubi.add(rule272)

    pattern273 = Pattern(Integral(x_**S(2)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons18, cons162, cons171)
    rule273 = ReplacementRule(pattern273, lambda a, x, d, b, c : -c*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/b + x*sqrt(a + b*x**S(2))/(b*sqrt(c + d*x**S(2))))
    rubi.add(rule273)


    def cons_f218(n):
        return Or(EqQ(n, S(2)), EqQ(n, S(4)))

    cons218 = CustomConstraint(cons_f218)

    def cons_f219(n, a, d, b, c):
        return Not(And(EqQ(n, S(2)), SimplerSqrtQ(-b/a, -d/c)))

    cons219 = CustomConstraint(cons_f219)
    pattern274 = Pattern(Integral(x_**n_/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons135, cons218, cons219)
    rule274 = ReplacementRule(pattern274, lambda n, a, x, d, b, c : -a*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b)
    rubi.add(rule274)


    def cons_f220(q, n, m, p):
        return IntegersQ(p + (m + S(1))/n, q)

    cons220 = CustomConstraint(cons_f220)
    pattern275 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons14, cons81, cons220, cons42, )
    def With275(q, n, m, a, x, p, d, b, c):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(c - x**k*(-a*d + b*c))**q*(-b*x**k + S(1))**(-p - q + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    rule275 = ReplacementRule(pattern275, lambda q, n, m, a, x, p, d, b, c : With275(q, n, m, a, x, p, d, b, c))
    rubi.add(rule275)

    pattern276 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons4, cons139, cons135, cons46, cons80)
    rule276 = ReplacementRule(pattern276, lambda q, n, m, a, x, p, d, b, c : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons4, cons139, cons46, cons114, )
    def With277(q, e, n, m, a, x, p, d, b, c):
        g = Denominator(m)
        return -g*Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(-g*n))**p*(c + d*e**(-n)*x**(-g*n))**q, x), x, (e*x)**(-S(1)/g))/e
    rule277 = ReplacementRule(pattern277, lambda q, e, n, m, a, x, p, d, b, c : With277(q, e, n, m, a, x, p, d, b, c))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons4, cons139, cons135, cons46, cons119)
    rule278 = ReplacementRule(pattern278, lambda q, e, n, m, a, x, p, d, b, c : -(e*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons68, cons4, cons139, cons135, cons47, )
    def With279(q, n, m, a, x, p, d, b, c):
        g = Denominator(n)
        return g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    rule279 = ReplacementRule(pattern279, lambda q, n, m, a, x, p, d, b, c : With279(q, n, m, a, x, p, d, b, c))
    rubi.add(rule279)

    pattern280 = Pattern(Integral((e_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons4, cons139, cons135, cons47)
    rule280 = ReplacementRule(pattern280, lambda q, e, n, m, a, x, p, d, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule280)

    pattern281 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons139, cons135, cons120, cons121)
    rule281 = ReplacementRule(pattern281, lambda q, n, m, a, x, p, d, b, c : Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule281)

    pattern282 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons135, cons120, cons121)
    rule282 = ReplacementRule(pattern282, lambda q, e, n, m, a, x, p, d, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons163, cons22, cons166, cons206)
    rule283 = ReplacementRule(pattern283, lambda q, e, n, m, a, x, p, d, b, c : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons163, cons22, cons164, cons206)
    rule284 = ReplacementRule(pattern284, lambda q, e, n, m, a, x, p, d, b, c : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons139, cons135, cons15, cons22, cons206)
    rule285 = ReplacementRule(pattern285, lambda q, e, n, m, a, x, p, d, b, c : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons163, cons144, cons16, cons206)
    rule286 = ReplacementRule(pattern286, lambda q, e, n, m, a, x, p, d, b, c : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons135, cons143, cons166, cons206)
    rule287 = ReplacementRule(pattern287, lambda q, e, n, m, a, x, p, d, b, c : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule287)


    def cons_f221(n, m):
        return Or(ZeroQ(m - n), ZeroQ(m - S(2)*n + S(1)))

    cons221 = CustomConstraint(cons_f221)
    pattern288 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons135, cons221)
    rule288 = ReplacementRule(pattern288, lambda n, m, a, x, d, b, c : -a*Int(x**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*Int(x**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons3, cons68, cons135)
    rule289 = ReplacementRule(pattern289, lambda e, n, m, a, x, d, b, c : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule289)


    def cons_f222(q, m, p):
        return IntegersQ(m, p, q)

    cons222 = CustomConstraint(cons_f222)

    def cons_f223(p):
        return GreaterEqual(p, S(-2))

    cons223 = CustomConstraint(cons_f223)

    def cons_f224(q, m):
        return Or(GreaterEqual(q, S(-2)), And(Equal(q, S(-3)), IntegerQ(m/S(2) + S(-1)/2)))

    cons224 = CustomConstraint(cons_f224)
    pattern290 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons135, cons222, cons223, cons224)
    rule290 = ReplacementRule(pattern290, lambda q, e, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons178, cons179, cons180)
    rule291 = ReplacementRule(pattern291, lambda q, mn, n, m, a, x, p, d, b, c : Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons139, cons178, cons181, cons66)
    rule292 = ReplacementRule(pattern292, lambda q, mn, n, m, a, x, p, d, b, c : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral((e_*x_)**m_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons178)
    rule293 = ReplacementRule(pattern293, lambda q, e, mn, n, m, a, x, p, d, b, c : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q, x))
    rubi.add(rule293)


    def cons_f225(n, m):
        return NonzeroQ(m - n + S(1))

    cons225 = CustomConstraint(cons_f225)
    pattern294 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons135, cons75, cons225, cons19, cons173)
    rule294 = ReplacementRule(pattern294, lambda q, e, n, m, a, x, p, d, b, c : a**p*c**q*(e*x)**(m + S(1))*AppellF1((m + S(1))/n, -p, -q, S(1) + (m + S(1))/n, -b*x**n/a, -d*x**n/c)/(e*(m + S(1))))
    rubi.add(rule294)

    pattern295 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_), cons6, cons2, cons67, cons137, cons186, cons68, cons3, cons4, cons139, cons135, cons75, cons225, cons20)
    rule295 = ReplacementRule(pattern295, lambda q, e, n, m, a, x, p, d, b, c : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((e*x)**m*(S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule295)

    pattern296 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons3, cons4, cons139, cons132, cons80, cons133)
    rule296 = ReplacementRule(pattern296, lambda q, n, m, a, v, p, d, b, c, x : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule296)

    pattern297 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons68, cons3, cons4, cons139, cons134)
    rule297 = ReplacementRule(pattern297, lambda q, u, n, m, a, v, p, d, b, c, x : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons3, cons4, cons139, cons189, cons56, cons57)
    rule298 = ReplacementRule(pattern298, lambda q, non2, a1, b2, a2, u, n, x, p, d, c, b1 : Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule298)


    def cons_f226(n2, n):
        return ZeroQ(-S(2)*n + n2)

    cons226 = CustomConstraint(cons_f226)
    pattern299 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons3, cons4, cons139, cons189, cons226, cons56, cons57)
    rule299 = ReplacementRule(pattern299, lambda q, non2, e, b2, a1, a2, u, n, x, p, d, n2, c, b1 : Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**p_*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons3, cons4, cons139, cons189, cons56)
    rule300 = ReplacementRule(pattern300, lambda q, non2, a1, a2, u, n, p, x, b2, d, c, b1 : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_), cons58, cons59, cons60, cons61, cons67, cons137, cons186, cons3, cons4, cons139, cons189, cons226, cons56)
    rule301 = ReplacementRule(pattern301, lambda q, non2, e, b2, a1, a2, u, n, x, p, d, n2, c, b1 : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule301)


    def cons_f227(q, r, p):
        return PositiveIntegerQ(p, q, r)

    cons227 = CustomConstraint(cons_f227)

    def cons_f228(f, x):
        return FreeQ(f, x)

    cons228 = CustomConstraint(cons_f228)
    pattern302 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons227)
    rule302 = ReplacementRule(pattern302, lambda q, f, e, r, n, a, x, p, d, b, c : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule302)


    def cons_f229(f, e, n, a, x, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, n), x)

    cons229 = CustomConstraint(cons_f229)
    pattern303 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons229)
    rule303 = ReplacementRule(pattern303, lambda f, e, n, a, x, d, b, c : (-a*f + b*e)*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - (-c*f + d*e)*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons229)
    rule304 = ReplacementRule(pattern304, lambda f, e, n, a, x, d, b, c : f*Int(S(1)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b)
    rubi.add(rule304)


    def cons_f230(n, a, d, b, c):
        return Not(And(ZeroQ(n + S(-2)), Or(And(PosQ(b/a), PosQ(d/c)), And(NegQ(b/a), Or(PosQ(d/c), And(PositiveQ(a), Or(Not(PositiveQ(c)), SimplerSqrtQ(-b/a, -d/c))))))))

    cons230 = CustomConstraint(cons_f230)
    pattern305 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons230)
    rule305 = ReplacementRule(pattern305, lambda f, e, n, a, x, d, b, c : f*Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b)
    rubi.add(rule305)

    pattern306 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons18, cons162)
    rule306 = ReplacementRule(pattern306, lambda f, e, a, x, d, b, c : (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons163, cons22, cons144)
    rule307 = ReplacementRule(pattern307, lambda q, f, e, n, a, x, p, d, b, c : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + S(1)) + b*e) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(n*q + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons139, cons15, cons22)
    rule308 = ReplacementRule(pattern308, lambda q, f, e, n, a, x, p, d, b, c : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule308)


    def cons_f231(q, n, p):
        return NonzeroQ(n*(p + q + S(1)) + S(1))

    cons231 = CustomConstraint(cons_f231)
    pattern309 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons4, cons143, cons144, cons231)
    rule309 = ReplacementRule(pattern309, lambda q, f, e, n, a, x, p, d, b, c : f*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(n*(p + q + S(1)) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + q + S(1)) + b*e) + x**n*(b*d*e*n*(p + q + S(1)) + d*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(n*(p + q + S(1)) + S(1))))
    rubi.add(rule309)


    def cons_f232(f, e, a, x, d, b, c):
        return FreeQ(List(a, b, c, d, e, f), x)

    cons232 = CustomConstraint(cons_f232)
    pattern310 = Pattern(Integral((e_ + x_**S(4)*WC('f', S(1)))/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule310 = ReplacementRule(pattern310, lambda f, e, a, x, d, b, c : (-a*f + b*e)*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - (-c*f + d*e)*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule310)


    def cons_f233(f, e, n, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, p, n), x)

    cons233 = CustomConstraint(cons_f233)
    pattern311 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons3, cons233)
    rule311 = ReplacementRule(pattern311, lambda f, e, n, a, x, p, d, b, c : f*Int((a + b*x**n)**p, x)/d + (-c*f + d*e)*Int((a + b*x**n)**p/(c + d*x**n), x)/d)
    rubi.add(rule311)


    def cons_f234(q, f, e, n, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, n, p, q), x)

    cons234 = CustomConstraint(cons_f234)
    pattern312 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons4, cons139, cons234)
    rule312 = ReplacementRule(pattern312, lambda q, f, e, n, a, x, p, d, b, c : e*Int((a + b*x**n)**p*(c + d*x**n)**q, x) + f*Int(x**n*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule312)

    pattern313 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule313 = ReplacementRule(pattern313, lambda f, e, a, x, d, b, c : b*Int(S(1)/((a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c) - d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c))
    rubi.add(rule313)


    def cons_f235(d, f, e, c):
        return NonzeroQ(-c*f + d*e)

    cons235 = CustomConstraint(cons_f235)
    pattern314 = Pattern(Integral(S(1)/(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons67, cons137, cons186, cons228, cons235)
    rule314 = ReplacementRule(pattern314, lambda f, e, x, d, c : -d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/c + Int(S(1)/(x**S(2)*sqrt(e + f*x**S(2))), x)/c)
    rubi.add(rule314)


    def cons_f236(d, c):
        return PositiveQ(d/c)

    cons236 = CustomConstraint(cons_f236)

    def cons_f237(f, e):
        return PositiveQ(f/e)

    cons237 = CustomConstraint(cons_f237)

    def cons_f238(d, f, c, e):
        return Not(SimplerSqrtQ(d/c, f/e))

    cons238 = CustomConstraint(cons_f238)
    pattern315 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons236, cons237, cons238)
    rule315 = ReplacementRule(pattern315, lambda f, e, a, x, d, b, c : d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule315)


    def cons_f239(d, f, e, c):
        return Not(SimplerSqrtQ(-f/e, -d/c))

    cons239 = CustomConstraint(cons_f239)
    pattern316 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons239)
    rule316 = ReplacementRule(pattern316, lambda f, e, a, x, d, b, c : d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule316)


    def cons_f240(f, e):
        return PosQ(f/e)

    cons240 = CustomConstraint(cons_f240)
    pattern317 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons162, cons240, cons238)
    rule317 = ReplacementRule(pattern317, lambda f, e, a, x, d, b, c : b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*f + b*e) - f*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*f + b*e))
    rubi.add(rule317)


    def cons_f241(e):
        return PositiveQ(e)

    cons241 = CustomConstraint(cons_f241)

    def cons_f242(d, f, e, c):
        return Not(And(NegQ(f/e), SimplerSqrtQ(-f/e, -d/c)))

    cons242 = CustomConstraint(cons_f242)
    pattern318 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons172, cons173, cons241, cons242)
    rule318 = ReplacementRule(pattern318, lambda f, e, a, x, d, b, c : EllipticPi(b*c/(a*d), asin(x*Rt(-d/c, S(2))), c*f/(d*e))/(a*sqrt(c)*sqrt(e)*Rt(-d/c, S(2))))
    rubi.add(rule318)

    pattern319 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons176)
    rule319 = ReplacementRule(pattern319, lambda f, e, a, x, d, b, c : sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule319)

    pattern320 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons162)
    rule320 = ReplacementRule(pattern320, lambda f, e, a, x, d, b, c : c*sqrt(e + f*x**S(2))*EllipticPi(S(1) - b*c/(a*d), ArcTan(x*Rt(d/c, S(2))), -c*f/(d*e) + S(1))/(a*e*sqrt(c*(e + f*x**S(2))/(e*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule320)

    pattern321 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons172)
    rule321 = ReplacementRule(pattern321, lambda f, e, a, x, d, b, c : d*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b + (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b)
    rubi.add(rule321)

    pattern322 = Pattern(Integral(sqrt(e_ + x_**S(2)*WC('f', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons162, cons240)
    rule322 = ReplacementRule(pattern322, lambda f, e, a, x, d, b, c : b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - d*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons162, cons240)
    rule323 = ReplacementRule(pattern323, lambda f, e, a, x, d, b, c : (-a*f + b*e)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons162, cons240)
    rule324 = ReplacementRule(pattern324, lambda f, e, a, x, d, b, c : d*Int(sqrt(e + f*x**S(2))*(-a*d + S(2)*b*c + b*d*x**S(2))/sqrt(c + d*x**S(2)), x)/b**S(2) + (-a*d + b*c)**S(2)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b**S(2))
    rubi.add(rule324)


    def cons_f243(q, r):
        return RationalQ(q, r)

    cons243 = CustomConstraint(cons_f243)

    def cons_f244(q):
        return Less(q, S(-1))

    cons244 = CustomConstraint(cons_f244)

    def cons_f245(r):
        return Greater(r, S(1))

    cons245 = CustomConstraint(cons_f245)
    pattern325 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons243, cons244, cons245)
    rule325 = ReplacementRule(pattern325, lambda q, f, e, a, x, r, d, b, c : b*(-a*f + b*e)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**(r + S(-1))/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - Int((c + d*x**S(2))**q*(e + f*x**S(2))**(r + S(-1))*(-a*d**S(2)*e - b*c**S(2)*f + S(2)*b*c*d*e + d**S(2)*x**S(2)*(-a*f + b*e)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule325)


    def cons_f246(r, x):
        return FreeQ(r, x)

    cons246 = CustomConstraint(cons_f246)
    pattern326 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons246, cons143, cons166)
    rule326 = ReplacementRule(pattern326, lambda q, f, e, a, x, r, d, b, c : d*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r, x)/b + (-a*d + b*c)*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/b)
    rubi.add(rule326)

    pattern327 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons246, cons143, cons244)
    rule327 = ReplacementRule(pattern327, lambda q, f, e, a, x, r, d, b, c : b**S(2)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r*(-a*d + S(2)*b*c + b*d*x**S(2)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule327)


    def cons_f247(q):
        return LessEqual(q, S(-1))

    cons247 = CustomConstraint(cons_f247)
    pattern328 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons246, cons143, cons247)
    rule328 = ReplacementRule(pattern328, lambda q, f, e, a, x, r, d, b, c : b*Int((c + d*x**S(2))**(q + S(1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r, x)/(-a*d + b*c))
    rubi.add(rule328)

    pattern329 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**S(2), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule329 = ReplacementRule(pattern329, lambda f, e, a, x, d, b, c : x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))) + d*f*Int((a - b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)) + (-a**S(2)*d*f + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)))
    rubi.add(rule329)

    pattern330 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**S(2)*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule330 = ReplacementRule(pattern330, lambda f, e, a, x, d, b, c : b**S(2)*x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))*(-a*d + b*c)*(-a*f + b*e)) - d*f*Int((a + b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)) + (S(3)*a**S(2)*d*f - S(2)*a*b*(c*f + d*e) + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule330)

    pattern331 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons246, cons145, cons143, cons144)
    rule331 = ReplacementRule(pattern331, lambda q, f, e, r, n, a, x, p, d, b, c : d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b + (-a*d + b*c)*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b)
    rubi.add(rule331)

    pattern332 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons139, cons145, cons143, cons247)
    rule332 = ReplacementRule(pattern332, lambda q, f, e, r, n, a, x, p, d, b, c : b*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(1))*(e + f*x**n)**r, x)/(-a*d + b*c) - d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(e + f*x**n)**r, x)/(-a*d + b*c))
    rubi.add(rule332)

    pattern333 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule333 = ReplacementRule(pattern333, lambda f, e, a, x, d, b, c : sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule333)

    pattern334 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule334 = ReplacementRule(pattern334, lambda f, e, a, x, d, b, c : a*sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)*(-b*x**S(2) + S(1))), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule334)

    pattern335 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule335 = ReplacementRule(pattern335, lambda f, e, a, x, d, b, c : sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)/sqrt(S(1) - x**S(2)*(-a*f + b*e)/e), x), x, x/sqrt(a + b*x**S(2)))/(a*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule335)


    def cons_f248(d, f, e, c):
        return PosQ((-c*f + d*e)/c)

    cons248 = CustomConstraint(cons_f248)
    pattern336 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons248)
    rule336 = ReplacementRule(pattern336, lambda f, e, a, x, d, b, c : b*c*(-c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f) - c*(-c*f + d*e)*Int(sqrt(a + b*x**S(2))/((c + d*x**S(2))**(S(3)/2)*sqrt(e + f*x**S(2))), x)/(S(2)*f) + d*x*sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*f*sqrt(c + d*x**S(2))) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f))
    rubi.add(rule336)


    def cons_f249(d, f, e, c):
        return NegQ((-c*f + d*e)/c)

    cons249 = CustomConstraint(cons_f249)
    pattern337 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons249)
    rule337 = ReplacementRule(pattern337, lambda f, e, a, x, d, b, c : e*(-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/(S(2)*f) + x*sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))/(S(2)*sqrt(e + f*x**S(2))) + (-a*f + b*e)*(-S(2)*c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*f**S(2)) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(e + f*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(S(2)*f**S(2)))
    rubi.add(rule337)

    pattern338 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/(e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons232)
    rule338 = ReplacementRule(pattern338, lambda f, e, a, x, d, b, c : b*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/f - (-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/f)
    rubi.add(rule338)

    pattern339 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons139, cons246, cons14, )
    def With339(q, f, e, r, n, a, x, p, d, b, c):
        u = ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
        if SumQ(u):
            return Int(u, x)
        print("Unable to Integrate")
    rule339 = ReplacementRule(pattern339, lambda q, f, e, r, n, a, x, p, d, b, c : With339(q, f, e, r, n, a, x, p, d, b, c))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons139, cons246, cons46)
    rule340 = ReplacementRule(pattern340, lambda q, f, e, r, n, a, x, p, d, b, c : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r/x**S(2), x), x, S(1)/x))
    rubi.add(rule340)


    def cons_f250(q, f, e, r, n, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, n, p, q, r), x)

    cons250 = CustomConstraint(cons_f250)
    pattern341 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons4, cons139, cons246, cons250)
    rule341 = ReplacementRule(pattern341, lambda q, f, e, r, n, a, x, p, d, b, c : Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule341)


    def cons_f251(u, v):
        return ZeroQ(u - v)

    cons251 = CustomConstraint(cons_f251)

    def cons_f252(u, w):
        return ZeroQ(u - w)

    cons252 = CustomConstraint(cons_f252)
    pattern342 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons3, cons139, cons246, cons251, cons252, cons54, cons55)
    rule342 = ReplacementRule(pattern342, lambda q, f, e, w, u, r, n, a, v, p, d, b, c, x : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule342)

    pattern343 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons4, cons246, cons178, cons179)
    rule343 = ReplacementRule(pattern343, lambda q, f, e, mn, r, n, a, x, p, d, b, c : Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule343)


    def cons_f253(r):
        return IntegerQ(r)

    cons253 = CustomConstraint(cons_f253)
    pattern344 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons139, cons178, cons12, cons253)
    rule344 = ReplacementRule(pattern344, lambda q, f, e, mn, r, n, a, x, p, d, b, c : Int(x**(n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule344)

    pattern345 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons3, cons4, cons139, cons246, cons178, cons181)
    rule345 = ReplacementRule(pattern345, lambda q, f, e, mn, r, n, a, x, p, d, b, c : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule345)


    def cons_f254(n2, n):
        return ZeroQ(-n/S(2) + n2)

    cons254 = CustomConstraint(cons_f254)

    def cons_f255(e2, f1, e1, f2):
        return ZeroQ(e1*f2 + e2*f1)

    cons255 = CustomConstraint(cons_f255)

    def cons_f256(e2, e1, r):
        return Or(IntegerQ(r), And(PositiveQ(e1), PositiveQ(e2)))

    cons256 = CustomConstraint(cons_f256)

    def cons_f257(e1, x):
        return FreeQ(e1, x)

    cons257 = CustomConstraint(cons_f257)

    def cons_f258(f1, x):
        return FreeQ(f1, x)

    cons258 = CustomConstraint(cons_f258)

    def cons_f259(e2, x):
        return FreeQ(e2, x)

    cons259 = CustomConstraint(cons_f259)

    def cons_f260(f2, x):
        return FreeQ(f2, x)

    cons260 = CustomConstraint(cons_f260)
    pattern346 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons257, cons258, cons259, cons260, cons3, cons4, cons139, cons246, cons254, cons255, cons256)
    rule346 = ReplacementRule(pattern346, lambda q, f2, e2, n, e1, f1, b, a, x, r, d, n2, c, p : Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule346)

    pattern347 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons257, cons258, cons259, cons260, cons3, cons4, cons139, cons246, cons254, cons255)
    rule347 = ReplacementRule(pattern347, lambda q, f2, e2, n, e1, f1, b, a, x, r, d, n2, c, p : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule347)


    def cons_f261(m, g):
        return Or(IntegerQ(m), PositiveQ(g))

    cons261 = CustomConstraint(cons_f261)

    def cons_f262(g, x):
        return FreeQ(g, x)

    cons262 = CustomConstraint(cons_f262)
    pattern348 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons261, cons70)
    rule348 = ReplacementRule(pattern348, lambda q, f, e, g, n, m, x, r, d, b, c, p : b**(S(1) - (m + S(1))/n)*g**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule348)

    pattern349 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons261, cons71)
    rule349 = ReplacementRule(pattern349, lambda q, f, e, g, n, m, x, r, d, b, c, p : b**IntPart(p)*g**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule349)

    pattern350 = Pattern(Integral((g_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons72)
    rule350 = ReplacementRule(pattern350, lambda q, f, e, g, n, m, x, r, d, b, c, p : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule350)


    def cons_f263(q, r, p):
        return PositiveIntegerQ(p + S(2), q, r)

    cons263 = CustomConstraint(cons_f263)
    pattern351 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons263)
    rule351 = ReplacementRule(pattern351, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons139, cons246, cons187)
    rule352 = ReplacementRule(pattern352, lambda q, f, e, r, n, m, a, x, p, d, b, c : Subst(Int((a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule352)


    def cons_f264(q, r, p):
        return IntegersQ(p, q, r)

    cons264 = CustomConstraint(cons_f264)
    pattern353 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons264, cons73)
    rule353 = ReplacementRule(pattern353, lambda q, f, e, r, n, m, a, x, p, d, b, c : Int(x**(m + n*(p + q + r))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q*(e*x**(-n) + f)**r, x))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons139, cons246, cons70)
    rule354 = ReplacementRule(pattern354, lambda q, f, e, r, n, m, a, x, p, d, b, c : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule354)

    pattern355 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons70)
    rule355 = ReplacementRule(pattern355, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons139, cons246, cons14, cons80, )
    def With356(q, f, e, r, n, m, a, x, p, d, b, c):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q*(e + f*x**(n/k))**r, x), x, x**k)/k
        print("Unable to Integrate")
    rule356 = ReplacementRule(pattern356, lambda q, f, e, r, n, m, a, x, p, d, b, c : With356(q, f, e, r, n, m, a, x, p, d, b, c))
    rubi.add(rule356)

    pattern357 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons4, cons139, cons246, cons14, cons114, )
    def With357(q, f, e, g, r, n, m, a, x, p, d, b, c):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(k*n))**p*(c + d*g**(-n)*x**(k*n))**q*(e + f*g**(-n)*x**(k*n))**r, x), x, (g*x)**(S(1)/k))/g
    rule357 = ReplacementRule(pattern357, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : With357(q, f, e, g, r, n, m, a, x, p, d, b, c))
    rubi.add(rule357)


    def cons_f265(q, f, e, a, d, b, c):
        return Not(And(Equal(q, S(1)), SimplerQ(-a*d + b*c, -a*f + b*e)))

    cons265 = CustomConstraint(cons_f265)
    pattern358 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons14, cons163, cons22, cons144, cons265)
    rule358 = ReplacementRule(pattern358, lambda q, f, e, g, n, m, a, x, p, d, b, c : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule358)

    pattern359 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons139, cons14, cons81, cons22, cons205)
    rule359 = ReplacementRule(pattern359, lambda q, f, e, g, n, m, a, x, p, d, b, c : -g**n*Int((g*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e)*(m - n + S(1)) + x**n*(-b*n*(p + S(1))*(c*f - d*e) + d*(-a*f + b*e)*(m + n*q + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)) + g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule359)

    pattern360 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons139, cons14, cons15, cons22)
    rule360 = ReplacementRule(pattern360, lambda q, f, e, g, n, m, a, x, p, d, b, c : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule360)


    def cons_f266(q, f, e, n, x, d, c):
        return Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n)))

    cons266 = CustomConstraint(cons_f266)
    pattern361 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons4, cons14, cons209, cons144, cons82, cons266)
    rule361 = ReplacementRule(pattern361, lambda q, f, e, g, n, m, a, x, p, d, b, c : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*g*(m + S(1))) - g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + e*n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))))
    rubi.add(rule361)

    pattern362 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons14, cons143, cons144, cons266)
    rule362 = ReplacementRule(pattern362, lambda q, f, e, g, n, m, a, x, p, d, b, c : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule362)

    pattern363 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons4, cons139, cons14, cons106, cons107)
    rule363 = ReplacementRule(pattern363, lambda q, f, e, g, n, m, a, x, p, d, b, c : f*g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q + S(1)) + S(1))) - g**n*Int((g*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m - n + S(1)) + x**n*(a*d*f*(m + n*q + S(1)) + b*(c*f*(m + n*p + S(1)) - d*e*(m + n*(p + q + S(1)) + S(1)))), x), x)/(b*d*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule363)

    pattern364 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons4, cons139, cons14, cons106, cons82)
    rule364 = ReplacementRule(pattern364, lambda q, f, e, g, n, m, a, x, p, d, b, c : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*g*(m + S(1))) + g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m + S(1)) - b*d*e*x**n*(m + n*(p + q + S(2)) + S(1)) - e*n*(a*d*q + b*c*p) - e*(a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))))
    rubi.add(rule364)

    pattern365 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons14)
    rule365 = ReplacementRule(pattern365, lambda f, e, g, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule365)

    pattern366 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons139, cons14)
    rule366 = ReplacementRule(pattern366, lambda q, f, e, g, n, m, a, x, p, d, b, c : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule366)


    def cons_f267(r):
        return PositiveIntegerQ(r)

    cons267 = CustomConstraint(cons_f267)
    pattern367 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons139, cons14, cons267)
    rule367 = ReplacementRule(pattern367, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x))
    rubi.add(rule367)

    pattern368 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons4, cons139, cons246, cons46, cons80)
    rule368 = ReplacementRule(pattern368, lambda q, f, e, r, n, m, a, x, p, d, b, c : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons4, cons139, cons246, cons46, cons114, )
    def With369(q, f, e, g, r, n, m, a, x, p, d, b, c):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(-k*n))**p*(c + d*g**(-n)*x**(-k*n))**q*(e + f*g**(-n)*x**(-k*n))**r, x), x, (g*x)**(-S(1)/k))/g
    rule369 = ReplacementRule(pattern369, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : With369(q, f, e, g, r, n, m, a, x, p, d, b, c))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons139, cons246, cons46, cons119)
    rule370 = ReplacementRule(pattern370, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : -(g*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x))
    rubi.add(rule370)

    pattern371 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons4, cons139, cons246, cons47, )
    def With371(q, f, e, r, n, m, a, x, p, d, b, c):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p*(c + d*x**(k*n))**q*(e + f*x**(k*n))**r, x), x, x**(S(1)/k))
    rule371 = ReplacementRule(pattern371, lambda q, f, e, r, n, m, a, x, p, d, b, c : With371(q, f, e, r, n, m, a, x, p, d, b, c))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((g_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons4, cons139, cons246, cons47)
    rule372 = ReplacementRule(pattern372, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule372)

    pattern373 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons139, cons246, cons120)
    rule373 = ReplacementRule(pattern373, lambda q, f, e, r, n, m, a, x, p, d, b, c : Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q*(e + f*x**(n/(m + S(1))))**r, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons120)
    rule374 = ReplacementRule(pattern374, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule374)

    pattern375 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons163, cons22, cons144, cons265)
    rule375 = ReplacementRule(pattern375, lambda q, f, e, g, n, m, a, x, p, d, b, c : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule375)

    pattern376 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons139, cons15, cons22)
    rule376 = ReplacementRule(pattern376, lambda q, f, e, g, n, m, a, x, p, d, b, c : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule376)

    pattern377 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons143, cons144, cons266)
    rule377 = ReplacementRule(pattern377, lambda q, f, e, g, n, m, a, x, p, d, b, c : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule377)


    def cons_f268(f, e, g, n, m, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, g, m, n, p), x)

    cons268 = CustomConstraint(cons_f268)
    pattern378 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons268)
    rule378 = ReplacementRule(pattern378, lambda f, e, g, n, m, a, x, p, d, b, c : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule378)


    def cons_f269(q, f, e, g, n, m, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q), x)

    cons269 = CustomConstraint(cons_f269)
    pattern379 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons269)
    rule379 = ReplacementRule(pattern379, lambda q, f, e, g, n, m, a, x, p, d, b, c : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + f*x**(-m)*(g*x)**m*Int(x**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule379)

    pattern380 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons246, cons178, cons179)
    rule380 = ReplacementRule(pattern380, lambda q, f, e, mn, r, n, m, a, x, p, d, b, c : Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule380)

    pattern381 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons139, cons178, cons12, cons253)
    rule381 = ReplacementRule(pattern381, lambda q, f, e, mn, r, n, m, a, x, p, d, b, c : Int(x**(m + n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule381)

    pattern382 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons139, cons246, cons178, cons181)
    rule382 = ReplacementRule(pattern382, lambda q, f, e, mn, r, n, m, a, x, p, d, b, c : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule382)

    pattern383 = Pattern(Integral((g_*x_)**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons178)
    rule383 = ReplacementRule(pattern383, lambda q, f, e, g, mn, r, n, m, a, x, p, d, b, c : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q*(e + f*x**n)**r, x))
    rubi.add(rule383)


    def cons_f270(q, f, e, g, r, n, m, a, x, p, d, b, c):
        return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q, r), x)

    cons270 = CustomConstraint(cons_f270)
    pattern384 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons262, cons68, cons3, cons4, cons139, cons246, cons270)
    rule384 = ReplacementRule(pattern384, lambda q, f, e, g, r, n, m, a, x, p, d, b, c : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule384)

    pattern385 = Pattern(Integral(u_**WC('m', S(1))*(e_ + v_**n_*WC('f', S(1)))**WC('r', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_), cons6, cons2, cons67, cons137, cons186, cons228, cons68, cons3, cons4, cons139, cons246, cons134)
    rule385 = ReplacementRule(pattern385, lambda q, f, e, u, r, n, m, a, v, p, d, b, c, x : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule385)

    pattern386 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons257, cons258, cons259, cons260, cons262, cons68, cons3, cons4, cons139, cons246, cons254, cons255, cons256)
    rule386 = ReplacementRule(pattern386, lambda q, f2, g, e2, n, e1, m, f1, b, a, x, r, d, n2, c, p : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule386)

    pattern387 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_), cons6, cons2, cons67, cons137, cons257, cons258, cons259, cons260, cons262, cons68, cons3, cons4, cons139, cons246, cons254, cons255)
    rule387 = ReplacementRule(pattern387, lambda q, f2, g, e2, n, e1, m, f1, b, a, x, r, d, n2, c, p : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule387)

    return rubi
