from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (sympy_op_factory, Int, Sum, Set, With, Module, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec, ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma, FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist)
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

def linear_products(rubi):
    pattern1 = Pattern(Integral(S(1)/x_, x_))
    rule1 = ReplacementRule(pattern1, lambda x : log(x))
    rubi.add(rule1)


    def cons_f1(m):
        return NonzeroQ(m + S(1))

    cons1 = CustomConstraint(cons_f1)

    def cons_f2(m, x):
        return FreeQ(m, x)

    cons2 = CustomConstraint(cons_f2)
    pattern2 = Pattern(Integral(x_**WC('m', S(1)), x_), cons2, cons1)
    rule2 = ReplacementRule(pattern2, lambda m, x : x**(m + S(1))/(m + S(1)))
    rubi.add(rule2)


    def cons_f3(x, a, b):
        return FreeQ(List(a, b), x)

    cons3 = CustomConstraint(cons_f3)

    def cons_f4(a, x):
        return FreeQ(a, x)

    cons4 = CustomConstraint(cons_f4)

    def cons_f5(b, x):
        return FreeQ(b, x)

    cons5 = CustomConstraint(cons_f5)
    pattern3 = Pattern(Integral(S(1)/(a_ + x_*WC('b', S(1))), x_), cons4, cons5, cons3)
    rule3 = ReplacementRule(pattern3, lambda x, a, b : log(RemoveContent(a + b*x, x))/b)
    rubi.add(rule3)

    pattern4 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons4, cons5, cons2, cons1)
    rule4 = ReplacementRule(pattern4, lambda m, x, a, b : (a + b*x)**(m + S(1))/(b*(m + S(1))))
    rubi.add(rule4)


    def cons_f6(u, x):
        return LinearQ(u, x)

    cons6 = CustomConstraint(cons_f6)

    def cons_f7(u, x):
        return NonzeroQ(u - x)

    cons7 = CustomConstraint(cons_f7)
    pattern5 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons4, cons5, cons2, cons6, cons7)
    rule5 = ReplacementRule(pattern5, lambda m, u, b, x, a : Subst(Int((a + b*x)**m, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule5)


    def cons_f8(a, d, b, c):
        return ZeroQ(a*d + b*c)

    cons8 = CustomConstraint(cons_f8)

    def cons_f9(c, x):
        return FreeQ(c, x)

    cons9 = CustomConstraint(cons_f9)

    def cons_f10(d, x):
        return FreeQ(d, x)

    cons10 = CustomConstraint(cons_f10)
    pattern6 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_), cons4, cons5, cons9, cons10, cons8)
    rule6 = ReplacementRule(pattern6, lambda c, d, b, x, a : Int(S(1)/(a*c + b*d*x**S(2)), x))
    rubi.add(rule6)


    def cons_f11(a, d, b, c):
        return NonzeroQ(-a*d + b*c)

    cons11 = CustomConstraint(cons_f11)
    pattern7 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons11)
    rule7 = ReplacementRule(pattern7, lambda d, c, b, x, a : b*Int(S(1)/(a + b*x), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule7)


    def cons_f12(n, m):
        return ZeroQ(m + n + S(2))

    cons12 = CustomConstraint(cons_f12)

    def cons_f13(n, x):
        return FreeQ(n, x)

    cons13 = CustomConstraint(cons_f13)
    pattern8 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_), cons4, cons5, cons9, cons10, cons2, cons13, cons11, cons12, cons1)
    rule8 = ReplacementRule(pattern8, lambda d, c, m, b, x, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule8)


    def cons_f14(m):
        return PositiveIntegerQ(m + S(1)/2)

    cons14 = CustomConstraint(cons_f14)
    pattern9 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons4, cons5, cons9, cons10, cons8, cons14)
    rule9 = ReplacementRule(pattern9, lambda d, c, m, b, x, a : S(2)*a*c*m*Int((a + b*x)**(m + S(-1))*(c + d*x)**(m + S(-1)), x)/(S(2)*m + S(1)) + x*(a + b*x)**m*(c + d*x)**m/(S(2)*m + S(1)))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(3)/2)*(c_ + x_*WC('d', S(1)))**(S(3)/2)), x_), cons4, cons5, cons9, cons10, cons8)
    rule10 = ReplacementRule(pattern10, lambda c, d, b, x, a : x/(a*c*sqrt(a + b*x)*sqrt(c + d*x)))
    rubi.add(rule10)


    def cons_f15(m):
        return NegativeIntegerQ(m + S(3)/2)

    cons15 = CustomConstraint(cons_f15)
    pattern11 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons4, cons5, cons9, cons10, cons8, cons15)
    rule11 = ReplacementRule(pattern11, lambda d, c, m, b, x, a : -x*(a + b*x)**(m + S(1))*(c + d*x)**(m + S(1))/(S(2)*a*c*(m + S(1))) + (S(2)*m + S(3))*Int((a + b*x)**(m + S(1))*(c + d*x)**(m + S(1)), x)/(S(2)*a*c*(m + S(1))))
    rubi.add(rule11)


    def cons_f16(m, c, a):
        return Or(IntegerQ(m), And(PositiveQ(a), PositiveQ(c)))

    cons16 = CustomConstraint(cons_f16)
    pattern12 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('m', S(1)), x_), cons4, cons5, cons9, cons10, cons2, cons8, cons16)
    rule12 = ReplacementRule(pattern12, lambda d, c, m, b, x, a : Int((a*c + b*d*x**S(2))**m, x))
    rubi.add(rule12)


    def cons_f17(a):
        return PositiveQ(a)

    cons17 = CustomConstraint(cons_f17)

    def cons_f18(a, c):
        return ZeroQ(a + c)

    cons18 = CustomConstraint(cons_f18)
    pattern13 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons4, cons5, cons9, cons10, cons8, cons17, cons18)
    rule13 = ReplacementRule(pattern13, lambda c, d, b, x, a : acosh(b*x/a)/b)
    rubi.add(rule13)

    pattern14 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons4, cons5, cons9, cons10, cons8)
    rule14 = ReplacementRule(pattern14, lambda c, d, b, x, a : S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)))
    rubi.add(rule14)


    def cons_f19(m):
        return Not(IntegerQ(S(2)*m))

    cons19 = CustomConstraint(cons_f19)
    pattern15 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_), cons4, cons5, cons9, cons10, cons2, cons8, cons19)
    rule15 = ReplacementRule(pattern15, lambda d, c, m, b, x, a : (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((a*c + b*d*x**S(2))**m, x))
    rubi.add(rule15)


    def cons_f20(a, c, b, d):
        return PosQ(b*d/(a*c))

    cons20 = CustomConstraint(cons_f20)
    pattern16 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(5)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons4, cons5, cons9, cons10, cons8, cons20)
    rule16 = ReplacementRule(pattern16, lambda c, d, b, x, a : (-a*d + b*c)*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(2)*b) - S(2)/(b*(a + b*x)**(S(1)/4)*(c + d*x)**(S(1)/4)))
    rubi.add(rule16)

    pattern17 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(9)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_), cons4, cons5, cons9, cons10, cons8, cons20)
    rule17 = ReplacementRule(pattern17, lambda c, d, b, x, a : -d*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(5)*b) - S(4)/(S(5)*b*(a + b*x)**(S(5)/4)*(c + d*x)**(S(1)/4)))
    rubi.add(rule17)


    def cons_f21(m):
        return IntegerQ(m + S(1)/2)

    cons21 = CustomConstraint(cons_f21)

    def cons_f22(n):
        return IntegerQ(n + S(1)/2)

    cons22 = CustomConstraint(cons_f22)

    def cons_f23(n, m):
        return Less(S(0), m, n)

    cons23 = CustomConstraint(cons_f23)
    pattern18 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons4, cons5, cons9, cons10, cons8, cons21, cons22, cons23)
    rule18 = ReplacementRule(pattern18, lambda d, c, m, b, x, n, a : S(2)*c*n*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(m + n + S(1)) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))))
    rubi.add(rule18)


    def cons_f24(n, m):
        return Less(m, n, S(0))

    cons24 = CustomConstraint(cons_f24)
    pattern19 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons4, cons5, cons9, cons10, cons8, cons21, cons22, cons24)
    rule19 = ReplacementRule(pattern19, lambda d, c, m, b, x, n, a : (m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(S(2)*a*(m + S(1))) - (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule19)


    def cons_f25(m):
        return PositiveIntegerQ(m)

    cons25 = CustomConstraint(cons_f25)

    def cons_f26(n, m, c):
        return Or(Not(IntegerQ(n)), And(ZeroQ(c), LessEqual(S(7)*m + S(4)*n, S(0))), Less(S(9)*m + S(5)*n + S(5), S(0)), Greater(m + n + S(2), S(0)))

    cons26 = CustomConstraint(cons_f26)
    pattern20 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons13, cons11, cons25, cons26)
    rule20 = ReplacementRule(pattern20, lambda d, c, m, b, x, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x))
    rubi.add(rule20)


    def cons_f27(m):
        return NegativeIntegerQ(m)

    cons27 = CustomConstraint(cons_f27)

    def cons_f28(n):
        return IntegerQ(n)

    cons28 = CustomConstraint(cons_f28)

    def cons_f29(n, m):
        return Not(And(PositiveIntegerQ(n), Less(m + n + S(2), S(0))))

    cons29 = CustomConstraint(cons_f29)
    pattern21 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons13, cons11, cons27, cons28, cons29)
    rule21 = ReplacementRule(pattern21, lambda d, c, m, b, x, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x))
    rubi.add(rule21)


    def cons_f30(n):
        return RationalQ(n)

    cons30 = CustomConstraint(cons_f30)

    def cons_f31(n):
        return Greater(n, S(0))

    cons31 = CustomConstraint(cons_f31)
    pattern22 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons11, cons30, cons31)
    rule22 = ReplacementRule(pattern22, lambda d, c, b, x, n, a : (-a*d + b*c)*Int((c + d*x)**(n + S(-1))/(a + b*x), x)/b + (c + d*x)**n/(b*n))
    rubi.add(rule22)


    def cons_f32(n):
        return Less(n, S(-1))

    cons32 = CustomConstraint(cons_f32)
    pattern23 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons11, cons30, cons32)
    rule23 = ReplacementRule(pattern23, lambda d, c, b, x, n, a : b*Int((c + d*x)**(n + S(1))/(a + b*x), x)/(-a*d + b*c) - (c + d*x)**(n + S(1))/((n + S(1))*(-a*d + b*c)))
    rubi.add(rule23)


    def cons_f33(a, d, b, c):
        return PosQ((-a*d + b*c)/b)

    cons33 = CustomConstraint(cons_f33)
    pattern24 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons33, )
    def With24(d, c, b, x, a):
        q = Rt((-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    rule24 = ReplacementRule(pattern24, lambda d, c, b, x, a : With24(d, c, b, x, a))
    rubi.add(rule24)


    def cons_f34(a, d, b, c):
        return NegQ((-a*d + b*c)/b)

    cons34 = CustomConstraint(cons_f34)
    pattern25 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons34, )
    def With25(d, c, b, x, a):
        q = Rt(-(-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    rule25 = ReplacementRule(pattern25, lambda d, c, b, x, a : With25(d, c, b, x, a))
    rubi.add(rule25)

    pattern26 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons4, cons5, cons9, cons10, cons33, )
    def With26(d, c, b, x, a):
        q = Rt((-a*d + b*c)/b, S(3))
        return -S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    rule26 = ReplacementRule(pattern26, lambda d, c, b, x, a : With26(d, c, b, x, a))
    rubi.add(rule26)

    pattern27 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons4, cons5, cons9, cons10, cons34, )
    def With27(d, c, b, x, a):
        q = Rt(-(-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    rule27 = ReplacementRule(pattern27, lambda d, c, b, x, a : With27(d, c, b, x, a))
    rubi.add(rule27)


    def cons_f35(n):
        return Less(S(-1), n, S(0))

    cons35 = CustomConstraint(cons_f35)
    pattern28 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons11, cons30, cons35, )
    def With28(d, c, b, x, n, a):
        p = Denominator(n)
        return p*Subst(Int(x**(p*(n + S(1)) + S(-1))/(a*d - b*c + b*x**p), x), x, (c + d*x)**(S(1)/p))
    rule28 = ReplacementRule(pattern28, lambda d, c, b, x, n, a : With28(d, c, b, x, n, a))
    rubi.add(rule28)


    def cons_f36(n):
        return Not(IntegerQ(n))

    cons36 = CustomConstraint(cons_f36)
    pattern29 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_/x_, x_), cons9, cons10, cons13, cons36)
    rule29 = ReplacementRule(pattern29, lambda n, c, x, d : -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), S(1) + d*x/c)/(c*(n + S(1))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons13, cons11, cons36)
    rule30 = ReplacementRule(pattern30, lambda d, c, b, x, n, a : -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), TogetherSimplify(b*(c + d*x)/(-a*d + b*c)))/((n + S(1))*(-a*d + b*c)))
    rubi.add(rule30)


    def cons_f37(n, m):
        return RationalQ(m, n)

    cons37 = CustomConstraint(cons_f37)

    def cons_f38(m):
        return Less(m, S(-1))

    cons38 = CustomConstraint(cons_f38)

    def cons_f39(n, m):
        return Not(And(IntegerQ(n), Not(IntegerQ(m))))

    cons39 = CustomConstraint(cons_f39)

    def cons_f40(n, m):
        return Not(And(IntegerQ(m + n), LessEqual(m + n + S(2), S(0)), Or(FractionQ(m), GreaterEqual(m + S(2)*n + S(1), S(0)))))

    cons40 = CustomConstraint(cons_f40)

    def cons_f41(d, c, m, b, x, n, a):
        return IntLinearcQ(a, b, c, d, m, n, x)

    cons41 = CustomConstraint(cons_f41)
    pattern31 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons11, cons37, cons38, cons31, cons39, cons40, cons41)
    rule31 = ReplacementRule(pattern31, lambda d, c, m, b, x, n, a : -d*n*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1)), x)/(b*(m + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + S(1))))
    rubi.add(rule31)


    def cons_f42(m, n, a, c):
        return Not(And(Less(n, S(-1)), Or(ZeroQ(a), And(NonzeroQ(c), Less(m, n), IntegerQ(n)))))

    cons42 = CustomConstraint(cons_f42)
    pattern32 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons11, cons37, cons38, cons42, cons41)
    rule32 = ReplacementRule(pattern32, lambda d, c, m, b, x, n, a : -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule32)


    def cons_f43(n, m):
        return Unequal(m + n + S(1), S(0))

    cons43 = CustomConstraint(cons_f43)

    def cons_f44(n, m):
        return Not(And(PositiveIntegerQ(m), Or(Not(IntegerQ(n)), Less(S(0), m, n))))

    cons44 = CustomConstraint(cons_f44)

    def cons_f45(n, m):
        return Not(And(IntegerQ(m + n), Less(m + n + S(2), S(0))))

    cons45 = CustomConstraint(cons_f45)
    pattern33 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons11, cons37, cons31, cons43, cons44, cons45, cons41)
    rule33 = ReplacementRule(pattern33, lambda d, c, m, b, x, n, a : n*(-a*d + b*c)*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(b*(m + n + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))))
    rubi.add(rule33)


    def cons_f46(b, d):
        return ZeroQ(b + d)

    cons46 = CustomConstraint(cons_f46)

    def cons_f47(a, c):
        return PositiveQ(a + c)

    cons47 = CustomConstraint(cons_f47)
    pattern34 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons46, cons47)
    rule34 = ReplacementRule(pattern34, lambda d, c, b, x, a : Int(S(1)/sqrt(a*c - b**S(2)*x**S(2) - b*x*(a - c)), x))
    rubi.add(rule34)


    def cons_f48(a, d, b, c):
        return PositiveQ(-a*d + b*c)

    cons48 = CustomConstraint(cons_f48)

    def cons_f49(b):
        return PositiveQ(b)

    cons49 = CustomConstraint(cons_f49)
    pattern35 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons48, cons49)
    rule35 = ReplacementRule(pattern35, lambda d, c, b, x, a : S(2)*Subst(Int(S(1)/sqrt(-a*d + b*c + d*x**S(2)), x), x, sqrt(a + b*x))/sqrt(b))
    rubi.add(rule35)


    def cons_f50(b, d):
        return ZeroQ(b - d)

    cons50 = CustomConstraint(cons_f50)
    pattern36 = Pattern(Integral(S(1)/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons11, cons50)
    rule36 = ReplacementRule(pattern36, lambda c, d, b, x, a : S(2)*Subst(Int(S(1)/sqrt(-a + c + x**S(2)), x), x, sqrt(a + b*x))/b)
    rubi.add(rule36)

    pattern37 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons11)
    rule37 = ReplacementRule(pattern37, lambda d, c, b, x, a : S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)))
    rubi.add(rule37)


    def cons_f51(m):
        return RationalQ(m)

    cons51 = CustomConstraint(cons_f51)

    def cons_f52(m):
        return Less(S(-1), m, S(0))

    cons52 = CustomConstraint(cons_f52)

    def cons_f53(m):
        return LessEqual(S(3), Denominator(m), S(4))

    cons53 = CustomConstraint(cons_f53)
    pattern38 = Pattern(Integral((c_ + x_*WC('d', S(1)))**m_*(x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_), cons4, cons5, cons9, cons10, cons11, cons51, cons52, cons53)
    rule38 = ReplacementRule(pattern38, lambda d, c, m, b, x, a : (a + b*x)**m*(c + d*x)**m*(a*c + b*d*x**S(2) + x*(a*d + b*c))**(-m)*Int((a*c + b*d*x**S(2) + x*(a*d + b*c))**m, x))
    rubi.add(rule38)


    def cons_f54(b, d):
        return PosQ(d/b)

    cons54 = CustomConstraint(cons_f54)
    pattern39 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons4, cons5, cons9, cons10, cons11, cons54, )
    def With39(d, c, b, x, a):
        q = Rt(d/b, S(3))
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d - q*log(c + d*x)/(S(2)*d) - S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(-1))/(S(2)*d)
    rule39 = ReplacementRule(pattern39, lambda d, c, b, x, a : With39(d, c, b, x, a))
    rubi.add(rule39)


    def cons_f55(b, d):
        return NegQ(d/b)

    cons55 = CustomConstraint(cons_f55)
    pattern40 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_), cons4, cons5, cons9, cons10, cons11, cons55, )
    def With40(d, c, b, x, a):
        q = Rt(-d/b, S(3))
        return sqrt(S(3))*q*ArcTan(-S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d + q*log(c + d*x)/(S(2)*d) + S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(1))/(S(2)*d)
    rule40 = ReplacementRule(pattern40, lambda d, c, b, x, a : With40(d, c, b, x, a))
    rubi.add(rule40)


    def cons_f56(n, m):
        return Equal(m + n + S(1), S(0))

    cons56 = CustomConstraint(cons_f56)
    pattern41 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons11, cons37, cons52, cons56, )
    def With41(d, c, m, b, x, n, a):
        p = Denominator(m)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))/(b - d*x**p), x), x, (a + b*x)**(S(1)/p)*(c + d*x)**(-S(1)/p))
    rule41 = ReplacementRule(pattern41, lambda d, c, m, b, x, n, a : With41(d, c, m, b, x, n, a))
    rubi.add(rule41)


    def cons_f57(n, m):
        return LessEqual(Denominator(n), Denominator(m))

    cons57 = CustomConstraint(cons_f57)
    pattern42 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons11, cons37, cons52, cons35, cons57, cons41, )
    def With42(d, c, m, b, x, n, a):
        p = Denominator(m)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))*(-a*d/b + c + d*x**p/b)**n, x), x, (a + b*x)**(S(1)/p))/b
    rule42 = ReplacementRule(pattern42, lambda d, c, m, b, x, n, a : With42(d, c, m, b, x, n, a))
    rubi.add(rule42)


    def cons_f58(n, m):
        return NegativeIntegerQ(m + n + S(2))

    cons58 = CustomConstraint(cons_f58)

    def cons_f59(n, m):
        return Or(SumSimplerQ(m, S(1)), Not(SumSimplerQ(n, S(1))))

    cons59 = CustomConstraint(cons_f59)
    pattern43 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_), cons4, cons5, cons9, cons10, cons2, cons13, cons11, cons58, cons1, cons59)
    rule43 = ReplacementRule(pattern43, lambda d, c, m, b, x, n, a : -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule43)


    def cons_f60(m):
        return Not(IntegerQ(m))

    cons60 = CustomConstraint(cons_f60)

    def cons_f61(n, d, b, c):
        return Or(IntegerQ(n), And(PositiveQ(c), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2)), PositiveQ(-d/(b*c))))))

    cons61 = CustomConstraint(cons_f61)
    pattern44 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons5, cons9, cons10, cons2, cons13, cons60, cons61)
    rule44 = ReplacementRule(pattern44, lambda c, d, m, b, x, n : c**n*(b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*x/c)/(b*(m + S(1))))
    rubi.add(rule44)


    def cons_f62(m, c, b, d):
        return Or(IntegerQ(m), PositiveQ(-d/(b*c)))

    cons62 = CustomConstraint(cons_f62)
    pattern45 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons5, cons9, cons10, cons2, cons13, cons36, cons62)
    rule45 = ReplacementRule(pattern45, lambda c, d, m, b, x, n : (-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*Hypergeometric2F1(-m, n + S(1), n + S(2), S(1) + d*x/c)/(d*(n + S(1))))
    rubi.add(rule45)


    def cons_f63(c):
        return Not(PositiveQ(c))

    cons63 = CustomConstraint(cons_f63)

    def cons_f64(c, b, d):
        return Not(PositiveQ(-d/(b*c)))

    cons64 = CustomConstraint(cons_f64)

    def cons_f65(n, m, d, c):
        return Or(And(RationalQ(m), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2))))), Not(RationalQ(n)))

    cons65 = CustomConstraint(cons_f65)
    pattern46 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons5, cons9, cons10, cons2, cons13, cons60, cons36, cons63, cons64, cons65)
    rule46 = ReplacementRule(pattern46, lambda c, d, m, b, x, n : c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n, x))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons5, cons9, cons10, cons2, cons13, cons60, cons36, cons63, cons64)
    rule47 = ReplacementRule(pattern47, lambda c, d, m, b, x, n : (b*x)**FracPart(m)*(-b*c/d)**IntPart(m)*(-d*x/c)**(-FracPart(m))*Int((-d*x/c)**m*(c + d*x)**n, x))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons4, cons5, cons9, cons10, cons2, cons11, cons60, cons28)
    rule48 = ReplacementRule(pattern48, lambda d, c, m, b, x, n, a : b**(-n + S(-1))*(a + b*x)**(m + S(1))*(-a*d + b*c)**n*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(m + S(1)))
    rubi.add(rule48)


    def cons_f66(a, d, b, c):
        return PositiveQ(b/(-a*d + b*c))

    cons66 = CustomConstraint(cons_f66)

    def cons_f67(c, d, m, b, n, a):
        return Or(RationalQ(m), Not(And(RationalQ(n), PositiveQ(-d/(-a*d + b*c)))))

    cons67 = CustomConstraint(cons_f67)
    pattern49 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons4, cons5, cons9, cons10, cons2, cons13, cons11, cons60, cons36, cons66, cons67)
    rule49 = ReplacementRule(pattern49, lambda d, c, m, b, x, n, a : (b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(b*(m + S(1))))
    rubi.add(rule49)


    def cons_f68(n, m):
        return Or(RationalQ(m), Not(SimplerQ(n + S(1), m + S(1))))

    cons68 = CustomConstraint(cons_f68)
    pattern50 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_), cons4, cons5, cons9, cons10, cons2, cons13, cons11, cons60, cons36, cons68)
    rule50 = ReplacementRule(pattern50, lambda d, c, m, b, x, n, a : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule50)


    def cons_f69(u, x):
        return NonzeroQ(Coefficient(u, x, S(0)))

    cons69 = CustomConstraint(cons_f69)
    pattern51 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons2, cons13, cons6, cons69)
    rule51 = ReplacementRule(pattern51, lambda d, c, m, u, b, x, n, a : Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule51)


    def cons_f70(n, m):
        return ZeroQ(m - n)

    cons70 = CustomConstraint(cons_f70)

    def cons_f71(m):
        return IntegerQ(m)

    cons71 = CustomConstraint(cons_f71)

    def cons_f72(e, x):
        return FreeQ(e, x)

    cons72 = CustomConstraint(cons_f72)

    def cons_f73(f, x):
        return FreeQ(f, x)

    cons73 = CustomConstraint(cons_f73)

    def cons_f74(p, x):
        return FreeQ(p, x)

    cons74 = CustomConstraint(cons_f74)
    pattern52 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons8, cons70, cons71)
    rule52 = ReplacementRule(pattern52, lambda e, d, f, c, m, b, x, p, n, a : Int((e + f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule52)


    def cons_f75(n, p):
        return NonzeroQ(n + p + S(2))

    cons75 = CustomConstraint(cons_f75)

    def cons_f76(e, f, d, c, b, p, n, a):
        return ZeroQ(a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))

    cons76 = CustomConstraint(cons_f76)
    pattern53 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons75, cons76)
    rule53 = ReplacementRule(pattern53, lambda e, d, c, f, b, x, p, n, a : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))))
    rubi.add(rule53)


    def cons_f77(p):
        return PositiveIntegerQ(p)

    cons77 = CustomConstraint(cons_f77)

    def cons_f78(e, b, f, a):
        return ZeroQ(a*f + b*e)

    cons78 = CustomConstraint(cons_f78)

    def cons_f79(n, p):
        return Not(And(NegativeIntegerQ(n + p + S(2)), Greater(n + S(2)*p, S(0))))

    cons79 = CustomConstraint(cons_f79)
    pattern54 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons10, cons72, cons73, cons13, cons77, cons78, cons79)
    rule54 = ReplacementRule(pattern54, lambda e, d, f, b, x, p, n, a : Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x))
    rubi.add(rule54)


    def cons_f80(n, p):
        return Or(NonzeroQ(n + S(1)), Equal(p, S(1)))

    cons80 = CustomConstraint(cons_f80)

    def cons_f81(e, b, f, a):
        return NonzeroQ(a*f + b*e)

    cons81 = CustomConstraint(cons_f81)

    def cons_f82(e, d, f, b, p, n, a):
        return Or(Not(IntegerQ(n)), Less(S(5)*n + S(9)*p, S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, d, e, f)))

    cons82 = CustomConstraint(cons_f82)
    pattern55 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons10, cons72, cons73, cons13, cons77, cons80, cons81, cons82)
    rule55 = ReplacementRule(pattern55, lambda e, d, f, b, x, p, n, a : Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x))
    rubi.add(rule55)


    def cons_f83(e, c, d, f, b, p, n, a):
        return Or(NegativeIntegerQ(n, p), ZeroQ(p + S(-1)), And(PositiveIntegerQ(p), Or(Not(IntegerQ(n)), LessEqual(S(5)*n + S(9)*p + S(10), S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, c, d, e, f)))))

    cons83 = CustomConstraint(cons_f83)
    pattern56 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons11, cons83)
    rule56 = ReplacementRule(pattern56, lambda e, d, f, c, b, x, p, n, a : Int(ExpandIntegrand((a + b*x)*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule56)


    def cons_f84(n, p):
        return ZeroQ(n + p + S(2))

    cons84 = CustomConstraint(cons_f84)

    def cons_f85(p):
        return NonzeroQ(p + S(1))

    cons85 = CustomConstraint(cons_f85)

    def cons_f86(n, p):
        return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1)))))

    cons86 = CustomConstraint(cons_f86)
    pattern57 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons84, cons85, cons86)
    rule57 = ReplacementRule(pattern57, lambda e, d, c, f, b, x, p, n, a : b*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/f - (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule57)


    def cons_f87(p):
        return RationalQ(p)

    cons87 = CustomConstraint(cons_f87)

    def cons_f88(p):
        return Less(p, S(-1))

    cons88 = CustomConstraint(cons_f88)

    def cons_f89(n, e, c, p):
        return Or(Not(And(RationalQ(n), Less(n, S(-1)))), IntegerQ(p), Not(Or(IntegerQ(n), Not(Or(ZeroQ(e), Not(Or(ZeroQ(c), Less(p, n))))))))

    cons89 = CustomConstraint(cons_f89)
    pattern58 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons75, cons87, cons88, cons89)
    rule58 = ReplacementRule(pattern58, lambda e, d, c, f, b, x, p, n, a : -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule58)


    def cons_f90(p):
        return Not(RationalQ(p))

    cons90 = CustomConstraint(cons_f90)

    def cons_f91(p):
        return SumSimplerQ(p, S(1))

    cons91 = CustomConstraint(cons_f91)
    pattern59 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons75, cons90, cons91)
    rule59 = ReplacementRule(pattern59, lambda e, d, c, f, b, x, p, n, a : -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule59)

    pattern60 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons75)
    rule60 = ReplacementRule(pattern60, lambda e, d, c, f, b, x, p, n, a : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))) + (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**p, x)/(d*f*(n + p + S(2))))
    rubi.add(rule60)


    def cons_f92(n, p):
        return NonzeroQ(n + p + S(3))

    cons92 = CustomConstraint(cons_f92)

    def cons_f93(e, f, d, c, b, p, n, a):
        return ZeroQ(-b*(c*f*(p + S(1)) + d*e*(n + S(1)))*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) + d*f*(a**S(2)*d*f*(n + p + S(3)) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e))*(n + p + S(2)))

    cons93 = CustomConstraint(cons_f93)
    pattern61 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons75, cons92, cons93)
    rule61 = ReplacementRule(pattern61, lambda e, d, c, f, b, x, p, n, a : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(S(2)*a*d*f*(n + p + S(3)) + b*d*f*x*(n + p + S(2)) - b*(c*f*(p + S(2)) + d*e*(n + S(2))))/(d**S(2)*f**S(2)*(n + p + S(2))*(n + p + S(3))))
    rubi.add(rule61)


    def cons_f94(n, m):
        return ZeroQ(m - n + S(-1))

    cons94 = CustomConstraint(cons_f94)

    def cons_f95(m):
        return Not(PositiveIntegerQ(m))

    cons95 = CustomConstraint(cons_f95)

    def cons_f96(n, m, p):
        return NonzeroQ(m + n + p + S(2))

    cons96 = CustomConstraint(cons_f96)
    pattern62 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons2, cons13, cons74, cons8, cons94, cons90, cons95, cons96)
    rule62 = ReplacementRule(pattern62, lambda f, c, d, m, b, x, p, n, a : a*Int((f*x)**p*(a + b*x)**n*(c + d*x)**n, x) + b*Int((f*x)**(p + S(1))*(a + b*x)**n*(c + d*x)**n, x)/f)
    rubi.add(rule62)


    def cons_f97(p):
        return IntegerQ(p)

    cons97 = CustomConstraint(cons_f97)
    pattern63 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons97)
    rule63 = ReplacementRule(pattern63, lambda e, d, c, f, b, x, p, a : Int(ExpandIntegrand((e + f*x)**p/((a + b*x)*(c + d*x)), x), x))
    rubi.add(rule63)


    def cons_f98(p):
        return Less(S(0), p, S(1))

    cons98 = CustomConstraint(cons_f98)
    pattern64 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons87, cons98)
    rule64 = ReplacementRule(pattern64, lambda e, d, c, f, b, x, p, a : (-a*f + b*e)*Int((e + f*x)**(p + S(-1))/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule64)


    def cons_f99(p):
        return Greater(p, S(1))

    cons99 = CustomConstraint(cons_f99)
    pattern65 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons87, cons99)
    rule65 = ReplacementRule(pattern65, lambda e, d, c, f, b, x, p, a : f*(e + f*x)**(p + S(-1))/(b*d*(p + S(-1))) + Int((e + f*x)**(p + S(-2))*(-a*c*f**S(2) + b*d*e**S(2) + f*x*(-a*d*f - b*c*f + S(2)*b*d*e))/((a + b*x)*(c + d*x)), x)/(b*d))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons87, cons88)
    rule66 = ReplacementRule(pattern66, lambda e, d, c, f, b, x, p, a : f*(e + f*x)**(p + S(1))/((p + S(1))*(-a*f + b*e)*(-c*f + d*e)) + Int((e + f*x)**(p + S(1))*(-a*d*f - b*c*f + b*d*e - b*d*f*x)/((a + b*x)*(c + d*x)), x)/((-a*f + b*e)*(-c*f + d*e)))
    rubi.add(rule66)


    def cons_f100(p):
        return Not(IntegerQ(p))

    cons100 = CustomConstraint(cons_f100)
    pattern67 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons100)
    rule67 = ReplacementRule(pattern67, lambda e, d, c, f, b, x, p, a : b*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - d*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule67)


    def cons_f101(n):
        return PositiveIntegerQ(n)

    cons101 = CustomConstraint(cons_f101)

    def cons_f102(p):
        return FractionQ(p)

    cons102 = CustomConstraint(cons_f102)
    pattern68 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons101, cons102, cons88)
    rule68 = ReplacementRule(pattern68, lambda e, d, c, f, b, x, p, n, a : Int(ExpandIntegrand((e + f*x)**FractionalPart(p), (c + d*x)**n*(e + f*x)**IntegerPart(p)/(a + b*x), x), x))
    rubi.add(rule68)


    def cons_f103(n, m):
        return IntegersQ(m, n)

    cons103 = CustomConstraint(cons_f103)

    def cons_f104(n, m, p):
        return Or(IntegerQ(p), And(Greater(m, S(0)), GreaterEqual(n, S(-1))))

    cons104 = CustomConstraint(cons_f104)
    pattern69 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons103, cons104)
    rule69 = ReplacementRule(pattern69, lambda e, d, c, f, m, b, x, p, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule69)


    def cons_f105(n, p):
        return Or(And(RationalQ(n), Less(n, S(-1))), And(ZeroQ(n + p + S(3)), NonzeroQ(n + S(1)), Or(SumSimplerQ(n, S(1)), Not(SumSimplerQ(p, S(1))))))

    cons105 = CustomConstraint(cons_f105)
    pattern70 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons105)
    rule70 = ReplacementRule(pattern70, lambda e, d, c, f, b, x, p, n, a : (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*d + b*c)**S(2)/(d**S(2)*(n + S(1))*(-c*f + d*e)) - Int((c + d*x)**(n + S(1))*(e + f*x)**p*Simp(a**S(2)*d**S(2)*f*(n + p + S(2)) - S(2)*a*b*d*(c*f*(p + S(1)) + d*e*(n + S(1))) + b**S(2)*c*(c*f*(p + S(1)) + d*e*(n + S(1))) - b**S(2)*d*x*(n + S(1))*(-c*f + d*e), x), x)/(d**S(2)*(n + S(1))*(-c*f + d*e)))
    rubi.add(rule70)

    pattern71 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons92)
    rule71 = ReplacementRule(pattern71, lambda e, d, c, f, b, x, p, n, a : b*(a + b*x)*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(3))) + Int((c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(n + p + S(3)) + b*x*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e), x), x)/(d*f*(n + p + S(3))))
    rubi.add(rule71)


    def cons_f106(e, d, c, f, b, x, a):
        return FreeQ(List(a, b, c, d, e, f), x)

    cons106 = CustomConstraint(cons_f106)
    pattern72 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons106, )
    def With72(e, d, c, f, b, x, a):
        q = Rt((-c*f + d*e)/(-a*f + b*e), S(3))
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(-c*f + d*e) + q*log(e + f*x)/(-S(2)*c*f + S(2)*d*e) - S(3)*q*log(q*(a + b*x)**(S(1)/3) - (c + d*x)**(S(1)/3))/(-S(2)*c*f + S(2)*d*e)
    rule72 = ReplacementRule(pattern72, lambda e, d, c, f, b, x, a : With72(e, d, c, f, b, x, a))
    rubi.add(rule72)


    def cons_f107(e, f, c, d, b, a):
        return ZeroQ(S(2)*b*d*e - f*(a*d + b*c))

    cons107 = CustomConstraint(cons_f107)
    pattern73 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons107)
    rule73 = ReplacementRule(pattern73, lambda e, d, c, f, b, x, a : b*f*Subst(Int(S(1)/(b*f**S(2)*x**S(2) + d*(-a*f + b*e)**S(2)), x), x, sqrt(a + b*x)*sqrt(c + d*x)))
    rubi.add(rule73)


    def cons_f108(n, m):
        return ZeroQ(m + n + S(1))

    cons108 = CustomConstraint(cons_f108)

    def cons_f109(d, c, b, x, a):
        return SimplerQ(a + b*x, c + d*x)

    cons109 = CustomConstraint(cons_f109)
    pattern74 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons108, cons37, cons52, cons109, )
    def With74(e, d, c, f, m, b, x, n, a):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))/(-a*f + b*e - x**q*(-c*f + d*e)), x), x, (a + b*x)**(S(1)/q)*(c + d*x)**(-S(1)/q))
    rule74 = ReplacementRule(pattern74, lambda e, d, c, f, m, b, x, n, a : With74(e, d, c, f, m, b, x, n, a))
    rubi.add(rule74)


    def cons_f110(n, m, p):
        return ZeroQ(m + n + p + S(2))

    cons110 = CustomConstraint(cons_f110)

    def cons_f111(m, p):
        return Not(And(SumSimplerQ(p, S(1)), Not(SumSimplerQ(m, S(1)))))

    cons111 = CustomConstraint(cons_f111)
    pattern75 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons74, cons110, cons30, cons31, cons111)
    rule75 = ReplacementRule(pattern75, lambda e, d, c, f, m, b, x, p, n, a : -n*(-c*f + d*e)*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p, x)/((m + S(1))*(-a*f + b*e)) + (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule75)


    def cons_f112(n, m, p):
        return ZeroQ(m + n + p + S(3))

    cons112 = CustomConstraint(cons_f112)

    def cons_f113(e, f, d, c, m, b, p, n, a):
        return ZeroQ(a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))

    cons113 = CustomConstraint(cons_f113)
    pattern76 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons112, cons113, cons1)
    rule76 = ReplacementRule(pattern76, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule76)


    def cons_f114(m):
        return Or(And(RationalQ(m), Less(m, S(-1))), SumSimplerQ(m, S(1)))

    cons114 = CustomConstraint(cons_f114)
    pattern77 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons112, cons114)
    rule77 = ReplacementRule(pattern77, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + (a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule77)


    def cons_f115(n, m, p):
        return RationalQ(m, n, p)

    cons115 = CustomConstraint(cons_f115)

    def cons_f116(p):
        return Greater(p, S(0))

    cons116 = CustomConstraint(cons_f116)

    def cons_f117(n, m, p):
        return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), IntegersQ(m, n + p), IntegersQ(p, m + n))

    cons117 = CustomConstraint(cons_f117)
    pattern78 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons115, cons38, cons31, cons116, cons117)
    rule78 = ReplacementRule(pattern78, lambda e, d, c, f, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p/(b*(m + S(1))) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(-1))*Simp(c*f*p + d*e*n + d*f*x*(n + p), x), x)/(b*(m + S(1))))
    rubi.add(rule78)


    def cons_f118(n):
        return Greater(n, S(1))

    cons118 = CustomConstraint(cons_f118)
    pattern79 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons115, cons38, cons118, cons117)
    rule79 = ReplacementRule(pattern79, lambda e, d, c, f, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(1))*(-a*d + b*c)/(b*(m + S(1))*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-2))*(e + f*x)**p*Simp(a*d*(c*f*(p + S(1)) + d*e*(n + S(-1))) + b*c*(-c*f*(m + p + S(2)) + d*e*(m - n + S(2))) + d*x*(a*d*f*(n + p) + b*(-c*f*(m + n + p + S(1)) + d*e*(m + S(1)))), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons115, cons38, cons31, cons117)
    rule80 = ReplacementRule(pattern80, lambda e, d, c, f, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(c*f*(m + p + S(2)) + d*e*n + d*f*x*(m + n + p + S(2)), x), x)/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule80)


    def cons_f119(m):
        return Greater(m, S(1))

    cons119 = CustomConstraint(cons_f119)

    def cons_f120(n, m, p):
        return NonzeroQ(m + n + p + S(1))

    cons120 = CustomConstraint(cons_f120)
    pattern81 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons51, cons119, cons120, cons71)
    rule81 = ReplacementRule(pattern81, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1))))
    rubi.add(rule81)


    def cons_f121(m):
        return Greater(m, S(0))

    cons121 = CustomConstraint(cons_f121)

    def cons_f122(n, m, p):
        return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), Or(IntegersQ(m, n + p), IntegersQ(p, m + n)))

    cons122 = CustomConstraint(cons_f122)
    pattern82 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons115, cons121, cons31, cons120, cons122)
    rule82 = ReplacementRule(pattern82, lambda e, d, c, f, m, b, x, p, n, a : (a + b*x)**m*(c + d*x)**n*(e + f*x)**(p + S(1))/(f*(m + n + p + S(1))) - Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(a*n*(-c*f + d*e) + c*m*(-a*f + b*e) + x*(b*n*(-c*f + d*e) + d*m*(-a*f + b*e)), x), x)/(f*(m + n + p + S(1))))
    rubi.add(rule82)


    def cons_f123(n, m, p):
        return IntegersQ(S(2)*m, S(2)*n, S(2)*p)

    cons123 = CustomConstraint(cons_f123)
    pattern83 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons51, cons119, cons120, cons123)
    rule83 = ReplacementRule(pattern83, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1))))
    rubi.add(rule83)


    def cons_f124(n, p):
        return Or(IntegerQ(n), IntegersQ(S(2)*n, S(2)*p))

    cons124 = CustomConstraint(cons_f124)
    pattern84 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons51, cons38, cons71, cons124)
    rule84 = ReplacementRule(pattern84, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons51, cons38, cons123)
    rule85 = ReplacementRule(pattern85, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule85)


    def cons_f125(n, m):
        return PositiveIntegerQ(m + n + S(1))

    cons125 = CustomConstraint(cons_f125)

    def cons_f126(n, m):
        return Or(And(RationalQ(m), Greater(m, S(0))), And(Not(RationalQ(m)), Or(SumSimplerQ(m, S(-1)), Not(SumSimplerQ(n, S(-1))))))

    cons126 = CustomConstraint(cons_f126)
    pattern86 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons125, cons126)
    rule86 = ReplacementRule(pattern86, lambda e, d, c, f, m, b, x, n, a : b*Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x)/f - (-a*f + b*e)*Int((a + b*x)**(m + S(-1))*(c + d*x)**n/(e + f*x), x)/f)
    rubi.add(rule86)


    def cons_f127(e, d, c, f):
        return PositiveQ(-f/(-c*f + d*e))

    cons127 = CustomConstraint(cons_f127)
    pattern87 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons127)
    rule87 = ReplacementRule(pattern87, lambda e, d, c, f, b, x, a : -S(4)*Subst(Int(x**S(2)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)))
    rubi.add(rule87)


    def cons_f128(e, d, c, f):
        return Not(PositiveQ(-f/(-c*f + d*e)))

    cons128 = CustomConstraint(cons_f128)
    pattern88 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons128)
    rule88 = ReplacementRule(pattern88, lambda e, d, c, f, b, x, a : sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(1)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons127)
    rule89 = ReplacementRule(pattern89, lambda e, d, c, f, b, x, a : -S(4)*Subst(Int(S(1)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)))
    rubi.add(rule89)

    pattern90 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons128)
    rule90 = ReplacementRule(pattern90, lambda e, d, c, f, b, x, a : sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(3)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x))
    rubi.add(rule90)


    def cons_f129(f, e, c, d):
        return NonzeroQ(-c*f + d*e)

    cons129 = CustomConstraint(cons_f129)

    def cons_f130(c):
        return PositiveQ(c)

    cons130 = CustomConstraint(cons_f130)

    def cons_f131(e):
        return PositiveQ(e)

    cons131 = CustomConstraint(cons_f131)

    def cons_f132(b, d):
        return Not(NegativeQ(-b/d))

    cons132 = CustomConstraint(cons_f132)
    pattern91 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons129, cons130, cons131, cons132)
    rule91 = ReplacementRule(pattern91, lambda e, f, c, d, b, x : S(2)*sqrt(e)*EllipticE(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/b)
    rubi.add(rule91)


    def cons_f133(b, d):
        return NegativeQ(-b/d)

    cons133 = CustomConstraint(cons_f133)
    pattern92 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons129, cons130, cons131, cons133)
    rule92 = ReplacementRule(pattern92, lambda e, f, c, d, b, x : sqrt(-b*x)*Int(sqrt(e + f*x)/(sqrt(-b*x)*sqrt(c + d*x)), x)/sqrt(b*x))
    rubi.add(rule92)


    def cons_f134(e, c):
        return Not(And(PositiveQ(c), PositiveQ(e)))

    cons134 = CustomConstraint(cons_f134)
    pattern93 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons129, cons134)
    rule93 = ReplacementRule(pattern93, lambda e, f, c, d, b, x : sqrt(S(1) + d*x/c)*sqrt(e + f*x)*Int(sqrt(S(1) + f*x/e)/(sqrt(b*x)*sqrt(S(1) + d*x/c)), x)/(sqrt(S(1) + f*x/e)*sqrt(c + d*x)))
    rubi.add(rule93)


    def cons_f135(e, b, f, a):
        return PositiveQ(b/(-a*f + b*e))

    cons135 = CustomConstraint(cons_f135)

    def cons_f136(a, d, b, c):
        return Not(NegativeQ(-(-a*d + b*c)/d))

    cons136 = CustomConstraint(cons_f136)

    def cons_f137(e, d, c, f, b, x, a):
        return Not(And(SimplerQ(c + d*x, a + b*x), PositiveQ(-d/(-a*d + b*c)), PositiveQ(d/(-c*f + d*e)), Not(NegativeQ((-a*d + b*c)/b))))

    cons137 = CustomConstraint(cons_f137)
    pattern94 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons66, cons135, cons136, cons137)
    rule94 = ReplacementRule(pattern94, lambda e, f, d, c, b, x, a : S(2)*EllipticE(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*f + b*e)/d, S(2))/b)
    rubi.add(rule94)


    def cons_f138(e, d, c, f, b, a):
        return Not(And(PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e))))

    cons138 = CustomConstraint(cons_f138)
    pattern95 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons138, cons136)
    rule95 = ReplacementRule(pattern95, lambda e, f, d, c, b, x, a : sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*Int(sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))), x)/(sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)))
    rubi.add(rule95)


    def cons_f139(f, b, d):
        return Or(PositiveQ(-b/d), NegativeQ(-b/f))

    cons139 = CustomConstraint(cons_f139)
    pattern96 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons130, cons131, cons139)
    rule96 = ReplacementRule(pattern96, lambda e, f, c, d, b, x : S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)))
    rubi.add(rule96)


    def cons_f140(f, b, d):
        return Or(PosQ(-b/d), NegQ(-b/f))

    cons140 = CustomConstraint(cons_f140)
    pattern97 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons130, cons131, cons140)
    rule97 = ReplacementRule(pattern97, lambda e, f, c, d, b, x : S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)))
    rubi.add(rule97)

    pattern98 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons5, cons9, cons10, cons72, cons73, cons134)
    rule98 = ReplacementRule(pattern98, lambda e, f, c, d, b, x : sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)*Int(S(1)/(sqrt(b*x)*sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)), x)/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule98)


    def cons_f141(e, f, b, x, a):
        return SimplerQ(a + b*x, e + f*x)

    cons141 = CustomConstraint(cons_f141)

    def cons_f142(e, d, c, f, b, a):
        return Or(PositiveQ(-(-a*d + b*c)/d), NegativeQ(-(-a*f + b*e)/f))

    cons142 = CustomConstraint(cons_f142)
    pattern99 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons66, cons135, cons109, cons141, cons142)
    rule99 = ReplacementRule(pattern99, lambda e, f, d, c, b, x, a : S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b)
    rubi.add(rule99)


    def cons_f143(e, d, c, f, b, a):
        return Or(PosQ(-(-a*d + b*c)/d), NegQ(-(-a*f + b*e)/f))

    cons143 = CustomConstraint(cons_f143)
    pattern100 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons66, cons135, cons109, cons141, cons143)
    rule100 = ReplacementRule(pattern100, lambda e, f, d, c, b, x, a : S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b)
    rubi.add(rule100)

    pattern101 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons138, cons109, cons141)
    rule101 = ReplacementRule(pattern101, lambda e, f, d, c, b, x, a : sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*Int(S(1)/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))*sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))), x)/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule101)


    def cons_f144(e, f, c, d, b, a):
        return ZeroQ(-a*d*f - b*c*f + S(2)*b*d*e)

    cons144 = CustomConstraint(cons_f144)
    pattern102 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons144, )
    def With102(e, d, c, f, b, x, a):
        q = Rt(b*(-a*f + b*e)/(-a*d + b*c)**S(2), S(3))
        return -sqrt(S(3))*ArcTan(S(2)*sqrt(S(3))*q*(c + d*x)**(S(2)/3)/(S(3)*(e + f*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*q*(-a*d + b*c)) - log(a + b*x)/(S(2)*q*(-a*d + b*c)) + S(3)*log(q*(c + d*x)**(S(2)/3) - (e + f*x)**(S(1)/3))/(S(4)*q*(-a*d + b*c))
    rule102 = ReplacementRule(pattern102, lambda e, d, c, f, b, x, a : With102(e, d, c, f, b, x, a))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_/((x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons144, cons71, cons38)
    rule103 = ReplacementRule(pattern103, lambda e, d, c, f, m, b, x, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(S(2)/3)*(e + f*x)**(S(2)/3)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + f*Int((a + b*x)**(m + S(1))*(a*d*(S(3)*m + S(1)) - S(3)*b*c*(S(3)*m + S(5)) - S(2)*b*d*x*(S(3)*m + S(7)))/((c + d*x)**(S(1)/3)*(e + f*x)**(S(1)/3)), x)/(S(6)*(m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule103)

    pattern104 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons2, cons13, cons74, cons8, cons70, cons17, cons130)
    rule104 = ReplacementRule(pattern104, lambda f, c, d, m, b, x, p, n, a : Int((f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons2, cons13, cons74, cons8, cons70)
    rule105 = ReplacementRule(pattern105, lambda f, c, d, m, b, x, p, n, a : (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule105)


    def cons_f145(n, m):
        return PositiveIntegerQ(m - n)

    cons145 = CustomConstraint(cons_f145)
    pattern106 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons2, cons13, cons74, cons8, cons145, cons96)
    rule106 = ReplacementRule(pattern106, lambda f, c, d, m, b, x, p, n, a : Int(ExpandIntegrand((f*x)**p*(a + b*x)**n*(c + d*x)**n, (a + b*x)**(m - n), x), x))
    rubi.add(rule106)


    def cons_f146(n, m):
        return Or(PositiveIntegerQ(m), NegativeIntegerQ(m, n))

    cons146 = CustomConstraint(cons_f146)
    pattern107 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons13, cons74, cons146)
    rule107 = ReplacementRule(pattern107, lambda e, d, c, f, m, b, x, p, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule107)


    def cons_f147(n, m, p):
        return NegativeIntegerQ(m + n + p + S(2))

    cons147 = CustomConstraint(cons_f147)

    def cons_f148(n, m, p):
        return Or(SumSimplerQ(m, S(1)), And(Not(And(NonzeroQ(n + S(1)), SumSimplerQ(n, S(1)))), Not(And(NonzeroQ(p + S(1)), SumSimplerQ(p, S(1))))))

    cons148 = CustomConstraint(cons_f148)
    pattern108 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons147, cons1, cons148)
    rule108 = ReplacementRule(pattern108, lambda e, d, c, f, m, b, x, p, n, a : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule108)


    def cons_f149(n):
        return NegativeIntegerQ(n)

    cons149 = CustomConstraint(cons_f149)
    pattern109 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons74, cons110, cons149)
    rule109 = ReplacementRule(pattern109, lambda e, d, c, f, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(e + f*x)**(-m + S(-1))*(-a*d + b*c)**n*(-a*f + b*e)**(-n + S(-1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/(m + S(1)))
    rubi.add(rule109)

    pattern110 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons110, cons36)
    rule110 = ReplacementRule(pattern110, lambda e, d, c, f, m, b, x, p, n, a : ((c + d*x)*(-a*f + b*e)/((e + f*x)*(-a*d + b*c)))**(-n)*(a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule110)


    def cons_f150(e, p):
        return Or(IntegerQ(p), PositiveQ(e))

    cons150 = CustomConstraint(cons_f150)
    pattern111 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons130, cons150)
    rule111 = ReplacementRule(pattern111, lambda e, f, d, c, m, b, x, p, n : c**n*e**p*(b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*x/c, -f*x/e)/(b*(m + S(1))))
    rubi.add(rule111)


    def cons_f151(c, b, d):
        return PositiveQ(-d/(b*c))

    cons151 = CustomConstraint(cons_f151)

    def cons_f152(e, f, c, d, p):
        return Or(IntegerQ(p), PositiveQ(d/(-c*f + d*e)))

    cons152 = CustomConstraint(cons_f152)
    pattern112 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons151, cons152)
    rule112 = ReplacementRule(pattern112, lambda e, f, d, c, m, b, x, p, n : (d/(-c*f + d*e))**(-p)*(-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*AppellF1(n + S(1), -m, -p, n + S(2), S(1) + d*x/c, -f*(c + d*x)/(-c*f + d*e))/(d*(n + S(1))))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_), cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons63)
    rule113 = ReplacementRule(pattern113, lambda e, f, d, c, m, b, x, p, n : c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n*(e + f*x)**p, x))
    rubi.add(rule113)


    def cons_f153(c, d, b, x, a):
        return Not(And(PositiveQ(d/(a*d - b*c)), SimplerQ(c + d*x, a + b*x)))

    cons153 = CustomConstraint(cons_f153)
    pattern114 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons60, cons36, cons97, cons66, cons153)
    rule114 = ReplacementRule(pattern114, lambda e, f, d, c, m, b, x, p, n, a : b**(-p + S(-1))*(b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*(-a*f + b*e)**p*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(m + S(1)))
    rubi.add(rule114)


    def cons_f154(a, d, b, c):
        return Not(PositiveQ(b/(-a*d + b*c)))

    cons154 = CustomConstraint(cons_f154)

    def cons_f155(d, c, b, x, a):
        return Not(SimplerQ(c + d*x, a + b*x))

    cons155 = CustomConstraint(cons_f155)
    pattern115 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons60, cons36, cons97, cons154, cons155)
    rule115 = ReplacementRule(pattern115, lambda e, f, d, c, m, b, x, p, n, a : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule115)


    def cons_f156(e, c, d, f, b, x, a):
        return Not(And(PositiveQ(d/(a*d - b*c)), PositiveQ(d/(-c*f + d*e)), SimplerQ(c + d*x, a + b*x)))

    cons156 = CustomConstraint(cons_f156)

    def cons_f157(e, c, f, d, b, x, a):
        return Not(And(PositiveQ(f/(a*f - b*e)), PositiveQ(f/(c*f - d*e)), SimplerQ(e + f*x, a + b*x)))

    cons157 = CustomConstraint(cons_f157)
    pattern116 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons100, cons66, cons135, cons156, cons157)
    rule116 = ReplacementRule(pattern116, lambda e, f, d, c, m, b, x, p, n, a : (b/(-a*d + b*c))**(-n)*(b/(-a*f + b*e))**(-p)*(a + b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(b*(m + S(1))))
    rubi.add(rule116)


    def cons_f158(e, b, f, a):
        return Not(PositiveQ(b/(-a*f + b*e)))

    cons158 = CustomConstraint(cons_f158)
    pattern117 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons100, cons66, cons158)
    rule117 = ReplacementRule(pattern117, lambda e, f, d, c, m, b, x, p, n, a : (b/(-a*f + b*e))**(-IntPart(p))*(b*(e + f*x)/(-a*f + b*e))**(-FracPart(p))*(e + f*x)**FracPart(p)*Int((a + b*x)**m*(c + d*x)**n*(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))**p, x))
    rubi.add(rule117)


    def cons_f159(e, f, b, x, a):
        return Not(SimplerQ(e + f*x, a + b*x))

    cons159 = CustomConstraint(cons_f159)
    pattern118 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons60, cons36, cons100, cons154, cons155, cons159)
    rule118 = ReplacementRule(pattern118, lambda e, f, d, c, m, b, x, p, n, a : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule118)

    pattern119 = Pattern(Integral((e_ + u_*WC('f', S(1)))**WC('p', S(1))*(u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons2, cons13, cons74, cons6, cons7)
    rule119 = ReplacementRule(pattern119, lambda e, d, c, f, m, u, b, x, p, n, a : Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule119)


    def cons_f160(n, m):
        return Or(PositiveIntegerQ(m), IntegersQ(m, n))

    cons160 = CustomConstraint(cons_f160)

    def cons_f161(g, x):
        return FreeQ(g, x)

    cons161 = CustomConstraint(cons_f161)

    def cons_f162(h, x):
        return FreeQ(h, x)

    cons162 = CustomConstraint(cons_f162)
    pattern120 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons160)
    rule120 = ReplacementRule(pattern120, lambda e, h, d, c, f, g, m, b, x, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)*(g + h*x), x), x))
    rubi.add(rule120)


    def cons_f163(n, m):
        return Not(And(SumSimplerQ(n, S(1)), Not(SumSimplerQ(m, S(1)))))

    cons163 = CustomConstraint(cons_f163)
    pattern121 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons12, cons1, cons163)
    rule121 = ReplacementRule(pattern121, lambda e, h, d, c, f, g, m, b, x, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(2)*d*f*h*m - a*b*(-c*f*h*(m + S(1)) + d*(e*h + f*g)) + b**S(2)*d*e*g + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)) + (a*d*f*h*m + b*(-c*f*h*(m + S(2)) + d*(e*h + f*g)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons37, cons38, cons32)
    rule122 = ReplacementRule(pattern122, lambda e, h, f, d, c, g, m, b, x, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*c*d*f*h*(n + S(1)) + a*b*(c**S(2)*f*h*(m + S(1)) - c*d*(e*h + f*g)*(m + n + S(2)) + d**S(2)*e*g*(m + S(1))) + b**S(2)*c*d*e*g*(n + S(1)) + x*(a**S(2)*d**S(2)*f*h*(n + S(1)) - a*b*d**S(2)*(n + S(1))*(e*h + f*g) + b**S(2)*(c**S(2)*f*h*(m + S(1)) - c*d*(m + S(1))*(e*h + f*g) + d**S(2)*e*g*(m + n + S(2)))))/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)) - (a**S(2)*d**S(2)*f*h*(n**S(2) + S(3)*n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m**S(2) + S(3)*m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m**S(2) + m*(S(2)*n + S(5)) + n**S(2) + S(5)*n + S(6))))*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1)), x)/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)))
    rubi.add(rule122)


    def cons_f164(n, m):
        return Or(And(RationalQ(m), Less(m, S(-2))), And(ZeroQ(m + n + S(3)), Not(And(RationalQ(n), Less(n, S(-2))))))

    cons164 = CustomConstraint(cons_f164)
    pattern123 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons164)
    rule123 = ReplacementRule(pattern123, lambda e, h, d, c, f, g, m, b, x, n, a : (-d*(m + n + S(3))*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)) + f*h/b**S(2))*Int((a + b*x)**(m + S(2))*(c + d*x)**n, x) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(3)*d*f*h*(n + S(2)) - a**S(2)*b*(c*f*h*m - d*(e*h + f*g)*(m + n + S(3))) - a*b**S(2)*(c*(e*h + f*g) + d*e*g*(S(2)*m + n + S(4))) + b**S(3)*c*e*g*(m + S(2)) + b*x*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2)))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)))
    rubi.add(rule123)


    def cons_f165(m):
        return Or(And(RationalQ(m), Inequality(S(-2), LessEqual, m, Less, S(-1))), SumSimplerQ(m, S(1)))

    cons165 = CustomConstraint(cons_f165)

    def cons_f166(n, m):
        return NonzeroQ(m + n + S(3))

    cons166 = CustomConstraint(cons_f166)
    pattern124 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons165, cons1, cons166)
    rule124 = ReplacementRule(pattern124, lambda e, h, d, c, f, g, m, b, x, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*d*f*h*(n + S(2)) + a*b*(c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*d*e*g*(m + n + S(3)) + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))) - (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))))
    rubi.add(rule124)


    def cons_f167(n, m):
        return NonzeroQ(m + n + S(2))

    cons167 = CustomConstraint(cons_f167)
    pattern125 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons167, cons166)
    rule125 = ReplacementRule(pattern125, lambda e, h, d, c, f, g, m, b, x, n, a : -(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a*d*f*h*(n + S(2)) + b*c*f*h*(m + S(2)) - b*d*f*h*x*(m + n + S(2)) - b*d*(e*h + f*g)*(m + n + S(3)))/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))) + (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**m*(c + d*x)**n, x)/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))))
    rubi.add(rule125)


    def cons_f168(n, m, p):
        return Or(IntegersQ(m, n, p), PositiveIntegerQ(n, p))

    cons168 = CustomConstraint(cons_f168)
    pattern126 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons168)
    rule126 = ReplacementRule(pattern126, lambda e, h, d, c, f, g, m, b, x, p, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x), x), x))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons37, cons38, cons31, cons71)
    rule127 = ReplacementRule(pattern127, lambda e, h, d, c, f, g, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule127)

    pattern128 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons37, cons38, cons31, cons123)
    rule128 = ReplacementRule(pattern128, lambda e, h, d, c, f, g, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule128)

    pattern129 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons51, cons38, cons71)
    rule129 = ReplacementRule(pattern129, lambda e, h, d, c, f, g, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons51, cons38, cons123)
    rule130 = ReplacementRule(pattern130, lambda e, h, d, c, f, g, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule130)

    pattern131 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons51, cons121, cons96, cons71)
    rule131 = ReplacementRule(pattern131, lambda e, h, d, c, f, g, m, b, x, p, n, a : h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2))))
    rubi.add(rule131)

    pattern132 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons51, cons121, cons96, cons123)
    rule132 = ReplacementRule(pattern132, lambda e, h, d, c, f, g, m, b, x, p, n, a : h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2))))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons147, cons1, cons148)
    rule133 = ReplacementRule(pattern133, lambda e, h, d, c, f, g, m, b, x, p, n, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule133)


    def cons_f169(e, h, d, c, f, g, b, x, a):
        return FreeQ(List(a, b, c, d, e, f, g, h), x)

    cons169 = CustomConstraint(cons_f169)
    pattern134 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule134 = ReplacementRule(pattern134, lambda e, h, d, c, f, g, b, x, p, a : (-a*h + b*g)*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - (-c*h + d*g)*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule134)


    def cons_f170(e, h, d, c, f, g, b, x, p, n, a):
        return FreeQ(List(a, b, c, d, e, f, g, h, n, p), x)

    cons170 = CustomConstraint(cons_f170)
    pattern135 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/(x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons13, cons74, cons170)
    rule135 = ReplacementRule(pattern135, lambda e, h, d, c, f, g, b, x, p, n, a : h*Int((c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((c + d*x)**n*(e + f*x)**p/(a + b*x), x)/b)
    rubi.add(rule135)


    def cons_f171(e, f, d, c, x):
        return SimplerQ(c + d*x, e + f*x)

    cons171 = CustomConstraint(cons_f171)
    pattern136 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons141, cons171)
    rule136 = ReplacementRule(pattern136, lambda e, h, f, d, g, c, b, x, a : h*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)), x)/f + (-e*h + f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)), x)/f)
    rubi.add(rule136)


    def cons_f172(n, m, p):
        return Or(SumSimplerQ(m, S(1)), And(Not(SumSimplerQ(n, S(1))), Not(SumSimplerQ(p, S(1)))))

    cons172 = CustomConstraint(cons_f172)
    pattern137 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons74, cons172)
    rule137 = ReplacementRule(pattern137, lambda e, h, d, c, f, g, m, b, x, p, n, a : h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x)/b)
    rubi.add(rule137)


    def cons_f173(q, x):
        return FreeQ(q, x)

    cons173 = CustomConstraint(cons_f173)
    pattern138 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons173, cons87, cons98)
    rule138 = ReplacementRule(pattern138, lambda e, h, d, c, f, g, q, b, x, p, a : (-a*f + b*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule138)

    pattern139 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule139 = ReplacementRule(pattern139, lambda e, h, d, c, f, g, b, x, a : -S(2)*sqrt(d*(e + f*x)/(-c*f + d*e))*sqrt(d*(g + h*x)/(-c*h + d*g))*EllipticPi(-b*(-c*f + d*e)/(f*(-a*d + b*c)), asin(sqrt(-f/(-c*f + d*e))*sqrt(c + d*x)), h*(-c*f + d*e)/(f*(-c*h + d*g)))/(sqrt(-f/(-c*f + d*e))*sqrt(e + f*x)*sqrt(g + h*x)*(-a*d + b*c)))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons22)
    rule140 = ReplacementRule(pattern140, lambda e, h, d, c, f, g, b, x, n, a : Int(ExpandIntegrand(S(1)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), (c + d*x)**(n + S(1)/2)/(a + b*x), x), x))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule141 = ReplacementRule(pattern141, lambda e, h, d, c, f, g, b, x, a : (-a*f + b*e)*(-a*h + b*g)*Int(S(1)/((a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2) + Int((-a*f*h + b*e*h + b*f*g + b*f*h*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule142 = ReplacementRule(pattern142, lambda e, h, d, c, f, g, b, x, a : -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule143 = ReplacementRule(pattern143, lambda e, h, d, c, f, g, b, x, a : -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*(-c*h + d*g)*Subst(Int(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))/sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1)), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)**S(2)))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule144 = ReplacementRule(pattern144, lambda e, h, d, c, f, g, b, x, a : S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/((-b*x**S(2) + h)*sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule145 = ReplacementRule(pattern145, lambda e, h, d, c, f, g, b, x, a : b*Int(sqrt(c + d*x)/((a + b*x)**(S(3)/2)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c) - d*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule146 = ReplacementRule(pattern146, lambda e, h, d, c, f, g, b, x, a : sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)/(h*sqrt(e + f*x)) - (-c*f + d*e)*(-e*h + f*g)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*(e + f*x)**(S(3)/2)*sqrt(g + h*x)), x)/(S(2)*f*h) + (-c*f + d*e)*(-S(2)*a*f*h + b*e*h + b*f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h) + (a*d*f*h - b*(-c*f*h + d*e*h + d*f*g))*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons169)
    rule147 = ReplacementRule(pattern147, lambda e, h, d, c, f, g, b, x, a : b*Int(sqrt(a + b*x)*sqrt(c + d*x)/(sqrt(e + f*x)*sqrt(g + h*x)), x)/d - (-a*d + b*c)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/d)
    rubi.add(rule147)


    def cons_f174(q, p):
        return IntegersQ(p, q)

    cons174 = CustomConstraint(cons_f174)
    pattern148 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons174)
    rule148 = ReplacementRule(pattern148, lambda e, h, d, c, f, g, q, m, b, x, p, n, a : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x))
    rubi.add(rule148)


    def cons_f175(q):
        return PositiveIntegerQ(q)

    cons175 = CustomConstraint(cons_f175)
    pattern149 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons74, cons175, cons172)
    rule149 = ReplacementRule(pattern149, lambda e, h, d, c, f, g, q, m, b, x, p, n, a : h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b)
    rubi.add(rule149)


    def cons_f176(e, h, d, c, f, g, q, m, b, x, p, n, a):
        return FreeQ(List(a, b, c, d, e, f, g, h, m, n, p, q), x)

    cons176 = CustomConstraint(cons_f176)
    pattern150 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons74, cons173, cons176)
    rule150 = ReplacementRule(pattern150, lambda e, h, d, c, f, g, q, m, b, x, p, n, a : Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(u_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(u_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons13, cons74, cons173, cons6, cons7)
    rule151 = ReplacementRule(pattern151, lambda e, h, d, c, f, g, q, m, u, b, x, p, n, a : Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule151)


    def cons_f177(e, h, d, c, f, g, q, m, b, r, i, x, p, n, a):
        return FreeQ(List(a, b, c, d, e, f, g, h, i, m, n, p, q, r), x)

    cons177 = CustomConstraint(cons_f177)

    def cons_f178(i, x):
        return FreeQ(i, x)

    cons178 = CustomConstraint(cons_f178)

    def cons_f179(r, x):
        return FreeQ(r, x)

    cons179 = CustomConstraint(cons_f179)
    pattern152 = Pattern(Integral(((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_*WC('i', S(1)))**r_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons178, cons2, cons13, cons74, cons173, cons179, cons177)
    rule152 = ReplacementRule(pattern152, lambda e, h, d, c, f, g, q, m, b, r, i, x, p, n, a : (i*(a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q)**r*(a + b*x)**(-m*r)*(c + d*x)**(-n*r)*(e + f*x)**(-p*r)*(g + h*x)**(-q*r)*Int((a + b*x)**(m*r)*(c + d*x)**(n*r)*(e + f*x)**(p*r)*(g + h*x)**(q*r), x))
    rubi.add(rule152)

    return rubi
