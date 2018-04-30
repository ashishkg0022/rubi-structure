
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

    from sympy.integrals.rubi.constraints import (constraint_freeq_a, constraint_freeq_b, constraint_freeq_c, constraint_freeq_d, constraint_freeq_e, constraint_freeq_f,
        constraint_freeq_g, constraint_freeq_h, constraint_freeq_m, constraint_freeq_n, constraint_freeq_p)

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False
    import functools, operator

def linear_products(rubi):
    def f1(x):
        return functools.reduce(operator.and_, [])
    pattern1 = Pattern(Integral(S(1)/x_, x_),CustomConstraint(f1))
    rule1 = ReplacementRule(pattern1, lambda x : log(x))
    rubi.add(rule1)

    def f2(m, x):
        return functools.reduce(operator.and_, [ NonzeroQ(m + 1), FreeQ(m, x)])
    pattern2 = Pattern(Integral(x_**WC('m', S(1)), x_),CustomConstraint(f2))
    rule2 = ReplacementRule(pattern2, lambda m, x : x**(m + S(1))/(m + S(1)))
    rubi.add(rule2)

    def f3(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern3 = Pattern(Integral(S(1)/(a_ + x_*WC('b', S(1))), x_),CustomConstraint(f3))
    rule3 = ReplacementRule(pattern3, lambda b, a, x : log(RemoveContent(a + b*x, x))/b)
    rubi.add(rule3)

    def f4(m, b, a, x):
        return functools.reduce(operator.and_, [ NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x)])
    pattern4 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_),CustomConstraint(f4))
    rule4 = ReplacementRule(pattern4, lambda m, b, a, x : (a + b*x)**(m + S(1))/(b*(m + S(1))))
    rubi.add(rule4)

    def f5(u, b, x, m, a):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x)])
    pattern5 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**m_, x_),CustomConstraint(f5))
    rule5 = ReplacementRule(pattern5, lambda u, b, x, m, a : Subst(Int((a + b*x)**m, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule5)

    def f6(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern6 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))*(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f6))
    rule6 = ReplacementRule(pattern6, lambda b, a, x, d, c : Int(S(1)/(a*c + b*d*x**S(2)), x))
    rubi.add(rule6)

    def f7(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern7 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f7))
    rule7 = ReplacementRule(pattern7, lambda b, a, x, d, c : b*Int(S(1)/(a + b*x), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule7)

    def f8(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(m + n + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern8 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1)), x_),CustomConstraint(f8))
    rule8 = ReplacementRule(pattern8, lambda c, b, n, x, d, m, a : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule8)

    def f9(c, b, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PositiveIntegerQ(m + 1/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern9 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_),CustomConstraint(f9))
    rule9 = ReplacementRule(pattern9, lambda c, b, x, d, m, a : S(2)*a*c*m*Int((a + b*x)**(m + S(-1))*(c + d*x)**(m + S(-1)), x)/(S(2)*m + S(1)) + x*(a + b*x)**m*(c + d*x)**m/(S(2)*m + S(1)))
    rubi.add(rule9)

    def f10(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern10 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(3)/2)*(c_ + x_*WC('d', S(1)))**(S(3)/2)), x_),CustomConstraint(f10))
    rule10 = ReplacementRule(pattern10, lambda b, a, x, d, c : x/(a*c*sqrt(a + b*x)*sqrt(c + d*x)))
    rubi.add(rule10)

    def f11(c, b, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), NegativeIntegerQ(m + 3/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern11 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_),CustomConstraint(f11))
    rule11 = ReplacementRule(pattern11, lambda c, b, x, d, m, a : -x*(a + b*x)**(m + S(1))*(c + d*x)**(m + S(1))/(S(2)*a*c*(m + S(1))) + (S(2)*m + S(3))*Int((a + b*x)**(m + S(1))*(c + d*x)**(m + S(1)), x)/(S(2)*a*c*(m + S(1))))
    rubi.add(rule11)

    def f12(c, b, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x)])
    pattern12 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('m', S(1)), x_),CustomConstraint(f12))
    rule12 = ReplacementRule(pattern12, lambda c, b, x, d, m, a : Int((a*c + b*d*x**S(2))**m, x))
    rubi.add(rule12)

    def f13(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PositiveQ(a), ZeroQ(a + c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern13 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f13))
    rule13 = ReplacementRule(pattern13, lambda b, a, x, d, c : acosh(b*x/a)/b)
    rubi.add(rule13)

    def f14(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern14 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f14))
    rule14 = ReplacementRule(pattern14, lambda b, a, x, d, c : S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)))
    rubi.add(rule14)

    def f15(c, b, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x)])
    pattern15 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**m_, x_),CustomConstraint(f15))
    rule15 = ReplacementRule(pattern15, lambda c, b, x, d, m, a : (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((a*c + b*d*x**S(2))**m, x))
    rubi.add(rule15)

    def f16(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PosQ(b*d/(a*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern16 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(5)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_),CustomConstraint(f16))
    rule16 = ReplacementRule(pattern16, lambda b, a, x, d, c : (-a*d + b*c)*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(2)*b) - S(2)/(b*(a + b*x)**(S(1)/4)*(c + d*x)**(S(1)/4)))
    rubi.add(rule16)

    def f17(b, a, x, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PosQ(b*d/(a*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern17 = Pattern(Integral(S(1)/((a_ + x_*WC('b', S(1)))**(S(9)/4)*(c_ + x_*WC('d', S(1)))**(S(1)/4)), x_),CustomConstraint(f17))
    rule17 = ReplacementRule(pattern17, lambda b, a, x, d, c : -d*Int(S(1)/((a + b*x)**(S(5)/4)*(c + d*x)**(S(5)/4)), x)/(S(5)*b) - S(4)/(S(5)*b*(a + b*x)**(S(5)/4)*(c + d*x)**(S(1)/4)))
    rubi.add(rule17)

    def f18(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), IntegerQ(m + 1/2), IntegerQ(n + 1/2), Less(0, m, n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern18 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f18))
    rule18 = ReplacementRule(pattern18, lambda c, b, n, x, d, m, a : S(2)*c*n*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(m + n + S(1)) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))))
    rubi.add(rule18)

    def f19(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), IntegerQ(m + 1/2), IntegerQ(n + 1/2), Less(m, n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern19 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f19))
    rule19 = ReplacementRule(pattern19, lambda c, b, n, x, d, m, a : (m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(S(2)*a*(m + S(1))) - (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/(S(2)*a*d*(m + S(1))))
    rubi.add(rule19)

    def f20(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern20 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f20))
    rule20 = ReplacementRule(pattern20, lambda b, n, a, x, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x))
    rubi.add(rule20)

    def f21(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegativeIntegerQ(m), IntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern21 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f21))
    rule21 = ReplacementRule(pattern21, lambda b, n, a, x, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x))
    rubi.add(rule21)

    def f22(b, n, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(n), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern22 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f22))
    rule22 = ReplacementRule(pattern22, lambda b, n, a, x, d, c : (-a*d + b*c)*Int((c + d*x)**(n + S(-1))/(a + b*x), x)/b + (c + d*x)**n/(b*n))
    rubi.add(rule22)

    def f23(b, n, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(n), Less(n, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern23 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f23))
    rule23 = ReplacementRule(pattern23, lambda b, n, a, x, d, c : b*Int((c + d*x)**(n + S(1))/(a + b*x), x)/(-a*d + b*c) - (c + d*x)**(n + S(1))/((n + S(1))*(-a*d + b*c)))
    rubi.add(rule23)

    def f24(b, a, x, d, c):
        return functools.reduce(operator.and_, [ PosQ((-a*d + b*c)/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern24 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_),CustomConstraint(f24), )
    def With24(b, a, x, d, c):
        q = Rt((-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    rule24 = ReplacementRule(pattern24, lambda b, a, x, d, c : With24(b, a, x, d, c))
    rubi.add(rule24)

    def f25(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NegQ((-a*d + b*c)/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern25 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)), x_),CustomConstraint(f25), )
    def With25(b, a, x, d, c):
        q = Rt(-(-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b) - S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + log(RemoveContent(a + b*x, x))/(S(2)*b*q)
    rule25 = ReplacementRule(pattern25, lambda b, a, x, d, c : With25(b, a, x, d, c))
    rubi.add(rule25)

    def f26(b, a, x, d, c):
        return functools.reduce(operator.and_, [ PosQ((-a*d + b*c)/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern26 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_),CustomConstraint(f26), )
    def With26(b, a, x, d, c):
        q = Rt((-a*d + b*c)/b, S(3))
        return -S(3)*Subst(Int(S(1)/(q**S(2) + q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) - S(3)*Subst(Int(S(1)/(q - x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    rule26 = ReplacementRule(pattern26, lambda b, a, x, d, c : With26(b, a, x, d, c))
    rubi.add(rule26)

    def f27(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NegQ((-a*d + b*c)/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern27 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_),CustomConstraint(f27), )
    def With27(b, a, x, d, c):
        q = Rt(-(-a*d + b*c)/b, S(3))
        return S(3)*Subst(Int(S(1)/(q**S(2) - q*x + x**S(2)), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q) + S(3)*Subst(Int(S(1)/(q + x), x), x, (c + d*x)**(S(1)/3))/(S(2)*b*q**S(2)) - log(RemoveContent(a + b*x, x))/(S(2)*b*q**S(2))
    rule27 = ReplacementRule(pattern27, lambda b, a, x, d, c : With27(b, a, x, d, c))
    rubi.add(rule27)

    def f28(b, n, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(n), Less(-1, n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern28 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f28), )
    def With28(b, n, a, x, d, c):
        p = Denominator(n)
        return p*Subst(Int(x**(p*(n + S(1)) + S(-1))/(a*d - b*c + b*x**p), x), x, (c + d*x)**(S(1)/p))
    rule28 = ReplacementRule(pattern28, lambda b, n, a, x, d, c : With28(b, n, a, x, d, c))
    rubi.add(rule28)

    def f29(d, n, c, x):
        return functools.reduce(operator.and_, [ FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern29 = Pattern(Integral((c_ + x_*WC('d', S(1)))**n_/x_, x_),CustomConstraint(f29))
    rule29 = ReplacementRule(pattern29, lambda d, n, c, x : -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), S(1) + d*x/c)/(c*(n + S(1))))
    rubi.add(rule29)

    def f30(c, b, n, x, d, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern30 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/(a_ + x_*WC('b', S(1))), x_),CustomConstraint(f30))
    rule30 = ReplacementRule(pattern30, lambda c, b, n, x, d, a : -(c + d*x)**(n + S(1))*Hypergeometric2F1(S(1), n + S(1), n + S(2), TogetherSimplify(b*(c + d*x)/(-a*d + b*c)))/((n + S(1))*(-a*d + b*c)))
    rubi.add(rule30)

    def f31(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m, n), Less(m, -1), Greater(n, 0), IntLinearcQ(a, b, c, d, m, n, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern31 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f31))
    rule31 = ReplacementRule(pattern31, lambda b, n, a, x, d, m, c : -d*n*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1)), x)/(b*(m + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + S(1))))
    rubi.add(rule31)

    def f32(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m, n), Less(m, -1), IntLinearcQ(a, b, c, d, m, n, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern32 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f32))
    rule32 = ReplacementRule(pattern32, lambda b, n, a, x, d, m, c : -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule32)

    def f33(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m, n), Greater(n, 0), Unequal(m + n + 1, 0), IntLinearcQ(a, b, c, d, m, n, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern33 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f33))
    rule33 = ReplacementRule(pattern33, lambda b, n, a, x, d, m, c : n*(-a*d + b*c)*Int((a + b*x)**m*(c + d*x)**(n + S(-1)), x)/(b*(m + n + S(1))) + (a + b*x)**(m + S(1))*(c + d*x)**n/(b*(m + n + S(1))))
    rubi.add(rule33)

    def f34(c, b, x, d, a):
        return functools.reduce(operator.and_, [ ZeroQ(b + d), PositiveQ(a + c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern34 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f34))
    rule34 = ReplacementRule(pattern34, lambda c, b, x, d, a : Int(S(1)/sqrt(a*c - b**S(2)*x**S(2) - b*x*(a - c)), x))
    rubi.add(rule34)

    def f35(b, a, x, d, c):
        return functools.reduce(operator.and_, [ PositiveQ(-a*d + b*c), PositiveQ(b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern35 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f35))
    rule35 = ReplacementRule(pattern35, lambda b, a, x, d, c : S(2)*Subst(Int(S(1)/sqrt(-a*d + b*c + d*x**S(2)), x), x, sqrt(a + b*x))/sqrt(b))
    rubi.add(rule35)

    def f36(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(b - d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern36 = Pattern(Integral(S(1)/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f36))
    rule36 = ReplacementRule(pattern36, lambda b, a, x, d, c : S(2)*Subst(Int(S(1)/sqrt(-a + c + x**S(2)), x), x, sqrt(a + b*x))/b)
    rubi.add(rule36)

    def f37(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern37 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f37))
    rule37 = ReplacementRule(pattern37, lambda b, a, x, d, c : S(2)*Subst(Int(S(1)/(b - d*x**S(2)), x), x, sqrt(a + b*x)/sqrt(c + d*x)))
    rubi.add(rule37)

    def f38(c, b, x, d, m, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m), Less(-1, m, 0), LessEqual(3, Denominator(m), 4), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern38 = Pattern(Integral((c_ + x_*WC('d', S(1)))**m_*(x_*WC('b', S(1)) + WC('a', S(0)))**m_, x_),CustomConstraint(f38))
    rule38 = ReplacementRule(pattern38, lambda c, b, x, d, m, a : (a + b*x)**m*(c + d*x)**m*(a*c + b*d*x**S(2) + x*(a*d + b*c))**(-m)*Int((a*c + b*d*x**S(2) + x*(a*d + b*c))**m, x))
    rubi.add(rule38)

    def f39(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PosQ(d/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern39 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_),CustomConstraint(f39), )
    def With39(b, a, x, d, c):
        q = Rt(d/b, S(3))
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d - q*log(c + d*x)/(S(2)*d) - S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(-1))/(S(2)*d)
    rule39 = ReplacementRule(pattern39, lambda b, a, x, d, c : With39(b, a, x, d, c))
    rubi.add(rule39)

    def f40(b, a, x, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegQ(d/b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern40 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)), x_),CustomConstraint(f40), )
    def With40(b, a, x, d, c):
        q = Rt(-d/b, S(3))
        return sqrt(S(3))*q*ArcTan(-S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/d + q*log(c + d*x)/(S(2)*d) + S(3)*q*log(q*(a + b*x)**(S(1)/3)/(c + d*x)**(S(1)/3) + S(1))/(S(2)*d)
    rule40 = ReplacementRule(pattern40, lambda b, a, x, d, c : With40(b, a, x, d, c))
    rubi.add(rule40)

    def f41(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m, n), Less(-1, m, 0), Equal(m + n + 1, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern41 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f41), )
    def With41(b, n, a, x, d, m, c):
        p = Denominator(m)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))/(b - d*x**p), x), x, (a + b*x)**(S(1)/p)*(c + d*x)**(-S(1)/p))
    rule41 = ReplacementRule(pattern41, lambda b, n, a, x, d, m, c : With41(b, n, a, x, d, m, c))
    rubi.add(rule41)

    def f42(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(m, n), Less(-1, m, 0), Less(-1, n, 0), LessEqual(Denominator(n), Denominator(m)), IntLinearcQ(a, b, c, d, m, n, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern42 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f42), )
    def With42(b, n, a, x, d, m, c):
        p = Denominator(m)
        return p*Subst(Int(x**(p*(m + S(1)) + S(-1))*(-a*d/b + c + d*x**p/b)**n, x), x, (a + b*x)**(S(1)/p))/b
    rule42 = ReplacementRule(pattern42, lambda b, n, a, x, d, m, c : With42(b, n, a, x, d, m, c))
    rubi.add(rule42)

    def f43(b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegativeIntegerQ(m + n + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern43 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_, x_),CustomConstraint(f43))
    rule43 = ReplacementRule(pattern43, lambda b, n, a, x, d, m, c : -d*(m + n + S(2))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/((m + S(1))*(-a*d + b*c)) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))/((m + S(1))*(-a*d + b*c)))
    rubi.add(rule43)

    def f44(b, n, x, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern44 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f44))
    rule44 = ReplacementRule(pattern44, lambda b, n, x, d, m, c : c**n*(b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*x/c)/(b*(m + S(1))))
    rubi.add(rule44)

    def f45(b, n, x, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern45 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f45))
    rule45 = ReplacementRule(pattern45, lambda b, n, x, d, m, c : (-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*Hypergeometric2F1(-m, n + S(1), n + S(2), S(1) + d*x/c)/(d*(n + S(1))))
    rubi.add(rule45)

    def f46(b, n, x, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern46 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f46))
    rule46 = ReplacementRule(pattern46, lambda b, n, x, d, m, c : c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n, x))
    rubi.add(rule46)

    def f47(b, n, x, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern47 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f47))
    rule47 = ReplacementRule(pattern47, lambda b, n, x, d, m, c : (b*x)**FracPart(m)*(-b*c/d)**IntPart(m)*(-d*x/c)**(-FracPart(m))*Int((-d*x/c)**m*(c + d*x)**n, x))
    rubi.add(rule47)

    def f48(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x)])
    pattern48 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f48))
    rule48 = ReplacementRule(pattern48, lambda c, b, n, x, d, m, a : b**(-n + S(-1))*(a + b*x)**(m + S(1))*(-a*d + b*c)**n*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(m + S(1)))
    rubi.add(rule48)

    def f49(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveQ(b/(-a*d + b*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern49 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f49))
    rule49 = ReplacementRule(pattern49, lambda c, b, n, x, d, m, a : (b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*Hypergeometric2F1(-n, m + S(1), m + S(2), -d*(a + b*x)/(-a*d + b*c))/(b*(m + S(1))))
    rubi.add(rule49)

    def f50(c, b, n, x, d, m, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern50 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_, x_),CustomConstraint(f50))
    rule50 = ReplacementRule(pattern50, lambda c, b, n, x, d, m, a : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule50)

    def f51(u, b, n, a, x, d, m, c):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(Coefficient(u, x, 0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern51 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f51))
    rule51 = ReplacementRule(pattern51, lambda u, b, n, a, x, d, m, c : Subst(Int((a + b*x)**m*(c + d*x)**n, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule51)

    def f52(c, b, n, p, f, x, e, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), ZeroQ(m - n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern52 = Pattern(Integral((a_ + x_*WC('b', S(1)))**WC('m', S(1))*(c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f52))
    rule52 = ReplacementRule(pattern52, lambda c, b, n, p, f, x, e, d, m, a : Int((e + f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule52)

    def f53(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 2), ZeroQ(a*d*f*(n + p + 2) - b*(c*f*(p + 1) + d*e*(n + 1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern53 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f53))
    rule53 = ReplacementRule(pattern53, lambda b, n, a, p, f, x, e, d, c : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))))
    rubi.add(rule53)

    def f54(b, n, p, f, x, e, d, a):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p), ZeroQ(a*f + b*e), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern54 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_),CustomConstraint(f54))
    rule54 = ReplacementRule(pattern54, lambda b, n, p, f, x, e, d, a : Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x))
    rubi.add(rule54)

    def f55(b, n, p, f, x, e, d, a):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p), NonzeroQ(a*f + b*e), FreeQ(a, x), FreeQ(b, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern55 = Pattern(Integral((x_*WC('d', S(1)))**WC('n', S(1))*(a_ + x_*WC('b', S(1)))*(e_ + x_*WC('f', S(1)))**WC('p', S(1)), x_),CustomConstraint(f55))
    rule55 = ReplacementRule(pattern55, lambda b, n, p, f, x, e, d, a : Int(ExpandIntegrand((d*x)**n*(a + b*x)*(e + f*x)**p, x), x))
    rubi.add(rule55)

    def f56(c, b, n, p, f, x, e, d, a):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern56 = Pattern(Integral((c_ + x_*WC('d', S(1)))**WC('n', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f56))
    rule56 = ReplacementRule(pattern56, lambda c, b, n, p, f, x, e, d, a : Int(ExpandIntegrand((a + b*x)*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule56)

    def f57(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(n + p + 2), NonzeroQ(p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern57 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f57))
    rule57 = ReplacementRule(pattern57, lambda b, n, a, p, f, x, e, d, c : b*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/f - (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule57)

    def f58(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern58 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f58))
    rule58 = ReplacementRule(pattern58, lambda b, n, a, p, f, x, e, d, c : -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule58)

    def f59(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 2), SumSimplerQ(p, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern59 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f59))
    rule59 = ReplacementRule(pattern59, lambda b, n, a, p, f, x, e, d, c : -(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*f + b*e)/(f*(p + S(1))*(c*f - d*e)) - (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**(p + S(1)), x)/(f*(p + S(1))*(c*f - d*e)))
    rubi.add(rule59)

    def f60(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern60 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f60))
    rule60 = ReplacementRule(pattern60, lambda b, n, a, p, f, x, e, d, c : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(2))) + (a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))*Int((c + d*x)**n*(e + f*x)**p, x)/(d*f*(n + p + S(2))))
    rubi.add(rule60)

    def f61(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 2), NonzeroQ(n + p + 3), ZeroQ(-b*(c*f*(p + 1) + d*e*(n + 1))*(a*d*f*(n + p + 4) - b*(c*f*(p + 2) + d*e*(n + 2))) + d*f*(a**2*d*f*(n + p + 3) - b*(a*(c*f*(p + 1) + d*e*(n + 1)) + b*c*e))*(n + p + 2)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern61 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f61))
    rule61 = ReplacementRule(pattern61, lambda b, n, a, p, f, x, e, d, c : b*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(S(2)*a*d*f*(n + p + S(3)) + b*d*f*x*(n + p + S(2)) - b*(c*f*(p + S(2)) + d*e*(n + S(2))))/(d**S(2)*f**S(2)*(n + p + S(2))*(n + p + S(3))))
    rubi.add(rule61)

    def f62(c, b, n, p, f, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), ZeroQ(m - n - 1), NonzeroQ(m + n + p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern62 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f62))
    rule62 = ReplacementRule(pattern62, lambda c, b, n, p, f, x, d, m, a : a*Int((f*x)**p*(a + b*x)**n*(c + d*x)**n, x) + b*Int((f*x)**(p + S(1))*(a + b*x)**n*(c + d*x)**n, x)/f)
    rubi.add(rule62)

    def f63(b, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern63 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f63))
    rule63 = ReplacementRule(pattern63, lambda b, a, p, f, x, e, d, c : Int(ExpandIntegrand((e + f*x)**p/((a + b*x)*(c + d*x)), x), x))
    rubi.add(rule63)

    def f64(b, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(0, p, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern64 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f64))
    rule64 = ReplacementRule(pattern64, lambda b, a, p, f, x, e, d, c : (-a*f + b*e)*Int((e + f*x)**(p + S(-1))/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule64)

    def f65(b, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ RationalQ(p), Greater(p, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern65 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f65))
    rule65 = ReplacementRule(pattern65, lambda b, a, p, f, x, e, d, c : f*(e + f*x)**(p + S(-1))/(b*d*(p + S(-1))) + Int((e + f*x)**(p + S(-2))*(-a*c*f**S(2) + b*d*e**S(2) + f*x*(-a*d*f - b*c*f + S(2)*b*d*e))/((a + b*x)*(c + d*x)), x)/(b*d))
    rubi.add(rule65)

    def f66(b, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern66 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f66))
    rule66 = ReplacementRule(pattern66, lambda b, a, p, f, x, e, d, c : f*(e + f*x)**(p + S(1))/((p + S(1))*(-a*f + b*e)*(-c*f + d*e)) + Int((e + f*x)**(p + S(1))*(-a*d*f - b*c*f + b*d*e - b*d*f*x)/((a + b*x)*(c + d*x)), x)/((-a*f + b*e)*(-c*f + d*e)))
    rubi.add(rule66)

    def f67(b, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x)])
    pattern67 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f67))
    rule67 = ReplacementRule(pattern67, lambda b, a, p, f, x, e, d, c : b*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - d*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule67)

    def f68(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FractionQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern68 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_/(x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f68))
    rule68 = ReplacementRule(pattern68, lambda b, n, a, p, f, x, e, d, c : Int(ExpandIntegrand((e + f*x)**FractionalPart(p), (c + d*x)**n*(e + f*x)**IntegerPart(p)/(a + b*x), x), x))
    rubi.add(rule68)

    def f69(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ IntegersQ(m, n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x)])
    pattern69 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f69))
    rule69 = ReplacementRule(pattern69, lambda b, n, a, p, f, x, e, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule69)

    def f70(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern70 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f70))
    rule70 = ReplacementRule(pattern70, lambda b, n, a, p, f, x, e, d, c : (c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*d + b*c)**S(2)/(d**S(2)*(n + S(1))*(-c*f + d*e)) - Int((c + d*x)**(n + S(1))*(e + f*x)**p*Simp(a**S(2)*d**S(2)*f*(n + p + S(2)) - S(2)*a*b*d*(c*f*(p + S(1)) + d*e*(n + S(1))) + b**S(2)*c*(c*f*(p + S(1)) + d*e*(n + S(1))) - b**S(2)*d*x*(n + S(1))*(-c*f + d*e), x), x)/(d**S(2)*(n + S(1))*(-c*f + d*e)))
    rubi.add(rule70)

    def f71(b, n, a, p, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(n + p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern71 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**S(2)*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f71))
    rule71 = ReplacementRule(pattern71, lambda b, n, a, p, f, x, e, d, c : b*(a + b*x)*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(n + p + S(3))) + Int((c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(n + p + S(3)) + b*x*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e), x), x)/(d*f*(n + p + S(3))))
    rubi.add(rule71)

    def f72(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern72 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(2)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))), x_),CustomConstraint(f72), )
    def With72(b, a, f, x, e, d, c):
        q = Rt((-c*f + d*e)/(-a*f + b*e), S(3))
        return -sqrt(S(3))*q*ArcTan(S(2)*sqrt(S(3))*q*(a + b*x)**(S(1)/3)/(S(3)*(c + d*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(-c*f + d*e) + q*log(e + f*x)/(-S(2)*c*f + S(2)*d*e) - S(3)*q*log(q*(a + b*x)**(S(1)/3) - (c + d*x)**(S(1)/3))/(-S(2)*c*f + S(2)*d*e)
    rule72 = ReplacementRule(pattern72, lambda b, a, f, x, e, d, c : With72(b, a, f, x, e, d, c))
    rubi.add(rule72)

    def f73(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(2*b*d*e - f*(a*d + b*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern73 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))), x_),CustomConstraint(f73))
    rule73 = ReplacementRule(pattern73, lambda b, a, f, x, e, d, c : b*f*Subst(Int(S(1)/(b*f**S(2)*x**S(2) + d*(-a*f + b*e)**S(2)), x), x, sqrt(a + b*x)*sqrt(c + d*x)))
    rubi.add(rule73)

    def f74(b, n, a, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + 1), RationalQ(m, n), Less(-1, m, 0), SimplerQ(a + b*x, c + d*x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern74 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_),CustomConstraint(f74), )
    def With74(b, n, a, f, x, e, d, m, c):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))/(-a*f + b*e - x**q*(-c*f + d*e)), x), x, (a + b*x)**(S(1)/q)*(c + d*x)**(-S(1)/q))
    rule74 = ReplacementRule(pattern74, lambda b, n, a, f, x, e, d, m, c : With74(b, n, a, f, x, e, d, m, c))
    rubi.add(rule74)

    def f75(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + p + 2), RationalQ(n), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(p, x)])
    pattern75 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f75))
    rule75 = ReplacementRule(pattern75, lambda b, n, a, p, f, x, e, d, m, c : -n*(-c*f + d*e)*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p, x)/((m + S(1))*(-a*f + b*e)) + (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule75)

    def f76(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + p + 3), ZeroQ(a*d*f*(m + 1) + b*c*f*(n + 1) + b*d*e*(p + 1)), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern76 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f76))
    rule76 = ReplacementRule(pattern76, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule76)

    def f77(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern77 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f77))
    rule77 = ReplacementRule(pattern77, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + (a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule77)

    def f78(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n, p), Less(m, -1), Greater(n, 0), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern78 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f78))
    rule78 = ReplacementRule(pattern78, lambda b, n, a, p, f, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p/(b*(m + S(1))) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(-1))*Simp(c*f*p + d*e*n + d*f*x*(n + p), x), x)/(b*(m + S(1))))
    rubi.add(rule78)

    def f79(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n, p), Less(m, -1), Greater(n, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x)])
    pattern79 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f79))
    rule79 = ReplacementRule(pattern79, lambda b, n, a, p, f, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**(p + S(1))*(-a*d + b*c)/(b*(m + S(1))*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-2))*(e + f*x)**p*Simp(a*d*(c*f*(p + S(1)) + d*e*(n + S(-1))) + b*c*(-c*f*(m + p + S(2)) + d*e*(m - n + S(2))) + d*x*(a*d*f*(n + p) + b*(-c*f*(m + n + p + S(1)) + d*e*(m + S(1)))), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule79)

    def f80(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n, p), Less(m, -1), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x)])
    pattern80 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f80))
    rule80 = ReplacementRule(pattern80, lambda b, n, a, p, f, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))/((m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(c*f*(m + p + S(2)) + d*e*n + d*f*x*(m + n + p + S(2)), x), x)/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule80)

    def f81(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Greater(m, 1), NonzeroQ(m + n + p + 1), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern81 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f81))
    rule81 = ReplacementRule(pattern81, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1))))
    rubi.add(rule81)

    def f82(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n, p), Greater(m, 0), Greater(n, 0), NonzeroQ(m + n + p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern82 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f82))
    rule82 = ReplacementRule(pattern82, lambda b, n, a, p, f, x, e, d, m, c : (a + b*x)**m*(c + d*x)**n*(e + f*x)**(p + S(1))/(f*(m + n + p + S(1))) - Int((a + b*x)**(m + S(-1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(a*n*(-c*f + d*e) + c*m*(-a*f + b*e) + x*(b*n*(-c*f + d*e) + d*m*(-a*f + b*e)), x), x)/(f*(m + n + p + S(1))))
    rubi.add(rule82)

    def f83(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Greater(m, 1), NonzeroQ(m + n + p + 1), IntegersQ(2*m, 2*n, 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern83 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f83))
    rule83 = ReplacementRule(pattern83, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(-1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(1))) + Int((a + b*x)**(m + S(-2))*(c + d*x)**n*(e + f*x)**p*Simp(a**S(2)*d*f*(m + n + p + S(1)) + b*x*(a*d*f*(S(2)*m + n + p) - b*(c*f*(m + p) + d*e*(m + n))) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*(m + S(-1))), x), x)/(d*f*(m + n + p + S(1))))
    rubi.add(rule83)

    def f84(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Less(m, -1), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern84 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f84))
    rule84 = ReplacementRule(pattern84, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule84)

    def f85(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Less(m, -1), IntegersQ(2*m, 2*n, 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern85 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f85))
    rule85 = ReplacementRule(pattern85, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule85)

    def f86(b, n, a, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m + n + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x)])
    pattern86 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_/(x_*WC('f', S(1)) + WC('e', S(0))), x_),CustomConstraint(f86))
    rule86 = ReplacementRule(pattern86, lambda b, n, a, f, x, e, d, m, c : b*Int((a + b*x)**(m + S(-1))*(c + d*x)**n, x)/f - (-a*f + b*e)*Int((a + b*x)**(m + S(-1))*(c + d*x)**n/(e + f*x), x)/f)
    rubi.add(rule86)

    def f87(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ PositiveQ(-f/(-c*f + d*e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern87 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_),CustomConstraint(f87))
    rule87 = ReplacementRule(pattern87, lambda b, a, f, x, e, d, c : -S(4)*Subst(Int(x**S(2)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)))
    rubi.add(rule87)

    def f88(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern88 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/4)), x_),CustomConstraint(f88))
    rule88 = ReplacementRule(pattern88, lambda b, a, f, x, e, d, c : sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(1)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x))
    rubi.add(rule88)

    def f89(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ PositiveQ(-f/(-c*f + d*e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern89 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_),CustomConstraint(f89))
    rule89 = ReplacementRule(pattern89, lambda b, a, f, x, e, d, c : -S(4)*Subst(Int(S(1)/(sqrt(c - d*e/f + d*x**S(4)/f)*(-a*f + b*e - b*x**S(4))), x), x, (e + f*x)**(S(1)/4)))
    rubi.add(rule89)

    def f90(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern90 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(3)/4)), x_),CustomConstraint(f90))
    rule90 = ReplacementRule(pattern90, lambda b, a, f, x, e, d, c : sqrt(-f*(c + d*x)/(-c*f + d*e))*Int(S(1)/((a + b*x)*(e + f*x)**(S(3)/4)*sqrt(-c*f/(-c*f + d*e) - d*f*x/(-c*f + d*e))), x)/sqrt(c + d*x))
    rubi.add(rule90)

    def f91(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-c*f + d*e), PositiveQ(c), PositiveQ(e), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern91 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f91))
    rule91 = ReplacementRule(pattern91, lambda b, f, x, e, d, c : S(2)*sqrt(e)*EllipticE(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/b)
    rubi.add(rule91)

    def f92(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-c*f + d*e), PositiveQ(c), PositiveQ(e), NegativeQ(-b/d), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern92 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f92))
    rule92 = ReplacementRule(pattern92, lambda b, f, x, e, d, c : sqrt(-b*x)*Int(sqrt(e + f*x)/(sqrt(-b*x)*sqrt(c + d*x)), x)/sqrt(b*x))
    rubi.add(rule92)

    def f93(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ NonzeroQ(-c*f + d*e), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern93 = Pattern(Integral(sqrt(e_ + x_*WC('f', S(1)))/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f93))
    rule93 = ReplacementRule(pattern93, lambda b, f, x, e, d, c : sqrt(S(1) + d*x/c)*sqrt(e + f*x)*Int(sqrt(S(1) + f*x/e)/(sqrt(b*x)*sqrt(S(1) + d*x/c)), x)/(sqrt(S(1) + f*x/e)*sqrt(c + d*x)))
    rubi.add(rule93)

    def f94(c, b, f, x, e, d, a):
        return functools.reduce(operator.and_, [ PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern94 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f94))
    rule94 = ReplacementRule(pattern94, lambda c, b, f, x, e, d, a : S(2)*EllipticE(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*f + b*e)/d, S(2))/b)
    rubi.add(rule94)

    def f95(c, b, f, x, e, d, a):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern95 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))), x_),CustomConstraint(f95))
    rule95 = ReplacementRule(pattern95, lambda c, b, f, x, e, d, a : sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(e + f*x)*Int(sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))), x)/(sqrt(b*(e + f*x)/(-a*f + b*e))*sqrt(c + d*x)))
    rubi.add(rule95)

    def f96(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ PositiveQ(c), PositiveQ(e), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern96 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f96))
    rule96 = ReplacementRule(pattern96, lambda b, f, x, e, d, c : S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)))
    rubi.add(rule96)

    def f97(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ PositiveQ(c), PositiveQ(e), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern97 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f97))
    rule97 = ReplacementRule(pattern97, lambda b, f, x, e, d, c : S(2)*EllipticF(asin(sqrt(b*x)/(sqrt(c)*Rt(-b/d, S(2)))), c*f/(d*e))*Rt(-b/d, S(2))/(b*sqrt(e)))
    rubi.add(rule97)

    def f98(b, f, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern98 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f98))
    rule98 = ReplacementRule(pattern98, lambda b, f, x, e, d, c : sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)*Int(S(1)/(sqrt(b*x)*sqrt(S(1) + d*x/c)*sqrt(S(1) + f*x/e)), x)/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule98)

    def f99(c, b, f, x, e, d, a):
        return functools.reduce(operator.and_, [ PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e)), SimplerQ(a + b*x, c + d*x), SimplerQ(a + b*x, e + f*x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern99 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f99))
    rule99 = ReplacementRule(pattern99, lambda c, b, f, x, e, d, a : S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b)
    rubi.add(rule99)

    def f100(c, b, f, x, e, d, a):
        return functools.reduce(operator.and_, [ PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e)), SimplerQ(a + b*x, c + d*x), SimplerQ(a + b*x, e + f*x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern100 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f100))
    rule100 = ReplacementRule(pattern100, lambda c, b, f, x, e, d, a : S(2)*sqrt(b**S(2)/((-a*d + b*c)*(-a*f + b*e)))*EllipticF(asin(sqrt(a + b*x)/Rt(-(-a*d + b*c)/d, S(2))), f*(-a*d + b*c)/(d*(-a*f + b*e)))*Rt(-(-a*d + b*c)/d, S(2))/b)
    rubi.add(rule100)

    def f101(c, b, f, x, e, d, a):
        return functools.reduce(operator.and_, [ SimplerQ(a + b*x, c + d*x), SimplerQ(a + b*x, e + f*x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern101 = Pattern(Integral(S(1)/(sqrt(a_ + x_*WC('b', S(1)))*sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))), x_),CustomConstraint(f101))
    rule101 = ReplacementRule(pattern101, lambda c, b, f, x, e, d, a : sqrt(b*(c + d*x)/(-a*d + b*c))*sqrt(b*(e + f*x)/(-a*f + b*e))*Int(S(1)/(sqrt(a + b*x)*sqrt(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))*sqrt(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))), x)/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule101)

    def f102(b, a, f, x, e, d, c):
        return functools.reduce(operator.and_, [ ZeroQ(-a*d*f - b*c*f + 2*b*d*e), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern102 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_),CustomConstraint(f102), )
    def With102(b, a, f, x, e, d, c):
        q = Rt(b*(-a*f + b*e)/(-a*d + b*c)**S(2), S(3))
        return -sqrt(S(3))*ArcTan(S(2)*sqrt(S(3))*q*(c + d*x)**(S(2)/3)/(S(3)*(e + f*x)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*q*(-a*d + b*c)) - log(a + b*x)/(S(2)*q*(-a*d + b*c)) + S(3)*log(q*(c + d*x)**(S(2)/3) - (e + f*x)**(S(1)/3))/(S(4)*q*(-a*d + b*c))
    rule102 = ReplacementRule(pattern102, lambda b, a, f, x, e, d, c : With102(b, a, f, x, e, d, c))
    rubi.add(rule102)

    def f103(b, a, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(-a*d*f - b*c*f + 2*b*d*e), IntegerQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern103 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_/((x_*WC('d', S(1)) + WC('c', S(0)))**(S(1)/3)*(x_*WC('f', S(1)) + WC('e', S(0)))**(S(1)/3)), x_),CustomConstraint(f103))
    rule103 = ReplacementRule(pattern103, lambda b, a, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(S(2)/3)*(e + f*x)**(S(2)/3)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + f*Int((a + b*x)**(m + S(1))*(a*d*(S(3)*m + S(1)) - S(3)*b*c*(S(3)*m + S(5)) - S(2)*b*d*x*(S(3)*m + S(7)))/((c + d*x)**(S(1)/3)*(e + f*x)**(S(1)/3)), x)/(S(6)*(m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule103)

    def f104(c, b, n, p, f, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), ZeroQ(m - n), PositiveQ(a), PositiveQ(c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern104 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f104))
    rule104 = ReplacementRule(pattern104, lambda c, b, n, p, f, x, d, m, a : Int((f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule104)

    def f105(c, b, n, p, f, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), ZeroQ(m - n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern105 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f105))
    rule105 = ReplacementRule(pattern105, lambda c, b, n, p, f, x, d, m, a : (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)*(a*c + b*d*x**S(2))**(-FracPart(m))*Int((f*x)**p*(a*c + b*d*x**S(2))**m, x))
    rubi.add(rule105)

    def f106(c, b, n, p, f, x, d, m, a):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PositiveIntegerQ(m - n), NonzeroQ(m + n + p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern106 = Pattern(Integral((x_*WC('f', S(1)))**WC('p', S(1))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f106))
    rule106 = ReplacementRule(pattern106, lambda c, b, n, p, f, x, d, m, a : Int(ExpandIntegrand((f*x)**p*(a + b*x)**n*(c + d*x)**n, (a + b*x)**(m - n), x), x))
    rubi.add(rule106)

    def f107(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern107 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f107))
    rule107 = ReplacementRule(pattern107, lambda b, n, a, p, f, x, e, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x))
    rubi.add(rule107)

    def f108(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(m + n + p + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern108 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1)), x_),CustomConstraint(f108))
    rule108 = ReplacementRule(pattern108, lambda b, n, a, p, f, x, e, d, m, c : b*(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*(m + S(1)) - b*d*f*x*(m + n + p + S(3)) - b*(c*f*(m + p + S(2)) + d*e*(m + n + S(2))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule108)

    def f109(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + p + 2), NegativeIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(p, x)])
    pattern109 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f109))
    rule109 = ReplacementRule(pattern109, lambda b, n, a, p, f, x, e, d, m, c : (a + b*x)**(m + S(1))*(e + f*x)**(-m + S(-1))*(-a*d + b*c)**n*(-a*f + b*e)**(-n + S(-1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/(m + S(1)))
    rubi.add(rule109)

    def f110(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern110 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f110))
    rule110 = ReplacementRule(pattern110, lambda b, n, a, p, f, x, e, d, m, c : ((c + d*x)*(-a*f + b*e)/((e + f*x)*(-a*d + b*c)))**(-n)*(a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*Hypergeometric2F1(m + S(1), -n, m + S(2), -(a + b*x)*(-c*f + d*e)/((e + f*x)*(-a*d + b*c)))/((m + S(1))*(-a*f + b*e)))
    rubi.add(rule110)

    def f111(b, n, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveQ(c), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern111 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_),CustomConstraint(f111))
    rule111 = ReplacementRule(pattern111, lambda b, n, p, f, x, e, d, m, c : c**n*e**p*(b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*x/c, -f*x/e)/(b*(m + S(1))))
    rubi.add(rule111)

    def f112(b, n, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveQ(-d/(b*c)), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern112 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_),CustomConstraint(f112))
    rule112 = ReplacementRule(pattern112, lambda b, n, p, f, x, e, d, m, c : (d/(-c*f + d*e))**(-p)*(-d/(b*c))**(-m)*(c + d*x)**(n + S(1))*AppellF1(n + S(1), -m, -p, n + S(2), S(1) + d*x/c, -f*(c + d*x)/(-c*f + d*e))/(d*(n + S(1))))
    rubi.add(rule112)

    def f113(b, n, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern113 = Pattern(Integral((x_*WC('b', S(1)))**m_*(c_ + x_*WC('d', S(1)))**n_*(e_ + x_*WC('f', S(1)))**p_, x_),CustomConstraint(f113))
    rule113 = ReplacementRule(pattern113, lambda b, n, p, f, x, e, d, m, c : c**IntPart(n)*(S(1) + d*x/c)**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((b*x)**m*(S(1) + d*x/c)**n*(e + f*x)**p, x))
    rubi.add(rule113)

    def f114(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ IntegerQ(p), PositiveQ(b/(-a*d + b*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x)])
    pattern114 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f114))
    rule114 = ReplacementRule(pattern114, lambda b, n, a, p, f, x, e, d, m, c : b**(-p + S(-1))*(b/(-a*d + b*c))**(-n)*(a + b*x)**(m + S(1))*(-a*f + b*e)**p*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(m + S(1)))
    rubi.add(rule114)

    def f115(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x)])
    pattern115 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f115))
    rule115 = ReplacementRule(pattern115, lambda b, n, a, p, f, x, e, d, m, c : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule115)

    def f116(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern116 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f116))
    rule116 = ReplacementRule(pattern116, lambda b, n, a, p, f, x, e, d, m, c : (b/(-a*d + b*c))**(-n)*(b/(-a*f + b*e))**(-p)*(a + b*x)**(m + S(1))*AppellF1(m + S(1), -n, -p, m + S(2), -d*(a + b*x)/(-a*d + b*c), -f*(a + b*x)/(-a*f + b*e))/(b*(m + S(1))))
    rubi.add(rule116)

    def f117(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveQ(b/(-a*d + b*c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern117 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f117))
    rule117 = ReplacementRule(pattern117, lambda b, n, a, p, f, x, e, d, m, c : (b/(-a*f + b*e))**(-IntPart(p))*(b*(e + f*x)/(-a*f + b*e))**(-FracPart(p))*(e + f*x)**FracPart(p)*Int((a + b*x)**m*(c + d*x)**n*(b*e/(-a*f + b*e) + b*f*x/(-a*f + b*e))**p, x))
    rubi.add(rule117)

    def f118(b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern118 = Pattern(Integral((a_ + x_*WC('b', S(1)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_, x_),CustomConstraint(f118))
    rule118 = ReplacementRule(pattern118, lambda b, n, a, p, f, x, e, d, m, c : (b/(-a*d + b*c))**(-IntPart(n))*(b*(c + d*x)/(-a*d + b*c))**(-FracPart(n))*(c + d*x)**FracPart(n)*Int((a + b*x)**m*(e + f*x)**p*(b*c/(-a*d + b*c) + b*d*x/(-a*d + b*c))**n, x))
    rubi.add(rule118)

    def f119(u, b, n, a, p, f, x, e, d, m, c):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern119 = Pattern(Integral((e_ + u_*WC('f', S(1)))**WC('p', S(1))*(u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1)), x_),CustomConstraint(f119))
    rule119 = ReplacementRule(pattern119, lambda u, b, n, a, p, f, x, e, d, m, c : Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule119)

    def f120(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern120 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f120))
    rule120 = ReplacementRule(pattern120, lambda g, b, n, a, f, h, x, e, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)*(g + h*x), x), x))
    rubi.add(rule120)

    def f121(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ ZeroQ(m + n + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x)])
    pattern121 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f121))
    rule121 = ReplacementRule(pattern121, lambda g, b, n, a, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(2)*d*f*h*m - a*b*(-c*f*h*(m + S(1)) + d*(e*h + f*g)) + b**S(2)*d*e*g + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)) + (a*d*f*h*m + b*(-c*f*h*(m + S(2)) + d*(e*h + f*g)))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d))
    rubi.add(rule121)

    def f122(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n), Less(m, -1), Less(n, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern122 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f122))
    rule122 = ReplacementRule(pattern122, lambda g, b, n, a, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*c*d*f*h*(n + S(1)) + a*b*(c**S(2)*f*h*(m + S(1)) - c*d*(e*h + f*g)*(m + n + S(2)) + d**S(2)*e*g*(m + S(1))) + b**S(2)*c*d*e*g*(n + S(1)) + x*(a**S(2)*d**S(2)*f*h*(n + S(1)) - a*b*d**S(2)*(n + S(1))*(e*h + f*g) + b**S(2)*(c**S(2)*f*h*(m + S(1)) - c*d*(m + S(1))*(e*h + f*g) + d**S(2)*e*g*(m + n + S(2)))))/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)) - (a**S(2)*d**S(2)*f*h*(n**S(2) + S(3)*n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m**S(2) + S(3)*m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m**S(2) + m*(S(2)*n + S(5)) + n**S(2) + S(5)*n + S(6))))*Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(1)), x)/(b*d*(m + S(1))*(n + S(1))*(-a*d + b*c)**S(2)))
    rubi.add(rule122)

    def f123(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x)])
    pattern123 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f123))
    rule123 = ReplacementRule(pattern123, lambda g, b, n, a, f, h, x, e, d, m, c : (-d*(m + n + S(3))*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)) + f*h/b**S(2))*Int((a + b*x)**(m + S(2))*(c + d*x)**n, x) + (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(-a**S(3)*d*f*h*(n + S(2)) - a**S(2)*b*(c*f*h*m - d*(e*h + f*g)*(m + n + S(3))) - a*b**S(2)*(c*(e*h + f*g) + d*e*g*(S(2)*m + n + S(4))) + b**S(3)*c*e*g*(m + S(2)) + b*x*(a**S(2)*d*f*h*(m - n) - a*b*(S(2)*c*f*h*(m + S(1)) - d*(n + S(1))*(e*h + f*g)) + b**S(2)*(c*(m + S(1))*(e*h + f*g) - d*e*g*(m + n + S(2)))))/(b**S(2)*(m + S(1))*(m + S(2))*(-a*d + b*c)**S(2)))
    rubi.add(rule123)

    def f124(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(m + 1), NonzeroQ(m + n + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x)])
    pattern124 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f124))
    rule124 = ReplacementRule(pattern124, lambda g, b, n, a, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a**S(2)*d*f*h*(n + S(2)) + a*b*(c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*d*e*g*(m + n + S(3)) + b*f*h*x*(m + S(1))*(-a*d + b*c))/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))) - (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**(m + S(1))*(c + d*x)**n, x)/(b**S(2)*d*(m + S(1))*(-a*d + b*c)*(m + n + S(3))))
    rubi.add(rule124)

    def f125(g, b, n, a, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ NonzeroQ(m + n + 2), NonzeroQ(m + n + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x)])
    pattern125 = Pattern(Integral((e_ + x_*WC('f', S(1)))*(x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f125))
    rule125 = ReplacementRule(pattern125, lambda g, b, n, a, f, h, x, e, d, m, c : -(a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(a*d*f*h*(n + S(2)) + b*c*f*h*(m + S(2)) - b*d*f*h*x*(m + n + S(2)) - b*d*(e*h + f*g)*(m + n + S(3)))/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))) + (a**S(2)*d**S(2)*f*h*(n + S(1))*(n + S(2)) + a*b*d*(n + S(1))*(S(2)*c*f*h*(m + S(1)) - d*(e*h + f*g)*(m + n + S(3))) + b**S(2)*(c**S(2)*f*h*(m + S(1))*(m + S(2)) - c*d*(m + S(1))*(e*h + f*g)*(m + n + S(3)) + d**S(2)*e*g*(m + n + S(2))*(m + n + S(3))))*Int((a + b*x)**m*(c + d*x)**n, x)/(b**S(2)*d**S(2)*(m + n + S(2))*(m + n + S(3))))
    rubi.add(rule125)

    def f126(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern126 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f126))
    rule126 = ReplacementRule(pattern126, lambda g, b, n, a, p, f, h, x, e, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x), x), x))
    rubi.add(rule126)

    def f127(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n), Less(m, -1), Greater(n, 0), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern127 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f127))
    rule127 = ReplacementRule(pattern127, lambda g, b, n, a, p, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule127)

    def f128(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m, n), Less(m, -1), Greater(n, 0), IntegersQ(2*m, 2*n, 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern128 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f128))
    rule128 = ReplacementRule(pattern128, lambda g, b, n, a, p, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**(p + S(1))*(-a*h + b*g)/(b*(m + S(1))*(-a*f + b*e)) - Int((a + b*x)**(m + S(1))*(c + d*x)**(n + S(-1))*(e + f*x)**p*Simp(b*c*(m + S(1))*(-e*h + f*g) + d*x*(b*(m + S(1))*(-e*h + f*g) + f*(-a*h + b*g)*(n + p + S(1))) + (-a*h + b*g)*(c*f*(p + S(1)) + d*e*n), x), x)/(b*(m + S(1))*(-a*f + b*e)))
    rubi.add(rule128)

    def f129(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Less(m, -1), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern129 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f129))
    rule129 = ReplacementRule(pattern129, lambda g, b, n, a, p, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule129)

    def f130(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Less(m, -1), IntegersQ(2*m, 2*n, 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern130 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f130))
    rule130 = ReplacementRule(pattern130, lambda g, b, n, a, p, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule130)

    def f131(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Greater(m, 0), NonzeroQ(m + n + p + 2), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern131 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f131))
    rule131 = ReplacementRule(pattern131, lambda g, b, n, a, p, f, h, x, e, d, m, c : h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2))))
    rubi.add(rule131)

    def f132(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ RationalQ(m), Greater(m, 0), NonzeroQ(m + n + p + 2), IntegersQ(2*m, 2*n, 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern132 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f132))
    rule132 = ReplacementRule(pattern132, lambda g, b, n, a, p, f, h, x, e, d, m, c : h*(a + b*x)**m*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))/(d*f*(m + n + p + S(2))) + Int((a + b*x)**(m + S(-1))*(c + d*x)**n*(e + f*x)**p*Simp(a*d*f*g*(m + n + p + S(2)) - h*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e*m) + x*(b*d*f*g*(m + n + p + S(2)) + h*(a*d*f*m - b*(c*f*(m + p + S(1)) + d*e*(m + n + S(1))))), x), x)/(d*f*(m + n + p + S(2))))
    rubi.add(rule132)

    def f133(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(m + n + p + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern133 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f133))
    rule133 = ReplacementRule(pattern133, lambda g, b, n, a, p, f, h, x, e, d, m, c : (a + b*x)**(m + S(1))*(c + d*x)**(n + S(1))*(e + f*x)**(p + S(1))*(-a*h + b*g)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)) + Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*Simp(-d*f*x*(-a*h + b*g)*(m + n + p + S(3)) + (m + S(1))*(a*d*f*g + b*c*e*h - b*g*(c*f + d*e)) - (-a*h + b*g)*(c*f*(p + S(1)) + d*e*(n + S(1))), x), x)/((m + S(1))*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule133)

    def f134(g, b, a, p, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern134 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f134))
    rule134 = ReplacementRule(pattern134, lambda g, b, a, p, f, h, x, e, d, c : (-a*h + b*g)*Int((e + f*x)**p/(a + b*x), x)/(-a*d + b*c) - (-c*h + d*g)*Int((e + f*x)**p/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule134)

    def f135(g, b, n, a, p, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h, n, p), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(n, x), FreeQ(p, x)])
    pattern135 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))/(x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f135))
    rule135 = ReplacementRule(pattern135, lambda g, b, n, a, p, f, h, x, e, d, c : h*Int((c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((c + d*x)**n*(e + f*x)**p/(a + b*x), x)/b)
    rubi.add(rule135)

    def f136(g, c, b, f, h, x, e, d, a):
        return functools.reduce(operator.and_, [ SimplerQ(a + b*x, e + f*x), SimplerQ(c + d*x, e + f*x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern136 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(c_ + x_*WC('d', S(1)))*sqrt(e_ + x_*WC('f', S(1)))*sqrt(x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f136))
    rule136 = ReplacementRule(pattern136, lambda g, c, b, f, h, x, e, d, a : h*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)), x)/f + (-e*h + f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)), x)/f)
    rubi.add(rule136)

    def f137(g, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern137 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0))), x_),CustomConstraint(f137))
    rule137 = ReplacementRule(pattern137, lambda g, b, n, a, p, f, h, x, e, d, m, c : h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p, x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p, x)/b)
    rubi.add(rule137)

    def f138(g, q, b, a, p, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(0, p, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern138 = Pattern(Integral((x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_/((x_*WC('b', S(1)) + WC('a', S(0)))*(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f138))
    rule138 = ReplacementRule(pattern138, lambda g, q, b, a, p, f, h, x, e, d, c : (-a*f + b*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(a + b*x), x)/(-a*d + b*c) - (-c*f + d*e)*Int((e + f*x)**(p + S(-1))*(g + h*x)**q/(c + d*x), x)/(-a*d + b*c))
    rubi.add(rule138)

    def f139(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern139 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f139))
    rule139 = ReplacementRule(pattern139, lambda g, b, a, f, h, x, e, d, c : -S(2)*sqrt(d*(e + f*x)/(-c*f + d*e))*sqrt(d*(g + h*x)/(-c*h + d*g))*EllipticPi(-b*(-c*f + d*e)/(f*(-a*d + b*c)), asin(sqrt(-f/(-c*f + d*e))*sqrt(c + d*x)), h*(-c*f + d*e)/(f*(-c*h + d*g)))/(sqrt(-f/(-c*f + d*e))*sqrt(e + f*x)*sqrt(g + h*x)*(-a*d + b*c)))
    rubi.add(rule139)

    def f140(g, b, n, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ IntegerQ(n + 1/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern140 = Pattern(Integral((x_*WC('d', S(1)) + WC('c', S(0)))**n_/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f140))
    rule140 = ReplacementRule(pattern140, lambda g, b, n, a, f, h, x, e, d, c : Int(ExpandIntegrand(S(1)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), (c + d*x)**(n + S(1)/2)/(a + b*x), x), x))
    rubi.add(rule140)

    def f141(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern141 = Pattern(Integral(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))), x_),CustomConstraint(f141))
    rule141 = ReplacementRule(pattern141, lambda g, b, a, f, h, x, e, d, c : (-a*f + b*e)*(-a*h + b*g)*Int(S(1)/((a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2) + Int((-a*f*h + b*e*h + b*f*g + b*f*h*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/b**S(2))
    rubi.add(rule141)

    def f142(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern142 = Pattern(Integral(S(1)/(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f142))
    rule142 = ReplacementRule(pattern142, lambda g, b, a, f, h, x, e, d, c : -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)))
    rubi.add(rule142)

    def f143(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern143 = Pattern(Integral(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f143))
    rule143 = ReplacementRule(pattern143, lambda g, b, a, f, h, x, e, d, c : -S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*(-c*h + d*g)*Subst(Int(sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))/sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1)), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)*(-a*h + b*g)**S(2)))
    rubi.add(rule143)

    def f144(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern144 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f144))
    rule144 = ReplacementRule(pattern144, lambda g, b, a, f, h, x, e, d, c : S(2)*sqrt((c + d*x)*(-a*h + b*g)/((a + b*x)*(-c*h + d*g)))*sqrt((e + f*x)*(-a*h + b*g)/((a + b*x)*(-e*h + f*g)))*(a + b*x)*Subst(Int(S(1)/((-b*x**S(2) + h)*sqrt(x**S(2)*(-a*d + b*c)/(-c*h + d*g) + S(1))*sqrt(x**S(2)*(-a*f + b*e)/(-e*h + f*g) + S(1))), x), x, sqrt(g + h*x)/sqrt(a + b*x))/(sqrt(c + d*x)*sqrt(e + f*x)))
    rubi.add(rule144)

    def f145(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern145 = Pattern(Integral(S(1)/((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f145))
    rule145 = ReplacementRule(pattern145, lambda g, b, a, f, h, x, e, d, c : b*Int(sqrt(c + d*x)/((a + b*x)**(S(3)/2)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c) - d*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(-a*d + b*c))
    rubi.add(rule145)

    def f146(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern146 = Pattern(Integral(sqrt(x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_*WC('d', S(1)) + WC('c', S(0)))/(sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f146))
    rule146 = ReplacementRule(pattern146, lambda g, b, a, f, h, x, e, d, c : sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)/(h*sqrt(e + f*x)) - (-c*f + d*e)*(-e*h + f*g)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*(e + f*x)**(S(3)/2)*sqrt(g + h*x)), x)/(S(2)*f*h) + (-c*f + d*e)*(-S(2)*a*f*h + b*e*h + b*f*g)*Int(S(1)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h) + (a*d*f*h - b*(-c*f*h + d*e*h + d*f*g))*Int(sqrt(e + f*x)/(sqrt(a + b*x)*sqrt(c + d*x)*sqrt(g + h*x)), x)/(S(2)*f**S(2)*h))
    rubi.add(rule146)

    def f147(g, b, a, f, h, x, e, d, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern147 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)/(sqrt(x_*WC('d', S(1)) + WC('c', S(0)))*sqrt(x_*WC('f', S(1)) + WC('e', S(0)))*sqrt(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f147))
    rule147 = ReplacementRule(pattern147, lambda g, b, a, f, h, x, e, d, c : b*Int(sqrt(a + b*x)*sqrt(c + d*x)/(sqrt(e + f*x)*sqrt(g + h*x)), x)/d - (-a*d + b*c)*Int(sqrt(a + b*x)/(sqrt(c + d*x)*sqrt(e + f*x)*sqrt(g + h*x)), x)/d)
    rubi.add(rule147)

    def f148(g, q, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ IntegersQ(p, q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x)])
    pattern148 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_),CustomConstraint(f148))
    rule148 = ReplacementRule(pattern148, lambda g, q, b, n, a, p, f, h, x, e, d, m, c : Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x))
    rubi.add(rule148)

    def f149(g, q, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern149 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_, x_),CustomConstraint(f149))
    rule149 = ReplacementRule(pattern149, lambda g, q, b, n, a, p, f, h, x, e, d, m, c : h*Int((a + b*x)**(m + S(1))*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b + (-a*h + b*g)*Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**(q + S(-1)), x)/b)
    rubi.add(rule149)

    def f150(g, q, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h, m, n, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern150 = Pattern(Integral((x_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(x_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(x_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_),CustomConstraint(f150))
    rule150 = ReplacementRule(pattern150, lambda g, q, b, n, a, p, f, h, x, e, d, m, c : Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x))
    rubi.add(rule150)

    def f151(g, u, q, b, n, a, p, f, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern151 = Pattern(Integral((u_*WC('b', S(1)) + WC('a', S(0)))**WC('m', S(1))*(u_*WC('d', S(1)) + WC('c', S(0)))**WC('n', S(1))*(u_*WC('f', S(1)) + WC('e', S(0)))**WC('p', S(1))*(u_*WC('h', S(1)) + WC('g', S(0)))**WC('q', S(1)), x_),CustomConstraint(f151))
    rule151 = ReplacementRule(pattern151, lambda g, u, q, b, n, a, p, f, h, x, e, d, m, c : Subst(Int((a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule151)

    def f152(g, r, q, b, n, a, p, f, i, h, x, e, d, m, c):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h, i, m, n, p, q, r), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(i, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern152 = Pattern(Integral(((x_*WC('b', S(1)) + WC('a', S(0)))**m_*(x_*WC('d', S(1)) + WC('c', S(0)))**n_*(x_*WC('f', S(1)) + WC('e', S(0)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**q_*WC('i', S(1)))**r_, x_),CustomConstraint(f152))
    rule152 = ReplacementRule(pattern152, lambda g, r, q, b, n, a, p, f, i, h, x, e, d, m, c : (i*(a + b*x)**m*(c + d*x)**n*(e + f*x)**p*(g + h*x)**q)**r*(a + b*x)**(-m*r)*(c + d*x)**(-n*r)*(e + f*x)**(-p*r)*(g + h*x)**(-q*r)*Int((a + b*x)**(m*r)*(c + d*x)**(n*r)*(e + f*x)**(p*r)*(g + h*x)**(q*r), x))
    rubi.add(rule152)

    return rubi
