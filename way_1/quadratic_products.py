
from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (Int, Set, With, Module, Scan, MapAnd, FalseQ, ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ, PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ, ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ, NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart, FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest, SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient, Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart, IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan, ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec, ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less, Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ, PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ, ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ, Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ, SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator, NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ, InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ, EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree, PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts, TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ, NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll, ExpandLinearProduct, GCD, ContentFactor, NumericFactor, NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst, ExpandExpression, Apart, SmartApart, MatchQ, PolynomialQuotientRemainder, FreeFactors, NonfreeFactors, RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms, ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup, AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor, RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon, MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ, GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList, PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ, RationalFunctionFactors, NonrationalFunctionFactors, Reverse, RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand, SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree, CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree, GeneralizedBinomialParts, GeneralizedTrinomialDegree, GeneralizedTrinomialParts, MonomialQ, MonomialSumQ, MinimumMonomialExponent, MonomialExponent, LinearMatchQ, PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ, TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ, QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms, NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial, PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD, AlgebraicFunctionFactors, NonalgebraicFunctionFactors, QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ, Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors, NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop, CombineExponents, FactorInteger, FactorAbsurdNumber, SubstForInverseFunction, SubstForFractionalPower, SubstForFractionalPowerOfQuotientOfLinears, FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ, SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ, FractionalPowerSubexpressionQ, Apply, FactorNumericGcd, MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ, TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest, OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors, PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn, PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree, FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify, FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand, NormalizeIntegrandAux, NormalizeIntegrandFactor, NormalizeIntegrandFactorBase, NormalizeTogether, NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors, SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm, TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum, UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear, PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ, IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor, FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ, FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator, SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand, SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM, SubstForFractionalPowerOfLinear, FractionalPowerOfLinear, InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig, FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ, PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ, KnownTangentIntegrandQ, KnownCotangentIntegrandQ, KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst, AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand, ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp, ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ, FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ, PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ, FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ, FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ, FunctionOfLog, PowerVariableExpn, PowerVariableDegree, PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic, SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ, Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ, SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2, ConstantFactor, SameQ, ReplacePart, CommonFactors, MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential, FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux, FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev, rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent, RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct, SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma, FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ, _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify, _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum, _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux, TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist)
    from sympy import Integral, S, sqrt
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)

    from sympy.integrals.rubi.constraints import (constraint_freeq_a, constraint_freeq_b, constraint_freeq_c, constraint_freeq_d, constraint_freeq_e, constraint_freeq_f,
        constraint_freeq_g, constraint_freeq_h, constraint_freeq_m, constraint_freeq_n, constraint_freeq_p)

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

    import functools, operator

def quadratic_products(rubi):
    def f1(x, a, c, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern1 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f1))
    rule1 = ReplacementRule(pattern1, lambda x, a, c, b : (b/S(2) + c*x)*Int(S(1)/(b/S(2) + c*x), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule1)

    def f2(a, c, x, p, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern2 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f2))
    rule2 = ReplacementRule(pattern2, lambda a, c, x, p, b : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule2)

    def f3(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), PositiveIntegerQ(p), PerfectSquareQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern3 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f3), )
    def With3(a, c, x, p, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c**(-p)*Int(Simp(b/S(2) + c*x - q/S(2), x)**p*Simp(b/S(2) + c*x + q/S(2), x)**p, x)
    rule3 = ReplacementRule(pattern3, lambda a, c, x, p, b : With3(a, c, x, p, b))
    rubi.add(rule3)

    def f4(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern4 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f4))
    rule4 = ReplacementRule(pattern4, lambda a, c, x, p, b : Int(ExpandIntegrand((a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule4)

    def f5(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Greater(p, 0), IntegerQ(4*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern5 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f5))
    rule5 = ReplacementRule(pattern5, lambda a, c, x, p, b : -p*(-S(4)*a*c + b**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*(S(2)*p + S(1))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule5)

    def f6(x, a, c, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern6 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_),CustomConstraint(f6))
    rule6 = ReplacementRule(pattern6, lambda x, a, c, b : -S(2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule6)

    def f7(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), Unequal(p, -3/2), IntegerQ(4*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern7 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f7))
    rule7 = ReplacementRule(pattern7, lambda a, c, x, p, b : -S(2)*c*(S(2)*p + S(3))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule7)

    def f8(x, a, c, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), PosQ(-4*a*c + b**2), PerfectSquareQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern8 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f8), )
    def With8(x, a, c, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c*Int(S(1)/Simp(b/S(2) + c*x - q/S(2), x), x)/q - c*Int(S(1)/Simp(b/S(2) + c*x + q/S(2), x), x)/q
    rule8 = ReplacementRule(pattern8, lambda x, a, c, b : With8(x, a, c, b))
    rubi.add(rule8)

    def f9(x, a, c, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern9 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f9))
    rule9 = ReplacementRule(pattern9, lambda x, a, c, b : -S(2)*Subst(Int(S(1)/Simp(-S(4)*a*c + b**S(2) - x**S(2), x), x), x, b + S(2)*c*x))
    rubi.add(rule9)

    def f10(a, c, x, p, b):
        return functools.reduce(operator.and_, [ PositiveQ(4*a - b**2/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern10 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f10))
    rule10 = ReplacementRule(pattern10, lambda a, c, x, p, b : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p, x), x, b + S(2)*c*x)/(S(2)*c))
    rubi.add(rule10)

    def f11(x, c, b):
        return functools.reduce(operator.and_, [ FreeQ(List(b, c), x), FreeQ(b, x), FreeQ(c, x)])
    pattern11 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f11))
    rule11 = ReplacementRule(pattern11, lambda x, c, b : S(2)*Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(b*x + c*x**S(2))))
    rubi.add(rule11)

    def f12(x, a, c, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern12 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f12))
    rule12 = ReplacementRule(pattern12, lambda x, a, c, b : S(2)*Subst(Int(S(1)/(S(4)*c - x**S(2)), x), x, (b + S(2)*c*x)/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule12)

    def f13(x, p, c, b):
        return functools.reduce(operator.and_, [ RationalQ(p), LessEqual(3, Denominator(p), 4), FreeQ(b, x), FreeQ(c, x)])
    pattern13 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f13))
    rule13 = ReplacementRule(pattern13, lambda x, p, c, b : (-c*(b*x + c*x**S(2))/b**S(2))**(-p)*(b*x + c*x**S(2))**p*Int((-c*x/b - c**S(2)*x**S(2)/b**S(2))**p, x))
    rubi.add(rule13)

    def f14(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern14 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f14), )
    def With14(a, c, x, p, b):
        d = Denominator(p)
        if LessEqual(S(3), d, S(4)):
            return d*sqrt((b + S(2)*c*x)**S(2))*Subst(Int(x**(d*(p + S(1)) + S(-1))/sqrt(-S(4)*a*c + b**S(2) + S(4)*c*x**d), x), x, (a + b*x + c*x**S(2))**(S(1)/d))/(b + S(2)*c*x)
        print("Unable to Integrate")
    rule14 = ReplacementRule(pattern14, lambda a, c, x, p, b : With14(a, c, x, p, b))
    rubi.add(rule14)

    def f15(a, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern15 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f15), )
    def With15(a, c, x, p, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -((-b - S(2)*c*x + q)/(S(2)*q))**(-p + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Hypergeometric2F1(-p, p + S(1), p + S(2), (b + S(2)*c*x + q)/(S(2)*q))/(q*(p + S(1)))
    rule15 = ReplacementRule(pattern15, lambda a, c, x, p, b : With15(a, c, x, p, b))
    rubi.add(rule15)

    def f16(u, a, c, x, p, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern16 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f16))
    rule16 = ReplacementRule(pattern16, lambda u, a, c, x, p, b : Subst(Int((a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule16)

    def f17(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), IntegerQ(m/2 + 1/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern17 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f17))
    rule17 = ReplacementRule(pattern17, lambda e, m, a, c, d, p, x, b : c**(-m/S(2) + S(-1)/2)*e**m*(a + b*x + c*x**S(2))**(m/S(2) + p + S(1)/2)/(m + S(2)*p + S(1)))
    rubi.add(rule17)

    def f18(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern18 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f18))
    rule18 = ReplacementRule(pattern18, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*log(RemoveContent(d + e*x, x))/e)
    rubi.add(rule18)

    def f19(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern19 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f19))
    rule19 = ReplacementRule(pattern19, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule19)

    def f20(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 2), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern20 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f20))
    rule20 = ReplacementRule(pattern20, lambda e, m, a, c, d, p, x, b : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule20)

    def f21(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern21 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0)))**S(2), x_),CustomConstraint(f21))
    rule21 = ReplacementRule(pattern21, lambda e, a, c, d, x, b : sqrt(a + b*x + c*x**S(2))*Int((b + S(2)*c*x)/(d + e*x)**S(2), x)/(b + S(2)*c*x))
    rubi.add(rule21)

    def f22(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), NonzeroQ(m + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern22 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f22))
    rule22 = ReplacementRule(pattern22, lambda e, m, a, c, d, x, b : (d + e*x)**(m + S(1))*sqrt(a + b*x + c*x**S(2))/(e*(m + S(2))) - (-b*e + S(2)*c*d)*sqrt(a + b*x + c*x**S(2))*Int((d + e*x)**m, x)/(e*(b + S(2)*c*x)*(m + S(2))))
    rubi.add(rule22)

    def f23(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern23 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f23))
    rule23 = ReplacementRule(pattern23, lambda e, a, c, d, x, b : -S(4)*c*e*sqrt(a + b*x + c*x**S(2))/((d + e*x)*(-b*e + S(2)*c*d)**S(2)) + S(2)*c*Int(S(1)/((d + e*x)*sqrt(a + b*x + c*x**S(2))), x)/(-b*e + S(2)*c*d))
    rubi.add(rule23)

    def f24(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 3), NonzeroQ(m + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern24 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f24))
    rule24 = ReplacementRule(pattern24, lambda e, m, a, c, d, p, x, b : -S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)**S(2)*(m*p + S(-1))) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(2))*(-b*e + S(2)*c*d)))
    rubi.add(rule24)

    def f25(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), NonzeroQ(p + 3/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern25 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f25))
    rule25 = ReplacementRule(pattern25, lambda e, a, c, d, p, x, b : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule25)

    def f26(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Greater(p, 1), Inequality(-2, LessEqual, m, Less, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern26 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f26))
    rule26 = ReplacementRule(pattern26, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))))
    rubi.add(rule26)

    def f27(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Greater(p, 1), Less(m, -2), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern27 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f27))
    rule27 = ReplacementRule(pattern27, lambda e, m, a, c, d, p, x, b : S(2)*c*p*(S(2)*p + S(-1))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2))))
    rubi.add(rule27)

    def f28(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Greater(p, 0), NonzeroQ(m + 2*p), NonzeroQ(m + 2*p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern28 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f28))
    rule28 = ReplacementRule(pattern28, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)**S(2)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))))
    rubi.add(rule28)

    def f29(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Less(p, -1), Inequality(0, Less, m, LessEqual, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern29 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f29))
    rule29 = ReplacementRule(pattern29, lambda e, m, a, c, d, p, x, b : e**S(2)*m*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule29)

    def f30(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern30 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f30))
    rule30 = ReplacementRule(pattern30, lambda e, m, a, c, d, p, x, b : e**S(2)*m*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) - e*m*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule30)

    def f31(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Less(p, -1), NonzeroQ(m + p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern31 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f31))
    rule31 = ReplacementRule(pattern31, lambda e, m, a, c, d, p, x, b : S(2)*c*e**S(2)*(m + S(2)*p + S(2))*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) - S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule31)

    def f32(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Greater(m, 0), NonzeroQ(m + 2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern32 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f32))
    rule32 = ReplacementRule(pattern32, lambda e, m, a, c, d, p, x, b : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*(m + S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(m + S(2)*p + S(1))))
    rubi.add(rule32)

    def f33(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Less(m, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern33 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f33))
    rule33 = ReplacementRule(pattern33, lambda e, m, a, c, d, p, x, b : S(2)*c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-b*e + S(2)*c*d)) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule33)

    def f34(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern34 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f34))
    rule34 = ReplacementRule(pattern34, lambda e, m, a, c, d, p, x, b : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m, x))
    rubi.add(rule34)

    def f35(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern35 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f35))
    rule35 = ReplacementRule(pattern35, lambda e, m, a, c, d, p, x, b : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule35)

    def f36(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern36 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f36))
    rule36 = ReplacementRule(pattern36, lambda e, m, a, c, d, p, x : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule36)

    def f37(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern37 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f37))
    rule37 = ReplacementRule(pattern37, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule37)

    def f38(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern38 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f38))
    rule38 = ReplacementRule(pattern38, lambda e, m, a, c, d, p, x : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule38)

    def f39(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern39 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f39))
    rule39 = ReplacementRule(pattern39, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule39)

    def f40(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern40 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f40))
    rule40 = ReplacementRule(pattern40, lambda e, m, a, c, d, p, x : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(p + S(1))))
    rubi.add(rule40)

    def f41(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern41 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f41))
    rule41 = ReplacementRule(pattern41, lambda e, a, c, d, p, x, b : -e**S(2)*(p + S(2))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule41)

    def f42(e, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern42 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**S(2), x_),CustomConstraint(f42))
    rule42 = ReplacementRule(pattern42, lambda e, a, c, d, p, x : -e**S(2)*(p + S(2))*Int((a + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)/(c*(p + S(1))))
    rubi.add(rule42)

    def f43(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), IntegerQ(m), RationalQ(p), Unequal(m, 2), Unequal(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern43 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f43))
    rule43 = ReplacementRule(pattern43, lambda e, m, a, c, d, p, x, b : Int((a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule43)

    def f44(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), IntegerQ(m), RationalQ(p), Unequal(m, 2), Unequal(m, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern44 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f44))
    rule44 = ReplacementRule(pattern44, lambda e, m, a, c, d, p, x : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m), x))
    rubi.add(rule44)

    def f45(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), PositiveIntegerQ(m + p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern45 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f45))
    rule45 = ReplacementRule(pattern45, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule45)

    def f46(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), PositiveIntegerQ(m + p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern46 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f46))
    rule46 = ReplacementRule(pattern46, lambda e, m, a, c, d, p, x : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule46)

    def f47(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), NegativeIntegerQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern47 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f47))
    rule47 = ReplacementRule(pattern47, lambda e, m, a, c, d, p, x, b : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule47)

    def f48(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), NegativeIntegerQ(m + 2*p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern48 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f48))
    rule48 = ReplacementRule(pattern48, lambda e, m, a, c, d, p, x : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule48)

    def f49(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern49 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f49))
    rule49 = ReplacementRule(pattern49, lambda e, a, c, d, x, b : S(2)*e*Subst(Int(S(1)/(-b*e + S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule49)

    def f50(e, a, c, d, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern50 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f50))
    rule50 = ReplacementRule(pattern50, lambda e, a, c, d, x : S(2)*e*Subst(Int(S(1)/(S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule50)

    def f51(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Greater(p, 0), NonzeroQ(m + p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern51 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f51))
    rule51 = ReplacementRule(pattern51, lambda e, m, a, c, d, p, x, b : -c*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + p + S(1))))
    rubi.add(rule51)

    def f52(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m, p), Greater(p, 0), NonzeroQ(m + p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern52 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f52))
    rule52 = ReplacementRule(pattern52, lambda e, m, a, c, d, p, x : -c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/(e**S(2)*(m + p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + p + S(1))))
    rubi.add(rule52)

    def f53(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Greater(p, 0), NonzeroQ(m + 2*p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern53 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f53))
    rule53 = ReplacementRule(pattern53, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule53)

    def f54(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m, p), Greater(p, 0), NonzeroQ(m + 2*p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern54 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f54))
    rule54 = ReplacementRule(pattern54, lambda e, m, a, c, d, p, x : -S(2)*c*d*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e**S(2)*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule54)

    def f55(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Less(p, -1), Inequality(0, Less, m, LessEqual, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern55 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f55))
    rule55 = ReplacementRule(pattern55, lambda e, m, a, c, d, p, x, b : -(-b*e + S(2)*c*d)*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule55)

    def f56(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Inequality(0, Less, m, LessEqual, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern56 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f56))
    rule56 = ReplacementRule(pattern56, lambda e, m, a, c, d, p, x : d*(m + S(2)*p + S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - d*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*e*(p + S(1))))
    rubi.add(rule56)

    def f57(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern57 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f57))
    rule57 = ReplacementRule(pattern57, lambda e, m, a, c, d, p, x, b : -e**S(2)*(m + p)*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule57)

    def f58(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern58 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f58))
    rule58 = ReplacementRule(pattern58, lambda e, m, a, c, d, p, x : -e**S(2)*(m + p)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule58)

    def f59(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m), GreaterEqual(m, 1), NonzeroQ(m + 2*p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern59 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f59))
    rule59 = ReplacementRule(pattern59, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule59)

    def f60(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m), GreaterEqual(m, 1), NonzeroQ(m + 2*p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern60 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f60))
    rule60 = ReplacementRule(pattern60, lambda e, m, a, c, d, p, x : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule60)

    def f61(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m), Less(m, 0), NonzeroQ(m + p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern61 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f61))
    rule61 = ReplacementRule(pattern61, lambda e, m, a, c, d, p, x, b : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule61)

    def f62(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), RationalQ(m), Less(m, 0), NonzeroQ(m + p + 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern62 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f62))
    rule62 = ReplacementRule(pattern62, lambda e, m, a, c, d, p, x : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule62)

    def f63(e, m, c, x, p, b):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(e, x), FreeQ(m, x)])
    pattern63 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f63))
    rule63 = ReplacementRule(pattern63, lambda e, m, c, x, p, b : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p, x))
    rubi.add(rule63)

    def f64(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), PositiveQ(a), PositiveQ(d), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern64 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f64))
    rule64 = ReplacementRule(pattern64, lambda e, m, a, c, d, p, x : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule64)

    def f65(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern65 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f65))
    rule65 = ReplacementRule(pattern65, lambda e, m, a, c, d, p, x, b : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule65)

    def f66(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern66 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f66))
    rule66 = ReplacementRule(pattern66, lambda e, m, a, c, d, p, x : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule66)

    def f67(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern67 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f67))
    rule67 = ReplacementRule(pattern67, lambda e, a, c, d, x, b : b**S(2)*Int((d + e*x)/(a + b*x + c*x**S(2)), x)/(d**S(2)*(-S(4)*a*c + b**S(2))) - S(4)*b*c*Int(S(1)/(b + S(2)*c*x), x)/(d*(-S(4)*a*c + b**S(2))))
    rubi.add(rule67)

    def f68(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 3), NonzeroQ(p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern68 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f68))
    rule68 = ReplacementRule(pattern68, lambda e, m, a, c, d, p, x, b : S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule68)

    def f69(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern69 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f69))
    rule69 = ReplacementRule(pattern69, lambda e, m, a, c, d, p, x, b : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule69)

    def f70(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(m, p), Greater(p, 0), Less(m, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern70 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f70))
    rule70 = ReplacementRule(pattern70, lambda e, m, a, c, d, p, x, b : -b*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(d*e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule70)

    def f71(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(p), Greater(p, 0), RationalQ(m), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern71 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f71))
    rule71 = ReplacementRule(pattern71, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - d*p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(b*e*(m + S(2)*p + S(1))))
    rubi.add(rule71)

    def f72(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern72 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f72))
    rule72 = ReplacementRule(pattern72, lambda e, m, a, c, d, p, x, b : -d*e*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))) + d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))))
    rubi.add(rule72)

    def f73(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(p), Less(p, -1), RationalQ(m), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern73 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f73))
    rule73 = ReplacementRule(pattern73, lambda e, m, a, c, d, p, x, b : -S(2)*c*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule73)

    def f74(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern74 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f74))
    rule74 = ReplacementRule(pattern74, lambda e, a, c, d, x, b : S(4)*c*Subst(Int(S(1)/(-S(4)*a*c*e + b**S(2)*e + S(4)*c*e*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule74)

    def f75(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NegativeQ(c/(-4*a*c + b**2)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern75 = Pattern(Integral(S(1)/(sqrt(d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f75))
    rule75 = ReplacementRule(pattern75, lambda e, a, c, d, x, b : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(S(1)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule75)

    def f76(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NegativeQ(c/(-4*a*c + b**2)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern76 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f76))
    rule76 = ReplacementRule(pattern76, lambda e, a, c, d, x, b : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(x**S(2)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule76)

    def f77(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), EqQ(m**2, 1/4), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern77 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f77))
    rule77 = ReplacementRule(pattern77, lambda e, m, a, c, d, x, b : sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*Int((d + e*x)**m/sqrt(-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2))), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule77)

    def f78(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(m), Greater(m, 1), NonzeroQ(m + 2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern78 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f78))
    rule78 = ReplacementRule(pattern78, lambda e, m, a, c, d, p, x, b : S(2)*d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(2)*p + S(1))) + d**S(2)*(m + S(-1))*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p, x)/(b**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule78)

    def f79(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 3), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern79 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f79))
    rule79 = ReplacementRule(pattern79, lambda e, m, a, c, d, p, x, b : b**S(2)*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/(d**S(2)*(m + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*b*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(d*(m + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule79)

    def f80(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern80 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f80))
    rule80 = ReplacementRule(pattern80, lambda e, m, a, c, d, p, x, b : Subst(Int(x**m*(a - b**S(2)/(S(4)*c) + c*x**S(2)/e**S(2))**p, x), x, d + e*x)/e)
    rubi.add(rule80)

    def f81(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern81 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f81))
    rule81 = ReplacementRule(pattern81, lambda e, m, a, c, d, p, x, b : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule81)

    def f82(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern82 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f82))
    rule82 = ReplacementRule(pattern82, lambda e, m, a, c, d, p, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m, x), x))
    rubi.add(rule82)

    def f83(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), NiceSqrtQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern83 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f83), )
    def With83(e, a, c, d, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (c*d - e*(b/S(2) - q/S(2)))*Int(S(1)/(b/S(2) + c*x - q/S(2)), x)/q - (c*d - e*(b/S(2) + q/S(2)))*Int(S(1)/(b/S(2) + c*x + q/S(2)), x)/q
    rule83 = ReplacementRule(pattern83, lambda e, a, c, d, x, b : With83(e, a, c, d, x, b))
    rubi.add(rule83)

    def f84(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), NiceSqrtQ(-a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern84 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f84), )
    def With84(e, a, c, d, x):
        q = Rt(-a*c, S(2))
        return (-c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x + q), x) + (c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x - q), x)
    rule84 = ReplacementRule(pattern84, lambda e, a, c, d, x : With84(e, a, c, d, x))
    rubi.add(rule84)

    def f85(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern85 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f85))
    rule85 = ReplacementRule(pattern85, lambda e, a, c, d, x, b : e*Int((b + S(2)*c*x)/(a + b*x + c*x**S(2)), x)/(S(2)*c) + (-b*e + S(2)*c*d)*Int(S(1)/(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule85)

    def f86(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern86 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f86))
    rule86 = ReplacementRule(pattern86, lambda e, a, c, d, x : d*Int(S(1)/(a + c*x**S(2)), x) + e*Int(x/(a + c*x**S(2)), x))
    rubi.add(rule86)

    def f87(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern87 = Pattern(Integral(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f87))
    rule87 = ReplacementRule(pattern87, lambda e, a, c, d, x, b : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule87)

    def f88(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern88 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f88))
    rule88 = ReplacementRule(pattern88, lambda e, a, c, d, x : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule88)

    def f89(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), IntegerQ(m), Greater(m, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern89 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f89))
    rule89 = ReplacementRule(pattern89, lambda e, m, a, c, d, x, b : Int(PolynomialDivide((d + e*x)**m, a + b*x + c*x**S(2), x), x))
    rubi.add(rule89)

    def f90(e, m, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), IntegerQ(m), Greater(m, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern90 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f90))
    rule90 = ReplacementRule(pattern90, lambda e, m, a, c, d, x : Int(PolynomialDivide((d + e*x)**m, a + c*x**S(2), x), x))
    rubi.add(rule90)

    def f91(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Greater(m, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern91 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f91))
    rule91 = ReplacementRule(pattern91, lambda e, m, a, c, d, x, b : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + e*x*(-b*e + S(2)*c*d), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule91)

    def f92(e, m, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m), Greater(m, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern92 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f92))
    rule92 = ReplacementRule(pattern92, lambda e, m, a, c, d, x : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + S(2)*c*d*e*x, x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule92)

    def f93(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern93 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f93))
    rule93 = ReplacementRule(pattern93, lambda e, a, c, d, x, b : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((-b*e + c*d - c*e*x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule93)

    def f94(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern94 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f94))
    rule94 = ReplacementRule(pattern94, lambda e, a, c, d, x : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((c*d - c*e*x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule94)

    def f95(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern95 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f95))
    rule95 = ReplacementRule(pattern95, lambda e, a, c, d, x, b : S(2)*e*Subst(Int(S(1)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule95)

    def f96(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern96 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f96))
    rule96 = ReplacementRule(pattern96, lambda e, a, c, d, x : S(2)*e*Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule96)

    def f97(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern97 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f97))
    rule97 = ReplacementRule(pattern97, lambda e, m, a, c, d, x, b : e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(-b*e + c*d - c*e*x, x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule97)

    def f98(e, m, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern98 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f98))
    rule98 = ReplacementRule(pattern98, lambda e, m, a, c, d, x : c*Int((d - e*x)*(d + e*x)**(m + S(1))/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)) + e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule98)

    def f99(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern99 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f99))
    rule99 = ReplacementRule(pattern99, lambda e, m, a, c, d, x, b : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule99)

    def f100(e, m, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern100 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f100))
    rule100 = ReplacementRule(pattern100, lambda e, m, a, c, d, x : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule100)

    def f101(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern101 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2), x_),CustomConstraint(f101))
    rule101 = ReplacementRule(pattern101, lambda e, a, c, d, x, b : -S(2)*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule101)

    def f102(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern102 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2), x_),CustomConstraint(f102))
    rule102 = ReplacementRule(pattern102, lambda e, a, c, d, x : (-a*e + c*d*x)/(a*c*sqrt(a + c*x**S(2))))
    rubi.add(rule102)

    def f103(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Less(p, -1), Unequal(p, -3/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern103 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f103))
    rule103 = ReplacementRule(pattern103, lambda e, a, c, d, p, x, b : -(S(2)*p + S(3))*(-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule103)

    def f104(e, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Less(p, -1), Unequal(p, -3/2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern104 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f104))
    rule104 = ReplacementRule(pattern104, lambda e, a, c, d, p, x : d*(S(2)*p + S(3))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule104)

    def f105(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern105 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f105))
    rule105 = ReplacementRule(pattern105, lambda e, a, c, d, p, x, b : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule105)

    def f106(e, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern106 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f106))
    rule106 = ReplacementRule(pattern106, lambda e, a, c, d, p, x : d*Int((a + c*x**S(2))**p, x) + e*(a + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule106)

    def f107(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(a*e + b*d), ZeroQ(b*e + c*d), PositiveIntegerQ(m - p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern107 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f107))
    rule107 = ReplacementRule(pattern107, lambda e, m, a, c, d, p, x, b : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m - p)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule107)

    def f108(e, m, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-b*e + c*d), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Equal(m**2, 1/4), NegativeQ(c), RationalQ(b), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern108 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f108))
    rule108 = ReplacementRule(pattern108, lambda e, m, c, d, x, b : Int((d + e*x)**m/(sqrt(b*x)*sqrt(S(1) + c*x/b)), x))
    rubi.add(rule108)

    def f109(e, m, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-b*e + c*d), NonzeroQ(-b*e + 2*c*d), RationalQ(m), Equal(m**2, 1/4), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern109 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f109))
    rule109 = ReplacementRule(pattern109, lambda e, m, c, d, x, b : sqrt(x)*sqrt(b + c*x)*Int((d + e*x)**m/(sqrt(x)*sqrt(b + c*x)), x)/sqrt(b*x + c*x**S(2)))
    rubi.add(rule109)

    def f110(m, a, c, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(m**2 - 1/4), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern110 = Pattern(Integral(x_**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f110))
    rule110 = ReplacementRule(pattern110, lambda m, a, c, x, b : S(2)*Subst(Int(x**(S(2)*m + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule110)

    def f111(e, m, a, c, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(m**2 - 1/4), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(e, x)])
    pattern111 = Pattern(Integral((e_*x_)**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f111))
    rule111 = ReplacementRule(pattern111, lambda e, m, a, c, x, b : x**(-m)*(e*x)**m*Int(x**m/sqrt(a + b*x + c*x**S(2)), x))
    rubi.add(rule111)

    def f112(e, m, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m**2 - 1/4), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern112 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f112))
    rule112 = ReplacementRule(pattern112, lambda e, m, a, c, d, x, b : S(2)*sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*(S(2)*c*(d + e*x)/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))))**(-m)*(d + e*x)**m*Rt(-S(4)*a*c + b**S(2), S(2))*Subst(Int((S(2)*e*x**S(2)*Rt(-S(4)*a*c + b**S(2), S(2))/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(S(2))*sqrt((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))/Rt(-S(4)*a*c + b**S(2), S(2)))/S(2))/(c*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule112)

    def f113(e, m, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), ZeroQ(m**2 - 1/4), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern113 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f113))
    rule113 = ReplacementRule(pattern113, lambda e, m, a, c, d, x : S(2)*a*(c*(d + e*x)/(-a*e*Rt(-c/a, S(2)) + c*d))**(-m)*sqrt(S(1) + c*x**S(2)/a)*(d + e*x)**m*Rt(-c/a, S(2))*Subst(Int((S(2)*a*e*x**S(2)*Rt(-c/a, S(2))/(-a*e*Rt(-c/a, S(2)) + c*d) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(-x*Rt(-c/a, S(2))/S(2) + S(1)/2))/(c*sqrt(a + c*x**S(2))))
    rubi.add(rule113)

    def f114(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Equal(m + 2*p + 2, 0), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern114 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f114))
    rule114 = ReplacementRule(pattern114, lambda e, m, a, c, d, p, x, b : p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule114)

    def f115(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Equal(m + 2*p + 2, 0), Greater(p, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern115 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f115))
    rule115 = ReplacementRule(pattern115, lambda e, m, a, c, d, p, x : -S(2)*a*c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*a*e + S(2)*c*d*x)/(S(2)*(m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule115)

    def f116(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Equal(m + 2*p + 2, 0), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern116 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f116))
    rule116 = ReplacementRule(pattern116, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*(S(2)*p + S(3))*(a*e**S(2) - b*d*e + c*d**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule116)

    def f117(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Equal(m + 2*p + 2, 0), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern117 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f117))
    rule117 = ReplacementRule(pattern117, lambda e, m, a, c, d, p, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) + (S(2)*p + S(3))*(a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule117)

    def f118(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern118 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f118))
    rule118 = ReplacementRule(pattern118, lambda e, a, c, d, x, b : -S(2)*Subst(Int(S(1)/(S(4)*a*e**S(2) - S(4)*b*d*e + S(4)*c*d**S(2) - x**S(2)), x), x, (S(2)*a*e - b*d - x*(-b*e + S(2)*c*d))/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule118)

    def f119(e, a, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern119 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f119))
    rule119 = ReplacementRule(pattern119, lambda e, a, c, d, x : -Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - x**S(2)), x), x, (a*e - c*d*x)/sqrt(a + c*x**S(2))))
    rubi.add(rule119)

    def f120(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern120 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f120))
    rule120 = ReplacementRule(pattern120, lambda e, m, a, c, d, p, x, b : -((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))**(-p)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), -S(4)*c*(d + e*x)*Rt(-S(4)*a*c + b**S(2), S(2))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))/((m + S(1))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))))
    rubi.add(rule120)

    def f121(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern121 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f121))
    rule121 = ReplacementRule(pattern121, lambda e, m, a, c, d, p, x : ((c*d + e*Rt(-a*c, S(2)))*(c*x + Rt(-a*c, S(2)))/((c*d - e*Rt(-a*c, S(2)))*(c*x - Rt(-a*c, S(2)))))**(-p)*(a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*x + Rt(-a*c, S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), S(2)*c*(d + e*x)*Rt(-a*c, S(2))/((c*d - e*Rt(-a*c, S(2)))*(-c*x + Rt(-a*c, S(2)))))/((m + S(1))*(c*d + e*Rt(-a*c, S(2)))))
    rubi.add(rule121)

    def f122(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 3), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern122 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f122))
    rule122 = ReplacementRule(pattern122, lambda e, m, a, c, d, p, x, b : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule122)

    def f123(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 3), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern123 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f123))
    rule123 = ReplacementRule(pattern123, lambda e, m, a, c, d, p, x : -d*m*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))))
    rubi.add(rule123)

    def f124(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern124 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f124))
    rule124 = ReplacementRule(pattern124, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + (-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule124)

    def f125(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern125 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f125))
    rule125 = ReplacementRule(pattern125, lambda e, m, a, c, d, p, x : c*d*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule125)

    def f126(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Greater(p, 0), NonzeroQ(m + 1), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern126 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f126))
    rule126 = ReplacementRule(pattern126, lambda e, m, a, c, d, p, x, b : -p*Int((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule126)

    def f127(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Greater(p, 0), NonzeroQ(m + 1), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern127 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f127))
    rule127 = ReplacementRule(pattern127, lambda e, m, a, c, d, p, x : -S(2)*c*p*Int(x*(a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e*(m + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule127)

    def f128(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Greater(p, 0), NonzeroQ(m + 2*p + 1), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern128 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f128))
    rule128 = ReplacementRule(pattern128, lambda e, m, a, c, d, p, x, b : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d), x), x)/(e*(m + S(2)*p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule128)

    def f129(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Greater(p, 0), NonzeroQ(m + 2*p + 1), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern129 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f129))
    rule129 = ReplacementRule(pattern129, lambda e, m, a, c, d, p, x : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*e - c*d*x, x), x)/(e*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule129)

    def f130(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Less(p, -1), Greater(m, 0), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern130 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f130))
    rule130 = ReplacementRule(pattern130, lambda e, m, a, c, d, p, x, b : (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(b*e*m + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(m + S(2)*p + S(3))), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule130)

    def f131(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 0), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern131 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f131))
    rule131 = ReplacementRule(pattern131, lambda e, m, a, c, d, p, x : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(d*(S(2)*p + S(3)) + e*x*(m + S(2)*p + S(3))), x)/(S(2)*a*(p + S(1))))
    rubi.add(rule131)

    def f132(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern132 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f132))
    rule132 = ReplacementRule(pattern132, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*c*d**S(2)*(S(2)*p + S(3)) + e*x*(b*e - S(2)*c*d)*(m + S(2)*p + S(2)) + e*(S(2)*a*e*(m + S(-1)) + b*d*(-m + S(2)*p + S(4))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule132)

    def f133(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern133 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f133))
    rule133 = ReplacementRule(pattern133, lambda e, m, a, c, d, p, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(a*e**S(2)*(m + S(-1)) - c*d**S(2)*(S(2)*p + S(3)) - c*d*e*x*(m + S(2)*p + S(2)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule133)

    def f134(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Less(p, -1), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern134 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f134))
    rule134 = ReplacementRule(pattern134, lambda e, m, a, c, d, p, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*e - b**S(2)*e + b*c*d + c*x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3)) - c*e*x*(-b*e + S(2)*c*d)*(m + S(2)*p + S(4)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule134)

    def f135(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Less(p, -1), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern135 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f135))
    rule135 = ReplacementRule(pattern135, lambda e, m, a, c, d, p, x : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(a*e + c*d*x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(a*e**S(2)*(m + S(2)*p + S(3)) + c*d**S(2)*(S(2)*p + S(3)) + c*d*e*x*(m + S(2)*p + S(4)), x), x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule135)

    def f136(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), If(RationalQ(m), Greater(m, 1), SumSimplerQ(m, -2)), NonzeroQ(m + 2*p + 1), IntQuadraticQ(a, b, c, d, e, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern136 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f136))
    rule136 = ReplacementRule(pattern136, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p*Simp(c*d**S(2)*(m + S(2)*p + S(1)) + e*x*(m + p)*(-b*e + S(2)*c*d) - e*(a*e*(m + S(-1)) + b*d*(p + S(1))), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule136)

    def f137(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), If(RationalQ(m), Greater(m, 1), SumSimplerQ(m, -2)), NonzeroQ(m + 2*p + 1), IntQuadraticQ(a, 0, c, d, e, m, p, x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern137 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f137))
    rule137 = ReplacementRule(pattern137, lambda e, m, a, c, d, p, x : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-2))*Simp(-a*e**S(2)*(m + S(-1)) + c*d**S(2)*(m + S(2)*p + S(1)) + S(2)*c*d*e*x*(m + p), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule137)

    def f138(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern138 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f138))
    rule138 = ReplacementRule(pattern138, lambda e, m, a, c, d, p, x, b : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*e*(m + p + S(2)) + c*d*(m + S(1)) - c*e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule138)

    def f139(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern139 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_),CustomConstraint(f139))
    rule139 = ReplacementRule(pattern139, lambda e, m, a, c, d, p, x : c*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(d*(m + S(1)) - e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule139)

    def f140(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-b*e + 2*c*d), ZeroQ(-3*a*c*e**2 + b**2*e**2 - b*c*d*e + c**2*d**2), PosQ(c*e**2*(-b*e + 2*c*d)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern140 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_),CustomConstraint(f140), )
    def With140(e, a, c, d, x, b):
        q = Rt(S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) + S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x - q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule140 = ReplacementRule(pattern140, lambda e, a, c, d, x, b : With140(e, a, c, d, x, b))
    rubi.add(rule140)

    def f141(e, a, c, d, x):
        return functools.reduce(operator.and_, [ ZeroQ(-3*a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern141 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f141), )
    def With141(e, a, c, d, x):
        q = Rt(S(6)*c**S(2)*e**S(2)/d**S(2), S(3))
        return -sqrt(S(3))*c*e*ArcTan(S(2)*sqrt(S(3))*c*(d - e*x)/(S(3)*d*q*(a + c*x**S(2))**(S(1)/3)) + sqrt(S(3))/S(3))/(d**S(2)*q**S(2)) - S(3)*c*e*log(d + e*x)/(S(2)*d**S(2)*q**S(2)) + S(3)*c*e*log(c*d - c*e*x - d*q*(a + c*x**S(2))**(S(1)/3))/(S(2)*d**S(2)*q**S(2))
    rule141 = ReplacementRule(pattern141, lambda e, a, c, d, x : With141(e, a, c, d, x))
    rubi.add(rule141)

    def f142(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-b*e + 2*c*d), ZeroQ(-3*a*c*e**2 + b**2*e**2 - b*c*d*e + c**2*d**2), NegQ(c*e**2*(-b*e + 2*c*d)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern142 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_),CustomConstraint(f142), )
    def With142(e, a, c, d, x, b):
        q = Rt(-S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) - S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x + q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule142 = ReplacementRule(pattern142, lambda e, a, c, d, x, b : With142(e, a, c, d, x, b))
    rubi.add(rule142)

    def f143(e, a, c, d, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(9*a*c*e**2 - 2*b**2*e**2 - b*c*d*e + c**2*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern143 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_),CustomConstraint(f143), )
    def With143(e, a, c, d, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)*Int(S(1)/((d + e*x)*(b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule143 = ReplacementRule(pattern143, lambda e, a, c, d, x, b : With143(e, a, c, d, x, b))
    rubi.add(rule143)

    def f144(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern144 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/4)*(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f144))
    rule144 = ReplacementRule(pattern144, lambda e, a, c, d, x : d*Int(S(1)/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule144)

    def f145(e, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern145 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/4)*(d_ + x_*WC('e', S(1)))), x_),CustomConstraint(f145))
    rule145 = ReplacementRule(pattern145, lambda e, a, c, d, x : d*Int(S(1)/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule145)

    def f146(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ PositiveQ(4*a - b**2/c), IntegerQ(4*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern146 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f146))
    rule146 = ReplacementRule(pattern146, lambda e, a, c, d, p, x, b : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p/Simp(-b*e + S(2)*c*d + e*x, x), x), x, b + S(2)*c*x))
    rubi.add(rule146)

    def f147(e, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ IntegerQ(4*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern147 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f147))
    rule147 = ReplacementRule(pattern147, lambda e, a, c, d, p, x, b : (-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))**(-p)*(a + b*x + c*x**S(2))**p*Int((-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2)))**p/(d + e*x), x))
    rubi.add(rule147)

    def f148(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), PositiveQ(a), NegativeQ(c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern148 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f148))
    rule148 = ReplacementRule(pattern148, lambda e, m, a, c, d, p, x : Int((d + e*x)**m*(-x*Rt(-c, S(2)) + Rt(a, S(2)))**p*(x*Rt(-c, S(2)) + Rt(a, S(2)))**p, x))
    rubi.add(rule148)

    def f149(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), NegativeIntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern149 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f149), )
    def With149(e, m, a, c, d, p, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -(e*(b + S(2)*c*x - q)/(S(2)*c*(d + e*x)))**(-p)*(e*(b + S(2)*c*x + q)/(S(2)*c*(d + e*x)))**(-p)*(a + b*x + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x*(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule149 = ReplacementRule(pattern149, lambda e, m, a, c, d, p, x, b : With149(e, m, a, c, d, p, x, b))
    rubi.add(rule149)

    def f150(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), NegativeIntegerQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern150 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f150), )
    def With150(e, m, a, c, d, p, x):
        q = Rt(-a*c, S(2))
        return -(e*(c*x + q)/(c*(d + e*x)))**(-p)*(-e*(-c*x + q)/(c*(d + e*x)))**(-p)*(a + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*q/c) + S(1), x)**p*Simp(-x*(d + e*q/c) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule150 = ReplacementRule(pattern150, lambda e, m, a, c, d, p, x : With150(e, m, a, c, d, p, x))
    rubi.add(rule150)

    def f151(e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(-b*e + 2*c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern151 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f151), )
    def With151(e, m, a, c, d, p, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (-(d + e*x)/(d - e*(b - q)/(S(2)*c)) + S(1))**(-p)*(-(d + e*x)/(d - e*(b + q)/(S(2)*c)) + S(1))**(-p)*(a + b*x + c*x**S(2))**p*Subst(Int(x**m*Simp(-x/(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x/(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, d + e*x)/e
    rule151 = ReplacementRule(pattern151, lambda e, m, a, c, d, p, x, b : With151(e, m, a, c, d, p, x, b))
    rubi.add(rule151)

    def f152(e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern152 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_),CustomConstraint(f152), )
    def With152(e, m, a, c, d, p, x):
        q = Rt(-a*c, S(2))
        return (a + c*x**S(2))**p*(-(d + e*x)/(d - e*q/c) + S(1))**(-p)*(-(d + e*x)/(d + e*q/c) + S(1))**(-p)*Subst(Int(x**m*Simp(-x/(d - e*q/c) + S(1), x)**p*Simp(-x/(d + e*q/c) + S(1), x)**p, x), x, d + e*x)/e
    rule152 = ReplacementRule(pattern152, lambda e, m, a, c, d, p, x : With152(e, m, a, c, d, p, x))
    rubi.add(rule152)

    def f153(u, e, m, a, c, d, p, x, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern153 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f153))
    rule153 = ReplacementRule(pattern153, lambda u, e, m, a, c, d, p, x, b : Subst(Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule153)

    def f154(u, e, m, a, c, d, p, x):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern154 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_),CustomConstraint(f154))
    rule154 = ReplacementRule(pattern154, lambda u, e, m, a, c, d, p, x : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule154)

    def f155(e, a, d, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(n), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern155 = Pattern(Integral(x_**WC('n', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f155))
    rule155 = ReplacementRule(pattern155, lambda e, a, d, c, n, p, x : d*Int(x**n*(a + c*x**S(2))**p, x) + e*Int(x**(n + S(1))*(a + c*x**S(2))**p, x))
    rubi.add(rule155)

    def f156(e, m, g, a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern156 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f156))
    rule156 = ReplacementRule(pattern156, lambda e, m, g, a, c, d, f, x, b : (f + g*x)*Int((d + e*x)**m, x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule156)

    def f157(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern157 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f157))
    rule157 = ReplacementRule(pattern157, lambda e, m, g, a, c, d, p, f, x, b : -f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule157)

    def f158(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), RationalQ(m, p), Less(p, -1), Greater(m, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern158 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f158))
    rule158 = ReplacementRule(pattern158, lambda e, m, g, a, c, d, p, f, x, b : -e*g*m*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule158)

    def f159(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern159 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f159))
    rule159 = ReplacementRule(pattern159, lambda e, m, g, a, c, d, p, f, x, b : e*f*g*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))*(-d*g + e*f)) - f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule159)

    def f160(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), RationalQ(m), Less(m, -1), NonzeroQ(2*p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern160 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f160))
    rule160 = ReplacementRule(pattern160, lambda e, m, g, a, c, d, p, f, x, b : -g*(S(2)*p + S(1))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule160)

    def f161(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), RationalQ(m), Less(m, -1), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern161 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f161))
    rule161 = ReplacementRule(pattern161, lambda e, m, g, a, c, d, p, f, x, b : -g*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-d*g + e*f)) + S(2)*f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(1))*(-d*g + e*f)))
    rubi.add(rule161)

    def f162(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), PositiveIntegerQ(m), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern162 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f162))
    rule162 = ReplacementRule(pattern162, lambda e, m, g, a, c, d, p, f, x, b : -b*m*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*f*(m + S(2)*p + S(2))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule162)

    def f163(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), ZeroQ(-b*g + 2*c*f), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern163 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f163))
    rule163 = ReplacementRule(pattern163, lambda e, m, g, a, c, d, p, f, x, b : (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(2))) + (S(2)*p + S(1))*(-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(2)*p + S(2))))
    rubi.add(rule163)

    def f164(e, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), NonzeroQ(-b*g + 2*c*f), NonzeroQ(-b*e + 2*c*d), RationalQ(p), Less(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern164 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f164))
    rule164 = ReplacementRule(pattern164, lambda e, g, a, c, d, p, f, x, b : (-b*g + S(2)*c*f)*Int((a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(d + e*x), x)/(-b*e + S(2)*c*d))
    rubi.add(rule164)

    def f165(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), NonzeroQ(-b*g + 2*c*f), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern165 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f165))
    rule165 = ReplacementRule(pattern165, lambda e, m, g, a, c, d, p, f, x, b : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule165)

    def f166(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), NonzeroQ(-b*g + 2*c*f), NonzeroQ(-b*e + 2*c*d), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern166 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f166))
    rule166 = ReplacementRule(pattern166, lambda e, m, g, a, c, d, p, f, x, b : (-b*g + S(2)*c*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d))
    rubi.add(rule166)

    def f167(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), NonzeroQ(-b*g + 2*c*f), NonzeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 2), NonzeroQ(m + 2*p + 3), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern167 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f167))
    rule167 = ReplacementRule(pattern167, lambda e, m, g, a, c, d, p, f, x, b : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))*(-b*e + S(2)*c*d)) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*d*(S(2)*p + S(1))))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule167)

    def f168(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), NonzeroQ(-b*g + 2*c*f), NonzeroQ(-b*e + 2*c*d), NonzeroQ(m + 2*p + 2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern168 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f168))
    rule168 = ReplacementRule(pattern168, lambda e, m, g, a, c, d, p, f, x, b : g*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*c*e*(m + S(2)*p + S(2))) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*(S(2)*d*p + d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*e*(m + S(2)*p + S(2))))
    rubi.add(rule168)

    def f169(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern169 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f169))
    rule169 = ReplacementRule(pattern169, lambda e, m, a, g, c, d, p, n, f, x, b : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule169)

    def f170(e, m, g, a, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern170 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f170))
    rule170 = ReplacementRule(pattern170, lambda e, m, g, a, d, c, n, p, f, x, b : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule170)

    def f171(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern171 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f171))
    rule171 = ReplacementRule(pattern171, lambda e, m, g, a, d, c, n, p, f, x : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule171)

    def f172(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 - b*d*e + c*d**2), NegativeIntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x), FreeQ(p, x)])
    pattern172 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f172))
    rule172 = ReplacementRule(pattern172, lambda e, m, a, g, d, c, n, p, f, x, b : d**m*e**m*Int((f + g*x)**n*(a*e + c*d*x)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule172)

    def f173(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), NegativeIntegerQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x), FreeQ(p, x)])
    pattern173 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f173))
    rule173 = ReplacementRule(pattern173, lambda e, m, g, a, d, c, n, p, f, x : d**m*e**m*Int((a + c*x**S(2))**(m + p)*(f + g*x)**n*(a*e + c*d*x)**(-m), x))
    rubi.add(rule173)

    def f174(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(e*(p + 1)*(-b*g + 2*c*f) + m*(c*e*f + g*(-b*e + c*d))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern174 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f174))
    rule174 = ReplacementRule(pattern174, lambda e, m, g, a, c, d, p, f, x, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule174)

    def f175(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(2*e*f*(p + 1) + m*(d*g + e*f)), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern175 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f175))
    rule175 = ReplacementRule(pattern175, lambda e, m, g, a, c, d, p, f, x : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule175)

    def f176(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern176 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f176))
    rule176 = ReplacementRule(pattern176, lambda e, m, g, a, c, d, p, f, x, b : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule176)

    def f177(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern177 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f177))
    rule177 = ReplacementRule(pattern177, lambda e, m, g, a, c, d, p, f, x : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule177)

    def f178(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), SumSimplerQ(p, 1), SumSimplerQ(m, -1), NonzeroQ(p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern178 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f178))
    rule178 = ReplacementRule(pattern178, lambda e, m, g, a, c, d, p, f, x, b : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule178)

    def f179(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), SumSimplerQ(p, 1), SumSimplerQ(m, -1), NonzeroQ(p + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern179 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f179))
    rule179 = ReplacementRule(pattern179, lambda e, m, g, a, c, d, p, f, x : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule179)

    def f180(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(m + p + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern180 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f180))
    rule180 = ReplacementRule(pattern180, lambda e, m, g, a, c, d, p, f, x, b : (d + e*x)**m*(d*g - e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule180)

    def f181(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), NonzeroQ(m + p + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern181 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f181))
    rule181 = ReplacementRule(pattern181, lambda e, m, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g - e*f)/(S(2)*c*d*(m + p + S(1))) + (S(2)*c*e*f*(p + S(1)) + m*(c*d*g + c*e*f))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*c*d*e*(m + p + S(1))))
    rubi.add(rule181)

    def f182(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern182 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f182))
    rule182 = ReplacementRule(pattern182, lambda e, m, g, a, c, d, p, f, x, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(c*e*(m + S(2)*p + S(2))))
    rubi.add(rule182)

    def f183(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern183 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f183))
    rule183 = ReplacementRule(pattern183, lambda e, m, g, a, c, d, p, f, x : (S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/(e*(m + S(2)*p + S(2))) + g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule183)

    def f184(a, g, c, x, p, f):
        return functools.reduce(operator.and_, [ ZeroQ(a*g**2 + c*f**2), RationalQ(p), Less(p, -2), FreeQ(a, x), FreeQ(c, x), FreeQ(f, x), FreeQ(g, x)])
    pattern184 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_),CustomConstraint(f184))
    rule184 = ReplacementRule(pattern184, lambda a, g, c, x, p, f : x**S(2)*(a + c*x**S(2))**(p + S(1))*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int(x*(a + c*x**S(2))**(p + S(1))*Simp(S(2)*a*g - c*f*x*(S(2)*p + S(5)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule184)

    def f185(a, g, c, x, p, f):
        return functools.reduce(operator.and_, [ ZeroQ(a*g**2 + c*f**2), FreeQ(a, x), FreeQ(c, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern185 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_),CustomConstraint(f185))
    rule185 = ReplacementRule(pattern185, lambda a, g, c, x, p, f : -f**S(2)*Int((a + c*x**S(2))**(p + S(1))/(f - g*x), x)/c + Int((a + c*x**S(2))**(p + S(1))*(f + g*x), x)/c)
    rubi.add(rule185)

    def f186(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), IntegersQ(m, n), RationalQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern186 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f186))
    rule186 = ReplacementRule(pattern186, lambda e, m, a, g, d, c, n, p, f, x, b : Int((f + g*x)**n*(a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule186)

    def f187(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), IntegersQ(m, n), RationalQ(p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern187 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f187))
    rule187 = ReplacementRule(pattern187, lambda e, m, g, a, d, c, n, p, f, x : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m)*(f + g*x)**n, x))
    rubi.add(rule187)

    def f188(e, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), PositiveIntegerQ(n), NegativeIntegerQ(n + 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern188 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f188))
    rule188 = ReplacementRule(pattern188, lambda e, a, g, d, c, n, p, f, x, b : -(f + g*x)**n*(a*(-b*e + S(2)*c*d) + c*x*(-S(2)*a*e + b*d))*(a + b*x + c*x**S(2))**p/(d*e*p*(-S(4)*a*c + b**S(2))) - Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*Simp(-S(2)*a*c*(d*g*n - e*f*(S(2)*p + S(1))) + b*(a*e*g*n - c*d*f*(S(2)*p + S(1))) - c*g*x*(-S(2)*a*e + b*d)*(n + S(2)*p + S(1)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))))
    rubi.add(rule188)

    def f189(e, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), PositiveIntegerQ(n), NegativeIntegerQ(n + 2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern189 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f189))
    rule189 = ReplacementRule(pattern189, lambda e, g, a, d, c, n, p, f, x : (a + c*x**S(2))**p*(d - e*x)*(f + g*x)**n/(S(2)*d*e*p) - Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*Simp(d*g*n - e*f*(S(2)*p + S(1)) - e*g*x*(n + S(2)*p + S(1)), x), x)/(S(2)*d*e*p))
    rubi.add(rule189)

    def f190(e, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), NegativeIntegerQ(n), NegativeIntegerQ(n + 2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern190 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f190))
    rule190 = ReplacementRule(pattern190, lambda e, a, g, d, c, n, p, f, x, b : -(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p*(a*c*d*(-b*g + S(2)*c*f) - a*e*(S(2)*a*c*g - b**S(2)*g + b*c*f) + c*x*(-a*e*(-b*g + S(2)*c*f) + c*d*(-S(2)*a*g + b*f)))/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**p*Simp(S(2)*a*c*(a*e*g**S(2)*(n + S(2)*p + S(1)) + c*f*(-d*g*n + S(2)*e*f*p + e*f)) + b**S(2)*g*(-a*e*g*(n + p + S(1)) + c*d*f*p) + b*c*(a*g*(d*g*(n + S(1)) + e*f*(n - S(2)*p)) - c*d*f**S(2)*(S(2)*p + S(1))) + c*g*x*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f))*(n + S(2)*p + S(2)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))))
    rubi.add(rule190)

    def f191(e, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), NegativeIntegerQ(n), NegativeIntegerQ(n + 2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern191 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f191))
    rule191 = ReplacementRule(pattern191, lambda e, g, a, d, c, n, p, f, x : (a + c*x**S(2))**p*(f + g*x)**(n + S(1))*(-a*e*g + c*d*f - c*x*(d*g + e*f))/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))) + Int((a + c*x**S(2))**p*(f + g*x)**n*Simp(a*e*g**S(2)*(n + S(2)*p + S(1)) - c*f*(d*g*n - e*(S(2)*f*p + f)) + c*g*x*(d*g + e*f)*(n + S(2)*p + S(2)), x), x)/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))))
    rubi.add(rule191)

    def f192(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), ZeroQ(-b*e*g + c*d*g + c*e*f), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern192 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f192))
    rule192 = ReplacementRule(pattern192, lambda e, m, a, g, d, c, n, p, f, x, b : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))))
    rubi.add(rule192)

    def f193(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), ZeroQ(d*g + e*f), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern193 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f193))
    rule193 = ReplacementRule(pattern193, lambda e, m, g, a, d, c, n, p, f, x : -e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule193)

    def f194(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), ZeroQ(m - n - 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern194 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f194))
    rule194 = ReplacementRule(pattern194, lambda e, m, a, g, c, d, p, n, f, x, b : -e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule194)

    def f195(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), ZeroQ(m - n - 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern195 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f195))
    rule195 = ReplacementRule(pattern195, lambda e, m, g, a, c, d, p, n, f, x : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(n + S(1))*(d*g + e*f)))
    rubi.add(rule195)

    def f196(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n, p), Greater(p, 0), Less(n, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern196 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f196))
    rule196 = ReplacementRule(pattern196, lambda e, m, a, g, c, d, p, n, f, x, b : c*m*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*g*(n + S(1))) + (d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(n + S(1))))
    rubi.add(rule196)

    def f197(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n, p), Greater(p, 0), Less(n, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern197 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f197))
    rule197 = ReplacementRule(pattern197, lambda e, m, g, a, c, d, p, n, f, x : c*m*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(1)), x)/(e*g*(n + S(1))) + (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(n + S(1))))
    rubi.add(rule197)

    def f198(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n, p), Greater(p, 0), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x)])
    pattern198 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f198))
    rule198 = ReplacementRule(pattern198, lambda e, m, a, g, d, c, n, p, f, x, b : -(d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(m - n + S(-1))) - m*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**(m + S(1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*g*(m - n + S(-1))))
    rubi.add(rule198)

    def f199(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n, p), Greater(p, 0), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x)])
    pattern199 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f199))
    rule199 = ReplacementRule(pattern199, lambda e, m, g, a, d, c, n, p, f, x : -c*m*(d*g + e*f)*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**n, x)/(e**S(2)*g*(m - n + S(-1))) - (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(m - n + S(-1))))
    rubi.add(rule199)

    def f200(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n, p), Less(p, -1), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern200 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f200))
    rule200 = ReplacementRule(pattern200, lambda e, m, a, g, d, c, n, p, f, x, b : -e*g*n*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule200)

    def f201(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n, p), Less(p, -1), Greater(n, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern201 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f201))
    rule201 = ReplacementRule(pattern201, lambda e, m, g, a, d, c, n, p, f, x : -e*g*n*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(p + S(1))))
    rubi.add(rule201)

    def f202(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n, p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x)])
    pattern202 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f202))
    rule202 = ReplacementRule(pattern202, lambda e, m, a, g, d, c, n, p, f, x, b : e**S(2)*g*(m - n + S(-2))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-b*e*g + c*d*g + c*e*f)) + e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule202)

    def f203(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n, p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x)])
    pattern203 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f203))
    rule203 = ReplacementRule(pattern203, lambda e, m, g, a, d, c, n, p, f, x : e**S(2)*g*(m - n + S(-2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(c*(p + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(p + S(1))*(d*g + e*f)))
    rubi.add(rule203)

    def f204(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n), Greater(n, 0), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern204 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f204))
    rule204 = ReplacementRule(pattern204, lambda e, m, a, g, d, c, n, p, f, x, b : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))) - n*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*e*(m - n + S(-1))))
    rubi.add(rule204)

    def f205(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n), Greater(n, 0), NonzeroQ(m - n - 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern205 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f205))
    rule205 = ReplacementRule(pattern205, lambda e, m, g, a, d, c, n, p, f, x : -n*(d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/(e*(m - n + S(-1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule205)

    def f206(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p), RationalQ(n), Less(n, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern206 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f206))
    rule206 = ReplacementRule(pattern206, lambda e, m, a, g, c, d, p, n, f, x, b : -c*e*(m - n + S(-2))*Int((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/((n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule206)

    def f207(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p), RationalQ(n), Less(n, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern207 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f207))
    rule207 = ReplacementRule(pattern207, lambda e, m, g, a, c, d, p, n, f, x : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/((n + S(1))*(c*d*g + c*e*f)) - e*(m - n + S(-2))*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1)), x)/((n + S(1))*(d*g + e*f)))
    rubi.add(rule207)

    def f208(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern208 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/((x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f208))
    rule208 = ReplacementRule(pattern208, lambda e, a, g, c, d, f, x, b : S(2)*e**S(2)*Subst(Int(S(1)/(-b*e*g + c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule208)

    def f209(e, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern209 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_),CustomConstraint(f209))
    rule209 = ReplacementRule(pattern209, lambda e, g, a, c, d, f, x : S(2)*e**S(2)*Subst(Int(S(1)/(c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule209)

    def f210(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p - 1), ZeroQ(b*e*g*(n + 1) - c*d*g*(2*n + p + 3) + c*e*f*(p + 1)), NonzeroQ(n + p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern210 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f210))
    rule210 = ReplacementRule(pattern210, lambda e, m, a, g, d, c, n, p, f, x, b : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule210)

    def f211(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p - 1), ZeroQ(-d*g*(2*n + p + 3) + e*f*(p + 1)), NonzeroQ(n + p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern211 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f211))
    rule211 = ReplacementRule(pattern211, lambda e, m, g, a, d, c, n, p, f, x : e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule211)

    def f212(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p - 1), RationalQ(n), Less(n, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern212 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f212))
    rule212 = ReplacementRule(pattern212, lambda e, m, a, g, c, d, p, n, f, x, b : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e*(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule212)

    def f213(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p - 1), RationalQ(n), Less(n, -1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern213 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f213))
    rule213 = ReplacementRule(pattern213, lambda e, m, g, a, c, d, p, n, f, x : -e*(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1)), x)/(g*(n + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)/(c*g*(n + S(1))*(d*g + e*f)))
    rubi.add(rule213)

    def f214(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + p - 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern214 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f214))
    rule214 = ReplacementRule(pattern214, lambda e, m, a, g, d, c, n, p, f, x, b : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))) - (b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)/(c*g*(n + p + S(2))))
    rubi.add(rule214)

    def f215(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), ZeroQ(m + p - 1), IntegerQ(2*p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern215 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f215))
    rule215 = ReplacementRule(pattern215, lambda e, m, g, a, d, c, n, p, f, x : -(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(g*(n + p + S(2))) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule215)

    def f216(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x), FreeQ(p, x)])
    pattern216 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f216))
    rule216 = ReplacementRule(pattern216, lambda e, m, a, g, d, c, n, p, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule216)

    def f217(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), IntegerQ(p - 1/2), IntegersQ(m, n), Unequal(p, 1/2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x), FreeQ(p, x)])
    pattern217 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f217))
    rule217 = ReplacementRule(pattern217, lambda e, m, g, a, d, c, n, p, f, x : Int(ExpandIntegrand(S(1)/sqrt(a + c*x**S(2)), (a + c*x**S(2))**(p + S(1)/2)*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule217)

    def f218(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(n, x), FreeQ(p, x)])
    pattern218 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f218))
    rule218 = ReplacementRule(pattern218, lambda e, m, g, a, d, c, n, p, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule218)

    def f219(e, a, d, c, x, p, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern219 = Pattern(Integral(x_**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f219))
    rule219 = ReplacementRule(pattern219, lambda e, a, d, c, x, p, b : d**S(2)*Int((a + b*x + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((d - e*x)*(a + b*x + c*x**S(2))**p, x)/e**S(2))
    rubi.add(rule219)

    def f220(e, a, d, c, x, p):
        return functools.reduce(operator.and_, [ ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern220 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_/(d_ + x_*WC('e', S(1))), x_),CustomConstraint(f220))
    rule220 = ReplacementRule(pattern220, lambda e, a, d, c, x, p : d**S(2)*Int((a + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((a + c*x**S(2))**p*(d - e*x), x)/e**S(2))
    rubi.add(rule220)

    def f221(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern221 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f221))
    rule221 = ReplacementRule(pattern221, lambda e, m, a, g, c, d, p, f, x, b : g*(d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(3))) - Int((d + e*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*e*g*(d*g + e*f*(m + p + S(1))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))) + e*g*x*(b*e*g*(m + p + S(2)) - c*(d*g*m + e*f*(m + S(2)*p + S(4)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule221)

    def f222(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern222 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_),CustomConstraint(f222))
    rule222 = ReplacementRule(pattern222, lambda e, m, g, a, c, d, p, f, x : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(f + g*x)/(c*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(d + e*x)**m*Simp(-c*e*g*x*(d*g*m + e*f*(m + S(2)*p + S(4))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule222)

    def f223(e, m, g, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(b, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern223 = Pattern(Integral((x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f223))
    rule223 = ReplacementRule(pattern223, lambda e, m, g, c, n, p, f, x, b : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p*(f + g*x)**n, x))
    rubi.add(rule223)

    def f224(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), PositiveQ(a), PositiveQ(d), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern224 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f224))
    rule224 = ReplacementRule(pattern224, lambda e, m, g, a, d, c, n, p, f, x : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule224)

    def f225(e, m, a, g, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), ZeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern225 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f225))
    rule225 = ReplacementRule(pattern225, lambda e, m, a, g, d, c, n, p, f, x, b : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule225)

    def f226(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), ZeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern226 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f226))
    rule226 = ReplacementRule(pattern226, lambda e, m, g, a, d, c, n, p, f, x : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule226)

    def f227(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern227 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f227))
    rule227 = ReplacementRule(pattern227, lambda e, m, g, a, c, d, p, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule227)

    def f228(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern228 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f228))
    rule228 = ReplacementRule(pattern228, lambda e, m, g, a, c, d, p, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x), x), x))
    rubi.add(rule228)

    def f229(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern229 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f229))
    rule229 = ReplacementRule(pattern229, lambda e, a, g, c, d, f, x, b : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int(Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule229)

    def f230(e, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern230 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f230))
    rule230 = ReplacementRule(pattern230, lambda e, g, a, c, d, f, x : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int(Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule230)

    def f231(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + 2*p + 3), ZeroQ(-2*a*e*g + b*(d*g + e*f) - 2*c*d*f), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern231 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f231))
    rule231 = ReplacementRule(pattern231, lambda e, m, a, g, c, d, p, f, x, b : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule231)

    def f232(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 3), ZeroQ(a*e*g + c*d*f), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern232 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f232))
    rule232 = ReplacementRule(pattern232, lambda e, m, g, a, c, d, p, f, x : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule232)

    def f233(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + 2*p + 3), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern233 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f233))
    rule233 = ReplacementRule(pattern233, lambda e, m, g, a, c, d, p, f, x, b : -m*(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule233)

    def f234(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 3), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern234 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f234))
    rule234 = ReplacementRule(pattern234, lambda e, m, g, a, c, d, p, f, x : -m*(a*e*g + c*d*f)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*c*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule234)

    def f235(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern235 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f235))
    rule235 = ReplacementRule(pattern235, lambda e, m, a, g, c, d, p, f, x, b : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule235)

    def f236(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern236 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f236))
    rule236 = ReplacementRule(pattern236, lambda e, m, g, a, c, d, p, f, x : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e*g + c*d*f)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule236)

    def f237(e, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), ZeroQ(-2*a*c*e*g + b**2*e*g*(p + 2) + c*(2*p + 3)*(-b*(d*g + e*f) + 2*c*d*f)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern237 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f237))
    rule237 = ReplacementRule(pattern237, lambda e, a, g, c, d, p, f, x, b : -(a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule237)

    def f238(e, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), ZeroQ(a*e*g - c*d*f*(2*p + 3)), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern238 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f238))
    rule238 = ReplacementRule(pattern238, lambda e, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule238)

    def f239(e, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern239 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f239))
    rule239 = ReplacementRule(pattern239, lambda e, a, g, c, d, p, f, x, b : -(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g - b*c*(d*g + e*f) + S(2)*c*(-a*e*g + c*d*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule239)

    def f240(e, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern240 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f240))
    rule240 = ReplacementRule(pattern240, lambda e, g, a, c, d, p, f, x : -(a + c*x**S(2))**(p + S(1))*(-a*(d*g + e*(f + g*x)) + c*d*f*x)/(S(2)*a*c*(p + S(1))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule240)

    def f241(e, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern241 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f241))
    rule241 = ReplacementRule(pattern241, lambda e, a, g, c, d, p, f, x, b : (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c**S(2)*(S(2)*p + S(3))) - (a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule241)

    def f242(e, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern242 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f242))
    rule242 = ReplacementRule(pattern242, lambda e, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**p, x)/(c*(S(2)*p + S(3))))
    rubi.add(rule242)

    def f243(e, m, g, a, c, x, p, f):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern243 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_),CustomConstraint(f243))
    rule243 = ReplacementRule(pattern243, lambda e, m, g, a, c, x, p, f : f*Int((e*x)**m*(a + c*x**S(2))**p, x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule243)

    def f244(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(m - p), ZeroQ(a*e + b*d), ZeroQ(b*e + c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern244 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f244))
    rule244 = ReplacementRule(pattern244, lambda e, m, g, a, c, d, p, f, x, b : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((f + g*x)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule244)

    def f245(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Greater(p, 0), Less(m, -2), Less(m + 2*p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern245 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f245))
    rule245 = ReplacementRule(pattern245, lambda e, m, a, g, c, d, p, f, x, b : -p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) + b**S(2)*e*(d*g*(p + S(1)) - e*f*(m + p + S(2))) + b*(a*e**S(2)*g*(m + S(1)) - c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))) - c*x*(S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2))) - e*(S(2)*a*e*g*(m + S(1)) - b*(d*g*(m - S(2)*p) + e*f*(m + S(2)*p + S(2))))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*p*(-b*e + S(2)*c*d)*(-d*g + e*f) - e*x*(g*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)) + p*(-b*e + S(2)*c*d)*(-d*g + e*f)) + (d*g - e*f*(m + S(2)))*(a*e**S(2) - b*d*e + c*d**S(2)))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule245)

    def f246(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Greater(p, 0), Less(m, -2), Less(m + 2*p, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern246 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f246))
    rule246 = ReplacementRule(pattern246, lambda e, m, g, a, c, d, p, f, x : -p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) - c*x*(-S(2)*a*e**S(2)*g*(m + S(1)) + S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*c*d**S(2)*p*(-d*g + e*f) - e*x*(S(2)*c*d*p*(-d*g + e*f) + g*(m + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e**S(2) + c*d**S(2))*(d*g - e*f*(m + S(2))))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule246)

    def f247(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(p), Greater(p, 0), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern247 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f247))
    rule247 = ReplacementRule(pattern247, lambda e, m, a, g, c, d, p, f, x, b : p*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-b*e*f*(m + S(2)*p + S(2)) + g*(S(2)*a*e*m + S(2)*a*e + S(2)*b*d*p + b*d) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(b*e*m + b*e + S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule247)

    def f248(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Greater(p, 0), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern248 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f248))
    rule248 = ReplacementRule(pattern248, lambda e, m, g, a, c, d, p, f, x : p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*Simp(g*(S(2)*a*e*m + S(2)*a*e) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule248)

    def f249(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(p), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern249 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f249))
    rule249 = ReplacementRule(pattern249, lambda e, m, a, g, c, d, p, f, x, b : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(c*e*f*(-S(2)*a*e + b*d)*(m + S(2)*p + S(2)) + g*(a*e*(b*e*m + b*e - S(2)*c*d*m) + b*d*(b*e*p - S(2)*c*d*p - c*d)) + x*(c*e*f*(-b*e + S(2)*c*d)*(m + S(2)*p + S(2)) + g*(b**S(2)*e**S(2)*(m + p + S(1)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(1)) - c*e*(S(2)*a*e*(m + S(2)*p + S(1)) + b*d*(m - S(2)*p)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)) - g*(-b*e*p + S(2)*c*d*p + c*d))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule249)

    def f250(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Greater(p, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern250 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f250))
    rule250 = ReplacementRule(pattern250, lambda e, m, g, a, c, d, p, f, x : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*c*d*e*g*m + a*c*e**S(2)*f*(m + S(2)*p + S(2)) - x*(c**S(2)*d*e*f*(m + S(2)*p + S(2)) - g*(a*c*e**S(2)*(m + S(2)*p + S(1)) + c**S(2)*d**S(2)*(S(2)*p + S(1)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*d*g*(S(2)*p + S(1)) + c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule250)

    def f251(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern251 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f251))
    rule251 = ReplacementRule(pattern251, lambda e, m, a, g, c, d, p, f, x, b : -(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g + S(2)*c**S(2)*d*f - c*(S(2)*a*e*g + b*d*g + b*e*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(b*e*g*(a*e*(m + S(-1)) + b*d*(p + S(2))) + S(2)*c**S(2)*d**S(2)*f*(S(2)*p + S(3)) - c*(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) + b*d*(d*g*(S(2)*p + S(3)) - e*f*(m - S(2)*p + S(-4)))) + e*x*(b**S(2)*e*g*(m + p + S(1)) + S(2)*c**S(2)*d*f*(m + S(2)*p + S(2)) - c*(S(2)*a*e*g*m + b*(d*g + e*f)*(m + S(2)*p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule251)

    def f252(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern252 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f252))
    rule252 = ReplacementRule(pattern252, lambda e, m, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(S(2)*a*(d*g + e*f) - x*(-S(2)*a*e*g + S(2)*c*d*f))/(S(4)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) - S(2)*c*d**S(2)*f*(S(2)*p + S(3)) + e*x*(S(2)*a*e*g*m - S(2)*c*d*f*(m + S(2)*p + S(2))), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule252)

    def f253(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern253 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f253))
    rule253 = ReplacementRule(pattern253, lambda e, m, g, a, c, d, p, f, x, b : (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-e*x*(-b*g + S(2)*c*f)*(m + S(2)*p + S(3)) - f*(b*e*m + S(2)*c*d*(S(2)*p + S(3))) + g*(S(2)*a*e*m + b*d*(S(2)*p + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule253)

    def f254(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(m, p), Less(p, -1), Greater(m, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern254 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f254))
    rule254 = ReplacementRule(pattern254, lambda e, m, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*Simp(a*e*g*m - c*d*f*(S(2)*p + S(3)) - c*e*f*x*(m + S(2)*p + S(3)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule254)

    def f255(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern255 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f255))
    rule255 = ReplacementRule(pattern255, lambda e, m, a, g, c, d, p, f, x, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-a*g*(-b*e + S(2)*c*d) + c*x*(f*(-b*e + S(2)*c*d) - g*(-S(2)*a*e + b*d)) + f*(S(2)*a*c*e - b**S(2)*e + b*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*e*x*(-f*(-b*e + S(2)*c*d) + g*(-S(2)*a*e + b*d))*(m + S(2)*p + S(4)) + f*(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3))) - g*(a*e*(b*e*m + b*e - S(2)*c*d*m) - b*d*(-b*e*p - b*e + S(2)*c*d*p + S(3)*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule255)

    def f256(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(p), Less(p, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern256 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f256))
    rule256 = ReplacementRule(pattern256, lambda e, m, g, a, c, d, p, f, x : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-a*c*d*g + a*c*e*f + c*x*(a*e*g + c*d*f))/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(-a*c*d*e*g*m + c*e*x*(a*e*g + c*d*f)*(m + S(2)*p + S(4)) + f*(a*c*e**S(2)*(m + S(2)*p + S(3)) + c**S(2)*d**S(2)*(S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule256)

    def f257(e, m, g, a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern257 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f257))
    rule257 = ReplacementRule(pattern257, lambda e, m, g, a, c, d, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule257)

    def f258(e, m, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), IntegerQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern258 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f258))
    rule258 = ReplacementRule(pattern258, lambda e, m, g, a, c, d, f, x : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule258)

    def f259(e, m, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FractionQ(m), Greater(m, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern259 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f259))
    rule259 = ReplacementRule(pattern259, lambda e, m, a, g, c, d, f, x, b : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule259)

    def f260(e, m, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FractionQ(m), Greater(m, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern260 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f260))
    rule260 = ReplacementRule(pattern260, lambda e, m, g, a, c, d, f, x : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule260)

    def f261(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern261 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f261))
    rule261 = ReplacementRule(pattern261, lambda e, a, g, c, d, f, x, b : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule261)

    def f262(e, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern262 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f262))
    rule262 = ReplacementRule(pattern262, lambda e, g, a, c, d, f, x : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule262)

    def f263(e, m, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FractionQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern263 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f263))
    rule263 = ReplacementRule(pattern263, lambda e, m, a, g, c, d, f, x, b : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule263)

    def f264(e, m, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FractionQ(m), Less(m, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern264 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f264))
    rule264 = ReplacementRule(pattern264, lambda e, m, g, a, c, d, f, x : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule264)

    def f265(e, m, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern265 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f265))
    rule265 = ReplacementRule(pattern265, lambda e, m, a, g, c, d, f, x, b : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule265)

    def f266(e, m, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern266 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f266))
    rule266 = ReplacementRule(pattern266, lambda e, m, g, a, c, d, f, x : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule266)

    def f267(e, m, g, a, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m), Greater(m, 0), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern267 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f267))
    rule267 = ReplacementRule(pattern267, lambda e, m, g, a, c, d, p, f, x, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p*Simp(d*(p + S(1))*(-b*g + S(2)*c*f) + m*(-a*e*g + c*d*f) + x*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(-b*e*g + c*d*g + c*e*f)), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule267)

    def f268(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(m), Greater(m, 0), NonzeroQ(m + 2*p + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern268 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f268))
    rule268 = ReplacementRule(pattern268, lambda e, m, g, a, c, d, p, f, x : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*Simp(-a*e*g*m + c*d*f*(m + S(2)*p + S(2)) + c*x*(d*g*m + e*f*(m + S(2)*p + S(2))), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule268)

    def f269(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern269 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f269))
    rule269 = ReplacementRule(pattern269, lambda e, m, a, g, c, d, p, f, x, b : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule269)

    def f270(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(m), Less(m, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern270 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f270))
    rule270 = ReplacementRule(pattern270, lambda e, m, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule270)

    def f271(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), NegativeIntegerQ(m + 2*p + 3), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern271 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f271))
    rule271 = ReplacementRule(pattern271, lambda e, m, a, g, c, d, p, f, x, b : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule271)

    def f272(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), NegativeIntegerQ(m + 2*p + 3), NonzeroQ(m + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern272 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f272))
    rule272 = ReplacementRule(pattern272, lambda e, m, g, a, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule272)

    def f273(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(4*c*(a - d) - (b - e)**2), ZeroQ(e*f*(b - e) - 2*g*(-a*e + b*d)), NonzeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern273 = Pattern(Integral((f_ + x_*WC('g', S(1)))/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f273))
    rule273 = ReplacementRule(pattern273, lambda e, a, g, c, d, f, x, b : S(4)*f*(a - d)*Subst(Int(S(1)/(S(4)*a - S(4)*d - x**S(2)), x), x, (S(2)*a - S(2)*d + x*(b - e))/sqrt(a + b*x + c*x**S(2)))/(-a*e + b*d))
    rubi.add(rule273)

    def f274(g, a, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(f, x), FreeQ(g, x)])
    pattern274 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f274))
    rule274 = ReplacementRule(pattern274, lambda g, a, c, x, f, b : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule274)

    def f275(a, g, c, x, f):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, f, g), x), FreeQ(a, x), FreeQ(c, x), FreeQ(f, x), FreeQ(g, x)])
    pattern275 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_),CustomConstraint(f275))
    rule275 = ReplacementRule(pattern275, lambda a, g, c, x, f : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule275)

    def f276(e, g, a, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern276 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f276))
    rule276 = ReplacementRule(pattern276, lambda e, g, a, c, x, f, b : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + b*x + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule276)

    def f277(e, g, a, c, x, f):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, e, f, g), x), FreeQ(a, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern277 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_),CustomConstraint(f277))
    rule277 = ReplacementRule(pattern277, lambda e, g, a, c, x, f : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule277)

    def f278(e, m, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern278 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f278))
    rule278 = ReplacementRule(pattern278, lambda e, m, a, g, c, d, p, f, x, b : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule278)

    def f279(e, m, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern279 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_),CustomConstraint(f279))
    rule279 = ReplacementRule(pattern279, lambda e, m, g, a, c, d, p, f, x : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/e)
    rubi.add(rule279)

    def f280(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), IntegersQ(m, n, p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern280 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f280))
    rule280 = ReplacementRule(pattern280, lambda e, m, a, g, c, d, p, n, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule280)

    def f281(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), IntegersQ(m, n, p), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern281 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f281))
    rule281 = ReplacementRule(pattern281, lambda e, m, g, a, c, d, p, n, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule281)

    def f282(e, a, g, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FractionQ(p), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern282 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_),CustomConstraint(f282))
    rule282 = ReplacementRule(pattern282, lambda e, a, g, c, d, p, f, x, b : (a*e**S(2) - b*d*e + c*d**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + b*x + c*x**S(2))**(p + S(-1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule282)

    def f283(e, g, a, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FractionQ(p), Greater(p, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern283 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_),CustomConstraint(f283))
    rule283 = ReplacementRule(pattern283, lambda e, g, a, c, d, p, f, x : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule283)

    def f284(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), IntegersQ(n, p), FractionQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern284 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f284), )
    def With284(e, m, a, g, c, d, p, n, f, x, b):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(c*x**(S(2)*q)/e**S(2) - x**q*(-b*e + S(2)*c*d)/e**S(2) + (a*e**S(2) - b*d*e + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule284 = ReplacementRule(pattern284, lambda e, m, a, g, c, d, p, n, f, x, b : With284(e, m, a, g, c, d, p, n, f, x, b))
    rubi.add(rule284)

    def f285(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), IntegersQ(n, p), FractionQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern285 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f285), )
    def With285(e, m, g, a, c, d, p, n, f, x):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(-S(2)*c*d*x**q/e**S(2) + c*x**(S(2)*q)/e**S(2) + (a*e**S(2) + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule285 = ReplacementRule(pattern285, lambda e, m, g, a, c, d, p, n, f, x : With285(e, m, g, a, c, d, p, n, f, x))
    rubi.add(rule285)

    def f286(e, m, g, a, d, c, n, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(m - n), ZeroQ(d*g + e*f), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern286 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f286))
    rule286 = ReplacementRule(pattern286, lambda e, m, g, a, d, c, n, p, f, x, b : Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule286)

    def f287(e, m, g, a, d, c, n, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(m - n), ZeroQ(d*g + e*f), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern287 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f287))
    rule287 = ReplacementRule(pattern287, lambda e, m, g, a, d, c, n, p, f, x : Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule287)

    def f288(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(m - n), ZeroQ(d*g + e*f), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern288 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f288))
    rule288 = ReplacementRule(pattern288, lambda e, m, a, g, c, d, p, n, f, x, b : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule288)

    def f289(e, m, a, g, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(m - n), ZeroQ(d*g + e*f), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern289 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f289))
    rule289 = ReplacementRule(pattern289, lambda e, m, a, g, c, d, p, n, f, x : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule289)

    def f290(e, m, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern290 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f290))
    rule290 = ReplacementRule(pattern290, lambda e, m, a, g, c, d, n, f, x, b : c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x) + Int((a + b*x)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule290)

    def f291(e, m, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, n), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern291 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f291))
    rule291 = ReplacementRule(pattern291, lambda e, m, g, a, c, d, n, f, x : a*Int((d + e*x)**m*(f + g*x)**n, x) + c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule291)

    def f292(e, m, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, n), Greater(m, 0), Greater(n, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern292 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f292))
    rule292 = ReplacementRule(pattern292, lambda e, m, a, g, c, d, n, f, x, b : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-b*e*g + c*d*g + S(2)*c*e*f + c*e*g*x, x), x)/c**S(2) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(a*b*e*g**S(2) - a*c*d*g**S(2) - S(2)*a*c*e*f*g + c**S(2)*d*f**S(2) + x*(-a*c*e*g**S(2) + b**S(2)*e*g**S(2) - b*c*d*g**S(2) - S(2)*b*c*e*f*g + S(2)*c**S(2)*d*f*g + c**S(2)*e*f**S(2)), x)/(a + b*x + c*x**S(2)), x)/c**S(2))
    rubi.add(rule292)

    def f293(e, m, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, n), Greater(m, 0), Greater(n, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern293 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f293))
    rule293 = ReplacementRule(pattern293, lambda e, m, g, a, c, d, n, f, x : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(d*g + S(2)*e*f + e*g*x, x), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-a*d*g**S(2) - S(2)*a*e*f*g + c*d*f**S(2) + x*(-a*e*g**S(2) + S(2)*c*d*f*g + c*e*f**S(2)), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule293)

    def f294(e, m, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, n), Greater(m, 0), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern294 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f294))
    rule294 = ReplacementRule(pattern294, lambda e, m, a, g, c, d, n, f, x, b : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule294)

    def f295(e, m, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, n), Greater(m, 0), Greater(n, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern295 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f295))
    rule295 = ReplacementRule(pattern295, lambda e, m, g, a, c, d, n, f, x : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule295)

    def f296(e, m, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(m, n), Greater(m, 0), Less(n, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern296 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f296))
    rule296 = ReplacementRule(pattern296, lambda e, m, a, g, c, d, n, f, x, b : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) - b*f*g + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g - b*d*g + c*d*f + c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*g**S(2) - b*f*g + c*f**S(2)))
    rubi.add(rule296)

    def f297(e, m, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), RationalQ(m, n), Greater(m, 0), Less(n, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern297 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f297))
    rule297 = ReplacementRule(pattern297, lambda e, m, g, a, c, d, n, f, x : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g + c*d*f + c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*g**S(2) + c*f**S(2)))
    rubi.add(rule297)

    def f298(e, m, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), PositiveIntegerQ(m + 1/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern298 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f298))
    rule298 = ReplacementRule(pattern298, lambda e, m, a, g, c, d, f, x, b : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule298)

    def f299(e, m, a, g, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), PositiveIntegerQ(m + 1/2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern299 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f299))
    rule299 = ReplacementRule(pattern299, lambda e, m, a, g, c, d, f, x : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + c*x**S(2)), x), x))
    rubi.add(rule299)

    def f300(e, m, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern300 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_),CustomConstraint(f300))
    rule300 = ReplacementRule(pattern300, lambda e, m, a, g, c, d, n, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule300)

    def f301(e, m, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern301 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_),CustomConstraint(f301))
    rule301 = ReplacementRule(pattern301, lambda e, m, g, a, c, d, n, f, x : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule301)

    def f302(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern302 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f302))
    rule302 = ReplacementRule(pattern302, lambda e, m, a, g, c, d, p, n, f, x, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule302)

    def f303(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern303 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f303))
    rule303 = ReplacementRule(pattern303, lambda e, m, g, a, c, d, p, n, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule303)

    def f304(e, m, g, a, d, c, n, p, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(m - p), ZeroQ(a*e + b*d), ZeroQ(b*e + c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern304 = Pattern(Integral((x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f304))
    rule304 = ReplacementRule(pattern304, lambda e, m, g, a, d, c, n, p, x, b : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((g*x)**n*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule304)

    def f305(e, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(n, p), Greater(p, 0), Less(n, -1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern305 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f305))
    rule305 = ReplacementRule(pattern305, lambda e, a, g, c, d, p, n, f, x, b : (a*e**S(2) - b*d*e + c*d**S(2))*Int((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1))*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule305)

    def f306(e, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(n, p), Greater(p, 0), Less(n, -1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern306 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f306))
    rule306 = ReplacementRule(pattern306, lambda e, g, a, c, d, p, n, f, x : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**(n + S(1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**n*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule306)

    def f307(e, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), RationalQ(n, p), Less(p, -1), Greater(n, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern307 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f307))
    rule307 = ReplacementRule(pattern307, lambda e, a, g, c, d, p, n, f, x, b : e*(-d*g + e*f)*Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule307)

    def f308(e, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), RationalQ(n, p), Less(p, -1), Greater(n, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern308 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f308))
    rule308 = ReplacementRule(pattern308, lambda e, g, a, c, d, p, n, f, x : e*(-d*g + e*f)*Int((a + c*x**S(2))**(p + S(1))*(f + g*x)**(n + S(-1))/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule308)

    def f309(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern309 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f309), )
    def With309(e, a, g, c, d, f, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -sqrt(S(2))*sqrt(-g*(b + S(2)*c*x - q)/(-b*g + S(2)*c*f + g*q))*sqrt(-g*(b + S(2)*c*x + q)/(-b*g + S(2)*c*f - g*q))*EllipticPi(e*(-b*g + S(2)*c*f + g*q)/(S(2)*c*(-d*g + e*f)), asin(sqrt(S(2))*sqrt(c/(-b*g + S(2)*c*f + g*q))*sqrt(f + g*x)), (-b*g + S(2)*c*f + g*q)/(-b*g + S(2)*c*f - g*q))/(sqrt(c/(-b*g + S(2)*c*f + g*q))*(-d*g + e*f)*sqrt(a + b*x + c*x**S(2)))
    rule309 = ReplacementRule(pattern309, lambda e, a, g, c, d, f, x, b : With309(e, a, g, c, d, f, x, b))
    rubi.add(rule309)

    def f310(e, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern310 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_),CustomConstraint(f310), )
    def With310(e, g, a, c, d, f, x):
        q = Rt(-a*c, S(2))
        return -S(2)*sqrt(g*(-c*x + q)/(c*f + g*q))*sqrt(-g*(c*x + q)/(c*f - g*q))*EllipticPi(e*(c*f + g*q)/(c*(-d*g + e*f)), asin(sqrt(c/(c*f + g*q))*sqrt(f + g*x)), (c*f + g*q)/(c*f - g*q))/(sqrt(c/(c*f + g*q))*sqrt(a + c*x**S(2))*(-d*g + e*f))
    rule310 = ReplacementRule(pattern310, lambda e, g, a, c, d, f, x : With310(e, g, a, c, d, f, x))
    rubi.add(rule310)

    def f311(e, a, g, c, d, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), IntegerQ(n + 1/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern311 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f311))
    rule311 = ReplacementRule(pattern311, lambda e, a, g, c, d, n, f, x, b : Int(ExpandIntegrand(S(1)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**S(2))), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule311)

    def f312(e, g, a, c, d, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), IntegerQ(n + 1/2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern312 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f312))
    rule312 = ReplacementRule(pattern312, lambda e, g, a, c, d, n, f, x : Int(ExpandIntegrand(S(1)/(sqrt(a + c*x**S(2))*sqrt(f + g*x)), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule312)

    def f313(e, a, g, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern313 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f313))
    rule313 = ReplacementRule(pattern313, lambda e, a, g, c, d, f, x, b : -S(2)*sqrt((-d*g + e*f)**S(2)*(a + b*x + c*x**S(2))/((d + e*x)**S(2)*(a*g**S(2) - b*f*g + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) - b*d*e + c*d**S(2))/(a*g**S(2) - b*f*g + c*f**S(2)) - x**S(2)*(S(2)*a*e*g - b*d*g - b*e*f + S(2)*c*d*f)/(a*g**S(2) - b*f*g + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/((-d*g + e*f)*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule313)

    def f314(e, g, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x)])
    pattern314 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_),CustomConstraint(f314))
    rule314 = ReplacementRule(pattern314, lambda e, g, a, c, d, f, x : -S(2)*sqrt((a + c*x**S(2))*(-d*g + e*f)**S(2)/((d + e*x)**S(2)*(a*g**S(2) + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) + c*d**S(2))/(a*g**S(2) + c*f**S(2)) - x**S(2)*(S(2)*a*e*g + S(2)*c*d*f)/(a*g**S(2) + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/(sqrt(a + c*x**S(2))*(-d*g + e*f)))
    rubi.add(rule314)

    def f315(e, m, g, a, c, x, p, f):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, e, f, g, m, p), x), FreeQ(a, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern315 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(2), x_),CustomConstraint(f315))
    rule315 = ReplacementRule(pattern315, lambda e, m, g, a, c, x, p, f : Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + g**S(2)*x**S(2)), x) + S(2)*f*g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule315)

    def f316(e, m, g, a, c, x, p, f):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, e, f, g, m, p), x), FreeQ(a, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern316 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(3), x_),CustomConstraint(f316))
    rule316 = ReplacementRule(pattern316, lambda e, m, g, a, c, x, p, f : f*Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + S(3)*g**S(2)*x**S(2)), x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p*(S(3)*f**S(2) + g**S(2)*x**S(2)), x)/e)
    rubi.add(rule316)

    def f317(e, m, a, g, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(-4*a*c + b**2), NonzeroQ(a*e**2 - b*d*e + c*d**2), PositiveIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern317 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f317))
    rule317 = ReplacementRule(pattern317, lambda e, m, a, g, c, d, p, n, f, x, b : g*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule317)

    def f318(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-d*g + e*f), NonzeroQ(a*e**2 + c*d**2), PositiveIntegerQ(n), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern318 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_),CustomConstraint(f318))
    rule318 = ReplacementRule(pattern318, lambda e, m, g, a, c, d, p, n, f, x : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/e)
    rubi.add(rule318)

    def f319(e, m, g, a, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, m, n, p), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern319 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f319))
    rule319 = ReplacementRule(pattern319, lambda e, m, g, a, c, d, p, n, f, x, b : Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule319)

    def f320(e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e, f, g, m, n, p), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern320 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f320))
    rule320 = ReplacementRule(pattern320, lambda e, m, g, a, c, d, p, n, f, x : Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule320)

    def f321(u, e, m, g, a, c, d, p, n, f, x, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern321 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f321))
    rule321 = ReplacementRule(pattern321, lambda u, e, m, g, a, c, d, p, n, f, x, b : Subst(Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule321)

    def f322(u, e, m, g, a, c, d, p, n, f, x):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern322 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_),CustomConstraint(f322))
    rule322 = ReplacementRule(pattern322, lambda u, e, m, g, a, c, d, p, n, f, x : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule322)

    def f323(e, a, d, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern323 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f323))
    rule323 = ReplacementRule(pattern323, lambda e, a, d, c, x, p, f, q, b : (c/f)**p*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule323)

    def f324(e, a, d, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern324 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f324))
    rule324 = ReplacementRule(pattern324, lambda e, a, d, c, x, p, f, q, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule324)

    def f325(e, a, d, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern325 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f325))
    rule325 = ReplacementRule(pattern325, lambda e, a, d, c, x, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule325)

    def f326(a, d, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern326 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f326))
    rule326 = ReplacementRule(pattern326, lambda a, d, c, x, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule326)

    def f327(e, a, q, c, d, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(c*(-2*d*f + e**2*(q + 2)) + f*(2*q + 3)*(2*a*f - b*e)), NonzeroQ(q + 1), NonzeroQ(2*q + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern327 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f327))
    rule327 = ReplacementRule(pattern327, lambda e, a, q, c, d, f, x, b : (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule327)

    def f328(e, a, q, c, d, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(2*a*f**2*(2*q + 3) + c*(-2*d*f + e**2*(q + 2))), NonzeroQ(q + 1), NonzeroQ(2*q + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern328 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f328))
    rule328 = ReplacementRule(pattern328, lambda e, a, q, c, d, f, x : (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule328)

    def f329(a, d, c, x, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(q + 1), ZeroQ(2*a*f*q + 3*a*f - c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(q, x)])
    pattern329 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f329))
    rule329 = ReplacementRule(pattern329, lambda a, d, c, x, f, q, b : (d + f*x**S(2))**(q + S(1))*(S(2)*a*f*x*(q + S(1)) + b*d)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule329)

    def f330(a, d, c, x, f, q, b):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(q + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(q, x)])
    pattern330 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f330))
    rule330 = ReplacementRule(pattern330, lambda a, d, c, x, f, q, b : b*Int(x*(d + f*x**S(2))**q, x) + Int((a + c*x**S(2))*(d + f*x**S(2))**q, x))
    rubi.add(rule330)

    def f331(e, a, x, c, d, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), PositiveIntegerQ(q + 2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern331 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f331))
    rule331 = ReplacementRule(pattern331, lambda e, a, x, c, d, f, q, b : Int(ExpandIntegrand((a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule331)

    def f332(e, a, x, c, d, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), PositiveIntegerQ(q + 2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern332 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f332))
    rule332 = ReplacementRule(pattern332, lambda e, a, x, c, d, f, q : Int(ExpandIntegrand((a + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule332)

    def f333(e, a, q, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(q), Less(q, -1), NonzeroQ(c*(-2*d*f + e**2*(q + 2)) + f*(2*q + 3)*(2*a*f - b*e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern333 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f333))
    rule333 = ReplacementRule(pattern333, lambda e, a, q, c, d, f, x, b : -(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f - S(2)*b*d*f + c*d*e + x*(c*(-S(2)*d*f + e**S(2)) + f*(S(2)*a*f - b*e)))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule333)

    def f334(e, a, q, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(q), Less(q, -1), NonzeroQ(2*a*f**2*(2*q + 3) + c*(-2*d*f + e**2*(q + 2))), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern334 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f334))
    rule334 = ReplacementRule(pattern334, lambda e, a, q, c, d, f, x : -(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f + c*d*e + x*(S(2)*a*f**S(2) + c*(-S(2)*d*f + e**S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule334)

    def f335(a, d, c, x, f, q, b):
        return functools.reduce(operator.and_, [ RationalQ(q), Less(q, -1), NonzeroQ(2*a*f*q + 3*a*f - c*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern335 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f335))
    rule335 = ReplacementRule(pattern335, lambda a, d, c, x, f, q, b : (d + f*x**S(2))**(q + S(1))*(b*d - x*(a*f - c*d))/(S(2)*d*f*(q + S(1))) + (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**(q + S(1)), x)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule335)

    def f336(e, a, q, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NonzeroQ(c*(-2*d*f + e**2*(q + 2)) + f*(2*q + 3)*(2*a*f - b*e)), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(q, x)])
    pattern336 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f336))
    rule336 = ReplacementRule(pattern336, lambda e, a, q, c, d, f, x, b : (c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule336)

    def f337(e, a, q, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NonzeroQ(2*a*f**2*(2*q + 3) + c*(-2*d*f + e**2*(q + 2))), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(a, x), FreeQ(c, x), FreeQ(q, x)])
    pattern337 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f337))
    rule337 = ReplacementRule(pattern337, lambda e, a, q, c, d, f, x : (S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule337)

    def f338(a, d, c, x, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(2*a*f*q + 3*a*f - c*d), FreeQ(d, x), FreeQ(f, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(q, x)])
    pattern338 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_),CustomConstraint(f338))
    rule338 = ReplacementRule(pattern338, lambda a, d, c, x, f, q, b : (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**q, x)/(f*(S(2)*q + S(3))) + (d + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule338)

    def f339(e, a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern339 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f339))
    rule339 = ReplacementRule(pattern339, lambda e, a, x, c, d, p, f, q, b : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(b*e*q + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(S(2)*b*f*q + S(2)*c*e*(S(2)*p + q + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule339)

    def f340(a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern340 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f340))
    rule340 = ReplacementRule(pattern340, lambda a, x, c, d, p, f, q, b : (b + S(2)*c*x)*(d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*b*f*q*x + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule340)

    def f341(e, a, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern341 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f341))
    rule341 = ReplacementRule(pattern341, lambda e, a, x, c, d, p, f, q : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(S(2)*p + q + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule341)

    def f342(e, a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern342 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f342))
    rule342 = ReplacementRule(pattern342, lambda e, a, x, c, d, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d) + c*x*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + S(2)*c*(p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) - e*(p + q + S(2))*(-S(2)*a*c**S(2)*e - b**S(3)*f + b**S(2)*c*e - b*c*(-S(3)*a*f + c*d)) + x*(S(2)*f*(p + q + S(2))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule342)

    def f343(a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), NonzeroQ(b**2*d*f + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(q, x)])
    pattern343 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f343))
    rule343 = ReplacementRule(pattern343, lambda a, x, c, d, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(b**S(3)*f + b*c*(-S(3)*a*f + c*d) + c*x*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*c*(p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + x*(-b*f*(p + S(1))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*f*(b**S(3)*f + b*c*(-S(3)*a*f + c*d))*(p + q + S(2))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule343)

    def f344(e, a, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(a*c*e**2 + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern344 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f344))
    rule344 = ReplacementRule(pattern344, lambda e, a, x, c, d, p, f, q : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c**S(2)*e + c*x*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(S(2)*a*c**S(2)*e**S(2)*(p + q + S(2)) + c*f*x**S(2)*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + S(2)*q + S(5)) + S(2)*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + x*(S(4)*a*c**S(2)*e*f*(p + q + S(2)) + c*e*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + q + S(4))) - (-S(2)*a*c*f + S(2)*c**S(2)*d)*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule344)

    def f345(e, a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 1), NonzeroQ(p + q), NonzeroQ(2*p + 2*q + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern345 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f345))
    rule345 = ReplacementRule(pattern345, lambda e, a, x, c, d, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(3)*p + S(2)*q) - c*e*(S(2)*p + q) + S(2)*c*f*x*(p + q))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(x**S(2)*(c*(p + q)*(-c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))) + f*(-S(2)*a*f + b*e)*(S(4)*p + S(2)*q + S(-1))) + p*(-p + S(1))*(-b*f + c*e)**S(2)) + x*(S(2)*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)*(-b*f + c*e) - (p + q)*(b*(c*(S(2)*p + q)*(-S(4)*d*f + e**S(2)) + f*(S(2)*p + S(2)*q + S(1))*(S(2)*a*f - b*e + S(2)*c*d)) + e*f*(-p + S(1))*(-S(4)*a*c + b**S(2)))) + (-p + S(1))*(S(2)*p + q)*(-a*e + b*d)*(-b*f + c*e) - (p + q)*(-a*(c*(S(2)*d*f - e**S(2)*(S(2)*p + q)) + f*(-S(2)*a*f + b*e)*(S(2)*p + S(2)*q + S(1))) + b**S(2)*d*f*(-p + S(1))), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule345)

    def f346(a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Greater(p, 1), NonzeroQ(p + q), NonzeroQ(2*p + 2*q + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(q, x)])
    pattern346 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f346))
    rule346 = ReplacementRule(pattern346, lambda a, x, c, d, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(b*(S(3)*p + S(2)*q) + S(2)*c*x*(p + q))*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-2))*Simp(b**S(2)*d*(p + S(-1))*(S(2)*p + q) + x**S(2)*(b**S(2)*f*p*(-p + S(1)) + S(2)*c*(p + q)*(-a*f*(S(4)*p + S(2)*q + S(-1)) + c*d*(S(2)*p + S(-1)))) - x*(S(2)*b*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d) - S(2)*b*(p + q)*(S(2)*c*d*(S(2)*p + q) - (a*f + c*d)*(S(2)*p + S(2)*q + S(1)))) - (p + q)*(-S(2)*a*(-a*f*(S(2)*p + S(2)*q + S(1)) + c*d) + b**S(2)*d*(-p + S(1))), x), x)/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule346)

    def f347(e, a, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 1), NonzeroQ(p + q), NonzeroQ(2*p + 2*q + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern347 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f347))
    rule347 = ReplacementRule(pattern347, lambda e, a, x, c, d, p, f, q : -c*(a + c*x**S(2))**(p + S(-1))*(e*(S(2)*p + q) - S(2)*f*x*(p + q))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(-a*c*e**S(2)*(-p + S(1))*(S(2)*p + q) + a*(p + q)*(-S(2)*a*f**S(2)*(S(2)*p + S(2)*q + S(1)) + c*(S(2)*d*f - e**S(2)*(S(2)*p + q))) + x**S(2)*(c**S(2)*e**S(2)*p*(-p + S(1)) - c*(p + q)*(S(2)*a*f**S(2)*(S(4)*p + S(2)*q + S(-1)) + c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))))) + x*(S(4)*a*c*e*f*(-p + S(1))*(p + q) + S(2)*c*e*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule347)

    def f348(e, a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern348 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f348), )
    def With348(e, a, d, c, x, f, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-a*c*f + b**S(2)*f - b*c*e + c**S(2)*d - x*(-b*c*f + c**S(2)*e))/(a + b*x + c*x**S(2)), x)/q + Int((a*f**S(2) - b*e*f - c*d*f + c*e**S(2) + x*(-b*f**S(2) + c*e*f))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule348 = ReplacementRule(pattern348, lambda e, a, d, c, x, f, b : With348(e, a, d, c, x, f, b))
    rubi.add(rule348)

    def f349(a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern349 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f349), )
    def With349(a, d, c, x, f, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return -Int((-a*f**S(2) + b*f**S(2)*x + c*d*f)/(d + f*x**S(2)), x)/q + Int((-a*c*f + b**S(2)*f + b*c*f*x + c**S(2)*d)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule349 = ReplacementRule(pattern349, lambda a, d, c, x, f, b : With349(a, d, c, x, f, b))
    rubi.add(rule349)

    def f350(e, a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-b*f + c*e), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern350 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f350))
    rule350 = ReplacementRule(pattern350, lambda e, a, c, d, f, x, b : -S(2)*e*Subst(Int(S(1)/(e*(-S(4)*a*f + b*e) - x**S(2)*(-a*e + b*d)), x), x, (e + S(2)*f*x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule350)

    def f351(e, a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(-b*f + c*e), PosQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern351 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f351), )
    def With351(e, a, c, d, f, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - S(2)*c*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule351 = ReplacementRule(pattern351, lambda e, a, c, d, f, x, b : With351(e, a, c, d, f, x, b))
    rubi.add(rule351)

    def f352(e, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), PosQ(-a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern352 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f352))
    rule352 = ReplacementRule(pattern352, lambda e, a, c, d, f, x : Int(S(1)/((a - x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2) + Int(S(1)/((a + x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2))
    rubi.add(rule352)

    def f353(a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), PosQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern353 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f353), )
    def With353(a, d, c, x, f, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule353 = ReplacementRule(pattern353, lambda a, d, c, x, f, b : With353(a, d, c, x, f, b))
    rubi.add(rule353)

    def f354(e, a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(-b*f + c*e), NegQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern354 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f354), )
    def With354(e, a, c, d, f, x, b):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d - q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule354 = ReplacementRule(pattern354, lambda e, a, c, d, f, x, b : With354(e, a, c, d, f, x, b))
    rubi.add(rule354)

    def f355(e, a, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NegQ(-a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern355 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f355), )
    def With355(e, a, c, d, f, x):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d + c*e*x - q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + c*e*x + q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule355 = ReplacementRule(pattern355, lambda e, a, c, d, f, x : With355(e, a, c, d, f, x))
    rubi.add(rule355)

    def f356(a, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NegQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern356 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f356), )
    def With356(a, c, d, f, x, b):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f - b*f*x + c*d - q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) + Int((-a*f - b*f*x + c*d + q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule356 = ReplacementRule(pattern356, lambda a, c, d, f, x, b : With356(a, c, d, f, x, b))
    rubi.add(rule356)

    def f357(e, a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern357 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_),CustomConstraint(f357))
    rule357 = ReplacementRule(pattern357, lambda e, a, d, c, x, f, b : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f + c*d + x*(-b*f + c*e))/(sqrt(a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule357)

    def f358(a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern358 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1))), x_),CustomConstraint(f358))
    rule358 = ReplacementRule(pattern358, lambda a, d, c, x, f, b : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f - b*f*x + c*d)/((d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), x)/f)
    rubi.add(rule358)

    def f359(e, a, d, c, x, f):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern359 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_),CustomConstraint(f359))
    rule359 = ReplacementRule(pattern359, lambda e, a, d, c, x, f : c*Int(S(1)/sqrt(a + c*x**S(2)), x)/f - Int((-a*f + c*d + c*e*x)/(sqrt(a + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule359)

    def f360(e, a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern360 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f360), )
    def With360(e, a, d, c, x, f, b):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*sqrt(d + e*x + f*x**S(2))), x)/sqrt(a + b*x + c*x**S(2))
    rule360 = ReplacementRule(pattern360, lambda e, a, d, c, x, f, b : With360(e, a, d, c, x, f, b))
    rubi.add(rule360)

    def f361(a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x)])
    pattern361 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f361), )
    def With361(a, d, c, x, f, b):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(d + f*x**S(2))*sqrt(b + S(2)*c*x + r)), x)/sqrt(a + b*x + c*x**S(2))
    rule361 = ReplacementRule(pattern361, lambda a, d, c, x, f, b : With361(a, d, c, x, f, b))
    rubi.add(rule361)

    def f362(e, a, q, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern362 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f362))
    rule362 = ReplacementRule(pattern362, lambda e, a, q, c, d, p, f, x, b : Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule362)

    def f363(e, a, q, c, d, p, f, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e, f, p, q), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern363 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f363))
    rule363 = ReplacementRule(pattern363, lambda e, a, q, c, d, p, f, x : Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule363)

    def f364(u, e, a, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern364 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f364))
    rule364 = ReplacementRule(pattern364, lambda u, e, a, x, c, d, p, f, q, b : Subst(Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule364)

    def f365(u, e, a, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x)])
    pattern365 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f365))
    rule365 = ReplacementRule(pattern365, lambda u, e, a, x, c, d, p, f, q : Subst(Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule365)

    def f366(e, m, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x), FreeQ(q, x)])
    pattern366 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f366))
    rule366 = ReplacementRule(pattern366, lambda e, m, g, a, d, h, c, x, p, f, q, b : (c/f)**p*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule366)

    def f367(e, m, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x), FreeQ(q, x)])
    pattern367 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f367))
    rule367 = ReplacementRule(pattern367, lambda e, m, g, a, d, h, c, x, p, f, q, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule367)

    def f368(e, m, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern368 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f368))
    rule368 = ReplacementRule(pattern368, lambda e, m, g, a, d, h, c, x, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule368)

    def f369(m, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern369 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_),CustomConstraint(f369))
    rule369 = ReplacementRule(pattern369, lambda m, g, a, d, h, c, x, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(g + h*x)**m, x))
    rubi.add(rule369)

    def f370(e, m, g, a, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(a*h**2 - b*g*h + c*g**2), ZeroQ(a**2*f*h**2 - a*c*e*g*h + c**2*d*g**2), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern370 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_),CustomConstraint(f370))
    rule370 = ReplacementRule(pattern370, lambda e, m, g, a, h, c, d, p, f, x, b : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule370)

    def f371(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*h**2 + c*g**2), ZeroQ(a**2*f*h**2 - a*c*e*g*h + c**2*d*g**2), IntegerQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern371 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_),CustomConstraint(f371))
    rule371 = ReplacementRule(pattern371, lambda e, m, g, a, h, c, d, p, f, x : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule371)

    def f372(m, g, a, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(a*h**2 - b*g*h + c*g**2), ZeroQ(a**2*f*h**2 + c**2*d*g**2), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern372 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f372))
    rule372 = ReplacementRule(pattern372, lambda m, g, a, h, c, d, p, f, x, b : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule372)

    def f373(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*h**2 + c*g**2), ZeroQ(a**2*f*h**2 + c**2*d*g**2), IntegerQ(m), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern373 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_),CustomConstraint(f373))
    rule373 = ReplacementRule(pattern373, lambda m, g, a, h, c, d, p, f, x : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule373)

    def f374(e, a, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(a*f**2 - b*e*f + c*e**2), IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern374 = Pattern(Integral(x_**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f374))
    rule374 = ReplacementRule(pattern374, lambda e, a, c, x, p, f, q, b : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule374)

    def f375(e, a, c, x, p, f, q):
        return functools.reduce(operator.and_, [ ZeroQ(a*f**2 + c*e**2), IntegerQ(p), FreeQ(a, x), FreeQ(c, x), FreeQ(e, x), FreeQ(f, x), FreeQ(q, x)])
    pattern375 = Pattern(Integral(x_**WC('p', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f375))
    rule375 = ReplacementRule(pattern375, lambda e, a, c, x, p, f, q : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule375)

    def f376(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(b*f*h*(m + p + 2) + c*(-e*h*(m + 2*p + 3) + 2*f*g*(p + 1))), ZeroQ(b*f*g*(p + 1) + h*(a*f*(m + 1) - c*d*(m + 2*p + 3))), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern376 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f376))
    rule376 = ReplacementRule(pattern376, lambda e, m, a, g, h, c, d, p, f, x, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule376)

    def f377(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(c*(-e*h*(m + 2*p + 3) + 2*f*g*(p + 1))), ZeroQ(h*(a*f*(m + 1) - c*d*(m + 2*p + 3))), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern377 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f377))
    rule377 = ReplacementRule(pattern377, lambda e, m, g, a, h, c, d, p, f, x : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule377)

    def f378(m, a, g, d, h, c, x, p, f, b):
        return functools.reduce(operator.and_, [ ZeroQ(b*f*h*(m + p + 2) + 2*c*f*g*(p + 1)), ZeroQ(b*f*g*(p + 1) + h*(a*f*(m + 1) - c*d*(m + 2*p + 3))), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern378 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f378))
    rule378 = ReplacementRule(pattern378, lambda m, a, g, d, h, c, x, p, f, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule378)

    def f379(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern379 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f379))
    rule379 = ReplacementRule(pattern379, lambda e, m, a, g, h, c, d, p, f, x, b : Int(ExpandIntegrand((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule379)

    def f380(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern380 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f380))
    rule380 = ReplacementRule(pattern380, lambda e, m, g, a, h, c, d, p, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule380)

    def f381(m, a, g, d, h, c, x, p, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern381 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f381))
    rule381 = ReplacementRule(pattern381, lambda m, a, g, d, h, c, x, p, f, b : Int(ExpandIntegrand((d + f*x**S(2))*(g + h*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule381)

    def f382(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern382 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_),CustomConstraint(f382))
    rule382 = ReplacementRule(pattern382, lambda m, g, a, h, c, d, p, f, x : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + f*x**S(2))*(g + h*x)**m, x), x))
    rubi.add(rule382)

    def f383(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(m), Less(m, -1), NonzeroQ(a*h**2 - b*g*h + c*g**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern383 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f383))
    rule383 = ReplacementRule(pattern383, lambda e, m, a, g, h, c, d, p, f, x, b : (g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(f*g**S(2)*(p + S(1)) - h*(-d*h*(m + p + S(2)) + e*g*(p + S(1)))) + h*(m + S(1))*(a*e*h - a*f*g + c*d*g) - x*(c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule383)

    def f384(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(m), Less(m, -1), NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern384 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f384))
    rule384 = ReplacementRule(pattern384, lambda e, m, g, a, h, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(a*e*h - a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule384)

    def f385(m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(m), Less(m, -1), NonzeroQ(a*h**2 - b*g*h + c*g**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern385 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f385))
    rule385 = ReplacementRule(pattern385, lambda m, a, g, h, c, d, p, f, x, b : (g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))*(a + b*x + c*x**S(2))**(p + S(1))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(d*h**S(2)*(m + p + S(2)) + f*g**S(2)*(p + S(1))) + h*(m + S(1))*(-a*f*g + c*d*g) - x*(c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule385)

    def f386(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ RationalQ(m), Less(m, -1), NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x)])
    pattern386 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_),CustomConstraint(f386))
    rule386 = ReplacementRule(pattern386, lambda m, g, a, h, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(-a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule386)

    def f387(e, a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(a*h**2 - b*g*h + c*g**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern387 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_),CustomConstraint(f387))
    rule387 = ReplacementRule(pattern387, lambda e, a, g, h, c, d, f, x, b : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((a*e*h - a*f*g - b*d*h + c*d*g + x*(a*f*h - b*f*g - c*d*h + c*e*g))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule387)

    def f388(e, g, a, h, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern388 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(x_*WC('h', S(1)) + WC('g', S(0)))), x_),CustomConstraint(f388))
    rule388 = ReplacementRule(pattern388, lambda e, g, a, h, c, d, f, x : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((a*e*h - a*f*g + c*d*g + x*(a*f*h - c*d*h + c*e*g))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule388)

    def f389(a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(a*h**2 - b*g*h + c*g**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern389 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_),CustomConstraint(f389))
    rule389 = ReplacementRule(pattern389, lambda a, g, h, c, d, f, x, b : (d*h**S(2) + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((-a*f*g - b*d*h + c*d*g - x*(-a*f*h + b*f*g + c*d*h))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule389)

    def f390(g, a, h, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern390 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(g_ + x_*WC('h', S(1)))), x_),CustomConstraint(f390))
    rule390 = ReplacementRule(pattern390, lambda g, a, h, c, d, f, x : (d*h**S(2) + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((-a*f*g + c*d*g - x*(-a*f*h + c*d*h))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule390)

    def f391(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern391 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f391))
    rule391 = ReplacementRule(pattern391, lambda e, m, a, g, h, c, d, p, f, x, b : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f - S(2)*a*c*e + b*c*d + x*(c*(-b*e + S(2)*c*d) + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(c*(S(2)*p + S(3))*(-b*e + S(2)*c*d) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f - S(2)*a*c*e + b*c*d) + h*x*(c*(-b*e + S(2)*c*d)*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule391)

    def f392(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern392 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f392))
    rule392 = ReplacementRule(pattern392, lambda e, m, g, a, h, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(a*e - x*(-a*f + c*d))/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*(e*h*m + f*g) - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule392)

    def f393(m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern393 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f393))
    rule393 = ReplacementRule(pattern393, lambda m, a, g, h, c, d, p, f, x, b : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f + b*c*d + x*(S(2)*c**S(2)*d + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(S(2)*c**S(2)*d*(S(2)*p + S(3)) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f + b*c*d) + h*x*(S(2)*c**S(2)*d*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule393)

    def f394(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ RationalQ(m, p), Less(p, -1), Greater(m, 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern394 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_),CustomConstraint(f394))
    rule394 = ReplacementRule(pattern394, lambda m, g, a, h, c, d, p, f, x : -x*(a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(-a*f + c*d)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*f*g - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule394)

    def f395(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(c*g**2 - h*(-a*h + b*g)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern395 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f395))
    rule395 = ReplacementRule(pattern395, lambda e, m, a, g, h, c, d, p, f, x, b : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(-S(2)*a*e*h + S(2)*a*f*g + b*d*h + b*e*g)) - (-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) - h*(-(-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + (p + S(1))*(c*g**S(2) - h*(-a*h + b*g))*(S(2)*a*f - b*e + S(2)*c*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule395)

    def f396(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern396 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f396))
    rule396 = ReplacementRule(pattern396, lambda e, m, g, a, h, c, d, p, f, x : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*c*e*g - a*h*(-a*f + c*d) - c*x*(a*e*h - a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(g*(p + S(2))*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g) + h*x*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g)*(m + S(2)*p + S(4)) - h*(a*c*e*g - a*h*(-a*f + c*d))*(m + p + S(2)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule396)

    def f397(m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), NonzeroQ(c*g**2 - h*(-a*h + b*g)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern397 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f397))
    rule397 = ReplacementRule(pattern397, lambda m, a, g, h, c, d, p, f, x, b : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*d*(-b*h + S(2)*c*g) - x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(S(2)*a*f*g + b*d*h)) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) - h*(-b*d*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + S(2)*(p + S(1))*(a*f + c*d)*(c*g**S(2) - h*(-a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule397)

    def f398(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(p, -1), NonzeroQ(a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x)])
    pattern398 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_),CustomConstraint(f398))
    rule398 = ReplacementRule(pattern398, lambda m, g, a, h, c, d, p, f, x : -(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*h*(-a*f + c*d) + c*x*(-a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(a*h**S(2)*(-a*f + c*d)*(m + p + S(2)) + g*(p + S(2))*(-a*c*f*g + c**S(2)*d*g) + h*x*(-a*c*f*g + c**S(2)*d*g)*(m + S(2)*p + S(4)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule398)

    def f399(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern399 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f399))
    rule399 = ReplacementRule(pattern399, lambda e, m, a, g, h, c, d, p, f, x, b : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule399)

    def f400(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern400 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f400))
    rule400 = ReplacementRule(pattern400, lambda e, m, g, a, h, c, d, p, f, x : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule400)

    def f401(m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern401 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f401))
    rule401 = ReplacementRule(pattern401, lambda m, a, g, h, c, d, p, f, x, b : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule401)

    def f402(m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern402 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_),CustomConstraint(f402))
    rule402 = ReplacementRule(pattern402, lambda m, g, a, h, c, d, p, f, x : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule402)

    def f403(e, m, a, g, h, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern403 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f403))
    rule403 = ReplacementRule(pattern403, lambda e, m, a, g, h, c, d, p, f, x, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))) + x*(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1)))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule403)

    def f404(e, m, g, a, h, c, d, p, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern404 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_),CustomConstraint(f404))
    rule404 = ReplacementRule(pattern404, lambda e, m, g, a, h, c, d, p, f, x : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(c*x*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule404)

    def f405(m, a, g, d, h, c, x, p, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern405 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f405))
    rule405 = ReplacementRule(pattern405, lambda m, a, g, d, h, c, x, p, f, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + f*x*(b*h*(m + p + S(2)) + S(2)*c*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule405)

    def f406(m, g, a, d, h, c, x, p, f):
        return functools.reduce(operator.and_, [ NonzeroQ(m + 2*p + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x)])
    pattern406 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)))*(g_ + x_*WC('h', S(1)))**WC('m', S(1)), x_),CustomConstraint(f406))
    rule406 = ReplacementRule(pattern406, lambda m, g, a, d, h, c, x, p, f : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(S(2)*c*f*g*x*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule406)

    def f407(e, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern407 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f407))
    rule407 = ReplacementRule(pattern407, lambda e, g, a, d, h, c, x, p, f, q, b : Int(ExpandIntegrand((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule407)

    def f408(e, g, a, d, h, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern408 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f408))
    rule408 = ReplacementRule(pattern408, lambda e, g, a, d, h, c, x, p, f, q : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule408)

    def f409(e, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern409 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f409))
    rule409 = ReplacementRule(pattern409, lambda e, g, a, d, h, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + e*q*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)) + x*(-e*(b*h - S(2)*c*g)*(S(2)*p + q + S(3)) + S(2)*f*q*(-S(2)*a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule409)

    def f410(e, g, a, d, h, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern410 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f410))
    rule410 = ReplacementRule(pattern410, lambda e, g, a, d, h, c, x, p, f, q : (a + c*x**S(2))**(p + S(1))*(a*h - c*g*x)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-a*e*h*q + c*d*g*(S(2)*p + S(3)) + c*f*g*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(-S(2)*a*f*h*q + c*e*g*(S(2)*p + q + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule410)

    def f411(g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern411 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f411))
    rule411 = ReplacementRule(pattern411, lambda g, a, d, h, c, x, p, f, q, b : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + S(2)*f*q*x*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule411)

    def f412(e, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern412 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f412))
    rule412 = ReplacementRule(pattern412, lambda e, g, a, d, h, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + c*x*(g*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - h*(a*b*f - S(2)*a*c*e + b*c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)) - e*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - x*(S(2)*f*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d))) + (p + S(1))*(b*h - S(2)*c*g)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule412)

    def f413(e, g, a, d, h, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(a*c*e**2 + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern413 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f413))
    rule413 = ReplacementRule(pattern413, lambda e, g, a, d, h, c, x, p, f, q : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d) + c*x*(S(2)*a*c*e*h + g*(-S(2)*a*c*f + S(2)*c**S(2)*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + S(2)*q + S(5)) - S(2)*c*g*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) - e*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2)) - x*(c*e*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + q + S(4)) + S(2)*f*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d)), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule413)

    def f414(g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), NonzeroQ(b**2*d*f + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern414 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f414))
    rule414 = ReplacementRule(pattern414, lambda g, a, d, h, c, x, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*c*g*(a*f + c*d) + c*x*(g*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - h*(a*b*f + b*c*d)) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) - x*(-b*f*(p + S(1))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) + S(2)*f*(-b*c*g*(a*f + c*d) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b*h - S(2)*c*g)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule414)

    def f415(e, g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern415 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f415))
    rule415 = ReplacementRule(pattern415, lambda e, g, a, d, h, c, x, p, f, q, b : h*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-a*e + b*d) + x**S(2)*(c*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-b*f + c*e)) + x*(b*(e*h - S(2)*f*g)*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule415)

    def f416(e, g, a, d, h, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern416 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f416))
    rule416 = ReplacementRule(pattern416, lambda e, g, a, d, h, c, x, p, f, q : h*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) + Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*e*h*p - a*(e*h - S(2)*f*g)*(p + q + S(1)) - S(2)*h*p*x*(-a*f + c*d) - x**S(2)*(c*e*h*p + c*(e*h - S(2)*f*g)*(p + q + S(1))), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule416)

    def f417(g, a, d, h, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(q, x)])
    pattern417 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f417))
    rule417 = ReplacementRule(pattern417, lambda g, a, d, h, c, x, p, f, q, b : h*(d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*f*(p + q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*f*g*(p + q + S(1)) + b*d*h*p + x**S(2)*(-b*f*h*p - S(2)*c*f*g*(p + q + S(1))) + x*(-S(2)*b*f*g*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule417)

    def f418(e, g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern418 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f418), )
    def With418(e, g, a, d, h, c, x, f, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(-a*b*f*h + a*c*e*h - a*c*f*g + b**S(2)*f*g - b*c*e*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(a + b*x + c*x**S(2)), x)/q + Int(Simp(a*f**S(2)*g + b*d*f*h - b*e*f*g - c*d*e*h - c*d*f*g + c*e**S(2)*g - f*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule418 = ReplacementRule(pattern418, lambda e, g, a, d, h, c, x, f, b : With418(e, g, a, d, h, c, x, f, b))
    rubi.add(rule418)

    def f419(g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern419 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f419), )
    def With419(g, a, d, h, c, x, f, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(a*f**S(2)*g + b*d*f*h - c*d*f*g - f*x*(-a*f*h + b*f*g + c*d*h), x)/(d + f*x**S(2)), x)/q + Int(Simp(-a*b*f*h - a*c*f*g + b**S(2)*f*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h), x)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule419 = ReplacementRule(pattern419, lambda g, a, d, h, c, x, f, b : With419(g, a, d, h, c, x, f, b))
    rubi.add(rule419)

    def f420(g, a, d, h, c, x, f):
        return functools.reduce(operator.and_, [ PositiveQ(a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern420 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f420))
    rule420 = ReplacementRule(pattern420, lambda g, a, d, h, c, x, f : g*Int(S(1)/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x) + h*Int(x/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x))
    rubi.add(rule420)

    def f421(g, a, d, h, c, x, f):
        return functools.reduce(operator.and_, [ FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern421 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f421), )
    def With421(g, a, d, h, c, x, f):
        q = Rt(-a*c, S(2))
        return -(c*g - h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(c*x + q)), x)/(S(2)*q) - (c*g + h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(-c*x + q)), x)/(S(2)*q)
    rule421 = ReplacementRule(pattern421, lambda g, a, d, h, c, x, f : With421(g, a, d, h, c, x, f))
    rubi.add(rule421)

    def f422(e, g, a, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-b*f + c*e), ZeroQ(e*h - 2*f*g), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern422 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f422))
    rule422 = ReplacementRule(pattern422, lambda e, g, a, h, c, d, f, x, b : -S(2)*g*Subst(Int(S(1)/(-a*e + b*d - b*x**S(2)), x), x, sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule422)

    def f423(e, g, a, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-b*f + c*e), NonzeroQ(e*h - 2*f*g), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern423 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f423))
    rule423 = ReplacementRule(pattern423, lambda e, g, a, h, c, d, f, x, b : h*Int((e + S(2)*f*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f) - (e*h - S(2)*f*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f))
    rubi.add(rule423)

    def f424(e, a, d, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern424 = Pattern(Integral(x_/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f424))
    rule424 = ReplacementRule(pattern424, lambda e, a, d, c, x, f, b : -S(2)*e*Subst(Int((-d*x**S(2) + S(1))/(-b*f + c*e + d**S(2)*x**S(4)*(-b*f + c*e) - e*x**S(2)*(S(2)*a*f - b*e + S(2)*c*d)), x), x, (S(1) + x*(e + sqrt(-S(4)*d*f + e**S(2)))/(S(2)*d))/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule424)

    def f425(e, g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-a*e + b*d), ZeroQ(2*d*h - e*g), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern425 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f425))
    rule425 = ReplacementRule(pattern425, lambda e, g, a, d, h, c, x, f, b : g*Subst(Int(S(1)/(a + x**S(2)*(-a*f + c*d)), x), x, x/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule425)

    def f426(e, g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), ZeroQ(-a*e + b*d), NonzeroQ(2*d*h - e*g), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern426 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f426))
    rule426 = ReplacementRule(pattern426, lambda e, g, a, d, h, c, x, f, b : h*Int((S(2)*d + e*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e - (S(2)*d*h - e*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e)
    rubi.add(rule426)

    def f427(e, a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(-a*e + b*d), ZeroQ(g**2*(-b*f + c*e) - 2*g*h*(-a*f + c*d) + h**2*(-a*e + b*d)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern427 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f427))
    rule427 = ReplacementRule(pattern427, lambda e, a, g, h, c, d, f, x, b : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g) - x**S(2)*(-a*e + b*d), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule427)

    def f428(e, g, a, h, c, d, f, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*e*h**2 - c*e*g**2 + 2*g*h*(-a*f + c*d)), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern428 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f428))
    rule428 = ReplacementRule(pattern428, lambda e, g, a, h, c, d, f, x : -S(2)*a*g*h*Subst(Int(S(1)/Simp(S(2)*a**S(2)*c*g*h + a*e*x**S(2), x), x), x, Simp(a*h - c*g*x, x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule428)

    def f429(a, g, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), ZeroQ(b*d*h**2 - b*f*g**2 - 2*g*h*(-a*f + c*d)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern429 = Pattern(Integral((g_ + x_*WC('h', S(1)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f429))
    rule429 = ReplacementRule(pattern429, lambda a, g, d, h, c, x, f, b : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(-b*d*x**S(2) + g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + f*x**S(2))))
    rubi.add(rule429)

    def f430(e, g, a, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), PosQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern430 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f430), )
    def With430(e, g, a, h, c, d, f, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule430 = ReplacementRule(pattern430, lambda e, g, a, h, c, d, f, x, b : With430(e, g, a, h, c, d, f, x, b))
    rubi.add(rule430)

    def f431(e, g, a, h, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), PosQ(-a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern431 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f431), )
    def With431(e, g, a, h, c, d, f, x):
        q = Rt(-a*c, S(2))
        return (-c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x + q)*sqrt(d + e*x + f*x**S(2))), x) + (c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x - q)*sqrt(d + e*x + f*x**S(2))), x)
    rule431 = ReplacementRule(pattern431, lambda e, g, a, h, c, d, f, x : With431(e, g, a, h, c, d, f, x))
    rubi.add(rule431)

    def f432(g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), PosQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern432 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f432), )
    def With432(g, a, d, h, c, x, f, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule432 = ReplacementRule(pattern432, lambda g, a, d, h, c, x, f, b : With432(g, a, d, h, c, x, f, b))
    rubi.add(rule432)

    def f433(e, a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), NonzeroQ(-a*e + b*d), NegQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern433 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f433), )
    def With433(e, a, g, h, c, d, f, x, b):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-g*(-a*f + c*d - q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d + q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-g*(-a*f + c*d + q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d - q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule433 = ReplacementRule(pattern433, lambda e, a, g, h, c, d, f, x, b : With433(e, a, g, h, c, d, f, x, b))
    rubi.add(rule433)

    def f434(e, g, a, h, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), NegQ(-a*c), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern434 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f434), )
    def With434(e, g, a, h, c, d, f, x):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-a*e*h - g*(-a*f + c*d - q) + x*(-c*e*g + h*(-a*f + c*d + q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-a*e*h - g*(-a*f + c*d + q) + x*(-c*e*g + h*(-a*f + c*d - q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule434 = ReplacementRule(pattern434, lambda e, g, a, h, c, d, f, x : With434(e, g, a, h, c, d, f, x))
    rubi.add(rule434)

    def f435(a, g, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NegQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern435 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_),CustomConstraint(f435), )
    def With435(a, g, d, h, c, x, f, b):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(b*d*h - g*(-a*f + c*d - q) + x*(b*f*g + h*(-a*f + c*d + q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) - Int(Simp(b*d*h - g*(-a*f + c*d + q) + x*(b*f*g + h*(-a*f + c*d - q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule435 = ReplacementRule(pattern435, lambda a, g, d, h, c, x, f, b : With435(a, g, d, h, c, x, f, b))
    rubi.add(rule435)

    def f436(e, g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern436 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f436), )
    def With436(e, g, a, d, h, c, x, f, b):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f + e**S(2), S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)), x)/(sqrt(a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2)))
    rule436 = ReplacementRule(pattern436, lambda e, g, a, d, h, c, x, f, b : With436(e, g, a, d, h, c, x, f, b))
    rubi.add(rule436)

    def f437(g, a, d, h, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern437 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f437), )
    def With437(g, a, d, h, c, x, f, b):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f, S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)), x)/(sqrt(d + f*x**S(2))*sqrt(a + b*x + c*x**S(2)))
    rule437 = ReplacementRule(pattern437, lambda g, a, d, h, c, x, f, b : With437(g, a, d, h, c, x, f, b))
    rubi.add(rule437)

    def f438(e, a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-b*f + c*e), ZeroQ(c**2*d - f*(-3*a*c + b**2)), ZeroQ(9*a*c*h**2 - 2*b**2*h**2 - b*c*g*h + c**2*g**2), PositiveQ(-9*c*h**2/(-b*h + 2*c*g)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern438 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f438), )
    def With438(e, a, g, h, c, d, f, x, b):
        q = S(3)**(S(2)/3)*(-c*h**S(2)/(-b*h + S(2)*c*g)**S(2))**(S(1)/3)
        return sqrt(S(3))*h*q*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3)/(S(3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3)) + sqrt(S(3))/S(3))/f - S(3)*h*q*log((-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3) + S(2)**(S(1)/3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3))/(S(2)*f) + h*q*log(d + e*x + f*x**S(2))/(S(2)*f)
    rule438 = ReplacementRule(pattern438, lambda e, a, g, h, c, d, f, x, b : With438(e, a, g, h, c, d, f, x, b))
    rubi.add(rule438)

    def f439(g, a, d, h, c, x, f):
        return functools.reduce(operator.and_, [ ZeroQ(3*a*f + c*d), ZeroQ(9*a*h**2 + c*g**2), PositiveQ(a), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern439 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f439))
    rule439 = ReplacementRule(pattern439, lambda g, a, d, h, c, x, f : S(2)**(S(1)/3)*sqrt(S(3))*h*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(S(1) - S(3)*h*x/g)**(S(2)/3)/(S(3)*(S(1) + S(3)*h*x/g)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*a**(S(1)/3)*f) + S(2)**(S(1)/3)*h*log(d + f*x**S(2))/(S(4)*a**(S(1)/3)*f) - S(3)*S(2)**(S(1)/3)*h*log((S(1) - S(3)*h*x/g)**(S(2)/3) + S(2)**(S(1)/3)*(S(1) + S(3)*h*x/g)**(S(1)/3))/(S(4)*a**(S(1)/3)*f))
    rubi.add(rule439)

    def f440(e, a, g, h, c, d, f, x, b):
        return functools.reduce(operator.and_, [ ZeroQ(-b*f + c*e), ZeroQ(c**2*d - f*(-3*a*c + b**2)), ZeroQ(9*a*c*h**2 - 2*b**2*h**2 - b*c*g*h + c**2*g**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern440 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f440), )
    def With440(e, a, g, h, c, d, f, x, b):
        q = -c/(-S(4)*a*c + b**S(2))
        return (q*(a + b*x + c*x**S(2)))**(S(1)/3)*Int((g + h*x)/((d + e*x + f*x**S(2))*(a*q + b*q*x + c*q*x**S(2))**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule440 = ReplacementRule(pattern440, lambda e, a, g, h, c, d, f, x, b : With440(e, a, g, h, c, d, f, x, b))
    rubi.add(rule440)

    def f441(g, a, d, h, c, x, f):
        return functools.reduce(operator.and_, [ ZeroQ(3*a*f + c*d), ZeroQ(9*a*h**2 + c*g**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x)])
    pattern441 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f441))
    rule441 = ReplacementRule(pattern441, lambda g, a, d, h, c, x, f : (S(1) + c*x**S(2)/a)**(S(1)/3)*Int((g + h*x)/((S(1) + c*x**S(2)/a)**(S(1)/3)*(d + f*x**S(2))), x)/(a + c*x**S(2))**(S(1)/3))
    rubi.add(rule441)

    def f442(e, a, g, h, q, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, h, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x), FreeQ(q, x)])
    pattern442 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f442))
    rule442 = ReplacementRule(pattern442, lambda e, a, g, h, q, c, d, p, f, x, b : Int((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule442)

    def f443(e, a, g, h, q, c, d, p, f, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e, f, g, h, p, q), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(p, x), FreeQ(q, x)])
    pattern443 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f443))
    rule443 = ReplacementRule(pattern443, lambda e, a, g, h, q, c, d, p, f, x : Int((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule443)

    def f444(u, e, m, a, g, h, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern444 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f444))
    rule444 = ReplacementRule(pattern444, lambda u, e, m, a, g, h, x, c, d, p, f, q, b : Subst(Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule444)

    def f445(u, e, m, a, g, h, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern445 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f445))
    rule445 = ReplacementRule(pattern445, lambda u, e, m, a, g, h, x, c, d, p, f, q : Subst(Int((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule445)

    def f446(u, m, v, z, x, p, q):
        return functools.reduce(operator.and_, [ LinearQ(z, x), QuadraticQ(List(u, v), x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern446 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*z_**WC('m', S(1)), x_),CustomConstraint(f446))
    rule446 = ReplacementRule(pattern446, lambda u, m, v, z, x, p, q : Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**m, x))
    rubi.add(rule446)

    def f447(e, m, g, i, a, d, h, x, c, n, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(d*g + e*f), ZeroQ(m - n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(i, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern447 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f447))
    rule447 = ReplacementRule(pattern447, lambda e, m, g, i, a, d, h, x, c, n, p, f, q, b : Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule447)

    def f448(e, m, g, i, a, h, x, c, d, p, n, f, q, b):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p), NegativeIntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(i, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern448 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f448))
    rule448 = ReplacementRule(pattern448, lambda e, m, g, i, a, h, x, c, d, p, n, f, q, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(h + i*x)**q*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule448)

    def f449(e, m, g, i, a, h, x, c, d, p, n, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(d*g + e*f), ZeroQ(m - n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(h, x), FreeQ(i, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern449 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f449))
    rule449 = ReplacementRule(pattern449, lambda e, m, g, i, a, h, x, c, d, p, n, f, q, b : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule449)

    def f450(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern450 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f450))
    rule450 = ReplacementRule(pattern450, lambda e, a, C, d, A, B, c, x, p, f, q, b : (c/f)**p*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule450)

    def f451(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern451 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f451))
    rule451 = ReplacementRule(pattern451, lambda e, a, C, d, A, c, x, p, f, q, b : (c/f)**p*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule451)

    def f452(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern452 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f452))
    rule452 = ReplacementRule(pattern452, lambda e, a, C, d, A, B, c, x, p, f, q, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule452)

    def f453(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-a*f + c*d), ZeroQ(-a*e + b*d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern453 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f453))
    rule453 = ReplacementRule(pattern453, lambda e, a, C, d, A, c, x, p, f, q, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule453)

    def f454(e, a, C, A, x, B, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern454 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f454))
    rule454 = ReplacementRule(pattern454, lambda e, a, C, A, x, B, c, d, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule454)

    def f455(e, a, C, A, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern455 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f455))
    rule455 = ReplacementRule(pattern455, lambda e, a, C, A, x, c, d, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule455)

    def f456(a, C, A, x, B, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern456 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f456))
    rule456 = ReplacementRule(pattern456, lambda a, C, A, x, B, c, d, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(A + B*x + C*x**S(2)), x))
    rubi.add(rule456)

    def f457(a, C, A, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ ZeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern457 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f457))
    rule457 = ReplacementRule(pattern457, lambda a, C, A, x, c, d, p, f, q, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule457)

    def f458(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern458 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f458))
    rule458 = ReplacementRule(pattern458, lambda e, a, C, d, A, B, c, x, p, f, q, b : Int(ExpandIntegrand((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule458)

    def f459(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), Greater(p, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern459 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f459))
    rule459 = ReplacementRule(pattern459, lambda e, a, C, d, A, c, x, p, f, q, b : Int(ExpandIntegrand((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule459)

    def f460(e, a, C, d, A, B, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern460 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f460))
    rule460 = ReplacementRule(pattern460, lambda e, a, C, d, A, B, c, x, p, f, q : Int(ExpandIntegrand((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule460)

    def f461(e, a, C, d, A, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), IntegersQ(p, q), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern461 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f461))
    rule461 = ReplacementRule(pattern461, lambda e, a, C, d, A, c, x, p, f, q : Int(ExpandIntegrand((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule461)

    def f462(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern462 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f462))
    rule462 = ReplacementRule(pattern462, lambda e, a, C, d, A, B, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + e*q*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))) + x*(-e*(C*(S(2)*a*c*(q + S(1)) - b**S(2)*(p + q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + q + S(3))) + S(2)*f*q*(A*b*c - S(2)*B*a*c + C*a*b)), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule462)

    def f463(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern463 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f463))
    rule463 = ReplacementRule(pattern463, lambda e, a, C, d, A, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*(b*e*q + S(2)*c*d*(S(2)*p + S(3))) - C*(-a*b*e*q + S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*c*(b*f*q + c*e*(S(2)*p + q + S(3))) + C*(S(2)*a*b*f*q - S(2)*a*c*e*(q + S(1)) + b**S(2)*e*(p + q + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule463)

    def f464(e, a, C, d, A, B, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern464 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f464))
    rule464 = ReplacementRule(pattern464, lambda e, a, C, d, A, B, c, x, p, f, q : (a + c*x**S(2))**(p + S(1))*(B*a - x*(A*c - C*a))*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - a*(B*e*q + C*d) - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - a*(S(2)*B*f*q + C*e*(q + S(1)))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule464)

    def f465(e, a, C, d, A, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern465 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f465))
    rule465 = ReplacementRule(pattern465, lambda e, a, C, d, A, c, x, p, f, q : -x*(a + c*x**S(2))**(p + S(1))*(A*c - C*a)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - C*a*d - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - C*a*e*(q + S(1))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule465)

    def f466(a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern466 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f466))
    rule466 = ReplacementRule(pattern466, lambda a, C, d, A, B, c, x, p, f, q, b : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + S(2)*f*q*x*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule466)

    def f467(a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p, q), Less(p, -1), Greater(q, 0), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern467 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f467))
    rule467 = ReplacementRule(pattern467, lambda a, C, d, A, c, x, p, f, q, b : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*A*c**S(2)*d*(S(2)*p + S(3)) - C*(S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*b*c*f*q + S(2)*C*a*b*f*q), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule467)

    def f468(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern468 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f468))
    rule468 = ReplacementRule(pattern468, lambda e, a, C, d, A, B, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - B*(a*b*f - S(2)*a*c*e + b*c*d) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)) - e*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e))) + (p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule468)

    def f469(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern469 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f469))
    rule469 = ReplacementRule(pattern469, lambda e, a, C, d, A, c, x, p, f, q, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)) - e*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule469)

    def f470(e, a, C, d, A, B, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(a*c*e**2 + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern470 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f470))
    rule470 = ReplacementRule(pattern470, lambda e, a, C, d, A, B, c, x, p, f, q : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*B*a*c*e - S(2)*C*a*(-a*f + c*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - e*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2)) - x*(c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + q + S(4)) + S(2)*f*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule470)

    def f471(e, a, C, d, A, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Less(p, -1), NonzeroQ(a*c*e**2 + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern471 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f471))
    rule471 = ReplacementRule(pattern471, lambda e, a, C, d, A, c, x, p, f, q : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) - S(2)*C*a*(-a*f + c*d)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-S(2)*a*c*e**S(2)*(A*c - C*a)*(p + q + S(2)) - c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - x*(S(4)*a*c*e*f*(A*c - C*a)*(p + q + S(2)) + c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + q + S(4))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule471)

    def f472(a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), NonzeroQ(b**2*d*f + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern472 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f472))
    rule472 = ReplacementRule(pattern472, lambda a, C, d, A, B, c, x, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - B*(a*b*f + b*c*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) + S(2)*f*(-b*(A*c - C*a)*(a*f + c*d) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule472)

    def f473(a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Less(p, -1), NonzeroQ(b**2*d*f + (-a*f + c*d)**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern473 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f473))
    rule473 = ReplacementRule(pattern473, lambda a, C, d, A, c, x, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) + S(2)*f*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule473)

    def f474(e, a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern474 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f474))
    rule474 = ReplacementRule(pattern474, lambda e, a, C, d, A, B, c, x, p, f, q, b : (a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) + S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(p*(-b*f + c*e)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule474)

    def f475(e, a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern475 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f475))
    rule475 = ReplacementRule(pattern475, lambda e, a, C, d, A, c, x, p, f, q, b : (S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + x**S(2)*(p*(-b*f + c*e)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2)))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule475)

    def f476(e, a, C, d, A, B, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern476 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f476))
    rule476 = ReplacementRule(pattern476, lambda e, a, C, d, A, B, c, x, p, f, q : (a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule476)

    def f477(e, a, C, d, A, c, x, p, f, q):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern477 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_),CustomConstraint(f477))
    rule477 = ReplacementRule(pattern477, lambda e, a, C, d, A, c, x, p, f, q : (a + c*x**S(2))**p*(-C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule477)

    def f478(a, C, d, A, B, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(q, x)])
    pattern478 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_),CustomConstraint(f478))
    rule478 = ReplacementRule(pattern478, lambda a, C, d, A, B, c, x, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p*(B*c*f*(S(2)*p + S(2)*q + S(3)) + C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(b*d*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + x**S(2)*(-b*f*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1)))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule478)

    def f479(a, C, d, A, c, x, p, f, q, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), RationalQ(p), Greater(p, 0), NonzeroQ(p + q + 1), NonzeroQ(2*p + 2*q + 3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(q, x)])
    pattern479 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f479))
    rule479 = ReplacementRule(pattern479, lambda a, C, d, A, c, x, p, f, q, b : (d + f*x**S(2))**(q + S(1))*(C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))*(a + b*x + c*x**S(2))**p/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-C*b**S(2)*d*f*p*(q + S(1)) + x**S(2)*(C*b**S(2)*f**S(2)*p*(q + S(1)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(2)*C*b*f*p*(q + S(1))*(-a*f + c*d) - b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule479)

    def f480(e, a, C, d, A, B, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern480 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f480), )
    def With480(e, a, C, d, A, B, c, x, f, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d - B*a*b*f + B*a*c*e + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) + B*b*d*f - B*c*d*e - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule480 = ReplacementRule(pattern480, lambda e, a, C, d, A, B, c, x, f, b : With480(e, a, C, d, A, B, c, x, f, b))
    rubi.add(rule480)

    def f481(e, a, C, d, A, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern481 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_),CustomConstraint(f481), )
    def With481(e, a, C, d, A, c, x, f, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule481 = ReplacementRule(pattern481, lambda e, a, C, d, A, c, x, f, b : With481(e, a, C, d, A, c, x, f, b))
    rubi.add(rule481)

    def f482(a, C, d, A, B, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern482 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f482), )
    def With482(a, C, d, A, B, c, x, f, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f + B*b*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d - B*a*b*f + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule482 = ReplacementRule(pattern482, lambda a, C, d, A, B, c, x, f, b : With482(a, C, d, A, B, c, x, f, b))
    rubi.add(rule482)

    def f483(a, C, d, A, c, x, f, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern483 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f483), )
    def With483(a, C, d, A, c, x, f, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule483 = ReplacementRule(pattern483, lambda a, C, d, A, c, x, f, b : With483(a, C, d, A, c, x, f, b))
    rubi.add(rule483)

    def f484(e, a, C, A, B, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern484 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f484))
    rule484 = ReplacementRule(pattern484, lambda e, a, C, A, B, c, d, f, x, b : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule484)

    def f485(e, a, C, A, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern485 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f485))
    rule485 = ReplacementRule(pattern485, lambda e, a, C, A, c, d, f, x, b : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule485)

    def f486(e, a, C, A, B, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern486 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f486))
    rule486 = ReplacementRule(pattern486, lambda e, a, C, A, B, c, d, f, x : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c + B*c*x - C*a)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule486)

    def f487(e, a, C, A, c, d, f, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*d*f + e**2), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern487 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_),CustomConstraint(f487))
    rule487 = ReplacementRule(pattern487, lambda e, a, C, A, c, d, f, x : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + (A*c - C*a)*Int(S(1)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule487)

    def f488(a, C, A, B, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x)])
    pattern488 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f488))
    rule488 = ReplacementRule(pattern488, lambda a, C, A, B, c, d, f, x, b : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule488)

    def f489(a, C, A, c, d, f, x, b):
        return functools.reduce(operator.and_, [ NonzeroQ(-4*a*c + b**2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x)])
    pattern489 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_),CustomConstraint(f489))
    rule489 = ReplacementRule(pattern489, lambda a, C, A, c, d, f, x, b : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule489)

    def f490(e, a, C, A, q, B, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, A, B, C, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern490 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f490))
    rule490 = ReplacementRule(pattern490, lambda e, a, C, A, q, B, c, d, p, f, x, b : Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule490)

    def f491(e, a, C, A, q, c, d, p, f, x, b):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, A, C, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern491 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f491))
    rule491 = ReplacementRule(pattern491, lambda e, a, C, A, q, c, d, p, f, x, b : Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule491)

    def f492(e, a, C, A, q, B, c, d, p, f, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e, f, A, B, C, p, q), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern492 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f492))
    rule492 = ReplacementRule(pattern492, lambda e, a, C, A, q, B, c, d, p, f, x : Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule492)

    def f493(e, a, C, A, q, c, d, p, f, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, c, d, e, f, A, C, p, q), x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern493 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_),CustomConstraint(f493))
    rule493 = ReplacementRule(pattern493, lambda e, a, C, A, q, c, d, p, f, x : Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule493)

    def f494(u, e, a, C, A, x, B, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern494 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f494))
    rule494 = ReplacementRule(pattern494, lambda u, e, a, C, A, x, B, c, d, p, f, q, b : Subst(Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule494)

    def f495(u, e, a, A, x, B, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern495 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f495))
    rule495 = ReplacementRule(pattern495, lambda u, e, a, A, x, B, c, d, p, f, q, b : Subst(Int((A + B*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule495)

    def f496(u, e, a, C, A, x, c, d, p, f, q, b):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern496 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f496))
    rule496 = ReplacementRule(pattern496, lambda u, e, a, C, A, x, c, d, p, f, q, b : Subst(Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule496)

    def f497(u, e, a, C, A, x, B, c, d, p, f, q):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern497 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f497))
    rule497 = ReplacementRule(pattern497, lambda u, e, a, C, A, x, B, c, d, p, f, q : Subst(Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule497)

    def f498(u, e, a, A, x, B, c, d, p, f, q):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(B, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern498 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f498))
    rule498 = ReplacementRule(pattern498, lambda u, e, a, A, x, B, c, d, p, f, q : Subst(Int((A + B*x)*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule498)

    def f499(u, e, a, C, A, x, c, d, p, f, q):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(A, x), FreeQ(C, x), FreeQ(p, x), FreeQ(q, x)])
    pattern499 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_),CustomConstraint(f499))
    rule499 = ReplacementRule(pattern499, lambda u, e, a, C, A, x, c, d, p, f, q : Subst(Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule499)

    return rubi
