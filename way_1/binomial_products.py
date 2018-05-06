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

    from sympy.integrals.rubi.constraints import (constraint_freeq_a, constraint_freeq_b, constraint_freeq_c, constraint_freeq_d, constraint_freeq_e, constraint_freeq_f,
        constraint_freeq_g, constraint_freeq_m, constraint_freeq_n, constraint_freeq_p)

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False


    import functools, operator
    from sympy import And, Or

def binomial_products(rubi):
    def f1(b, p, x, n):
        return functools.reduce(operator.and_, [ FreeQ(List(b, n, p), x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern1 = Pattern(Integral((x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f1))
    rule1 = ReplacementRule(pattern1, lambda b, p, x, n : b**IntPart(p)*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(n*p), x))
    rubi.add(rule1)

    def f2(b, a, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(p + S(1) + S(1)/n), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern2 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f2))
    rule2 = ReplacementRule(pattern2, lambda b, a, n, p, x : x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule2)

    def f3(b, a, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(p + S(1) + S(1)/n), NonzeroQ(p + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern3 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f3))
    rule3 = ReplacementRule(pattern3, lambda b, a, n, p, x : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule3)

    def f4(b, a, x, n):
        return functools.reduce(operator.and_, [ NonzeroQ(S(3)*n + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x)])
    pattern4 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**S(2), x_),CustomConstraint(f4))
    rule4 = ReplacementRule(pattern4, lambda b, a, x, n : Int(a**S(2) + S(2)*a*b*x**n + b**S(2)*x**(S(2)*n), x))
    rubi.add(rule4)

    def f5(b, a, n, p, x):
        return functools.reduce(operator.and_, [ RationalQ(n), Less(n, S(0)), IntegerQ(p), FreeQ(a, x), FreeQ(b, x)])
    pattern5 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f5))
    rule5 = ReplacementRule(pattern5, lambda b, a, n, p, x : Int(x**(n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule5)

    def f6(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n, p), FreeQ(a, x), FreeQ(b, x)])
    pattern6 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f6))
    rule6 = ReplacementRule(pattern6, lambda b, a, n, p, x : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule6)

    def f7(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Greater(p, S(0)), Or(IntegerQ(S(2)*p), And(Equal(n, S(2)), IntegerQ(S(4)*p)), And(Equal(n, S(2)), IntegerQ(S(3)*p)), Less(Denominator(p + S(1)/n), Denominator(p))), FreeQ(a, x), FreeQ(b, x)])
    pattern7 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f7))
    rule7 = ReplacementRule(pattern7, lambda b, a, n, p, x : a*n*p*Int((a + b*x**n)**(p + S(-1)), x)/(n*p + S(1)) + x*(a + b*x**n)**p/(n*p + S(1)))
    rubi.add(rule7)

    def f8(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), PositiveQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern8 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_),CustomConstraint(f8))
    rule8 = ReplacementRule(pattern8, lambda b, a, x : S(2)*EllipticE(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(5)/4)*Rt(b/a, S(2))))
    rubi.add(rule8)

    def f9(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x)])
    pattern9 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-5)/4), x_),CustomConstraint(f9))
    rule9 = ReplacementRule(pattern9, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-5)/4), x)/(a*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule9)

    def f10(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern10 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-7)/6), x_),CustomConstraint(f10))
    rule10 = ReplacementRule(pattern10, lambda b, a, x : Subst(Int((-b*x**S(2) + S(1))**(S(-1)/3), x), x, x/sqrt(a + b*x**S(2)))/((a/(a + b*x**S(2)))**(S(2)/3)*(a + b*x**S(2))**(S(2)/3)))
    rubi.add(rule10)

    def f11(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(p, S(-1)), Or(IntegerQ(S(2)*p), And(Equal(n, S(2)), IntegerQ(S(4)*p)), And(Equal(n, S(2)), IntegerQ(S(3)*p)), Less(Denominator(p + S(1)/n), Denominator(p))), FreeQ(a, x), FreeQ(b, x)])
    pattern11 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f11))
    rule11 = ReplacementRule(pattern11, lambda b, a, n, p, x : -x*(a + b*x**n)**(p + S(1))/(a*n*(p + S(1))) + (n*(p + S(1)) + S(1))*Int((a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))))
    rubi.add(rule11)

    def f12(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern12 = Pattern(Integral(S(1)/(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f12))
    rule12 = ReplacementRule(pattern12, lambda b, a, x : Int((-x*Rt(b, S(3)) + S(2)*Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))**S(2)) + Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))**S(2)))
    rubi.add(rule12)

    def f13(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(2) + S(-3)/2), PosQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern13 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f13), )
    def With13(b, a, x, n):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r + s*x), x)/(a*n)
    rule13 = ReplacementRule(pattern13, lambda b, a, x, n : With13(b, a, x, n))
    rubi.add(rule13)

    def f14(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(2) + S(-3)/2), NegQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern14 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f14), )
    def With14(b, a, x, n):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r*Int(S(1)/(r - s*x), x)/(a*n)
    rule14 = ReplacementRule(pattern14, lambda b, a, x, n : With14(b, a, x, n))
    rubi.add(rule14)

    def f15(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(a/b), Or(PositiveQ(a), PositiveQ(b)), FreeQ(a, x), FreeQ(b, x)])
    pattern15 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f15))
    rule15 = ReplacementRule(pattern15, lambda b, a, x : ArcTan(x*Rt(b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(b, S(2))))
    rubi.add(rule15)

    def f16(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(a/b), Or(NegativeQ(a), NegativeQ(b)), FreeQ(a, x), FreeQ(b, x)])
    pattern16 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f16))
    rule16 = ReplacementRule(pattern16, lambda b, a, x : -ArcTan(x*Rt(-b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(-b, S(2))))
    rubi.add(rule16)

    def f17(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern17 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f17))
    rule17 = ReplacementRule(pattern17, lambda b, a, x : ArcTan(x/Rt(a/b, S(2)))*Rt(a/b, S(2))/a)
    rubi.add(rule17)

    def f18(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(a/b), Or(PositiveQ(a), NegativeQ(b)), FreeQ(a, x), FreeQ(b, x)])
    pattern18 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f18))
    rule18 = ReplacementRule(pattern18, lambda b, a, x : atanh(x*Rt(-b, S(2))/Rt(a, S(2)))/(Rt(a, S(2))*Rt(-b, S(2))))
    rubi.add(rule18)

    def f19(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(a/b), Or(NegativeQ(a), PositiveQ(b)), FreeQ(a, x), FreeQ(b, x)])
    pattern19 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f19))
    rule19 = ReplacementRule(pattern19, lambda b, a, x : -atanh(x*Rt(b, S(2))/Rt(-a, S(2)))/(Rt(-a, S(2))*Rt(b, S(2))))
    rubi.add(rule19)

    def f20(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern20 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f20))
    rule20 = ReplacementRule(pattern20, lambda b, a, x : Rt(-a/b, S(2))*atanh(x/Rt(-a/b, S(2)))/a)
    rubi.add(rule20)

    def f21(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(4) + S(-1)/2), PosQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern21 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f21), )
    def With21(b, a, x, n):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        v = Symbol('v')
        u = Int((r - s*x*cos(Pi*(2*k - 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r + s*x*cos(Pi*(2*k - 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n)
    rule21 = ReplacementRule(pattern21, lambda b, a, x, n : With21(b, a, x, n))
    rubi.add(rule21)

    def f22(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(4) + S(-1)/2), NegQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern22 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f22), )
    def With22(b, a, x, n):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r - s*x*cos(2*Pi*k/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r + s*x*cos(2*Pi*k/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**S(2)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    rule22 = ReplacementRule(pattern22, lambda b, a, x, n : With22(b, a, x, n))
    rubi.add(rule22)

    def f23(b, a, x):
        return functools.reduce(operator.and_, [ Or(PositiveQ(a/b), And(PosQ(a/b), AtomQ(SplitProduct(SumBaseQ, a)), AtomQ(SplitProduct(SumBaseQ, b)))), FreeQ(a, x), FreeQ(b, x)])
    pattern23 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f23), )
    def With23(b, a, x):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*r)
    rule23 = ReplacementRule(pattern23, lambda b, a, x : With23(b, a, x))
    rubi.add(rule23)

    def f24(b, a, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(a/b)), FreeQ(a, x), FreeQ(b, x)])
    pattern24 = Pattern(Integral(S(1)/(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f24), )
    def With24(b, a, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*a)
    rule24 = ReplacementRule(pattern24, lambda b, a, x : With24(b, a, x))
    rubi.add(rule24)

    def f25(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(4) + S(-1)), PositiveQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern25 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f25), )
    def With25(b, a, x, n):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*r*Int((sqrt(S(2))*r - s*x**(n/S(4)))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a) + sqrt(S(2))*r*Int((sqrt(S(2))*r + s*x**(n/S(4)))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*a)
    rule25 = ReplacementRule(pattern25, lambda b, a, x, n : With25(b, a, x, n))
    rubi.add(rule25)

    def f26(b, a, x, n):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(4) + S(-1)), Not(PositiveQ(a/b)), FreeQ(a, x), FreeQ(b, x)])
    pattern26 = Pattern(Integral(S(1)/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f26), )
    def With26(b, a, x, n):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(S(1)/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(S(1)/(r + s*x**(n/S(2))), x)/(S(2)*a)
    rule26 = ReplacementRule(pattern26, lambda b, a, x, n : With26(b, a, x, n))
    rubi.add(rule26)

    def f27(b, a, x):
        return functools.reduce(operator.and_, [ PositiveQ(a), PosQ(b), FreeQ(a, x), FreeQ(b, x)])
    pattern27 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f27))
    rule27 = ReplacementRule(pattern27, lambda b, a, x : asinh(x*Rt(b, S(2))/sqrt(a))/Rt(b, S(2)))
    rubi.add(rule27)

    def f28(b, a, x):
        return functools.reduce(operator.and_, [ PositiveQ(a), NegQ(b), FreeQ(a, x), FreeQ(b, x)])
    pattern28 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f28))
    rule28 = ReplacementRule(pattern28, lambda b, a, x : asin(x*Rt(-b, S(2))/sqrt(a))/Rt(-b, S(2)))
    rubi.add(rule28)

    def f29(b, a, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x)])
    pattern29 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f29))
    rule29 = ReplacementRule(pattern29, lambda b, a, x : Subst(Int(S(1)/(-b*x**S(2) + S(1)), x), x, x/sqrt(a + b*x**S(2))))
    rubi.add(rule29)

    def f30(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern30 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f30), )
    def With30(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(-sqrt(S(3)) + S(1)))/(r*x + s*(S(1) + sqrt(S(3))))), S(-7) - S(4)*sqrt(S(3)))/(S(3)*r*sqrt(s*(r*x + s)/(r*x + s*(S(1) + sqrt(S(3))))**S(2))*sqrt(a + b*x**S(3)))
    rule30 = ReplacementRule(pattern30, lambda b, a, x : With30(b, a, x))
    rubi.add(rule30)

    def f31(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern31 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f31), )
    def With31(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(2)*S(3)**(S(3)/4)*sqrt((r**S(2)*x**S(2) - r*s*x + s**S(2))/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(-sqrt(S(3)) + S(2))*(r*x + s)*EllipticF(asin((r*x + s*(S(1) + sqrt(S(3))))/(r*x + s*(-sqrt(S(3)) + S(1)))), S(-7) + S(4)*sqrt(S(3)))/(S(3)*r*sqrt(-s*(r*x + s)/(r*x + s*(-sqrt(S(3)) + S(1)))**S(2))*sqrt(a + b*x**S(3)))
    rule31 = ReplacementRule(pattern31, lambda b, a, x : With31(b, a, x))
    rubi.add(rule31)

    def f32(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern32 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f32), )
    def With32(b, a, x):
        q = Rt(b/a, S(4))
        return sqrt((a + b*x**S(4))/(a*(q**S(2)*x**S(2) + S(1))**S(2)))*(q**S(2)*x**S(2) + S(1))*EllipticF(S(2)*ArcTan(q*x), S(1)/2)/(S(2)*q*sqrt(a + b*x**S(4)))
    rule32 = ReplacementRule(pattern32, lambda b, a, x : With32(b, a, x))
    rubi.add(rule32)

    def f33(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), PositiveQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern33 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f33))
    rule33 = ReplacementRule(pattern33, lambda b, a, x : EllipticF(asin(x*Rt(-b, S(4))/Rt(a, S(4))), S(-1))/(Rt(a, S(4))*Rt(-b, S(4))))
    rubi.add(rule33)

    def f34(b, a, x):
        return functools.reduce(operator.and_, [ NegativeQ(a), PositiveQ(b), FreeQ(a, x), FreeQ(b, x)])
    pattern34 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f34), )
    def With34(b, a, x):
        q = Rt(-a*b, S(2))
        if IntegerQ(q):
            return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt(-a + q*x**S(2))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(-a)*sqrt(a + b*x**S(4)))
        print("Unable to Integrate")
    rule34 = ReplacementRule(pattern34, lambda b, a, x : With34(b, a, x))
    rubi.add(rule34)

    def f35(b, a, x):
        return functools.reduce(operator.and_, [ NegativeQ(a), PositiveQ(b), FreeQ(a, x), FreeQ(b, x)])
    pattern35 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f35), )
    def With35(b, a, x):
        q = Rt(-a*b, S(2))
        return sqrt(S(2))*sqrt((a + q*x**S(2))/q)*sqrt((a - q*x**S(2))/(a + q*x**S(2)))*EllipticF(asin(sqrt(S(2))*x/sqrt((a + q*x**S(2))/q)), S(1)/2)/(S(2)*sqrt(a/(a + q*x**S(2)))*sqrt(a + b*x**S(4)))
    rule35 = ReplacementRule(pattern35, lambda b, a, x : With35(b, a, x))
    rubi.add(rule35)

    def f36(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x)])
    pattern36 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f36))
    rule36 = ReplacementRule(pattern36, lambda b, a, x : sqrt(S(1) + b*x**S(4)/a)*Int(S(1)/sqrt(S(1) + b*x**S(4)/a), x)/sqrt(a + b*x**S(4)))
    rubi.add(rule36)

    def f37(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern37 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_),CustomConstraint(f37), )
    def With37(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return S(3)**(S(3)/4)*x*sqrt((r**S(2)*x**S(4) - r*s*x**S(2) + s**S(2))/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*(r*x**S(2) + s)*EllipticF(acos((r*x**S(2)*(-sqrt(S(3)) + S(1)) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)), sqrt(S(3))/S(4) + S(1)/2)/(S(6)*s*sqrt(r*x**S(2)*(r*x**S(2) + s)/(r*x**S(2)*(S(1) + sqrt(S(3))) + s)**S(2))*sqrt(a + b*x**S(6)))
    rule37 = ReplacementRule(pattern37, lambda b, a, x : With37(b, a, x))
    rubi.add(rule37)

    def f38(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern38 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_),CustomConstraint(f38))
    rule38 = ReplacementRule(pattern38, lambda b, a, x : Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/S(2))
    rubi.add(rule38)

    def f39(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern39 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_),CustomConstraint(f39))
    rule39 = ReplacementRule(pattern39, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-5)/4), x) + S(2)*x/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule39)

    def f40(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), PositiveQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern40 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_),CustomConstraint(f40))
    rule40 = ReplacementRule(pattern40, lambda b, a, x : S(2)*EllipticE(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(1)/4)*Rt(-b/a, S(2))))
    rubi.add(rule40)

    def f41(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x)])
    pattern41 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/4), x_),CustomConstraint(f41))
    rule41 = ReplacementRule(pattern41, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(1)/4)*Int((S(1) + b*x**S(2)/a)**(S(-1)/4), x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule41)

    def f42(b, a, x):
        return functools.reduce(operator.and_, [ PositiveQ(a), PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern42 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_),CustomConstraint(f42))
    rule42 = ReplacementRule(pattern42, lambda b, a, x : S(2)*EllipticF(ArcTan(x*Rt(b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(b/a, S(2))))
    rubi.add(rule42)

    def f43(b, a, x):
        return functools.reduce(operator.and_, [ PositiveQ(a), NegQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern43 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_),CustomConstraint(f43))
    rule43 = ReplacementRule(pattern43, lambda b, a, x : S(2)*EllipticF(asin(x*Rt(-b/a, S(2)))/S(2), S(2))/(a**(S(3)/4)*Rt(-b/a, S(2))))
    rubi.add(rule43)

    def f44(b, a, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x)])
    pattern44 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-3)/4), x_),CustomConstraint(f44))
    rule44 = ReplacementRule(pattern44, lambda b, a, x : (S(1) + b*x**S(2)/a)**(S(3)/4)*Int((S(1) + b*x**S(2)/a)**(S(-3)/4), x)/(a + b*x**S(2))**(S(3)/4))
    rubi.add(rule44)

    def f45(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern45 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/3), x_),CustomConstraint(f45))
    rule45 = ReplacementRule(pattern45, lambda b, a, x : S(3)*sqrt(b*x**S(2))*Subst(Int(x/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule45)

    def f46(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern46 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-2)/3), x_),CustomConstraint(f46))
    rule46 = ReplacementRule(pattern46, lambda b, a, x : S(3)*sqrt(b*x**S(2))*Subst(Int(S(1)/sqrt(-a + x**S(3)), x), x, (a + b*x**S(2))**(S(1)/3))/(S(2)*b*x))
    rubi.add(rule46)

    def f47(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern47 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(-3)/4), x_),CustomConstraint(f47))
    rule47 = ReplacementRule(pattern47, lambda b, a, x : x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(3)/4)), x)/(a + b*x**S(4))**(S(3)/4))
    rubi.add(rule47)

    def f48(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern48 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(-1)/6), x_),CustomConstraint(f48))
    rule48 = ReplacementRule(pattern48, lambda b, a, x : -a*Int((a + b*x**S(2))**(S(-7)/6), x)/S(2) + S(3)*x/(S(2)*(a + b*x**S(2))**(S(1)/6)))
    rubi.add(rule48)

    def f49(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), IntegerQ(p + S(1)/n), FreeQ(a, x), FreeQ(b, x)])
    pattern49 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f49))
    rule49 = ReplacementRule(pattern49, lambda b, a, n, p, x : a**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule49)

    def f50(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), Less(Denominator(p + S(1)/n), Denominator(p)), FreeQ(a, x), FreeQ(b, x)])
    pattern50 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f50))
    rule50 = ReplacementRule(pattern50, lambda b, a, n, p, x : (a/(a + b*x**n))**(p + S(1)/n)*(a + b*x**n)**(p + S(1)/n)*Subst(Int((-b*x**n + S(1))**(-p + S(-1) - S(1)/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule50)

    def f51(b, a, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(p, x)])
    pattern51 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f51))
    rule51 = ReplacementRule(pattern51, lambda b, a, n, p, x : -Subst(Int((a + b*x**(-n))**p/x**S(2), x), x, S(1)/x))
    rubi.add(rule51)

    def f52(b, a, n, p, x):
        return functools.reduce(operator.and_, [ FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(p, x)])
    pattern52 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f52), )
    def With52(b, a, n, p, x):
        k = Denominator(n)
        return k*Subst(Int(x**(k + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    rule52 = ReplacementRule(pattern52, lambda b, a, n, p, x : With52(b, a, n, p, x))
    rubi.add(rule52)

    def f53(b, a, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x)])
    pattern53 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f53))
    rule53 = ReplacementRule(pattern53, lambda b, a, n, p, x : Int(ExpandIntegrand((a + b*x**n)**p, x), x))
    rubi.add(rule53)

    def f54(b, a, n, p, x):
        return functools.reduce(operator.and_, [ Not(PositiveIntegerQ(p)), Not(IntegerQ(S(1)/n)), Not(NegativeIntegerQ(p + S(1)/n)), Or(IntegerQ(p), PositiveQ(a)), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern54 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f54))
    rule54 = ReplacementRule(pattern54, lambda b, a, n, p, x : a**p*x*Hypergeometric2F1(-p, S(1)/n, S(1) + S(1)/n, -b*x**n/a))
    rubi.add(rule54)

    def f55(b, a, n, p, x):
        return functools.reduce(operator.and_, [ Not(PositiveIntegerQ(p)), Not(IntegerQ(S(1)/n)), Not(NegativeIntegerQ(p + S(1)/n)), Not(Or(IntegerQ(p), PositiveQ(a))), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern55 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f55))
    rule55 = ReplacementRule(pattern55, lambda b, a, n, p, x : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p, x))
    rubi.add(rule55)

    def f56(b, a, n, u, p, x):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern56 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f56))
    rule56 = ReplacementRule(pattern56, lambda b, a, n, u, p, x : Subst(Int((a + b*x**n)**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule56)

    def f57(b1, b2, a2, a1, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(n, x), FreeQ(p, x)])
    pattern57 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**WC('p', S(1))*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**WC('p', S(1)), x_),CustomConstraint(f57))
    rule57 = ReplacementRule(pattern57, lambda b1, b2, a2, a1, n, p, x : Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule57)

    def f58(b1, b2, a1, a2, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(p), Greater(p, S(0)), Or(IntegerQ(S(2)*p), Less(Denominator(p + S(1)/n), Denominator(p))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x)])
    pattern58 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**WC('p', S(1)), x_),CustomConstraint(f58))
    rule58 = ReplacementRule(pattern58, lambda b1, b2, a1, a2, n, p, x : S(2)*a1*a2*n*p*Int((a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(S(2)*n*p + S(1)) + x*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(S(2)*n*p + S(1)))
    rubi.add(rule58)

    def f59(b1, b2, a1, a2, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(p), Less(p, S(-1)), Or(IntegerQ(S(2)*p), Less(Denominator(p + S(1)/n), Denominator(p))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x)])
    pattern59 = Pattern(Integral((a1_ + x_**WC('n', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('n', S(1))*WC('b2', S(1)))**p_, x_),CustomConstraint(f59))
    rule59 = ReplacementRule(pattern59, lambda b1, b2, a1, a2, n, p, x : -x*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*n*(p + S(1))) + (S(2)*n*(p + S(1)) + S(1))*Int((a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))))
    rubi.add(rule59)

    def f60(b1, b2, a1, a2, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ(S(2)*n), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(p, x)])
    pattern60 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f60))
    rule60 = ReplacementRule(pattern60, lambda b1, b2, a1, a2, n, p, x : -Subst(Int((a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p/x**S(2), x), x, S(1)/x))
    rubi.add(rule60)

    def f61(b1, b2, a1, a2, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), FractionQ(S(2)*n), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(p, x)])
    pattern61 = Pattern(Integral((a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f61), )
    def With61(b1, b2, a1, a2, n, p, x):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    rule61 = ReplacementRule(pattern61, lambda b1, b2, a1, a2, n, p, x : With61(b1, b2, a1, a2, n, p, x))
    rubi.add(rule61)

    def f62(b1, b2, a2, a1, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), Not(IntegerQ(p)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(n, x), FreeQ(p, x)])
    pattern62 = Pattern(Integral((x_**n_*WC('b1', S(1)) + WC('a1', S(0)))**p_*(x_**n_*WC('b2', S(1)) + WC('a2', S(0)))**p_, x_),CustomConstraint(f62))
    rule62 = ReplacementRule(pattern62, lambda b1, b2, a2, a1, n, p, x : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule62)

    def f63(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern63 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f63))
    rule63 = ReplacementRule(pattern63, lambda b1, b2, a1, a2, m, c, n, p, x : Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule63)

    def f64(b, m, c, n, p, x):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(c)), IntegerQ((m + S(1))/n), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern64 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f64))
    rule64 = ReplacementRule(pattern64, lambda b, m, c, n, p, x : b**(S(1) - (m + S(1))/n)*c**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n), x), x, x**n)/n)
    rubi.add(rule64)

    def f65(b, m, c, n, p, x):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(c)), Not(IntegerQ((m + S(1))/n)), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern65 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_),CustomConstraint(f65))
    rule65 = ReplacementRule(pattern65, lambda b, m, c, n, p, x : b**IntPart(p)*c**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p), x))
    rubi.add(rule65)

    def f66(b, m, c, n, p, x):
        return functools.reduce(operator.and_, [ Not(IntegerQ(m)), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern66 = Pattern(Integral((c_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_, x_),CustomConstraint(f66))
    rule66 = ReplacementRule(pattern66, lambda b, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(b*x**n)**p, x))
    rubi.add(rule66)

    def f67(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(p), NegQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x)])
    pattern67 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f67))
    rule67 = ReplacementRule(pattern67, lambda b, a, m, n, p, x : Int(x**(m + n*p)*(a*x**(-n) + b)**p, x))
    rubi.add(rule67)

    def f68(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(p + S(1) + (m + S(1))/n), NonzeroQ(m + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern68 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f68))
    rule68 = ReplacementRule(pattern68, lambda b, a, m, c, n, p, x : (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule68)

    def f69(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), ZeroQ(p + S(1) + (m + S(1))/(S(2)*n)), NonzeroQ(m + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern69 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f69))
    rule69 = ReplacementRule(pattern69, lambda b1, b2, a1, a2, m, c, n, p, x : (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule69)

    def f70(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern70 = Pattern(Integral(x_**WC('m', S(1))*(x_**n_*WC('b', S(1)) + WC('a', S(0)))**p_, x_),CustomConstraint(f70))
    rule70 = ReplacementRule(pattern70, lambda b, a, m, n, p, x : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p, x), x, x**n)/n)
    rubi.add(rule70)

    def f71(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ((m + S(1))/(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern71 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f71))
    rule71 = ReplacementRule(pattern71, lambda b1, b2, a1, a2, m, n, p, x : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a1 + b1*x)**p*(a2 + b2*x)**p, x), x, x**n)/n)
    rubi.add(rule71)

    def f72(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern72 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f72))
    rule72 = ReplacementRule(pattern72, lambda b, a, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule72)

    def f73(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ((m + S(1))/(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern73 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f73))
    rule73 = ReplacementRule(pattern73, lambda b1, b2, a1, a2, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule73)

    def f74(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern74 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f74))
    rule74 = ReplacementRule(pattern74, lambda b, a, m, c, n, p, x : Int(ExpandIntegrand((c*x)**m*(a + b*x**n)**p, x), x))
    rubi.add(rule74)

    def f75(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n), NonzeroQ(m + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern75 = Pattern(Integral(x_**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f75))
    rule75 = ReplacementRule(pattern75, lambda b, a, m, n, p, x : -b*(m + n*(p + S(1)) + S(1))*Int(x**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + x**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*(m + S(1))))
    rubi.add(rule75)

    def f76(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)), NonzeroQ(m + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern76 = Pattern(Integral(x_**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f76))
    rule76 = ReplacementRule(pattern76, lambda b1, b2, a1, a2, m, n, p, x : -b1*b2*(m + S(2)*n*(p + S(1)) + S(1))*Int(x**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + x**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*(m + S(1))))
    rubi.add(rule76)

    def f77(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n), NonzeroQ(p + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern77 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f77))
    rule77 = ReplacementRule(pattern77, lambda b, a, m, c, n, p, x : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule77)

    def f78(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)), NonzeroQ(p + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern78 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f78))
    rule78 = ReplacementRule(pattern78, lambda b1, b2, a1, a2, m, c, n, p, x : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule78)

    def f79(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(p, x)])
    pattern79 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f79), )
    def With79(b, a, m, n, p, x):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p, x), x, x**k)/k
        print("Unable to Integrate")
    rule79 = ReplacementRule(pattern79, lambda b, a, m, n, p, x : With79(b, a, m, n, p, x))
    rubi.add(rule79)

    def f80(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), IntegerQ(m), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(p, x)])
    pattern80 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f80), )
    def With80(b1, b2, a1, a2, m, n, p, x):
        k = GCD(m + S(1), S(2)*n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a1 + b1*x**(n/k))**p*(a2 + b2*x**(n/k))**p, x), x, x**k)/k
        print("Unable to Integrate")
    rule80 = ReplacementRule(pattern80, lambda b1, b2, a1, a2, m, n, p, x : With80(b1, b2, a1, a2, m, n, p, x))
    rubi.add(rule80)

    def f81(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), Greater(p, S(0)), Less(m, S(-1)), Not(NegativeIntegerQ((m + n*p + n + S(1))/n)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern81 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f81))
    rule81 = ReplacementRule(pattern81, lambda b, a, m, c, n, p, x : -b*c**(-n)*n*p*Int((c*x)**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + S(1))))
    rubi.add(rule81)

    def f82(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(m, p), Greater(p, S(0)), NonzeroQ(m + S(2)*n*p + S(1)), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x)])
    pattern82 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f82))
    rule82 = ReplacementRule(pattern82, lambda b1, b2, a1, a2, m, c, n, p, x : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule82)

    def f83(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), Greater(p, S(0)), NonzeroQ(m + n*p + S(1)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x)])
    pattern83 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f83))
    rule83 = ReplacementRule(pattern83, lambda b, a, m, c, n, p, x : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule83)

    def f84(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern84 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f84))
    rule84 = ReplacementRule(pattern84, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule84)

    def f85(b, a, x, m):
        return functools.reduce(operator.and_, [ PosQ(b/a), PositiveIntegerQ(m/S(4) + S(-1)/2), FreeQ(a, x), FreeQ(b, x)])
    pattern85 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f85))
    rule85 = ReplacementRule(pattern85, lambda b, a, x, m : -a*(m + S(-3))*Int(x**(m + S(-4))/(a + b*x**S(4))**(S(5)/4), x)/(b*(m + S(-4))) + x**(m + S(-3))/(b*(a + b*x**S(4))**(S(1)/4)*(m + S(-4))))
    rubi.add(rule85)

    def f86(b, a, x, m):
        return functools.reduce(operator.and_, [ PosQ(b/a), NegativeIntegerQ(m/S(4) + S(-1)/2), FreeQ(a, x), FreeQ(b, x)])
    pattern86 = Pattern(Integral(x_**m_/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f86))
    rule86 = ReplacementRule(pattern86, lambda b, a, x, m : -b*m*Int(x**(m + S(4))/(a + b*x**S(4))**(S(5)/4), x)/(a*(m + S(1))) + x**(m + S(1))/(a*(a + b*x**S(4))**(S(1)/4)*(m + S(1))))
    rubi.add(rule86)

    def f87(b, a, x, c):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern87 = Pattern(Integral(sqrt(x_*WC('c', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f87))
    rule87 = ReplacementRule(pattern87, lambda b, a, x, c : sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(5)/4)), x)/(b*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule87)

    def f88(b, a, m, c, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), IntegerQ(S(2)*m), Greater(m, S(3)/2), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern88 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f88))
    rule88 = ReplacementRule(pattern88, lambda b, a, m, c, x : -S(2)*a*c**S(2)*(m + S(-1))*Int((c*x)**(m + S(-2))/(a + b*x**S(2))**(S(5)/4), x)/(b*(S(2)*m + S(-3))) + S(2)*c*(c*x)**(m + S(-1))/(b*(a + b*x**S(2))**(S(1)/4)*(S(2)*m + S(-3))))
    rubi.add(rule88)

    def f89(b, a, m, c, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), IntegerQ(S(2)*m), Less(m, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern89 = Pattern(Integral((x_*WC('c', S(1)))**m_/(a_ + x_**S(2)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f89))
    rule89 = ReplacementRule(pattern89, lambda b, a, m, c, x : -b*(S(2)*m + S(1))*Int((c*x)**(m + S(2))/(a + b*x**S(2))**(S(5)/4), x)/(S(2)*a*c**S(2)*(m + S(1))) + (c*x)**(m + S(1))/(a*c*(a + b*x**S(2))**(S(1)/4)*(m + S(1))))
    rubi.add(rule89)

    def f90(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern90 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(5)/4), x_),CustomConstraint(f90))
    rule90 = ReplacementRule(pattern90, lambda b, a, x : -Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/b - S(1)/(b*x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule90)

    def f91(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), Less(p, S(-1)), Greater(m + S(1), n), Not(NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern91 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f91))
    rule91 = ReplacementRule(pattern91, lambda b, a, m, c, n, p, x : -c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**(p + S(1)), x)/(b*n*(p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*n*(p + S(1))))
    rubi.add(rule91)

    def f92(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(m, p), Less(p, S(-1)), Greater(m + S(1), S(2)*n), Not(NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n))), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x)])
    pattern92 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f92))
    rule92 = ReplacementRule(pattern92, lambda b1, b2, a1, a2, m, c, n, p, x : -c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*b1*b2*n*(p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*b1*b2*n*(p + S(1))))
    rubi.add(rule92)

    def f93(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), Less(p, S(-1)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x)])
    pattern93 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f93))
    rule93 = ReplacementRule(pattern93, lambda b, a, m, c, n, p, x : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule93)

    def f94(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(m, p), Less(p, S(-1)), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x)])
    pattern94 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f94))
    rule94 = ReplacementRule(pattern94, lambda b1, b2, a1, a2, m, c, n, p, x : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule94)

    def f95(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern95 = Pattern(Integral(x_/(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f95))
    rule95 = ReplacementRule(pattern95, lambda b, a, x : Int((x*Rt(b, S(3)) + Rt(a, S(3)))/(x**S(2)*Rt(b, S(3))**S(2) - x*Rt(a, S(3))*Rt(b, S(3)) + Rt(a, S(3))**S(2)), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))) - Int(S(1)/(x*Rt(b, S(3)) + Rt(a, S(3))), x)/(S(3)*Rt(a, S(3))*Rt(b, S(3))))
    rubi.add(rule95)

    def f96(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n/S(2) + S(-1)/2), PositiveIntegerQ(m), Less(m, n + S(-1)), PosQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern96 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f96), )
    def With96(b, a, m, n, x):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) - s**(-m)*(-r)**(m + S(1))*Int(S(1)/(r + s*x), x)/(a*n)
    rule96 = ReplacementRule(pattern96, lambda b, a, m, n, x : With96(b, a, m, n, x))
    rubi.add(rule96)

    def f97(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(2) + S(-1)/2), PositiveIntegerQ(m), Less(m, n + S(-1)), NegQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern97 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f97), )
    def With97(b, a, m, n, x):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return -Dist(S(2)*s**(-m)*(-r)**(m + S(1))/(a*n), Sum(u, List(k, S(1), n/S(2) + S(-1)/2)), x) + r**(m + S(1))*s**(-m)*Int(S(1)/(r - s*x), x)/(a*n)
    rule97 = ReplacementRule(pattern97, lambda b, a, m, n, x : With97(b, a, m, n, x))
    rubi.add(rule97)

    def f98(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(4) + S(-1)/2), PositiveIntegerQ(m), Less(m, n + S(-1)), PosQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern98 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f98), )
    def With98(b, a, m, n, x):
        r = Numerator(Rt(a/b, n))
        s = Denominator(Rt(a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(Pi*m*(2*k - 1)/n) - s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 - 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x) + Int((r*cos(Pi*m*(2*k - 1)/n) + s*x*cos(Pi*(2*k - 1)*(m + 1)/n))/(r**2 + 2*r*s*x*cos(Pi*(2*k - 1)/n) + s**2*x**2), x)
        return S(2)*(S(-1))**(m/S(2))*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) + s**S(2)*x**S(2)), x)/(a*n) + Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x)
    rule98 = ReplacementRule(pattern98, lambda b, a, m, n, x : With98(b, a, m, n, x))
    rubi.add(rule98)

    def f99(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(4) + S(-1)/2), PositiveIntegerQ(m), Less(m, n + S(-1)), NegQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern99 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f99), )
    def With99(b, a, m, n, x):
        r = Numerator(Rt(-a/b, n))
        s = Denominator(Rt(-a/b, n))
        k = Symbol('k')
        u = Symbol('u')
        u = Int((r*cos(2*Pi*k*m/n) - s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 - 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x) + Int((r*cos(2*Pi*k*m/n) + s*x*cos(2*Pi*k*(m + 1)/n))/(r**2 + 2*r*s*x*cos(2*Pi*k/n) + s**2*x**2), x)
        return Dist(S(2)*r**(m + S(1))*s**(-m)/(a*n), Sum(u, List(k, S(1), n/S(4) + S(-1)/2)), x) + S(2)*r**(m + S(2))*s**(-m)*Int(S(1)/(r**S(2) - s**S(2)*x**S(2)), x)/(a*n)
    rule99 = ReplacementRule(pattern99, lambda b, a, m, n, x : With99(b, a, m, n, x))
    rubi.add(rule99)

    def f100(b, a, x):
        return functools.reduce(operator.and_, [ Or(PositiveQ(a/b), And(PosQ(a/b), AtomQ(SplitProduct(SumBaseQ, a)), AtomQ(SplitProduct(SumBaseQ, b)))), FreeQ(a, x), FreeQ(b, x)])
    pattern100 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f100), )
    def With100(b, a, x):
        r = Numerator(Rt(a/b, S(2)))
        s = Denominator(Rt(a/b, S(2)))
        return -Int((r - s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s) + Int((r + s*x**S(2))/(a + b*x**S(4)), x)/(S(2)*s)
    rule100 = ReplacementRule(pattern100, lambda b, a, x : With100(b, a, x))
    rubi.add(rule100)

    def f101(b, a, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(a/b)), FreeQ(a, x), FreeQ(b, x)])
    pattern101 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f101), )
    def With101(b, a, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(r - s*x**S(2)), x)/(S(2)*b) + s*Int(S(1)/(r + s*x**S(2)), x)/(S(2)*b)
    rule101 = ReplacementRule(pattern101, lambda b, a, x : With101(b, a, x))
    rubi.add(rule101)

    def f102(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(4)), PositiveIntegerQ(m), Less(m, n + S(-1)), PositiveQ(a/b), FreeQ(a, x), FreeQ(b, x)])
    pattern102 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f102), )
    def With102(b, a, m, n, x):
        r = Numerator(Rt(a/b, S(4)))
        s = Denominator(Rt(a/b, S(4)))
        return sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) - sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r) - sqrt(S(2))*s**S(3)*Int(x**(m - n/S(4))/(r**S(2) + sqrt(S(2))*r*s*x**(n/S(4)) + s**S(2)*x**(n/S(2))), x)/(S(4)*b*r)
    rule102 = ReplacementRule(pattern102, lambda b, a, m, n, x : With102(b, a, m, n, x))
    rubi.add(rule102)

    def f103(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(4)), PositiveIntegerQ(m), Less(m, n/S(2)), Not(PositiveQ(a/b)), FreeQ(a, x), FreeQ(b, x)])
    pattern103 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f103), )
    def With103(b, a, m, n, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return r*Int(x**m/(r - s*x**(n/S(2))), x)/(S(2)*a) + r*Int(x**m/(r + s*x**(n/S(2))), x)/(S(2)*a)
    rule103 = ReplacementRule(pattern103, lambda b, a, m, n, x : With103(b, a, m, n, x))
    rubi.add(rule103)

    def f104(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n/S(4)), PositiveIntegerQ(m), Inequality(n/S(2), LessEqual, m, Less, n), Not(PositiveQ(a/b)), FreeQ(a, x), FreeQ(b, x)])
    pattern104 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f104), )
    def With104(b, a, m, n, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(x**(m - n/S(2))/(r - s*x**(n/S(2))), x)/(S(2)*b) + s*Int(x**(m - n/S(2))/(r + s*x**(n/S(2))), x)/(S(2)*b)
    rule104 = ReplacementRule(pattern104, lambda b, a, m, n, x : With104(b, a, m, n, x))
    rubi.add(rule104)

    def f105(b, a, m, n, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(m, n), Greater(m, S(2)*n + S(-1)), FreeQ(a, x), FreeQ(b, x)])
    pattern105 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f105))
    rule105 = ReplacementRule(pattern105, lambda b, a, m, n, x : Int(PolynomialDivide(x**m, a + b*x**n, x), x))
    rubi.add(rule105)

    def f106(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern106 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f106), )
    def With106(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(sqrt(S(3)) + S(2))) + Int((r*x + s*(-sqrt(S(3)) + S(1)))/sqrt(a + b*x**S(3)), x)/r
    rule106 = ReplacementRule(pattern106, lambda b, a, x : With106(b, a, x))
    rubi.add(rule106)

    def f107(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(a), FreeQ(a, x), FreeQ(b, x)])
    pattern107 = Pattern(Integral(x_/sqrt(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f107), )
    def With107(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return -sqrt(S(2))*s*Int(S(1)/sqrt(a + b*x**S(3)), x)/(r*sqrt(-sqrt(S(3)) + S(2))) + Int((r*x + s*(S(1) + sqrt(S(3))))/sqrt(a + b*x**S(3)), x)/r
    rule107 = ReplacementRule(pattern107, lambda b, a, x : With107(b, a, x))
    rubi.add(rule107)

    def f108(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern108 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f108), )
    def With108(b, a, x):
        q = Rt(b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule108 = ReplacementRule(pattern108, lambda b, a, x : With108(b, a, x))
    rubi.add(rule108)

    def f109(b, a, x):
        return functools.reduce(operator.and_, [ NegativeQ(a), PositiveQ(b), FreeQ(a, x), FreeQ(b, x)])
    pattern109 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f109), )
    def With109(b, a, x):
        q = Rt(-b/a, S(2))
        return -Int((-q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q + Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule109 = ReplacementRule(pattern109, lambda b, a, x : With109(b, a, x))
    rubi.add(rule109)

    def f110(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern110 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f110), )
    def With110(b, a, x):
        q = Rt(-b/a, S(2))
        return Int((q*x**S(2) + S(1))/sqrt(a + b*x**S(4)), x)/q - Int(S(1)/sqrt(a + b*x**S(4)), x)/q
    rule110 = ReplacementRule(pattern110, lambda b, a, x : With110(b, a, x))
    rubi.add(rule110)

    def f111(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern111 = Pattern(Integral(x_**S(4)/sqrt(a_ + x_**S(6)*WC('b', S(1))), x_),CustomConstraint(f111), )
    def With111(b, a, x):
        r = Numer(Rt(b/a, S(3)))
        s = Denom(Rt(b/a, S(3)))
        return s**S(2)*(S(-1) + sqrt(S(3)))*Int(S(1)/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2)) - Int((-S(2)*r**S(2)*x**S(4) + s**S(2)*(S(-1) + sqrt(S(3))))/sqrt(a + b*x**S(6)), x)/(S(2)*r**S(2))
    rule111 = ReplacementRule(pattern111, lambda b, a, x : With111(b, a, x))
    rubi.add(rule111)

    def f112(b, a, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b), x), FreeQ(a, x), FreeQ(b, x)])
    pattern112 = Pattern(Integral(x_**S(2)/sqrt(a_ + x_**S(8)*WC('b', S(1))), x_),CustomConstraint(f112))
    rule112 = ReplacementRule(pattern112, lambda b, a, x : -Int((-x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))) + Int((x**S(2)*Rt(b/a, S(4)) + S(1))/sqrt(a + b*x**S(8)), x)/(S(2)*Rt(b/a, S(4))))
    rubi.add(rule112)

    def f113(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern113 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_),CustomConstraint(f113))
    rule113 = ReplacementRule(pattern113, lambda b, a, x : -a*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x)/S(2) + x**S(3)/(S(2)*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule113)

    def f114(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern114 = Pattern(Integral(x_**S(2)/(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4), x_),CustomConstraint(f114))
    rule114 = ReplacementRule(pattern114, lambda b, a, x : a*Int(S(1)/(x**S(2)*(a + b*x**S(4))**(S(1)/4)), x)/(S(2)*b) + (a + b*x**S(4))**(S(3)/4)/(S(2)*b*x))
    rubi.add(rule114)

    def f115(b, a, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern115 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_),CustomConstraint(f115))
    rule115 = ReplacementRule(pattern115, lambda b, a, x : -b*Int(x**S(2)/(a + b*x**S(4))**(S(5)/4), x) - S(1)/(x*(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule115)

    def f116(b, a, x):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x)])
    pattern116 = Pattern(Integral(S(1)/(x_**S(2)*(a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)), x_),CustomConstraint(f116))
    rule116 = ReplacementRule(pattern116, lambda b, a, x : x*(a/(b*x**S(4)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(3)*(a/(b*x**S(4)) + S(1))**(S(1)/4)), x)/(a + b*x**S(4))**(S(1)/4))
    rubi.add(rule116)

    def f117(b, a, x, c):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern117 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_),CustomConstraint(f117))
    rule117 = ReplacementRule(pattern117, lambda b, a, x, c : -a*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/S(2) + x*sqrt(c*x)/(a + b*x**S(2))**(S(1)/4))
    rubi.add(rule117)

    def f118(b, a, x, c):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern118 = Pattern(Integral(sqrt(c_*x_)/(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4), x_),CustomConstraint(f118))
    rule118 = ReplacementRule(pattern118, lambda b, a, x, c : a*c**S(2)*Int(S(1)/((c*x)**(S(3)/2)*(a + b*x**S(2))**(S(1)/4)), x)/(S(2)*b) + c*(a + b*x**S(2))**(S(3)/4)/(b*sqrt(c*x)))
    rubi.add(rule118)

    def f119(b, a, x, c):
        return functools.reduce(operator.and_, [ PosQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern119 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_),CustomConstraint(f119))
    rule119 = ReplacementRule(pattern119, lambda b, a, x, c : -b*Int(sqrt(c*x)/(a + b*x**S(2))**(S(5)/4), x)/c**S(2) - S(2)/(c*sqrt(c*x)*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule119)

    def f120(b, a, x, c):
        return functools.reduce(operator.and_, [ NegQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x)])
    pattern120 = Pattern(Integral(S(1)/((x_*WC('c', S(1)))**(S(3)/2)*(a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)), x_),CustomConstraint(f120))
    rule120 = ReplacementRule(pattern120, lambda b, a, x, c : sqrt(c*x)*(a/(b*x**S(2)) + S(1))**(S(1)/4)*Int(S(1)/(x**S(2)*(a/(b*x**S(2)) + S(1))**(S(1)/4)), x)/(c**S(2)*(a + b*x**S(2))**(S(1)/4)))
    rubi.add(rule120)

    def f121(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m), Greater(m, n + S(-1)), NonzeroQ(m + n*p + S(1)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern121 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f121))
    rule121 = ReplacementRule(pattern121, lambda b, a, m, c, n, p, x : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule121)

    def f122(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), SumSimplerQ(m, -n), NonzeroQ(m + n*p + S(1)), NegativeIntegerQ((m + n*p + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern122 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f122))
    rule122 = ReplacementRule(pattern122, lambda b, a, m, c, n, p, x : -a*c**n*(m - n + S(1))*Int((c*x)**(m - n)*(a + b*x**n)**p, x)/(b*(m + n*p + S(1))) + c**(n + S(-1))*(c*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*p + S(1))))
    rubi.add(rule122)

    def f123(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(m), Greater(m, S(2)*n + S(-1)), NonzeroQ(m + S(2)*n*p + S(1)), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(p, x)])
    pattern123 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f123))
    rule123 = ReplacementRule(pattern123, lambda b1, b2, a1, a2, m, c, n, p, x : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule123)

    def f124(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), SumSimplerQ(m, -S(2)*n), NonzeroQ(m + S(2)*n*p + S(1)), NegativeIntegerQ((m + S(2)*n*p + S(1))/(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern124 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f124))
    rule124 = ReplacementRule(pattern124, lambda b1, b2, a1, a2, m, c, n, p, x : -a1*a2*c**(S(2)*n)*(m - S(2)*n + S(1))*Int((c*x)**(m - S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(b1*b2*(m + S(2)*n*p + S(1))) + c**(S(2)*n + S(-1))*(c*x)**(m - S(2)*n + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(b1*b2*(m + S(2)*n*p + S(1))))
    rubi.add(rule124)

    def f125(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m), Less(m, S(-1)), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern125 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f125))
    rule125 = ReplacementRule(pattern125, lambda b, a, m, c, n, p, x : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule125)

    def f126(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), SumSimplerQ(m, n), NegativeIntegerQ((m + n*p + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern126 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f126))
    rule126 = ReplacementRule(pattern126, lambda b, a, m, c, n, p, x : -b*c**(-n)*(m + n*(p + S(1)) + S(1))*Int((c*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))) + (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*(m + S(1))))
    rubi.add(rule126)

    def f127(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(m), Less(m, S(-1)), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(p, x)])
    pattern127 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f127))
    rule127 = ReplacementRule(pattern127, lambda b1, b2, a1, a2, m, c, n, p, x : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule127)

    def f128(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), SumSimplerQ(m, S(2)*n), NegativeIntegerQ((m + S(2)*n*p + S(1))/(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern128 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f128))
    rule128 = ReplacementRule(pattern128, lambda b1, b2, a1, a2, m, c, n, p, x : -b1*b2*c**(-S(2)*n)*(m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**(m + S(2)*n)*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x)/(a1*a2*(m + S(1))) + (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(a1*a2*c*(m + S(1))))
    rubi.add(rule128)

    def f129(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FractionQ(m), IntBinomialQ(a, b, c, n, m, p, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern129 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f129), )
    def With129(b, a, m, c, n, p, x):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    rule129 = ReplacementRule(pattern129, lambda b, a, m, c, n, p, x : With129(b, a, m, c, n, p, x))
    rubi.add(rule129)

    def f130(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), FractionQ(m), IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(p, x)])
    pattern130 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f130), )
    def With130(b1, b2, a1, a2, m, c, n, p, x):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(k*n))**p*(a2 + b2*c**(-n)*x**(k*n))**p, x), x, (c*x)**(S(1)/k))/c
    rule130 = ReplacementRule(pattern130, lambda b1, b2, a1, a2, m, c, n, p, x : With130(b1, b2, a1, a2, m, c, n, p, x))
    rubi.add(rule130)

    def f131(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), IntegersQ(m, p + (m + S(1))/n), FreeQ(a, x), FreeQ(b, x)])
    pattern131 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f131))
    rule131 = ReplacementRule(pattern131, lambda b, a, m, n, p, x : a**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule131)

    def f132(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), IntegersQ(m, p + (m + S(1))/(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x)])
    pattern132 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f132))
    rule132 = ReplacementRule(pattern132, lambda b1, b2, a1, a2, m, n, p, x : (a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule132)

    def f133(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), IntegerQ(m), Less(Denominator(p + (m + S(1))/n), Denominator(p)), FreeQ(a, x), FreeQ(b, x)])
    pattern133 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f133))
    rule133 = ReplacementRule(pattern133, lambda b, a, m, n, p, x : (a/(a + b*x**n))**(p + (m + S(1))/n)*(a + b*x**n)**(p + (m + S(1))/n)*Subst(Int(x**m*(-b*x**n + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule133)

    def f134(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), PositiveIntegerQ(S(2)*n), RationalQ(p), Less(S(-1), p, S(0)), Unequal(p, S(-1)/2), IntegerQ(m), Less(Denominator(p + (m + S(1))/(S(2)*n)), Denominator(p)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x)])
    pattern134 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f134))
    rule134 = ReplacementRule(pattern134, lambda b1, b2, a1, a2, m, n, p, x : (a1/(a1 + b1*x**n))**(p + (m + S(1))/(S(2)*n))*(a2/(a2 + b2*x**n))**(p + (m + S(1))/(S(2)*n))*(a1 + b1*x**n)**(p + (m + S(1))/(S(2)*n))*(a2 + b2*x**n)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**m*(-b1*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**n + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x*(a1 + b1*x**n)**(-S(1)/(S(2)*n))*(a2 + b2*x**n)**(-S(1)/(S(2)*n))))
    rubi.add(rule134)

    def f135(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(p, x)])
    pattern135 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f135))
    rule135 = ReplacementRule(pattern135, lambda b, a, m, n, p, x : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule135)

    def f136(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ(S(2)*n), IntegerQ(m), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(p, x)])
    pattern136 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f136))
    rule136 = ReplacementRule(pattern136, lambda b1, b2, a1, a2, m, n, p, x : -Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule136)

    def f137(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), FractionQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(p, x)])
    pattern137 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f137), )
    def With137(b, a, m, c, n, p, x):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    rule137 = ReplacementRule(pattern137, lambda b, a, m, c, n, p, x : With137(b, a, m, c, n, p, x))
    rubi.add(rule137)

    def f138(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ(S(2)*n), FractionQ(m), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(p, x)])
    pattern138 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f138), )
    def With138(b1, b2, a1, a2, m, c, n, p, x):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a1 + b1*c**(-n)*x**(-k*n))**p*(a2 + b2*c**(-n)*x**(-k*n))**p, x), x, (c*x)**(-S(1)/k))/c
    rule138 = ReplacementRule(pattern138, lambda b1, b2, a1, a2, m, c, n, p, x : With138(b1, b2, a1, a2, m, c, n, p, x))
    rubi.add(rule138)

    def f139(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), Not(RationalQ(m)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern139 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f139))
    rule139 = ReplacementRule(pattern139, lambda b, a, m, c, n, p, x : -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule139)

    def f140(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), NegativeIntegerQ(S(2)*n), Not(RationalQ(m)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern140 = Pattern(Integral((x_*WC('c', S(1)))**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f140))
    rule140 = ReplacementRule(pattern140, lambda b1, b2, a1, a2, m, c, n, p, x : -(c*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a1 + b1*x**(-n))**p*(a2 + b2*x**(-n))**p, x), x, S(1)/x))
    rubi.add(rule140)

    def f141(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(p, x)])
    pattern141 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f141), )
    def With141(b, a, m, n, p, x):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p, x), x, x**(S(1)/k))
    rule141 = ReplacementRule(pattern141, lambda b, a, m, n, p, x : With141(b, a, m, n, p, x))
    rubi.add(rule141)

    def f142(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), FractionQ(S(2)*n), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(p, x)])
    pattern142 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f142), )
    def With142(b1, b2, a1, a2, m, n, p, x):
        k = Denominator(S(2)*n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a1 + b1*x**(k*n))**p*(a2 + b2*x**(k*n))**p, x), x, x**(S(1)/k))
    rule142 = ReplacementRule(pattern142, lambda b1, b2, a1, a2, m, n, p, x : With142(b1, b2, a1, a2, m, n, p, x))
    rubi.add(rule142)

    def f143(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern143 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f143))
    rule143 = ReplacementRule(pattern143, lambda b, a, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule143)

    def f144(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), FractionQ(S(2)*n), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(p, x)])
    pattern144 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f144))
    rule144 = ReplacementRule(pattern144, lambda b1, b2, a1, a2, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule144)

    def f145(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(n/(m + S(1))), Not(IntegerQ(n)), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern145 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f145))
    rule145 = ReplacementRule(pattern145, lambda b, a, m, n, p, x : Subst(Int((a + b*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule145)

    def f146(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(S(2)*n/(m + S(1))), Not(IntegerQ(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern146 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f146))
    rule146 = ReplacementRule(pattern146, lambda b1, b2, a1, a2, m, n, p, x : Subst(Int((a1 + b1*x**(n/(m + S(1))))**p*(a2 + b2*x**(n/(m + S(1))))**p, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule146)

    def f147(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(n/(m + S(1))), Not(IntegerQ(n)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern147 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f147))
    rule147 = ReplacementRule(pattern147, lambda b, a, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule147)

    def f148(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(S(2)*n/(m + S(1))), Not(IntegerQ(S(2)*n)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern148 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f148))
    rule148 = ReplacementRule(pattern148, lambda b1, b2, a1, a2, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule148)

    def f149(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(p + (m + S(1))/n), RationalQ(p), Greater(p, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x)])
    pattern149 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f149))
    rule149 = ReplacementRule(pattern149, lambda b, a, m, n, p, x : -b*n*p*Int(x**(m + n)*(a + b*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a + b*x**n)**p/(m + S(1)))
    rubi.add(rule149)

    def f150(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), ZeroQ(p + (m + S(1))/(S(2)*n)), RationalQ(p), Greater(p, S(0)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(n, x)])
    pattern150 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f150))
    rule150 = ReplacementRule(pattern150, lambda b1, b2, a1, a2, m, n, p, x : -S(2)*b1*b2*n*p*Int(x**(m + n)*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(1)) + x**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(m + S(1)))
    rubi.add(rule150)

    def f151(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(p + (m + S(1))/n), RationalQ(p), Greater(p, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern151 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f151))
    rule151 = ReplacementRule(pattern151, lambda b, a, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule151)

    def f152(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), ZeroQ(p + (m + S(1))/(S(2)*n)), RationalQ(p), Greater(p, S(0)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern152 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f152))
    rule152 = ReplacementRule(pattern152, lambda b1, b2, a1, a2, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule152)

    def f153(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(p + (m + S(1))/n), RationalQ(p), Greater(p, S(0)), NonzeroQ(m + n*p + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern153 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f153))
    rule153 = ReplacementRule(pattern153, lambda b, a, m, c, n, p, x : a*n*p*Int((c*x)**m*(a + b*x**n)**(p + S(-1)), x)/(m + n*p + S(1)) + (c*x)**(m + S(1))*(a + b*x**n)**p/(c*(m + n*p + S(1))))
    rubi.add(rule153)

    def f154(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(p + (m + S(1))/(S(2)*n)), RationalQ(p), Greater(p, S(0)), NonzeroQ(m + S(2)*n*p + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern154 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f154))
    rule154 = ReplacementRule(pattern154, lambda b1, b2, a1, a2, m, c, n, p, x : S(2)*a1*a2*n*p*Int((c*x)**m*(a1 + b1*x**n)**(p + S(-1))*(a2 + b2*x**n)**(p + S(-1)), x)/(m + S(2)*n*p + S(1)) + (c*x)**(m + S(1))*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p/(c*(m + S(2)*n*p + S(1))))
    rubi.add(rule154)

    def f155(b, a, m, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(p + (m + S(1))/n), RationalQ(p), Less(S(-1), p, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x)])
    pattern155 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f155), )
    def With155(b, a, m, n, p, x):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(-b*x**k + S(1))**(-p + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    rule155 = ReplacementRule(pattern155, lambda b, a, m, n, p, x : With155(b, a, m, n, p, x))
    rubi.add(rule155)

    def f156(b1, b2, a1, a2, m, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(p + (m + S(1))/(S(2)*n)), RationalQ(p), Less(S(-1), p, S(0)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(m, x), FreeQ(n, x)])
    pattern156 = Pattern(Integral(x_**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f156), )
    def With156(b1, b2, a1, a2, m, n, p, x):
        k = Denominator(p)
        return k*(a1*a2)**(p + (m + S(1))/(S(2)*n))*Subst(Int(x**(k*(m + S(1))/(S(2)*n) + S(-1))*(-b1*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n))*(-b2*x**k + S(1))**(-p + S(-1) - (m + S(1))/(S(2)*n)), x), x, x**(S(2)*n/k)*(a1 + b1*x**n)**(-S(1)/k)*(a2 + b2*x**n)**(-S(1)/k))/(S(2)*n)
    rule156 = ReplacementRule(pattern156, lambda b1, b2, a1, a2, m, n, p, x : With156(b1, b2, a1, a2, m, n, p, x))
    rubi.add(rule156)

    def f157(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(p + (m + S(1))/n), RationalQ(p), Less(S(-1), p, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern157 = Pattern(Integral((c_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f157))
    rule157 = ReplacementRule(pattern157, lambda b, a, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p, x))
    rubi.add(rule157)

    def f158(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(p + (m + S(1))/(S(2)*n)), RationalQ(p), Less(S(-1), p, S(0)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern158 = Pattern(Integral((c_*x_)**m_*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f158))
    rule158 = ReplacementRule(pattern158, lambda b1, b2, a1, a2, m, c, n, p, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m*(a1 + b1*x**n)**p*(a2 + b2*x**n)**p, x))
    rubi.add(rule158)

    def f159(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ IntegerQ(p + (m + S(1))/n), RationalQ(p), Less(p, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern159 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f159))
    rule159 = ReplacementRule(pattern159, lambda b, a, m, c, n, p, x : (m + n*(p + S(1)) + S(1))*Int((c*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*n*(p + S(1))) - (c*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*c*n*(p + S(1))))
    rubi.add(rule159)

    def f160(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), IntegerQ(p + (m + S(1))/n), RationalQ(p), Less(p, S(-1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern160 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f160))
    rule160 = ReplacementRule(pattern160, lambda b1, b2, a1, a2, m, c, n, p, x : (m + S(2)*n*(p + S(1)) + S(1))*Int((c*x)**m*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1)), x)/(S(2)*a1*a2*n*(p + S(1))) - (c*x)**(m + S(1))*(a1 + b1*x**n)**(p + S(1))*(a2 + b2*x**n)**(p + S(1))/(S(2)*a1*a2*c*n*(p + S(1))))
    rubi.add(rule160)

    def f161(b, a, m, n, x):
        return functools.reduce(operator.and_, [ FractionQ((m + S(1))/n), SumSimplerQ(m, -n), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x)])
    pattern161 = Pattern(Integral(x_**WC('m', S(1))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f161), )
    def With161(b, a, m, n, x):
        mn = m - n
        return -a*Int(x**mn/(a + b*x**n), x)/b + x**(mn + S(1))/(b*(mn + S(1)))
    rule161 = ReplacementRule(pattern161, lambda b, a, m, n, x : With161(b, a, m, n, x))
    rubi.add(rule161)

    def f162(b, a, m, n, x):
        return functools.reduce(operator.and_, [ FractionQ((m + S(1))/n), SumSimplerQ(m, n), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x)])
    pattern162 = Pattern(Integral(x_**m_/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f162))
    rule162 = ReplacementRule(pattern162, lambda b, a, m, n, x : -b*Int(x**(m + n)/(a + b*x**n), x)/a + x**(m + S(1))/(a*(m + S(1))))
    rubi.add(rule162)

    def f163(b, a, m, c, n, x):
        return functools.reduce(operator.and_, [ FractionQ((m + S(1))/n), Or(SumSimplerQ(m, n), SumSimplerQ(m, -n)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x)])
    pattern163 = Pattern(Integral((c_*x_)**m_/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f163))
    rule163 = ReplacementRule(pattern163, lambda b, a, m, c, n, x : c**IntPart(m)*x**(-FracPart(m))*(c*x)**FracPart(m)*Int(x**m/(a + b*x**n), x))
    rubi.add(rule163)

    def f164(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ Not(PositiveIntegerQ(p)), Or(NegativeIntegerQ(p), PositiveQ(a)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern164 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f164))
    rule164 = ReplacementRule(pattern164, lambda b, a, m, c, n, p, x : a**p*(c*x)**(m + S(1))*Hypergeometric2F1(-p, (m + S(1))/n, S(1) + (m + S(1))/n, -b*x**n/a)/(c*(m + S(1))))
    rubi.add(rule164)

    def f165(b, a, m, c, n, p, x):
        return functools.reduce(operator.and_, [ Not(PositiveIntegerQ(p)), Not(Or(NegativeIntegerQ(p), PositiveQ(a))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern165 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_, x_),CustomConstraint(f165))
    rule165 = ReplacementRule(pattern165, lambda b, a, m, c, n, p, x : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((c*x)**m*(S(1) + b*x**n/a)**p, x))
    rubi.add(rule165)

    def f166(b, a, m, n, p, x, v):
        return functools.reduce(operator.and_, [ LinearQ(v, x), IntegerQ(m), NonzeroQ(v - x), FreeQ(a, x), FreeQ(b, x), FreeQ(n, x), FreeQ(p, x)])
    pattern166 = Pattern(Integral(x_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f166))
    rule166 = ReplacementRule(pattern166, lambda b, a, m, n, p, x, v : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule166)

    def f167(b, a, m, n, u, p, x, v):
        return functools.reduce(operator.and_, [ LinearPairQ(u, v, x), FreeQ(a, x), FreeQ(b, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern167 = Pattern(Integral(u_**WC('m', S(1))*(a_ + v_**n_*WC('b', S(1)))**WC('p', S(1)), x_),CustomConstraint(f167))
    rule167 = ReplacementRule(pattern167, lambda b, a, m, n, u, p, x, v : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule167)

    def f168(b1, b2, a1, a2, m, c, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(a1*b2 + a2*b1), Not(IntegerQ(p)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern168 = Pattern(Integral((x_*WC('c', S(1)))**WC('m', S(1))*(a1_ + x_**n_*WC('b1', S(1)))**p_*(a2_ + x_**n_*WC('b2', S(1)))**p_, x_),CustomConstraint(f168))
    rule168 = ReplacementRule(pattern168, lambda b1, b2, a1, a2, m, c, n, p, x : (a1 + b1*x**n)**FracPart(p)*(a2 + b2*x**n)**FracPart(p)*(a1*a2 + b1*b2*x**(S(2)*n))**(-FracPart(p))*Int((c*x)**m*(a1*a2 + b1*b2*x**(S(2)*n))**p, x))
    rubi.add(rule168)

    def f169(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(p, q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern169 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f169))
    rule169 = ReplacementRule(pattern169, lambda b, a, c, d, n, p, q, x : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule169)

    def f170(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegersQ(p, q), NegQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern170 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f170))
    rule170 = ReplacementRule(pattern170, lambda b, a, c, d, n, p, q, x : Int(x**(n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule170)

    def f171(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegativeIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(p, x), FreeQ(q, x)])
    pattern171 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f171))
    rule171 = ReplacementRule(pattern171, lambda b, a, c, d, n, p, q, x : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q/x**S(2), x), x, S(1)/x))
    rubi.add(rule171)

    def f172(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(p, x), FreeQ(q, x)])
    pattern172 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f172), )
    def With172(b, a, c, d, n, p, q, x):
        g = Denominator(n)
        return g*Subst(Int(x**(g + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    rule172 = ReplacementRule(pattern172, lambda b, a, c, d, n, p, q, x : With172(b, a, c, d, n, p, q, x))
    rubi.add(rule172)

    def f173(b, a, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*p + S(1)), IntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern173 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f173))
    rule173 = ReplacementRule(pattern173, lambda b, a, c, d, n, p, x : Subst(Int(S(1)/(c - x**n*(-a*d + b*c)), x), x, x*(a + b*x**n)**(-S(1)/n)))
    rubi.add(rule173)

    def f174(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*(p + q + S(1)) + S(1)), RationalQ(q), Greater(q, S(0)), NonzeroQ(p + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x)])
    pattern174 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f174))
    rule174 = ReplacementRule(pattern174, lambda b, a, c, d, n, p, q, x : -c*q*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1)), x)/(a*(p + S(1))) - x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))))
    rubi.add(rule174)

    def f175(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*(p + q + S(1)) + S(1)), NegativeIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(q, x)])
    pattern175 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f175))
    rule175 = ReplacementRule(pattern175, lambda b, a, c, d, n, p, q, x : a**p*c**(-p + S(-1))*x*(c + d*x**n)**(-S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n))))
    rubi.add(rule175)

    def f176(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*(p + q + S(1)) + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern176 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f176))
    rule176 = ReplacementRule(pattern176, lambda b, a, c, d, n, p, q, x : x*(c*(a + b*x**n)/(a*(c + d*x**n)))**(-p)*(a + b*x**n)**p*(c + d*x**n)**(-p - S(1)/n)*Hypergeometric2F1(S(1)/n, -p, S(1) + S(1)/n, -x**n*(-a*d + b*c)/(a*(c + d*x**n)))/c)
    rubi.add(rule176)

    def f177(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*(p + q + S(2)) + S(1)), ZeroQ(a*d*(p + S(1)) + b*c*(q + S(1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern177 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f177))
    rule177 = ReplacementRule(pattern177, lambda b, a, c, d, n, p, q, x : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c))
    rubi.add(rule177)

    def f178(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(n*(p + q + S(2)) + S(1)), Or(And(RationalQ(p), Less(p, S(-1))), Not(And(RationalQ(q), Less(q, S(-1))))), NonzeroQ(p + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(q, x)])
    pattern178 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f178))
    rule178 = ReplacementRule(pattern178, lambda b, a, c, d, n, p, q, x : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + (b*c + n*(p + S(1))*(-a*d + b*c))*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q, x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule178)

    def f179(b, a, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(a*d - b*c*(n*(p + S(1)) + S(1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x)])
    pattern179 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f179))
    rule179 = ReplacementRule(pattern179, lambda b, a, c, d, n, p, x : c*x*(a + b*x**n)**(p + S(1))/a)
    rubi.add(rule179)

    def f180(b, a, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), Or(And(RationalQ(p), Less(p, S(-1))), NegativeIntegerQ(p + S(1)/n)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x)])
    pattern180 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f180))
    rule180 = ReplacementRule(pattern180, lambda b, a, c, d, n, p, x : -x*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*n*(p + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))))
    rubi.add(rule180)

    def f181(b, a, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(n), Less(n, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern181 = Pattern(Integral((c_ + x_**n_*WC('d', S(1)))/(a_ + x_**n_*WC('b', S(1))), x_),CustomConstraint(f181))
    rule181 = ReplacementRule(pattern181, lambda b, a, c, d, n, x : c*x/a - (-a*d + b*c)*Int(S(1)/(a*x**(-n) + b), x)/a)
    rubi.add(rule181)

    def f182(b, a, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(n*(p + S(1)) + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern182 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f182))
    rule182 = ReplacementRule(pattern182, lambda b, a, c, d, n, p, x : d*x*(a + b*x**n)**(p + S(1))/(b*(n*(p + S(1)) + S(1))) - (a*d - b*c*(n*(p + S(1)) + S(1)))*Int((a + b*x**n)**p, x)/(b*(n*(p + S(1)) + S(1))))
    rubi.add(rule182)

    def f183(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n, p), NegativeIntegerQ(q), GreaterEqual(p, -q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern183 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f183))
    rule183 = ReplacementRule(pattern183, lambda b, a, c, d, n, p, q, x : Int(PolynomialDivide((a + b*x**n)**p, (c + d*x**n)**(-q), x), x))
    rubi.add(rule183)

    def f184(b, a, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern184 = Pattern(Integral(S(1)/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f184))
    rule184 = ReplacementRule(pattern184, lambda b, a, c, d, n, x : b*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - d*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule184)

    def f185(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(S(3)*a*d + b*c), PosQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern185 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f185))
    rule185 = ReplacementRule(pattern185, lambda b, a, c, d, x : sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(-x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c) + sqrt(S(3))*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(x*Rt(b/a, S(2)) + sqrt(S(3)))), x)/(S(2)*c))
    rubi.add(rule185)

    def f186(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(S(3)*a*d + b*c), NegQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern186 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/3)*(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f186))
    rule186 = ReplacementRule(pattern186, lambda b, a, c, d, x : Int((-x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6) + Int((x*Rt(-b/a, S(2)) + S(3))/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/S(6))
    rubi.add(rule186)

    def f187(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(S(3)*a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern187 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**(S(2)/3)/(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f187))
    rule187 = ReplacementRule(pattern187, lambda b, a, c, d, x : b*Int((a + b*x**S(2))**(S(-1)/3), x)/d - (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))**(S(1)/3)*(c + d*x**S(2))), x)/d)
    rubi.add(rule187)

    def f188(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern188 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(1)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f188))
    rule188 = ReplacementRule(pattern188, lambda b, a, c, d, x : sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(1)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule188)

    def f189(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern189 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f189))
    rule189 = ReplacementRule(pattern189, lambda b, a, c, d, x : sqrt(-b*x**S(2)/a)*Subst(Int(S(1)/(sqrt(-b*x/a)*(a + b*x)**(S(3)/4)*(c + d*x)), x), x, x**S(2))/(S(2)*x))
    rubi.add(rule189)

    def f190(b, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Greater(p, S(0)), Or(Equal(p, S(1)/2), Equal(Denominator(p), S(4))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**WC('p', S(1))/(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f190))
    rule190 = ReplacementRule(pattern190, lambda b, a, c, d, p, x : b*Int((a + b*x**S(2))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(2))**(p + S(-1))/(c + d*x**S(2)), x)/d)
    rubi.add(rule190)

    def f191(b, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), Equal(Denominator(p), S(4)), Or(Equal(p, S(-5)/4), Equal(p, S(-7)/4)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern191 = Pattern(Integral((a_ + x_**S(2)*WC('b', S(1)))**p_/(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f191))
    rule191 = ReplacementRule(pattern191, lambda b, a, c, d, p, x : b*Int((a + b*x**S(2))**p, x)/(-a*d + b*c) - d*Int((a + b*x**S(2))**(p + S(1))/(c + d*x**S(2)), x)/(-a*d + b*c))
    rubi.add(rule191)

    def f192(b, a, c, d, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), PosQ(a*b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern192 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_),CustomConstraint(f192))
    rule192 = ReplacementRule(pattern192, lambda b, a, c, d, x : a*Subst(Int(S(1)/(-S(4)*a*b*x**S(4) + S(1)), x), x, x/sqrt(a + b*x**S(4)))/c)
    rubi.add(rule192)

    def f193(b, a, c, d, x):
        return functools.reduce(operator.and_, [ ZeroQ(a*d + b*c), NegQ(a*b), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern193 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_),CustomConstraint(f193), )
    def With193(b, a, c, d, x):
        q = Rt(-a*b, S(4))
        return a*ArcTan(q*x*(a + q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q) + a*atanh(q*x*(a - q**S(2)*x**S(2))/(a*sqrt(a + b*x**S(4))))/(S(2)*c*q)
    rule193 = ReplacementRule(pattern193, lambda b, a, c, d, x : With193(b, a, c, d, x))
    rubi.add(rule193)

    def f194(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern194 = Pattern(Integral(sqrt(a_ + x_**S(4)*WC('b', S(1)))/(c_ + x_**S(4)*WC('d', S(1))), x_),CustomConstraint(f194))
    rule194 = ReplacementRule(pattern194, lambda b, a, c, d, x : b*Int(S(1)/sqrt(a + b*x**S(4)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(4))*(c + d*x**S(4))), x)/d)
    rubi.add(rule194)

    def f195(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern195 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**(S(1)/4)/(c_ + x_**S(4)*WC('d', S(1))), x_),CustomConstraint(f195))
    rule195 = ReplacementRule(pattern195, lambda b, a, c, d, x : sqrt(a/(a + b*x**S(4)))*sqrt(a + b*x**S(4))*Subst(Int(S(1)/((c - x**S(4)*(-a*d + b*c))*sqrt(-b*x**S(4) + S(1))), x), x, x/(a + b*x**S(4))**(S(1)/4)))
    rubi.add(rule195)

    def f196(b, a, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Or(Equal(p, S(3)/4), Equal(p, S(5)/4)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern196 = Pattern(Integral((a_ + x_**S(4)*WC('b', S(1)))**p_/(c_ + x_**S(4)*WC('d', S(1))), x_),CustomConstraint(f196))
    rule196 = ReplacementRule(pattern196, lambda b, a, c, d, p, x : b*Int((a + b*x**S(4))**(p + S(-1)), x)/d - (-a*d + b*c)*Int((a + b*x**S(4))**(p + S(-1))/(c + d*x**S(4)), x)/d)
    rubi.add(rule196)

    def f197(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern197 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(4)*WC('b', S(1)))*(c_ + x_**S(4)*WC('d', S(1)))), x_),CustomConstraint(f197))
    rule197 = ReplacementRule(pattern197, lambda b, a, c, d, x : Int(S(1)/(sqrt(a + b*x**S(4))*(-x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c) + Int(S(1)/(sqrt(a + b*x**S(4))*(x**S(2)*Rt(-d/c, S(2)) + S(1))), x)/(S(2)*c))
    rubi.add(rule197)

    def f198(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern198 = Pattern(Integral(S(1)/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_),CustomConstraint(f198))
    rule198 = ReplacementRule(pattern198, lambda b, a, c, d, x : b*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - d*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule198)

    def f199(b, a, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), PosQ(d/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern199 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2), x_),CustomConstraint(f199))
    rule199 = ReplacementRule(pattern199, lambda b, a, c, d, x : sqrt(a + b*x**S(2))*EllipticE(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(c*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule199)

    def f200(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Less(p, S(-1)), Less(S(0), q, S(1)), IntBinomialQ(a, b, c, d, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern200 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f200))
    rule200 = ReplacementRule(pattern200, lambda b, a, c, d, n, p, q, x : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(n*(p + S(1)) + S(1)) + d*x**n*(n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))))
    rubi.add(rule200)

    def f201(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Less(p, S(-1)), Greater(q, S(1)), IntBinomialQ(a, b, c, d, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern201 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f201))
    rule201 = ReplacementRule(pattern201, lambda b, a, c, d, n, p, q, x : x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(a*d - b*c)/(a*b*n*(p + S(1))) - Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(a*d - b*c*(n*(p + S(1)) + S(1))) + d*x**n*(a*d*(n*(q + S(-1)) + S(1)) - b*c*(n*(p + q) + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule201)

    def f202(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), Not(And(Not(IntegerQ(p)), IntegerQ(q), Less(q, S(-1)))), IntBinomialQ(a, b, c, d, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(q, x)])
    pattern202 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f202))
    rule202 = ReplacementRule(pattern202, lambda b, a, c, d, n, p, q, x : -b*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c + b*d*x**n*(n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule202)

    def f203(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), IntegersQ(p, q), Greater(p + q, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern203 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f203))
    rule203 = ReplacementRule(pattern203, lambda b, a, c, d, n, p, q, x : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule203)

    def f204(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(q), Greater(q, S(1)), NonzeroQ(n*(p + q) + S(1)), Not(And(IntegerQ(p), Greater(p, S(1)))), IntBinomialQ(a, b, c, d, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x)])
    pattern204 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f204))
    rule204 = ReplacementRule(pattern204, lambda b, a, c, d, n, p, q, x : d*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*(n*(p + q) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(-a*d + b*c*(n*(p + q) + S(1))) + d*x**n*(-a*d*(n*(q + S(-1)) + S(1)) + b*c*(n*(p + S(2)*q + S(-1)) + S(1))), x), x)/(b*(n*(p + q) + S(1))))
    rubi.add(rule204)

    def f205(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Greater(q, S(0)), Greater(p, S(0)), IntBinomialQ(a, b, c, d, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x)])
    pattern205 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f205))
    rule205 = ReplacementRule(pattern205, lambda b, a, c, d, n, p, q, x : n*Int((a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(n*(p + q) + S(1)) + x*(a + b*x**n)**p*(c + d*x**n)**q/(n*(p + q) + S(1)))
    rubi.add(rule205)

    def f206(b, a, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(b/a), Not(SimplerSqrtQ(b/a, d/c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern206 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f206))
    rule206 = ReplacementRule(pattern206, lambda b, a, c, d, x : sqrt(a + b*x**S(2))*EllipticF(ArcTan(x*Rt(d/c, S(2))), S(1) - b*c/(a*d))/(a*sqrt(c*(a + b*x**S(2))/(a*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule206)

    def f207(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), PositiveQ(a), Not(And(NegQ(b/a), SimplerSqrtQ(-b/a, -d/c))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern207 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f207))
    rule207 = ReplacementRule(pattern207, lambda b, a, c, d, x : EllipticF(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(a)*sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule207)

    def f208(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), PositiveQ(a - b*c/d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern208 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f208))
    rule208 = ReplacementRule(pattern208, lambda b, a, c, d, x : -EllipticF(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*sqrt(a - b*c/d)*Rt(-d/c, S(2))))
    rubi.add(rule208)

    def f209(b, a, c, d, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern209 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f209))
    rule209 = ReplacementRule(pattern209, lambda b, a, c, d, x : sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*sqrt(a + b*x**S(2))), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule209)

    def f210(b, a, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern210 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f210))
    rule210 = ReplacementRule(pattern210, lambda b, a, c, d, x : a*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x) + b*Int(x**S(2)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x))
    rubi.add(rule210)

    def f211(b, a, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), NegQ(b/a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern211 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f211))
    rule211 = ReplacementRule(pattern211, lambda b, a, c, d, x : b*Int(sqrt(c + d*x**S(2))/sqrt(a + b*x**S(2)), x)/d - (-a*d + b*c)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/d)
    rubi.add(rule211)

    def f212(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), PositiveQ(a), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern212 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f212))
    rule212 = ReplacementRule(pattern212, lambda b, a, c, d, x : sqrt(a)*EllipticE(asin(x*Rt(-d/c, S(2))), b*c/(a*d))/(sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule212)

    def f213(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), PositiveQ(a - b*c/d), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern213 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f213))
    rule213 = ReplacementRule(pattern213, lambda b, a, c, d, x : -sqrt(a - b*c/d)*EllipticE(acos(x*Rt(-d/c, S(2))), b*c/(-a*d + b*c))/(sqrt(c)*Rt(-d/c, S(2))))
    rubi.add(rule213)

    def f214(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern214 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f214))
    rule214 = ReplacementRule(pattern214, lambda b, a, c, d, x : sqrt(a + b*x**S(2))*Int(sqrt(S(1) + b*x**S(2)/a)/sqrt(c + d*x**S(2)), x)/sqrt(S(1) + b*x**S(2)/a))
    rubi.add(rule214)

    def f215(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), Not(PositiveQ(c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern215 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/sqrt(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f215))
    rule215 = ReplacementRule(pattern215, lambda b, a, c, d, x : sqrt(S(1) + d*x**S(2)/c)*Int(sqrt(a + b*x**S(2))/sqrt(S(1) + d*x**S(2)/c), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule215)

    def f216(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(q, x)])
    pattern216 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f216))
    rule216 = ReplacementRule(pattern216, lambda b, a, c, d, n, p, q, x : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule216)

    def f217(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(n + S(1)), PositiveQ(a), PositiveQ(c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern217 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f217))
    rule217 = ReplacementRule(pattern217, lambda b, a, c, d, n, p, q, x : a**p*c**q*x*AppellF1(S(1)/n, -p, -q, S(1) + S(1)/n, -b*x**n/a, -d*x**n/c))
    rubi.add(rule217)

    def f218(b, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(n + S(1)), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern218 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f218))
    rule218 = ReplacementRule(pattern218, lambda b, a, c, d, n, p, q, x : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule218)

    def f219(b, mn, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(q), Or(PosQ(n), Not(IntegerQ(p))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x)])
    pattern219 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f219))
    rule219 = ReplacementRule(pattern219, lambda b, mn, a, c, d, n, p, q, x : Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule219)

    def f220(b, mn, a, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), Not(IntegerQ(q)), Not(IntegerQ(p)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern220 = Pattern(Integral((a_ + x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_),CustomConstraint(f220))
    rule220 = ReplacementRule(pattern220, lambda b, mn, a, c, d, n, p, q, x : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule220)

    def f221(b, a, c, d, n, u, p, q, x):
        return functools.reduce(operator.and_, [ LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern221 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_),CustomConstraint(f221))
    rule221 = ReplacementRule(pattern221, lambda b, a, c, d, n, u, p, q, x : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule221)

    def f222(v, u, p, q, x):
        return functools.reduce(operator.and_, [ PseudoBinomialPairQ(u, v, x), FreeQ(p, x), FreeQ(q, x)])
    pattern222 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1)), x_),CustomConstraint(f222))
    rule222 = ReplacementRule(pattern222, lambda v, u, p, q, x : Int(NormalizePseudoBinomial(u, x)**p*NormalizePseudoBinomial(v, x)**q, x))
    rubi.add(rule222)

    def f223(x, m, u, p, q, v):
        return functools.reduce(operator.and_, [ IntegersQ(p, m/p), PseudoBinomialPairQ(u*x**(m/p), v, x), FreeQ(p, x), FreeQ(q, x)])
    pattern223 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*x_**WC('m', S(1)), x_),CustomConstraint(f223))
    rule223 = ReplacementRule(pattern223, lambda x, m, u, p, q, v : Int(NormalizePseudoBinomial(v, x)**q*NormalizePseudoBinomial(u*x**(m/p), x)**p, x))
    rubi.add(rule223)

    def f224(b, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(e)), IntegerQ((m + S(1))/n), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern224 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f224))
    rule224 = ReplacementRule(pattern224, lambda b, e, m, c, d, n, p, q, x : b**(S(1) - (m + S(1))/n)*e**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule224)

    def f225(b, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(e)), Not(IntegerQ((m + S(1))/n)), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern225 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f225))
    rule225 = ReplacementRule(pattern225, lambda b, e, m, c, d, n, p, q, x : b**IntPart(p)*e**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q, x))
    rubi.add(rule225)

    def f226(b, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ Not(IntegerQ(m)), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern226 = Pattern(Integral((e_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f226))
    rule226 = ReplacementRule(pattern226, lambda b, e, m, c, d, n, p, q, x : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule226)

    def f227(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(m - n + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern227 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f227))
    rule227 = ReplacementRule(pattern227, lambda b, a, m, c, d, n, p, q, x : Subst(Int((a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule227)

    def f228(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegersQ(p, q), NegQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern228 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f228))
    rule228 = ReplacementRule(pattern228, lambda b, a, m, c, d, n, p, q, x : Int(x**(m + n*(p + q))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q, x))
    rubi.add(rule228)

    def f229(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern229 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f229))
    rule229 = ReplacementRule(pattern229, lambda b, a, m, c, d, n, p, q, x : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q, x), x, x**n)/n)
    rubi.add(rule229)

    def f230(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern230 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f230))
    rule230 = ReplacementRule(pattern230, lambda b, e, a, m, c, d, n, p, q, x : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule230)

    def f231(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(p, q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern231 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f231))
    rule231 = ReplacementRule(pattern231, lambda b, e, a, m, c, d, n, p, q, x : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule231)

    def f232(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1))), NonzeroQ(m + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern232 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f232))
    rule232 = ReplacementRule(pattern232, lambda b, e, a, m, c, d, n, p, x : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule232)

    def f233(non2, e, b1, b2, a1, a2, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), ZeroQ(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1))), NonzeroQ(m + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern233 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f233))
    rule233 = ReplacementRule(pattern233, lambda non2, e, b1, b2, a1, a2, m, c, d, n, p, x : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))))
    rubi.add(rule233)

    def f234(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(m + n*(p + S(1)) + S(1)), Or(IntegerQ(n), PositiveQ(e)), RationalQ(m, n), Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1)))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern234 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f234))
    rule234 = ReplacementRule(pattern234, lambda b, e, a, m, c, d, n, p, x : d*e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p, x) + c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))))
    rubi.add(rule234)

    def f235(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(m + n*(p + S(1)) + S(1)), NonzeroQ(m + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern235 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f235))
    rule235 = ReplacementRule(pattern235, lambda b, e, a, m, c, d, n, p, x : d*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/b + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*(m + S(1))))
    rubi.add(rule235)

    def f236(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), Or(IntegerQ(n), PositiveQ(e)), RationalQ(m, n), Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1)))), Not(And(IntegerQ(p), Less(p, S(-1)))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern236 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f236))
    rule236 = ReplacementRule(pattern236, lambda b, e, a, m, c, d, n, p, x : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) + e**(-n)*(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a + b*x**n)**p, x)/(a*(m + S(1))))
    rubi.add(rule236)

    def f237(non2, e, b1, b2, a1, a2, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), Or(IntegerQ(n), PositiveQ(e)), RationalQ(m, n), Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1)))), Not(And(IntegerQ(p), Less(p, S(-1)))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern237 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f237))
    rule237 = ReplacementRule(pattern237, lambda non2, e, b1, b2, a1, a2, m, c, d, n, p, x : c*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(a1*a2*e*(m + S(1))) + e**(-n)*(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**(m + n)*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(a1*a2*(m + S(1))))
    rubi.add(rule237)

    def f238(b, a, m, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), PositiveIntegerQ(m/S(2)), Or(IntegerQ(p), Equal(m + S(2)*p + S(1), S(0))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern238 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f238))
    rule238 = ReplacementRule(pattern238, lambda b, a, m, c, d, p, x : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int((a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*x**S(2)*(p + S(1))*(b**(m/S(2))*x**(m + S(-2))*(c + d*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - (-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule238)

    def f239(b, a, m, c, d, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), NegativeIntegerQ(m/S(2)), Or(IntegerQ(p), Equal(m + S(2)*p + S(1), S(0))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern239 = Pattern(Integral(x_**m_*(a_ + x_**S(2)*WC('b', S(1)))**p_*(c_ + x_**S(2)*WC('d', S(1))), x_),CustomConstraint(f239))
    rule239 = ReplacementRule(pattern239, lambda b, a, m, c, d, p, x : b**(-m/S(2) + S(-1))*x*(-a)**(m/S(2) + S(-1))*(a + b*x**S(2))**(p + S(1))*(-a*d + b*c)/(S(2)*(p + S(1))) + b**(-m/S(2) + S(-1))*Int(x**m*(a + b*x**S(2))**(p + S(1))*ExpandToSum(S(2)*b*(p + S(1))*(b**(m/S(2))*(c + d*x**S(2)) - x**(-m + S(2))*(-a)**(m/S(2) + S(-1))*(-a*d + b*c))/(a + b*x**S(2)) - x**(-m)*(-a)**(m/S(2) + S(-1))*(-a*d + b*c), x), x)/(S(2)*(p + S(1))))
    rubi.add(rule239)

    def f240(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), Or(IntegerQ(p), Not(RationalQ(m)), And(PositiveIntegerQ(n), NegativeIntegerQ(p + S(1)/2), LessEqual(S(-1), m, -n*(p + S(1))))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern240 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f240))
    rule240 = ReplacementRule(pattern240, lambda b, e, a, m, c, d, n, p, x : -(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**(p + S(1)), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule240)

    def f241(non2, e, b1, b2, a1, a2, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), RationalQ(p), Less(p, S(-1)), Or(IntegerQ(p), Not(RationalQ(m)), And(PositiveIntegerQ(n), NegativeIntegerQ(p + S(1)/2), LessEqual(S(-1), m, -n*(p + S(1))))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern241 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f241))
    rule241 = ReplacementRule(pattern241, lambda non2, e, b1, b2, a1, a2, m, c, d, n, p, x : -(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1)), x)/(a1*a2*b1*b2*n*(p + S(1))) - (e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))*(-a1*a2*d + b1*b2*c)/(a1*a2*b1*b2*e*n*(p + S(1))))
    rubi.add(rule241)

    def f242(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(m + n*(p + S(1)) + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern242 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f242))
    rule242 = ReplacementRule(pattern242, lambda b, e, a, m, c, d, n, p, x : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(b*e*(m + n*(p + S(1)) + S(1))) - (a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a + b*x**n)**p, x)/(b*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule242)

    def f243(non2, e, b1, b2, a1, a2, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), NonzeroQ(m + n*(p + S(1)) + S(1)), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern243 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f243))
    rule243 = ReplacementRule(pattern243, lambda non2, e, b1, b2, a1, a2, m, c, d, n, p, x : d*(e*x)**(m + S(1))*(a1 + b1*x**(n/S(2)))**(p + S(1))*(a2 + b2*x**(n/S(2)))**(p + S(1))/(b1*b2*e*(m + n*(p + S(1)) + S(1))) - (a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))*Int((e*x)**m*(a1 + b1*x**(n/S(2)))**p*(a2 + b2*x**(n/S(2)))**p, x)/(b1*b2*(m + n*(p + S(1)) + S(1))))
    rubi.add(rule243)

    def f244(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), PositiveIntegerQ(p), Or(IntegerQ(m), PositiveIntegerQ(S(2)*m + S(2)), Not(RationalQ(m))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern244 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_/(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f244))
    rule244 = ReplacementRule(pattern244, lambda b, e, a, m, c, d, n, p, x : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p/(c + d*x**n), x), x))
    rubi.add(rule244)

    def f245(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, n), Less(m, S(-1)), Greater(n, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern245 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_),CustomConstraint(f245))
    rule245 = ReplacementRule(pattern245, lambda b, e, a, m, c, d, n, p, x : c**S(2)*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*Simp(-a*d**S(2)*x**n*(m + S(1)) + b*c**S(2)*n*(p + S(1)) + c*(m + S(1))*(-S(2)*a*d + b*c), x), x)/(a*(m + S(1))))
    rubi.add(rule245)

    def f246(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(p), Less(p, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern246 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_),CustomConstraint(f246))
    rule246 = ReplacementRule(pattern246, lambda b, e, a, m, c, d, n, p, x : Int((e*x)**m*(a + b*x**n)**(p + S(1))*Simp(a*b*d**S(2)*n*x**n*(p + S(1)) + b**S(2)*c**S(2)*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)**S(2), x), x)/(a*b**S(2)*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(-a*d + b*c)**S(2)/(a*b**S(2)*e*n*(p + S(1))))
    rubi.add(rule246)

    def f247(b, e, a, m, c, d, n, p, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), NonzeroQ(m + n*(p + S(2)) + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern247 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**S(2), x_),CustomConstraint(f247))
    rule247 = ReplacementRule(pattern247, lambda b, e, a, m, c, d, n, p, x : d**S(2)*e**(-n + S(-1))*(e*x)**(m + n + S(1))*(a + b*x**n)**(p + S(1))/(b*(m + n*(p + S(2)) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*Simp(b*c**S(2)*(m + n*(p + S(2)) + S(1)) + d*x**n*(S(2)*b*c*n*(p + S(1)) + (-a*d + S(2)*b*c)*(m + n + S(1))), x), x)/(b*(m + n*(p + S(2)) + S(1))))
    rubi.add(rule247)

    def f248(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(p, x), FreeQ(q, x)])
    pattern248 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f248), )
    def With248(b, a, m, c, d, n, p, q, x):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q, x), x, x**k)/k
        print("Unable to Integrate")
    rule248 = ReplacementRule(pattern248, lambda b, a, m, c, d, n, p, q, x : With248(b, a, m, c, d, n, p, q, x))
    rubi.add(rule248)

    def f249(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), FractionQ(m), IntegerQ(p), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x), FreeQ(q, x)])
    pattern249 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f249), )
    def With249(b, e, a, m, c, d, n, p, q, x):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(k*n))**p*(c + d*e**(-n)*x**(k*n))**q, x), x, (e*x)**(S(1)/k))/e
    rule249 = ReplacementRule(pattern249, lambda b, e, a, m, c, d, n, p, q, x : With249(b, e, a, m, c, d, n, p, q, x))
    rubi.add(rule249)

    def f250(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, p, q), Less(p, S(-1)), Greater(q, S(0)), Greater(m - n + S(1), S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern250 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f250))
    rule250 = ReplacementRule(pattern250, lambda b, e, a, m, c, d, n, p, q, x : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(q + S(-1)) + S(1)), x), x)/(b*n*(p + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*n*(p + S(1))))
    rubi.add(rule250)

    def f251(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(p, q), Less(p, S(-1)), Greater(q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern251 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f251))
    rule251 = ReplacementRule(pattern251, lambda b, e, a, m, c, d, n, p, q, x : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule251)

    def f252(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(p, q), Less(p, S(-1)), Less(S(0), q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern252 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f252))
    rule252 = ReplacementRule(pattern252, lambda b, e, a, m, c, d, n, p, q, x : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule252)

    def f253(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, p), Less(p, S(-1)), Greater(m - n + S(1), n), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(q, x)])
    pattern253 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f253))
    rule253 = ReplacementRule(pattern253, lambda b, e, a, m, c, d, n, p, q, x : -a*e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*n*(p + S(1))*(-a*d + b*c)) + e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*q - n + S(1)) + b*c*n*(p + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule253)

    def f254(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, p), Less(p, S(-1)), Inequality(n, GreaterEqual, m - n + S(1), Greater, S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(q, x)])
    pattern254 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f254))
    rule254 = ReplacementRule(pattern254, lambda b, e, a, m, c, d, n, p, q, x : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m - n + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(n*(p + S(1))*(-a*d + b*c)) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule254)

    def f255(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(p), Less(p, S(-1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(q, x)])
    pattern255 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f255))
    rule255 = ReplacementRule(pattern255, lambda b, e, a, m, c, d, n, p, q, x : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule255)

    def f256(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, p, q), Greater(q, S(0)), Less(m, S(-1)), Greater(p, S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x)])
    pattern256 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f256))
    rule256 = ReplacementRule(pattern256, lambda b, e, a, m, c, d, n, p, q, x : -e**(-n)*n*Int((e*x)**(m + n)*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*d*q + b*c*p + b*d*x**n*(p + q), x), x)/(m + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + S(1))))
    rubi.add(rule256)

    def f257(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, q), Greater(q, S(1)), Less(m, S(-1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern257 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f257))
    rule257 = ReplacementRule(pattern257, lambda b, e, a, m, c, d, n, p, q, x : c*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(a*e*(m + S(1))) - e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*n*(a*d*(q + S(-1)) + b*c*(p + S(1))) + c*(m + S(1))*(-a*d + b*c) + d*x**n*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)), x), x)/(a*(m + S(1))))
    rubi.add(rule257)

    def f258(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, q), Less(S(0), q, S(1)), Less(m, S(-1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern258 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f258))
    rule258 = ReplacementRule(pattern258, lambda b, e, a, m, c, d, n, p, q, x : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(b*c*(m + S(1)) + d*x**n*(b*n*(p + q + S(1)) + b*(m + S(1))) + n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*(m + S(1))))
    rubi.add(rule258)

    def f259(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(p, q), Greater(q, S(0)), Greater(p, S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern259 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f259))
    rule259 = ReplacementRule(pattern259, lambda b, e, a, m, c, d, n, p, q, x : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule259)

    def f260(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(q), Greater(q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x)])
    pattern260 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f260))
    rule260 = ReplacementRule(pattern260, lambda b, e, a, m, c, d, n, p, q, x : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule260)

    def f261(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m, q), Greater(q, S(0)), Greater(m - n + S(1), S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x)])
    pattern261 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f261))
    rule261 = ReplacementRule(pattern261, lambda b, e, a, m, c, d, n, p, q, x : -e**n*Int((e*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(a*c*(m - n + S(1)) + x**n*(a*d*(m - n + S(1)) - n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))) + e**(n + S(-1))*(e*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule261)

    def f262(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m), Greater(m - n + S(1), n), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x), FreeQ(q, x)])
    pattern262 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f262))
    rule262 = ReplacementRule(pattern262, lambda b, e, a, m, c, d, n, p, q, x : -e**(S(2)*n)*Int((e*x)**(m - S(2)*n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*(m - S(2)*n + S(1)) + x**n*(a*d*(m + n*(q + S(-1)) + S(1)) + b*c*(m + n*(p + S(-1)) + S(1))), x), x)/(b*d*(m + n*(p + q) + S(1))) + e**(S(2)*n + S(-1))*(e*x)**(m - S(2)*n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q) + S(1))))
    rubi.add(rule262)

    def f263(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m), Less(m, S(-1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x), FreeQ(q, x)])
    pattern263 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f263))
    rule263 = ReplacementRule(pattern263, lambda b, e, a, m, c, d, n, p, q, x : -e**(-n)*Int((e*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(a*d*q + b*c*p) + (a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))) + (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*e*(m + S(1))))
    rubi.add(rule263)

    def f264(b, e, a, m, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), RationalQ(m), LessEqual(n, m, S(2)*n + S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern264 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f264))
    rule264 = ReplacementRule(pattern264, lambda b, e, a, m, c, d, n, x : -a*e**n*Int((e*x)**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*e**n*Int((e*x)**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule264)

    def f265(b, e, a, m, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PositiveIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x)])
    pattern265 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f265))
    rule265 = ReplacementRule(pattern265, lambda b, e, a, m, c, d, n, x : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule265)

    def f266(b, a, m, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegersQ(m/S(2), n/S(2)), Less(S(0), m - n + S(1), n), LessEqual(n, S(4)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern266 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f266))
    rule266 = ReplacementRule(pattern266, lambda b, a, m, c, d, n, x : -a*Int(x**(m - n)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(x**(m - n)/sqrt(c + d*x**n), x)/b)
    rubi.add(rule266)

    def f267(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern267 = Pattern(Integral(x_**S(2)/((a_ + x_**S(4)*WC('b', S(1)))*sqrt(c_ + x_**S(4)*WC('d', S(1)))), x_),CustomConstraint(f267), )
    def With267(b, a, c, d, x):
        r = Numerator(Rt(-a/b, S(2)))
        s = Denominator(Rt(-a/b, S(2)))
        return -s*Int(S(1)/(sqrt(c + d*x**S(4))*(r - s*x**S(2))), x)/(S(2)*b) + s*Int(S(1)/(sqrt(c + d*x**S(4))*(r + s*x**S(2))), x)/(S(2)*b)
    rule267 = ReplacementRule(pattern267, lambda b, a, c, d, x : With267(b, a, c, d, x))
    rubi.add(rule267)

    def f268(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(-a*d + S(4)*b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern268 = Pattern(Integral(x_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_),CustomConstraint(f268), )
    def With268(b, a, c, d, x):
        q = Rt(d/c, S(3))
        return -S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) - sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*sqrt(S(3))*q*ArcTan(sqrt(S(3))/S(3) + S(2)**(S(2)/3)*sqrt(S(3))*(sqrt(c) + sqrt(c + d*x**S(3)))/(S(3)*sqrt(c)*q*x))/(S(18)*b*sqrt(c)) + S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) - sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) - S(2)**(S(1)/3)*q*log(-S(2)**(S(1)/3)*q*x + S(1) + sqrt(c + d*x**S(3))/sqrt(c))/(S(12)*b*sqrt(c)) + S(2)**(S(1)/3)*q*atanh(sqrt(c + d*x**S(3))/sqrt(c))/(S(18)*b*sqrt(c))
    rule268 = ReplacementRule(pattern268, lambda b, a, c, d, x : With268(b, a, c, d, x))
    rubi.add(rule268)

    def f269(b, a, m, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(-a*d + S(4)*b*c), PositiveIntegerQ(m/S(3) + S(-1)/3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern269 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_),CustomConstraint(f269))
    rule269 = ReplacementRule(pattern269, lambda b, a, m, c, d, x : -a*Int(x**(m + S(-3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b + Int(x**(m + S(-3))/sqrt(c + d*x**S(3)), x)/b)
    rubi.add(rule269)

    def f270(b, a, m, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(-a*d + S(4)*b*c), NegativeIntegerQ(m/S(3) + S(-1)/3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern270 = Pattern(Integral(x_**m_/((a_ + x_**S(3)*WC('b', S(1)))*sqrt(c_ + x_**S(3)*WC('d', S(1)))), x_),CustomConstraint(f270))
    rule270 = ReplacementRule(pattern270, lambda b, a, m, c, d, x : -b*Int(x**(m + S(3))/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/a + Int(x**m/sqrt(c + d*x**S(3)), x)/a)
    rubi.add(rule270)

    def f271(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern271 = Pattern(Integral(x_**S(2)*sqrt(c_ + x_**S(4)*WC('d', S(1)))/(a_ + x_**S(4)*WC('b', S(1))), x_),CustomConstraint(f271))
    rule271 = ReplacementRule(pattern271, lambda b, a, c, d, x : d*Int(x**S(2)/sqrt(c + d*x**S(4)), x)/b + (-a*d + b*c)*Int(x**S(2)/((a + b*x**S(4))*sqrt(c + d*x**S(4))), x)/b)
    rubi.add(rule271)

    def f272(b, a, m, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), ZeroQ(-a*d + S(4)*b*c), IntegerQ(m/S(3) + S(-1)/3), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern272 = Pattern(Integral(x_**WC('m', S(1))*sqrt(c_ + x_**S(3)*WC('d', S(1)))/(a_ + x_**S(3)*WC('b', S(1))), x_),CustomConstraint(f272))
    rule272 = ReplacementRule(pattern272, lambda b, a, m, c, d, x : d*Int(x**m/sqrt(c + d*x**S(3)), x)/b + (-a*d + b*c)*Int(x**m/((a + b*x**S(3))*sqrt(c + d*x**S(3))), x)/b)
    rubi.add(rule272)

    def f273(b, a, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), PosQ(b/a), PosQ(d/c), Not(SimplerSqrtQ(b/a, d/c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern273 = Pattern(Integral(x_**S(2)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))), x_),CustomConstraint(f273))
    rule273 = ReplacementRule(pattern273, lambda b, a, c, d, x : -c*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/b + x*sqrt(a + b*x**S(2))/(b*sqrt(c + d*x**S(2))))
    rubi.add(rule273)

    def f274(b, a, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), Or(EqQ(n, S(2)), EqQ(n, S(4))), Not(And(EqQ(n, S(2)), SimplerSqrtQ(-b/a, -d/c))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern274 = Pattern(Integral(x_**n_/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f274))
    rule274 = ReplacementRule(pattern274, lambda b, a, c, d, n, x : -a*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b + Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b)
    rubi.add(rule274)

    def f275(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), IntegersQ(p + (m + S(1))/n, q), Less(S(-1), p, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x)])
    pattern275 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f275), )
    def With275(b, a, m, c, d, n, p, q, x):
        k = Denominator(p)
        return a**(p + (m + S(1))/n)*k*Subst(Int(x**(k*(m + S(1))/n + S(-1))*(c - x**k*(-a*d + b*c))**q*(-b*x**k + S(1))**(-p - q + S(-1) - (m + S(1))/n), x), x, x**(n/k)*(a + b*x**n)**(-S(1)/k))/n
    rule275 = ReplacementRule(pattern275, lambda b, a, m, c, d, n, p, q, x : With275(b, a, m, c, d, n, p, q, x))
    rubi.add(rule275)

    def f276(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegativeIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(p, x), FreeQ(q, x)])
    pattern276 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f276))
    rule276 = ReplacementRule(pattern276, lambda b, a, m, c, d, n, p, q, x : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x))
    rubi.add(rule276)

    def f277(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), FractionQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(p, x), FreeQ(q, x)])
    pattern277 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f277), )
    def With277(b, e, a, m, c, d, n, p, q, x):
        g = Denominator(m)
        return -g*Subst(Int(x**(-g*(m + S(1)) + S(-1))*(a + b*e**(-n)*x**(-g*n))**p*(c + d*e**(-n)*x**(-g*n))**q, x), x, (e*x)**(-S(1)/g))/e
    rule277 = ReplacementRule(pattern277, lambda b, e, a, m, c, d, n, p, q, x : With277(b, e, a, m, c, d, n, p, q, x))
    rubi.add(rule277)

    def f278(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NegativeIntegerQ(n), Not(RationalQ(m)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern278 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f278))
    rule278 = ReplacementRule(pattern278, lambda b, e, a, m, c, d, n, p, q, x : -(e*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q, x), x, S(1)/x))
    rubi.add(rule278)

    def f279(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern279 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f279), )
    def With279(b, a, m, c, d, n, p, q, x):
        g = Denominator(n)
        return g*Subst(Int(x**(g*(m + S(1)) + S(-1))*(a + b*x**(g*n))**p*(c + d*x**(g*n))**q, x), x, x**(S(1)/g))
    rule279 = ReplacementRule(pattern279, lambda b, a, m, c, d, n, p, q, x : With279(b, a, m, c, d, n, p, q, x))
    rubi.add(rule279)

    def f280(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern280 = Pattern(Integral((e_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f280))
    rule280 = ReplacementRule(pattern280, lambda b, e, a, m, c, d, n, p, q, x : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule280)

    def f281(b, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegerQ(n/(m + S(1))), Not(IntegerQ(n)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern281 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f281))
    rule281 = ReplacementRule(pattern281, lambda b, a, m, c, d, n, p, q, x : Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule281)

    def f282(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegerQ(n/(m + S(1))), Not(IntegerQ(n)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern282 = Pattern(Integral((e_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f282))
    rule282 = ReplacementRule(pattern282, lambda b, e, a, m, c, d, n, p, q, x : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule282)

    def f283(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Less(p, S(-1)), Greater(q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern283 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f283))
    rule283 = ReplacementRule(pattern283, lambda b, e, a, m, c, d, n, p, q, x : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + S(1)) + (m + S(1))*(-a*d + b*c)) + d*x**n*(b*c*n*(p + S(1)) + (-a*d + b*c)*(m + n*(q + S(-1)) + S(1))), x), x)/(a*b*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(-a*d + b*c)/(a*b*e*n*(p + S(1))))
    rubi.add(rule283)

    def f284(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Less(p, S(-1)), Less(S(0), q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern284 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f284))
    rule284 = ReplacementRule(pattern284, lambda b, e, a, m, c, d, n, p, q, x : Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(m + n*(p + S(1)) + S(1)) + d*x**n*(m + n*(p + q + S(1)) + S(1)), x), x)/(a*n*(p + S(1))) - (e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*e*n*(p + S(1))))
    rubi.add(rule284)

    def f285(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p), Less(p, S(-1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(q, x)])
    pattern285 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f285))
    rule285 = ReplacementRule(pattern285, lambda b, e, a, m, c, d, n, p, q, x : -b*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*e*n*(p + S(1))*(-a*d + b*c)) + Int((e*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(b*c*(m + S(1)) + b*d*x**n*(m + n*(p + q + S(2)) + S(1)) + n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule285)

    def f286(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(p, q), Greater(q, S(0)), Greater(p, S(0)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern286 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f286))
    rule286 = ReplacementRule(pattern286, lambda b, e, a, m, c, d, n, p, q, x : n*Int((e*x)**m*(a + b*x**n)**(p + S(-1))*(c + d*x**n)**(q + S(-1))*Simp(a*c*(p + q) + x**n*(a*d*(p + q) + q*(-a*d + b*c)), x), x)/(m + n*(p + q) + S(1)) + (e*x)**(m + S(1))*(a + b*x**n)**p*(c + d*x**n)**q/(e*(m + n*(p + q) + S(1))))
    rubi.add(rule286)

    def f287(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), RationalQ(q), Greater(q, S(1)), IntBinomialQ(a, b, c, d, e, m, n, p, q, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern287 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f287))
    rule287 = ReplacementRule(pattern287, lambda b, e, a, m, c, d, n, p, q, x : d*(e*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))/(b*e*(m + n*(p + q) + S(1))) + Int((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-2))*Simp(c*(b*c*n*(p + q) + (m + S(1))*(-a*d + b*c)) + x**n*(b*c*d*n*(p + q) + d*n*(q + S(-1))*(-a*d + b*c) + d*(m + S(1))*(-a*d + b*c)), x), x)/(b*(m + n*(p + q) + S(1))))
    rubi.add(rule287)

    def f288(b, a, m, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), Or(ZeroQ(m - n), ZeroQ(m - S(2)*n + S(1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x)])
    pattern288 = Pattern(Integral(x_**m_/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f288))
    rule288 = ReplacementRule(pattern288, lambda b, a, m, c, d, n, x : -a*Int(x**(m - n)/(a + b*x**n), x)/(-a*d + b*c) + c*Int(x**(m - n)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule288)

    def f289(b, e, a, m, c, d, n, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(n, x), FreeQ(m, x)])
    pattern289 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f289))
    rule289 = ReplacementRule(pattern289, lambda b, e, a, m, c, d, n, x : b*Int((e*x)**m/(a + b*x**n), x)/(-a*d + b*c) - d*Int((e*x)**m/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule289)

    def f290(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), IntegersQ(m, p, q), GreaterEqual(p, S(-2)), Or(GreaterEqual(q, S(-2)), And(Equal(q, S(-3)), IntegerQ(m/S(2) + S(-1)/2))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x)])
    pattern290 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f290))
    rule290 = ReplacementRule(pattern290, lambda b, e, a, m, c, d, n, p, q, x : Int(ExpandIntegrand((e*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x))
    rubi.add(rule290)

    def f291(b, mn, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(q), Or(PosQ(n), Not(IntegerQ(p))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern291 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1)), x_),CustomConstraint(f291))
    rule291 = ReplacementRule(pattern291, lambda b, mn, a, m, c, d, n, p, q, x : Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule291)

    def f292(b, mn, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), Not(IntegerQ(q)), Not(IntegerQ(p)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern292 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_, x_),CustomConstraint(f292))
    rule292 = ReplacementRule(pattern292, lambda b, mn, a, m, c, d, n, p, q, x : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(c*x**n + d)**q, x))
    rubi.add(rule292)

    def f293(b, a, mn, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern293 = Pattern(Integral((e_*x_)**m_*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f293))
    rule293 = ReplacementRule(pattern293, lambda b, a, mn, e, m, c, d, n, p, q, x : e**IntPart(m)*x**(-FracPart(m))*(e*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q, x))
    rubi.add(rule293)

    def f294(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(m + S(1)), NonzeroQ(m - n + S(1)), PositiveQ(a), PositiveQ(c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern294 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f294))
    rule294 = ReplacementRule(pattern294, lambda b, e, a, m, c, d, n, p, q, x : a**p*c**q*(e*x)**(m + S(1))*AppellF1((m + S(1))/n, -p, -q, S(1) + (m + S(1))/n, -b*x**n/a, -d*x**n/c)/(e*(m + S(1))))
    rubi.add(rule294)

    def f295(b, e, a, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-a*d + b*c), NonzeroQ(m + S(1)), NonzeroQ(m - n + S(1)), Not(PositiveQ(a)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern295 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_, x_),CustomConstraint(f295))
    rule295 = ReplacementRule(pattern295, lambda b, e, a, m, c, d, n, p, q, x : a**IntPart(p)*(S(1) + b*x**n/a)**(-FracPart(p))*(a + b*x**n)**FracPart(p)*Int((e*x)**m*(S(1) + b*x**n/a)**p*(c + d*x**n)**q, x))
    rubi.add(rule295)

    def f296(b, a, x, c, m, d, n, p, q, v):
        return functools.reduce(operator.and_, [ LinearQ(v, x), IntegerQ(m), NonzeroQ(v - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern296 = Pattern(Integral(x_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_),CustomConstraint(f296))
    rule296 = ReplacementRule(pattern296, lambda b, a, x, c, m, d, n, p, q, v : Coefficient(v, x, S(1))**(-m + S(-1))*Subst(Int(SimplifyIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(x - Coefficient(v, x, S(0)))**m, x), x), x, v))
    rubi.add(rule296)

    def f297(b, a, x, c, m, d, n, u, p, q, v):
        return functools.reduce(operator.and_, [ LinearPairQ(u, v, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern297 = Pattern(Integral(u_**WC('m', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_),CustomConstraint(f297))
    rule297 = ReplacementRule(pattern297, lambda b, a, x, c, m, d, n, u, p, q, v : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule297)

    def f298(non2, b1, b2, a1, a2, c, d, n, u, p, q, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern298 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_),CustomConstraint(f298))
    rule298 = ReplacementRule(pattern298, lambda non2, b1, b2, a1, a2, c, d, n, u, p, q, x : Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule298)

    def f299(non2, e, b1, b2, n2, a1, a2, c, d, n, u, p, q, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(-S(2)*n + n2), ZeroQ(a1*b2 + a2*b1), Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2))), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern299 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_),CustomConstraint(f299))
    rule299 = ReplacementRule(pattern299, lambda non2, e, b1, b2, n2, a1, a2, c, d, n, u, p, q, x : Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule299)

    def f300(non2, b1, b2, a1, a2, c, d, n, u, p, q, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(a1*b2 + a2*b1), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern300 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**p_*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**p_*(c_ + x_**WC('n', S(1))*WC('d', S(1)))**WC('q', S(1))*WC('u', S(1)), x_),CustomConstraint(f300))
    rule300 = ReplacementRule(pattern300, lambda non2, b1, b2, a1, a2, c, d, n, u, p, q, x : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(c + d*x**n)**q*(a1*a2 + b1*b2*x**n)**p, x))
    rubi.add(rule300)

    def f301(non2, e, b1, b2, n2, a1, a2, c, d, n, u, p, q, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + non2), ZeroQ(-S(2)*n + n2), ZeroQ(a1*b2 + a2*b1), FreeQ(a1, x), FreeQ(b1, x), FreeQ(a2, x), FreeQ(b2, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern301 = Pattern(Integral((a1_ + x_**WC('non2', S(1))*WC('b1', S(1)))**WC('p', S(1))*(a2_ + x_**WC('non2', S(1))*WC('b2', S(1)))**WC('p', S(1))*(c_ + x_**WC('n', S(1))*WC('d', S(1)) + x_**WC('n2', S(1))*WC('e', S(1)))**WC('q', S(1))*WC('u', S(1)), x_),CustomConstraint(f301))
    rule301 = ReplacementRule(pattern301, lambda non2, e, b1, b2, n2, a1, a2, c, d, n, u, p, q, x : (a1 + b1*x**(n/S(2)))**FracPart(p)*(a2 + b2*x**(n/S(2)))**FracPart(p)*(a1*a2 + b1*b2*x**n)**(-FracPart(p))*Int(u*(a1*a2 + b1*b2*x**n)**p*(c + d*x**n + e*x**(S(2)*n))**q, x))
    rubi.add(rule301)

    def f302(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p, q, r), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern302 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f302))
    rule302 = ReplacementRule(pattern302, lambda b, a, e, f, r, c, d, n, p, q, x : Int(ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule302)

    def f303(b, a, e, f, c, d, n, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, n), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern303 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f303))
    rule303 = ReplacementRule(pattern303, lambda b, a, e, f, c, d, n, x : (-a*f + b*e)*Int(S(1)/(a + b*x**n), x)/(-a*d + b*c) - (-c*f + d*e)*Int(S(1)/(c + d*x**n), x)/(-a*d + b*c))
    rubi.add(rule303)

    def f304(b, a, e, f, c, d, n, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, n), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern304 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/((a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f304))
    rule304 = ReplacementRule(pattern304, lambda b, a, e, f, c, d, n, x : f*Int(S(1)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/((a + b*x**n)*sqrt(c + d*x**n)), x)/b)
    rubi.add(rule304)

    def f305(b, a, e, f, c, d, n, x):
        return functools.reduce(operator.and_, [ Not(And(ZeroQ(n + S(-2)), Or(And(PosQ(b/a), PosQ(d/c)), And(NegQ(b/a), Or(PosQ(d/c), And(PositiveQ(a), Or(Not(PositiveQ(c)), SimplerSqrtQ(-b/a, -d/c)))))))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern305 = Pattern(Integral((e_ + x_**n_*WC('f', S(1)))/(sqrt(a_ + x_**n_*WC('b', S(1)))*sqrt(c_ + x_**n_*WC('d', S(1)))), x_),CustomConstraint(f305))
    rule305 = ReplacementRule(pattern305, lambda b, a, e, f, c, d, n, x : f*Int(sqrt(a + b*x**n)/sqrt(c + d*x**n), x)/b + (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**n)*sqrt(c + d*x**n)), x)/b)
    rubi.add(rule305)

    def f306(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(b/a), PosQ(d/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern306 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_),CustomConstraint(f306))
    rule306 = ReplacementRule(pattern306, lambda b, a, e, f, c, d, x : (-a*f + b*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(a + b*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule306)

    def f307(b, a, e, f, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ RationalQ(p, q), Less(p, S(-1)), Greater(q, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x)])
    pattern307 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f307))
    rule307 = ReplacementRule(pattern307, lambda b, a, e, f, c, d, n, p, q, x : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*n*(p + S(1))) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + S(1)) + b*e) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(n*q + S(1))), x), x)/(a*b*n*(p + S(1))))
    rubi.add(rule307)

    def f308(b, a, e, f, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(p, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(q, x)])
    pattern308 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f308))
    rule308 = ReplacementRule(pattern308, lambda b, a, e, f, c, d, n, p, q, x : -x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*n*(p + S(1))*(-a*d + b*c)) + Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule308)

    def f309(b, a, e, f, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ RationalQ(q), Greater(q, S(0)), NonzeroQ(n*(p + q + S(1)) + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x)])
    pattern309 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f309))
    rule309 = ReplacementRule(pattern309, lambda b, a, e, f, c, d, n, p, q, x : f*x*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*(n*(p + q + S(1)) + S(1))) + Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(-a*f + b*e*n*(p + q + S(1)) + b*e) + x**n*(b*d*e*n*(p + q + S(1)) + d*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(n*(p + q + S(1)) + S(1))))
    rubi.add(rule309)

    def f310(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern310 = Pattern(Integral((e_ + x_**S(4)*WC('f', S(1)))/((a_ + x_**S(4)*WC('b', S(1)))**(S(3)/4)*(c_ + x_**S(4)*WC('d', S(1)))), x_),CustomConstraint(f310))
    rule310 = ReplacementRule(pattern310, lambda b, a, e, f, c, d, x : (-a*f + b*e)*Int((a + b*x**S(4))**(S(-3)/4), x)/(-a*d + b*c) - (-c*f + d*e)*Int((a + b*x**S(4))**(S(1)/4)/(c + d*x**S(4)), x)/(-a*d + b*c))
    rubi.add(rule310)

    def f311(b, a, e, f, c, d, n, p, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, p, n), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(n, x)])
    pattern311 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f311))
    rule311 = ReplacementRule(pattern311, lambda b, a, e, f, c, d, n, p, x : f*Int((a + b*x**n)**p, x)/d + (-c*f + d*e)*Int((a + b*x**n)**p/(c + d*x**n), x)/d)
    rubi.add(rule311)

    def f312(b, a, e, f, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, n, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern312 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f312))
    rule312 = ReplacementRule(pattern312, lambda b, a, e, f, c, d, n, p, q, x : e*Int((a + b*x**n)**p*(c + d*x**n)**q, x) + f*Int(x**n*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule312)

    def f313(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern313 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f313))
    rule313 = ReplacementRule(pattern313, lambda b, a, e, f, c, d, x : b*Int(S(1)/((a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c) - d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*d + b*c))
    rubi.add(rule313)

    def f314(e, f, c, d, x):
        return functools.reduce(operator.and_, [ NonzeroQ(-c*f + d*e), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern314 = Pattern(Integral(S(1)/(x_**S(2)*(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f314))
    rule314 = ReplacementRule(pattern314, lambda e, f, c, d, x : -d*Int(S(1)/((c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/c + Int(S(1)/(x**S(2)*sqrt(e + f*x**S(2))), x)/c)
    rubi.add(rule314)

    def f315(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PositiveQ(d/c), PositiveQ(f/e), Not(SimplerSqrtQ(d/c, f/e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern315 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f315))
    rule315 = ReplacementRule(pattern315, lambda b, a, e, f, c, d, x : d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule315)

    def f316(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ Not(SimplerSqrtQ(-f/e, -d/c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern316 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f316))
    rule316 = ReplacementRule(pattern316, lambda b, a, e, f, c, d, x : d*Int(sqrt(e + f*x**S(2))/sqrt(c + d*x**S(2)), x)/b + (-a*d + b*c)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b)
    rubi.add(rule316)

    def f317(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(f/e), Not(SimplerSqrtQ(d/c, f/e)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern317 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f317))
    rule317 = ReplacementRule(pattern317, lambda b, a, e, f, c, d, x : b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*f + b*e) - f*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(-a*f + b*e))
    rubi.add(rule317)

    def f318(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), PositiveQ(c), PositiveQ(e), Not(And(NegQ(f/e), SimplerSqrtQ(-f/e, -d/c))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern318 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f318))
    rule318 = ReplacementRule(pattern318, lambda b, a, e, f, c, d, x : EllipticPi(b*c/(a*d), asin(x*Rt(-d/c, S(2))), c*f/(d*e))/(a*sqrt(c)*sqrt(e)*Rt(-d/c, S(2))))
    rubi.add(rule318)

    def f319(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ Not(PositiveQ(c)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern319 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f319))
    rule319 = ReplacementRule(pattern319, lambda b, a, e, f, c, d, x : sqrt(S(1) + d*x**S(2)/c)*Int(S(1)/(sqrt(S(1) + d*x**S(2)/c)*(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/sqrt(c + d*x**S(2)))
    rubi.add(rule319)

    def f320(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern320 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f320))
    rule320 = ReplacementRule(pattern320, lambda b, a, e, f, c, d, x : c*sqrt(e + f*x**S(2))*EllipticPi(S(1) - b*c/(a*d), ArcTan(x*Rt(d/c, S(2))), -c*f/(d*e) + S(1))/(a*e*sqrt(c*(e + f*x**S(2))/(e*(c + d*x**S(2))))*sqrt(c + d*x**S(2))*Rt(d/c, S(2))))
    rubi.add(rule320)

    def f321(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ NegQ(d/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern321 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f321))
    rule321 = ReplacementRule(pattern321, lambda b, a, e, f, c, d, x : d*Int(S(1)/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b + (-a*d + b*c)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/b)
    rubi.add(rule321)

    def f322(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(f/e), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern322 = Pattern(Integral(sqrt(e_ + x_**S(2)*WC('f', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_),CustomConstraint(f322))
    rule322 = ReplacementRule(pattern322, lambda b, a, e, f, c, d, x : b*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - d*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule322)

    def f323(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(f/e), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern323 = Pattern(Integral((e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2)/((a_ + x_**S(2)*WC('b', S(1)))*(c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)), x_),CustomConstraint(f323))
    rule323 = ReplacementRule(pattern323, lambda b, a, e, f, c, d, x : (-a*f + b*e)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(-a*d + b*c) - (-c*f + d*e)*Int(sqrt(e + f*x**S(2))/(c + d*x**S(2))**(S(3)/2), x)/(-a*d + b*c))
    rubi.add(rule323)

    def f324(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ(d/c), PosQ(f/e), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern324 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f324))
    rule324 = ReplacementRule(pattern324, lambda b, a, e, f, c, d, x : d*Int(sqrt(e + f*x**S(2))*(-a*d + S(2)*b*c + b*d*x**S(2))/sqrt(c + d*x**S(2)), x)/b**S(2) + (-a*d + b*c)**S(2)*Int(sqrt(e + f*x**S(2))/((a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/b**S(2))
    rubi.add(rule324)

    def f325(b, a, e, f, r, c, d, q, x):
        return functools.reduce(operator.and_, [ RationalQ(q, r), Less(q, S(-1)), Greater(r, S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern325 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f325))
    rule325 = ReplacementRule(pattern325, lambda b, a, e, f, r, c, d, q, x : b*(-a*f + b*e)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**(r + S(-1))/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - Int((c + d*x**S(2))**q*(e + f*x**S(2))**(r + S(-1))*(-a*d**S(2)*e - b*c**S(2)*f + S(2)*b*c*d*e + d**S(2)*x**S(2)*(-a*f + b*e)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule325)

    def f326(b, a, e, f, r, c, d, q, x):
        return functools.reduce(operator.and_, [ RationalQ(q), Greater(q, S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(r, x)])
    pattern326 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f326))
    rule326 = ReplacementRule(pattern326, lambda b, a, e, f, r, c, d, q, x : d*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r, x)/b + (-a*d + b*c)*Int((c + d*x**S(2))**(q + S(-1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/b)
    rubi.add(rule326)

    def f327(b, a, e, f, r, c, d, q, x):
        return functools.reduce(operator.and_, [ RationalQ(q), Less(q, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(r, x)])
    pattern327 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f327))
    rule327 = ReplacementRule(pattern327, lambda b, a, e, f, r, c, d, q, x : b**S(2)*Int((c + d*x**S(2))**(q + S(2))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c)**S(2) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r*(-a*d + S(2)*b*c + b*d*x**S(2)), x)/(-a*d + b*c)**S(2))
    rubi.add(rule327)

    def f328(b, a, e, f, r, c, d, q, x):
        return functools.reduce(operator.and_, [ RationalQ(q), LessEqual(q, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(r, x)])
    pattern328 = Pattern(Integral((c_ + x_**S(2)*WC('d', S(1)))**q_*(e_ + x_**S(2)*WC('f', S(1)))**r_/(a_ + x_**S(2)*WC('b', S(1))), x_),CustomConstraint(f328))
    rule328 = ReplacementRule(pattern328, lambda b, a, e, f, r, c, d, q, x : b*Int((c + d*x**S(2))**(q + S(1))*(e + f*x**S(2))**r/(a + b*x**S(2)), x)/(-a*d + b*c) - d*Int((c + d*x**S(2))**q*(e + f*x**S(2))**r, x)/(-a*d + b*c))
    rubi.add(rule328)

    def f329(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern329 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))/(a_ + x_**S(2)*WC('b', S(1)))**S(2), x_),CustomConstraint(f329))
    rule329 = ReplacementRule(pattern329, lambda b, a, e, f, c, d, x : x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))) + d*f*Int((a - b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)) + (-a**S(2)*d*f + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*b**S(2)))
    rubi.add(rule329)

    def f330(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern330 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('b', S(1)))**S(2)*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f330))
    rule330 = ReplacementRule(pattern330, lambda b, a, e, f, c, d, x : b**S(2)*x*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*a*(a + b*x**S(2))*(-a*d + b*c)*(-a*f + b*e)) - d*f*Int((a + b*x**S(2))/(sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)) + (S(3)*a**S(2)*d*f - S(2)*a*b*(c*f + d*e) + b**S(2)*c*e)*Int(S(1)/((a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*a*(-a*d + b*c)*(-a*f + b*e)))
    rubi.add(rule330)

    def f331(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(p), RationalQ(q), Greater(q, S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(r, x)])
    pattern331 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_),CustomConstraint(f331))
    rule331 = ReplacementRule(pattern331, lambda b, a, e, f, r, c, d, n, p, q, x : d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b + (-a*d + b*c)*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*(e + f*x**n)**r, x)/b)
    rubi.add(rule331)

    def f332(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(p), RationalQ(q), LessEqual(q, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(q, x)])
    pattern332 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_),CustomConstraint(f332))
    rule332 = ReplacementRule(pattern332, lambda b, a, e, f, r, c, d, n, p, q, x : b*Int((a + b*x**n)**p*(c + d*x**n)**(q + S(1))*(e + f*x**n)**r, x)/(-a*d + b*c) - d*Int((a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(e + f*x**n)**r, x)/(-a*d + b*c))
    rubi.add(rule332)

    def f333(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern333 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f333))
    rule333 = ReplacementRule(pattern333, lambda b, a, e, f, c, d, x : sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule333)

    def f334(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern334 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))/(sqrt(c_ + x_**S(2)*WC('d', S(1)))*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f334))
    rule334 = ReplacementRule(pattern334, lambda b, a, e, f, c, d, x : a*sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(S(1)/(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)*sqrt(S(1) - x**S(2)*(-a*f + b*e)/e)*(-b*x**S(2) + S(1))), x), x, x/sqrt(a + b*x**S(2)))/(c*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule334)

    def f335(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern335 = Pattern(Integral(sqrt(c_ + x_**S(2)*WC('d', S(1)))/((a_ + x_**S(2)*WC('b', S(1)))**(S(3)/2)*sqrt(e_ + x_**S(2)*WC('f', S(1)))), x_),CustomConstraint(f335))
    rule335 = ReplacementRule(pattern335, lambda b, a, e, f, c, d, x : sqrt(a*(e + f*x**S(2))/(e*(a + b*x**S(2))))*sqrt(c + d*x**S(2))*Subst(Int(sqrt(S(1) - x**S(2)*(-a*d + b*c)/c)/sqrt(S(1) - x**S(2)*(-a*f + b*e)/e), x), x, x/sqrt(a + b*x**S(2)))/(a*sqrt(a*(c + d*x**S(2))/(c*(a + b*x**S(2))))*sqrt(e + f*x**S(2))))
    rubi.add(rule335)

    def f336(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ PosQ((-c*f + d*e)/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern336 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_),CustomConstraint(f336))
    rule336 = ReplacementRule(pattern336, lambda b, a, e, f, c, d, x : b*c*(-c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f) - c*(-c*f + d*e)*Int(sqrt(a + b*x**S(2))/((c + d*x**S(2))**(S(3)/2)*sqrt(e + f*x**S(2))), x)/(S(2)*f) + d*x*sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))/(S(2)*f*sqrt(c + d*x**S(2))) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*d*f))
    rubi.add(rule336)

    def f337(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ NegQ((-c*f + d*e)/c), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern337 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/sqrt(e_ + x_**S(2)*WC('f', S(1))), x_),CustomConstraint(f337))
    rule337 = ReplacementRule(pattern337, lambda b, a, e, f, c, d, x : e*(-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/(S(2)*f) + x*sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))/(S(2)*sqrt(e + f*x**S(2))) + (-a*f + b*e)*(-S(2)*c*f + d*e)*Int(S(1)/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))*sqrt(e + f*x**S(2))), x)/(S(2)*f**S(2)) - (-a*d*f - b*c*f + b*d*e)*Int(sqrt(e + f*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(c + d*x**S(2))), x)/(S(2)*f**S(2)))
    rubi.add(rule337)

    def f338(b, a, e, f, c, d, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x)])
    pattern338 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('b', S(1)))*sqrt(c_ + x_**S(2)*WC('d', S(1)))/(e_ + x_**S(2)*WC('f', S(1)))**(S(3)/2), x_),CustomConstraint(f338))
    rule338 = ReplacementRule(pattern338, lambda b, a, e, f, c, d, x : b*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*sqrt(e + f*x**S(2))), x)/f - (-a*f + b*e)*Int(sqrt(c + d*x**S(2))/(sqrt(a + b*x**S(2))*(e + f*x**S(2))**(S(3)/2)), x)/f)
    rubi.add(rule338)

    def f339(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern339 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_),CustomConstraint(f339), )
    def With339(b, a, e, f, r, c, d, n, p, q, x):
        u = ExpandIntegrand((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x)
        if SumQ(u):
            return Int(u, x)
        print("Unable to Integrate")
    rule339 = ReplacementRule(pattern339, lambda b, a, e, f, r, c, d, n, p, q, x : With339(b, a, e, f, r, c, d, n, p, q, x))
    rubi.add(rule339)

    def f340(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern340 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_),CustomConstraint(f340))
    rule340 = ReplacementRule(pattern340, lambda b, a, e, f, r, c, d, n, p, q, x : -Subst(Int((a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r/x**S(2), x), x, S(1)/x))
    rubi.add(rule340)

    def f341(b, a, e, f, r, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, n, p, q, r), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern341 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f341))
    rule341 = ReplacementRule(pattern341, lambda b, a, e, f, r, c, d, n, p, q, x : Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule341)

    def f342(b, a, e, f, r, x, c, d, n, u, p, w, q, v):
        return functools.reduce(operator.and_, [ ZeroQ(u - v), ZeroQ(u - w), LinearQ(u, x), NonzeroQ(u - x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(n, x), FreeQ(q, x), FreeQ(r, x)])
    pattern342 = Pattern(Integral((u_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1))*(w_**n_*WC('f', S(1)) + WC('e', S(0)))**WC('r', S(1)), x_),CustomConstraint(f342))
    rule342 = ReplacementRule(pattern342, lambda b, a, e, f, r, x, c, d, n, u, p, w, q, v : Subst(Int((a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule342)

    def f343(b, a, mn, f, r, e, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x), FreeQ(r, x)])
    pattern343 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f343))
    rule343 = ReplacementRule(pattern343, lambda b, a, mn, f, r, e, c, d, n, p, q, x : Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule343)

    def f344(b, a, mn, f, r, e, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(p), IntegerQ(r), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(q, x)])
    pattern344 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f344))
    rule344 = ReplacementRule(pattern344, lambda b, a, mn, f, r, e, c, d, n, p, q, x : Int(x**(n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule344)

    def f345(b, a, mn, f, r, e, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), Not(IntegerQ(q)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern345 = Pattern(Integral((c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f345))
    rule345 = ReplacementRule(pattern345, lambda b, a, mn, f, r, e, c, d, n, p, q, x : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(-n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule345)

    def f346(b, a, r, f1, e2, n2, e1, c, d, n, q, p, f2, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + n2), ZeroQ(e1*f2 + e2*f1), Or(IntegerQ(r), And(PositiveQ(e1), PositiveQ(e2))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e1, x), FreeQ(f1, x), FreeQ(e2, x), FreeQ(f2, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern346 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_),CustomConstraint(f346))
    rule346 = ReplacementRule(pattern346, lambda b, a, r, f1, e2, n2, e1, c, d, n, q, p, f2, x : Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule346)

    def f347(b, a, r, f1, e2, n2, e1, c, d, n, q, p, f2, x):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + n2), ZeroQ(e1*f2 + e2*f1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e1, x), FreeQ(f1, x), FreeQ(e2, x), FreeQ(f2, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern347 = Pattern(Integral((a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_),CustomConstraint(f347))
    rule347 = ReplacementRule(pattern347, lambda b, a, r, f1, e2, n2, e1, c, d, n, q, p, f2, x : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule347)

    def f348(b, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(g)), IntegerQ((m + S(1))/n), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern348 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f348))
    rule348 = ReplacementRule(pattern348, lambda b, e, f, r, x, m, c, d, n, p, q, g : b**(S(1) - (m + S(1))/n)*g**m*Subst(Int((b*x)**(p + S(-1) + (m + S(1))/n)*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule348)

    def f349(b, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ Or(IntegerQ(m), PositiveQ(g)), Not(IntegerQ((m + S(1))/n)), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern349 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f349))
    rule349 = ReplacementRule(pattern349, lambda b, e, f, r, x, m, c, d, n, p, q, g : b**IntPart(p)*g**m*x**(-n*FracPart(p))*(b*x**n)**FracPart(p)*Int(x**(m + n*p)*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule349)

    def f350(b, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ Not(IntegerQ(m)), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern350 = Pattern(Integral((g_*x_)**m_*(x_**WC('n', S(1))*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f350))
    rule350 = ReplacementRule(pattern350, lambda b, e, f, r, x, m, c, d, n, p, q, g : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule350)

    def f351(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(p + S(2), q, r), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern351 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f351))
    rule351 = ReplacementRule(pattern351, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x))
    rubi.add(rule351)

    def f352(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ ZeroQ(m - n + S(1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern352 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f352))
    rule352 = ReplacementRule(pattern352, lambda b, a, e, f, r, m, c, d, n, p, q, x : Subst(Int((a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule352)

    def f353(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ IntegersQ(p, q, r), NegQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x)])
    pattern353 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f353))
    rule353 = ReplacementRule(pattern353, lambda b, a, e, f, r, m, c, d, n, p, q, x : Int(x**(m + n*(p + q + r))*(a*x**(-n) + b)**p*(c*x**(-n) + d)**q*(e*x**(-n) + f)**r, x))
    rubi.add(rule353)

    def f354(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern354 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f354))
    rule354 = ReplacementRule(pattern354, lambda b, a, e, f, r, m, c, d, n, p, q, x : Subst(Int(x**(S(-1) + (m + S(1))/n)*(a + b*x)**p*(c + d*x)**q*(e + f*x)**r, x), x, x**n)/n)
    rubi.add(rule354)

    def f355(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ IntegerQ((m + S(1))/n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern355 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f355))
    rule355 = ReplacementRule(pattern355, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule355)

    def f356(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern356 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f356), )
    def With356(b, a, e, f, r, m, c, d, n, p, q, x):
        k = GCD(m + S(1), n)
        if Unequal(k, S(1)):
            return Subst(Int(x**(S(-1) + (m + S(1))/k)*(a + b*x**(n/k))**p*(c + d*x**(n/k))**q*(e + f*x**(n/k))**r, x), x, x**k)/k
        print("Unable to Integrate")
    rule356 = ReplacementRule(pattern356, lambda b, a, e, f, r, m, c, d, n, p, q, x : With356(b, a, e, f, r, m, c, d, n, p, q, x))
    rubi.add(rule356)

    def f357(b, a, e, f, r, g, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FractionQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern357 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1)))**r_, x_),CustomConstraint(f357), )
    def With357(b, a, e, f, r, g, m, c, d, n, p, q, x):
        k = Denominator(m)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(k*n))**p*(c + d*g**(-n)*x**(k*n))**q*(e + f*g**(-n)*x**(k*n))**r, x), x, (g*x)**(S(1)/k))/g
    rule357 = ReplacementRule(pattern357, lambda b, a, e, f, r, g, m, c, d, n, p, q, x : With357(b, a, e, f, r, g, m, c, d, n, p, q, x))
    rubi.add(rule357)

    def f358(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p, q), Less(p, S(-1)), Greater(q, S(0)), Not(And(Equal(q, S(1)), SimplerQ(-a*d + b*c, -a*f + b*e))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x)])
    pattern358 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f358))
    rule358 = ReplacementRule(pattern358, lambda b, a, e, f, x, m, c, d, n, p, q, g : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule358)

    def f359(b, a, e, f, m, c, d, n, q, p, x, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, p), Less(p, S(-1)), Greater(m - n + S(1), S(0)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(q, x)])
    pattern359 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f359))
    rule359 = ReplacementRule(pattern359, lambda b, a, e, f, m, c, d, n, q, p, x, g : -g**n*Int((g*x)**(m - n)*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(-a*f + b*e)*(m - n + S(1)) + x**n*(-b*n*(p + S(1))*(c*f - d*e) + d*(-a*f + b*e)*(m + n*q + S(1))), x), x)/(b*n*(p + S(1))*(-a*d + b*c)) + g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(b*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule359)

    def f360(b, a, e, f, m, c, d, n, q, p, x, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(p), Less(p, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(q, x)])
    pattern360 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f360))
    rule360 = ReplacementRule(pattern360, lambda b, a, e, f, m, c, d, n, q, p, x, g : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule360)

    def f361(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m, q), Greater(q, S(0)), Less(m, S(-1)), Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x)])
    pattern361 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f361))
    rule361 = ReplacementRule(pattern361, lambda b, a, e, f, x, m, c, d, n, p, q, g : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(a*g*(m + S(1))) - g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + e*n*(a*d*q + b*c*(p + S(1))), x), x)/(a*(m + S(1))))
    rubi.add(rule361)

    def f362(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(q), Greater(q, S(0)), Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern362 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f362))
    rule362 = ReplacementRule(pattern362, lambda b, a, e, f, x, m, c, d, n, p, q, g : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule362)

    def f363(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m), Greater(m, n + S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x), FreeQ(q, x)])
    pattern363 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f363))
    rule363 = ReplacementRule(pattern363, lambda b, a, e, f, x, m, c, d, n, p, q, g : f*g**(n + S(-1))*(g*x)**(m - n + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(b*d*(m + n*(p + q + S(1)) + S(1))) - g**n*Int((g*x)**(m - n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m - n + S(1)) + x**n*(a*d*f*(m + n*q + S(1)) + b*(c*f*(m + n*p + S(1)) - d*e*(m + n*(p + q + S(1)) + S(1)))), x), x)/(b*d*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule363)

    def f364(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), RationalQ(m), Less(m, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x), FreeQ(q, x)])
    pattern364 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f364))
    rule364 = ReplacementRule(pattern364, lambda b, a, e, f, x, m, c, d, n, p, q, g : e*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))/(a*c*g*(m + S(1))) + g**(-n)*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*Simp(a*c*f*(m + S(1)) - b*d*e*x**n*(m + n*(p + q + S(2)) + S(1)) - e*n*(a*d*q + b*c*p) - e*(a*d + b*c)*(m + n + S(1)), x), x)/(a*c*(m + S(1))))
    rubi.add(rule364)

    def f365(b, a, e, f, m, c, d, n, p, x, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x)])
    pattern365 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f365))
    rule365 = ReplacementRule(pattern365, lambda b, a, e, f, m, c, d, n, p, x, g : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule365)

    def f366(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern366 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f366))
    rule366 = ReplacementRule(pattern366, lambda b, a, e, f, x, m, c, d, n, p, q, g : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule366)

    def f367(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ PositiveIntegerQ(n), PositiveIntegerQ(r), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x)])
    pattern367 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f367))
    rule367 = ReplacementRule(pattern367, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x) + e**(-n)*f*Int((g*x)**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**(r + S(-1)), x))
    rubi.add(rule367)

    def f368(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), IntegerQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern368 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f368))
    rule368 = ReplacementRule(pattern368, lambda b, a, e, f, r, m, c, d, n, p, q, x : -Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x))
    rubi.add(rule368)

    def f369(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), FractionQ(m), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern369 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f369), )
    def With369(b, a, e, f, r, x, m, c, d, n, p, q, g):
        k = Denominator(m)
        return -k*Subst(Int(x**(-k*(m + S(1)) + S(-1))*(a + b*g**(-n)*x**(-k*n))**p*(c + d*g**(-n)*x**(-k*n))**q*(e + f*g**(-n)*x**(-k*n))**r, x), x, (g*x)**(-S(1)/k))/g
    rule369 = ReplacementRule(pattern369, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : With369(b, a, e, f, r, x, m, c, d, n, p, q, g))
    rubi.add(rule369)

    def f370(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ NegativeIntegerQ(n), Not(RationalQ(m)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern370 = Pattern(Integral((x_*WC('g', S(1)))**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f370))
    rule370 = ReplacementRule(pattern370, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : -(g*x)**m*(S(1)/x)**m*Subst(Int(x**(-m + S(-2))*(a + b*x**(-n))**p*(c + d*x**(-n))**q*(e + f*x**(-n))**r, x), x, S(1)/x))
    rubi.add(rule370)

    def f371(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern371 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f371), )
    def With371(b, a, e, f, r, m, c, d, n, p, q, x):
        k = Denominator(n)
        return k*Subst(Int(x**(k*(m + S(1)) + S(-1))*(a + b*x**(k*n))**p*(c + d*x**(k*n))**q*(e + f*x**(k*n))**r, x), x, x**(S(1)/k))
    rule371 = ReplacementRule(pattern371, lambda b, a, e, f, r, m, c, d, n, p, q, x : With371(b, a, e, f, r, m, c, d, n, p, q, x))
    rubi.add(rule371)

    def f372(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ FractionQ(n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern372 = Pattern(Integral((g_*x_)**m_*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f372))
    rule372 = ReplacementRule(pattern372, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule372)

    def f373(b, a, e, f, r, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ IntegerQ(n/(m + S(1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern373 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f373))
    rule373 = ReplacementRule(pattern373, lambda b, a, e, f, r, m, c, d, n, p, q, x : Subst(Int((a + b*x**(n/(m + S(1))))**p*(c + d*x**(n/(m + S(1))))**q*(e + f*x**(n/(m + S(1))))**r, x), x, x**(m + S(1)))/(m + S(1)))
    rubi.add(rule373)

    def f374(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ IntegerQ(n/(m + S(1))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern374 = Pattern(Integral((g_*x_)**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f374))
    rule374 = ReplacementRule(pattern374, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule374)

    def f375(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ RationalQ(p, q), Less(p, S(-1)), Greater(q, S(0)), Not(And(Equal(q, S(1)), SimplerQ(-a*d + b*c, -a*f + b*e))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x)])
    pattern375 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f375))
    rule375 = ReplacementRule(pattern375, lambda b, a, e, f, x, m, c, d, n, p, q, g : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + S(1)) + (m + S(1))*(-a*f + b*e)) + d*x**n*(b*e*n*(p + S(1)) + (-a*f + b*e)*(m + n*q + S(1))), x), x)/(a*b*n*(p + S(1))) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*(-a*f + b*e)/(a*b*g*n*(p + S(1))))
    rubi.add(rule375)

    def f376(b, a, e, f, m, c, d, n, q, p, x, g):
        return functools.reduce(operator.and_, [ RationalQ(p), Less(p, S(-1)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(q, x)])
    pattern376 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f376))
    rule376 = ReplacementRule(pattern376, lambda b, a, e, f, m, c, d, n, q, p, x, g : Int((g*x)**m*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q*Simp(c*(m + S(1))*(-a*f + b*e) + d*x**n*(-a*f + b*e)*(m + n*(p + q + S(2)) + S(1)) + e*n*(p + S(1))*(-a*d + b*c), x), x)/(a*n*(p + S(1))*(-a*d + b*c)) - (g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**(q + S(1))*(-a*f + b*e)/(a*g*n*(p + S(1))*(-a*d + b*c)))
    rubi.add(rule376)

    def f377(b, a, e, f, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ RationalQ(q), Greater(q, S(0)), Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern377 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f377))
    rule377 = ReplacementRule(pattern377, lambda b, a, e, f, x, m, c, d, n, p, q, g : f*(g*x)**(m + S(1))*(a + b*x**n)**(p + S(1))*(c + d*x**n)**q/(b*g*(m + n*(p + q + S(1)) + S(1))) + Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**(q + S(-1))*Simp(c*(b*e*n*(p + q + S(1)) + (m + S(1))*(-a*f + b*e)) + x**n*(b*d*e*n*(p + q + S(1)) + d*(m + S(1))*(-a*f + b*e) + f*n*q*(-a*d + b*c)), x), x)/(b*(m + n*(p + q + S(1)) + S(1))))
    rubi.add(rule377)

    def f378(b, a, e, f, m, c, d, n, p, x, g):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, m, n, p), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x)])
    pattern378 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(e_ + x_**n_*WC('f', S(1)))/(c_ + x_**n_*WC('d', S(1))), x_),CustomConstraint(f378))
    rule378 = ReplacementRule(pattern378, lambda b, a, e, f, m, c, d, n, p, x, g : Int(ExpandIntegrand((g*x)**m*(a + b*x**n)**p*(e + f*x**n)/(c + d*x**n), x), x))
    rubi.add(rule378)

    def f379(b, a, e, f, m, c, d, n, q, p, x, g):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, m, n, p, q), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x)])
    pattern379 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**p_*(c_ + x_**n_*WC('d', S(1)))**q_*(e_ + x_**n_*WC('f', S(1))), x_),CustomConstraint(f379))
    rule379 = ReplacementRule(pattern379, lambda b, a, e, f, m, c, d, n, q, p, x, g : e*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q, x) + f*x**(-m)*(g*x)**m*Int(x**(m + n)*(a + b*x**n)**p*(c + d*x**n)**q, x))
    rubi.add(rule379)

    def f380(b, mn, a, f, r, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(q), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(r, x)])
    pattern380 = Pattern(Integral(x_**WC('m', S(1))*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f380))
    rule380 = ReplacementRule(pattern380, lambda b, mn, a, f, r, e, m, c, d, n, p, q, x : Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule380)

    def f381(b, a, mn, f, r, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), IntegerQ(p), IntegerQ(r), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(q, x)])
    pattern381 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f381))
    rule381 = ReplacementRule(pattern381, lambda b, a, mn, f, r, e, m, c, d, n, p, q, x : Int(x**(m + n*(p + r))*(c + d*x**(-n))**q*(a*x**(-n) + b)**p*(e*x**(-n) + f)**r, x))
    rubi.add(rule381)

    def f382(b, a, mn, f, r, e, m, c, d, n, p, q, x):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), Not(IntegerQ(q)), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern382 = Pattern(Integral(x_**WC('m', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**q_*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1))*(x_**WC('n', S(1))*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_),CustomConstraint(f382))
    rule382 = ReplacementRule(pattern382, lambda b, a, mn, f, r, e, m, c, d, n, p, q, x : x**(n*FracPart(q))*(c + d*x**(-n))**FracPart(q)*(c*x**n + d)**(-FracPart(q))*Int(x**(m - n*q)*(a + b*x**n)**p*(e + f*x**n)**r*(c*x**n + d)**q, x))
    rubi.add(rule382)

    def f383(b, mn, a, f, r, e, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ EqQ(mn, -n), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern383 = Pattern(Integral((g_*x_)**m_*(a_ + x_**WC('n', S(1))*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**WC('mn', S(1))*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**WC('n', S(1))*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f383))
    rule383 = ReplacementRule(pattern383, lambda b, mn, a, f, r, e, x, m, c, d, n, p, q, g : g**IntPart(m)*x**(-FracPart(m))*(g*x)**FracPart(m)*Int(x**m*(a + b*x**n)**p*(c + d*x**(-n))**q*(e + f*x**n)**r, x))
    rubi.add(rule383)

    def f384(b, a, e, f, r, x, m, c, d, n, p, q, g):
        return functools.reduce(operator.and_, [ FreeQ(List(a, b, c, d, e, f, g, m, n, p, q, r), x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern384 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e_ + x_**n_*WC('f', S(1)))**WC('r', S(1)), x_),CustomConstraint(f384))
    rule384 = ReplacementRule(pattern384, lambda b, a, e, f, r, x, m, c, d, n, p, q, g : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x))
    rubi.add(rule384)

    def f385(b, a, e, f, r, x, c, m, d, n, u, p, q, v):
        return functools.reduce(operator.and_, [ LinearPairQ(u, v, x), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e, x), FreeQ(f, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern385 = Pattern(Integral(u_**WC('m', S(1))*(e_ + v_**n_*WC('f', S(1)))**WC('r', S(1))*(v_**n_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(v_**n_*WC('d', S(1)) + WC('c', S(0)))**WC('q', S(1)), x_),CustomConstraint(f385))
    rule385 = ReplacementRule(pattern385, lambda b, a, e, f, r, x, c, m, d, n, u, p, q, v : u**m*v**(-m)*Subst(Int(x**m*(a + b*x**n)**p*(c + d*x**n)**q*(e + f*x**n)**r, x), x, v)/Coefficient(v, x, S(1)))
    rubi.add(rule385)

    def f386(b, a, r, f1, x, e2, n2, e1, m, c, d, n, q, p, f2, g):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + n2), ZeroQ(e1*f2 + e2*f1), Or(IntegerQ(r), And(PositiveQ(e1), PositiveQ(e2))), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e1, x), FreeQ(f1, x), FreeQ(e2, x), FreeQ(f2, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern386 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_),CustomConstraint(f386))
    rule386 = ReplacementRule(pattern386, lambda b, a, r, f1, x, e2, n2, e1, m, c, d, n, q, p, f2, g : Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule386)

    def f387(b, a, r, f1, x, e2, n2, e1, m, c, d, n, q, p, f2, g):
        return functools.reduce(operator.and_, [ ZeroQ(-n/S(2) + n2), ZeroQ(e1*f2 + e2*f1), FreeQ(a, x), FreeQ(b, x), FreeQ(c, x), FreeQ(d, x), FreeQ(e1, x), FreeQ(f1, x), FreeQ(e2, x), FreeQ(f2, x), FreeQ(g, x), FreeQ(m, x), FreeQ(n, x), FreeQ(p, x), FreeQ(q, x), FreeQ(r, x)])
    pattern387 = Pattern(Integral((x_*WC('g', S(1)))**WC('m', S(1))*(a_ + x_**n_*WC('b', S(1)))**WC('p', S(1))*(c_ + x_**n_*WC('d', S(1)))**WC('q', S(1))*(e1_ + x_**WC('n2', S(1))*WC('f1', S(1)))**WC('r', S(1))*(e2_ + x_**WC('n2', S(1))*WC('f2', S(1)))**WC('r', S(1)), x_),CustomConstraint(f387))
    rule387 = ReplacementRule(pattern387, lambda b, a, r, f1, x, e2, n2, e1, m, c, d, n, q, p, f2, g : (e1 + f1*x**(n/S(2)))**FracPart(r)*(e2 + f2*x**(n/S(2)))**FracPart(r)*(e1*e2 + f1*f2*x**n)**(-FracPart(r))*Int((g*x)**m*(a + b*x**n)**p*(c + d*x**n)**q*(e1*e2 + f1*f2*x**n)**r, x))
    rubi.add(rule387)

    return rubi
