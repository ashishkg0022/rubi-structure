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

    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False
    from sympy import And, Or

def quadratic_products(rubi):

    def cons_f1(c, a, b):
        return ZeroQ(-S(4)*a*c + b**S(2))

    cons1 = CustomConstraint(cons_f1)

    def cons_f2(a, x):
        return FreeQ(a, x)

    cons2 = CustomConstraint(cons_f2)

    def cons_f3(b, x):
        return FreeQ(b, x)

    cons3 = CustomConstraint(cons_f3)

    def cons_f4(c, x):
        return FreeQ(c, x)

    cons4 = CustomConstraint(cons_f4)
    pattern1 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons1)
    rule1 = ReplacementRule(pattern1, lambda c, a, x, b : (b/S(2) + c*x)*Int(S(1)/(b/S(2) + c*x), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule1)


    def cons_f5(p):
        return NonzeroQ(S(2)*p + S(1))

    cons5 = CustomConstraint(cons_f5)

    def cons_f6(p, x):
        return FreeQ(p, x)

    cons6 = CustomConstraint(cons_f6)
    pattern2 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons6, cons1, cons5)
    rule2 = ReplacementRule(pattern2, lambda x, c, p, a, b : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule2)


    def cons_f7(c, a, b):
        return NonzeroQ(-S(4)*a*c + b**S(2))

    cons7 = CustomConstraint(cons_f7)

    def cons_f8(p):
        return PositiveIntegerQ(p)

    cons8 = CustomConstraint(cons_f8)

    def cons_f9(c, a, b):
        return PerfectSquareQ(-S(4)*a*c + b**S(2))

    cons9 = CustomConstraint(cons_f9)
    pattern3 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons7, cons8, cons9, )
    def With3(x, c, p, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c**(-p)*Int(Simp(b/S(2) + c*x - q/S(2), x)**p*Simp(b/S(2) + c*x + q/S(2), x)**p, x)
    rule3 = ReplacementRule(pattern3, lambda x, c, p, a, b : With3(x, c, p, a, b))
    rubi.add(rule3)


    def cons_f10(c, a, b):
        return Not(PerfectSquareQ(-S(4)*a*c + b**S(2)))

    cons10 = CustomConstraint(cons_f10)
    pattern4 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons7, cons8, cons10)
    rule4 = ReplacementRule(pattern4, lambda x, c, p, a, b : Int(ExpandIntegrand((a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule4)


    def cons_f11(p):
        return RationalQ(p)

    cons11 = CustomConstraint(cons_f11)

    def cons_f12(p):
        return Greater(p, S(0))

    cons12 = CustomConstraint(cons_f12)

    def cons_f13(p):
        return IntegerQ(S(4)*p)

    cons13 = CustomConstraint(cons_f13)
    pattern5 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons7, cons11, cons12, cons13)
    rule5 = ReplacementRule(pattern5, lambda x, c, p, a, b : -p*(-S(4)*a*c + b**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*(S(2)*p + S(1))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), cons2, cons3, cons4, cons7)
    rule6 = ReplacementRule(pattern6, lambda c, a, x, b : -S(2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule6)


    def cons_f14(p):
        return Less(p, S(-1))

    cons14 = CustomConstraint(cons_f14)

    def cons_f15(p):
        return Unequal(p, S(-3)/2)

    cons15 = CustomConstraint(cons_f15)
    pattern7 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons7, cons11, cons14, cons15, cons13)
    rule7 = ReplacementRule(pattern7, lambda x, c, p, a, b : -S(2)*c*(S(2)*p + S(3))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule7)


    def cons_f16(c, a, b):
        return PosQ(-S(4)*a*c + b**S(2))

    cons16 = CustomConstraint(cons_f16)
    pattern8 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons7, cons16, cons9, )
    def With8(c, a, x, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c*Int(S(1)/Simp(b/S(2) + c*x - q/S(2), x), x)/q - c*Int(S(1)/Simp(b/S(2) + c*x + q/S(2), x), x)/q
    rule8 = ReplacementRule(pattern8, lambda c, a, x, b : With8(c, a, x, b))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons7, )
    def With9(c, a, x, b):
        q = -S(4)*a*c/b**S(2) + S(1)
        if And(RationalQ(q), Or(EqQ(q**S(2), S(1)), Not(RationalQ(-S(4)*a*c + b**S(2))))):
            return -S(2)*Subst(Int(S(1)/(q - x**S(2)), x), x, S(1) + S(2)*c*x/b)/b
        print("Unable to Integrate")
    rule9 = ReplacementRule(pattern9, lambda c, a, x, b : With9(c, a, x, b))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons7)
    rule10 = ReplacementRule(pattern10, lambda c, a, x, b : -S(2)*Subst(Int(S(1)/Simp(-S(4)*a*c + b**S(2) - x**S(2), x), x), x, b + S(2)*c*x))
    rubi.add(rule10)


    def cons_f17(c, a, b):
        return PositiveQ(S(4)*a - b**S(2)/c)

    cons17 = CustomConstraint(cons_f17)
    pattern11 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons6, cons17)
    rule11 = ReplacementRule(pattern11, lambda x, c, p, a, b : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p, x), x, b + S(2)*c*x)/(S(2)*c))
    rubi.add(rule11)


    def cons_f18(c, x, b):
        return FreeQ(List(b, c), x)

    cons18 = CustomConstraint(cons_f18)
    pattern12 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons4, cons18)
    rule12 = ReplacementRule(pattern12, lambda c, x, b : S(2)*Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(b*x + c*x**S(2))))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons7)
    rule13 = ReplacementRule(pattern13, lambda c, a, x, b : S(2)*Subst(Int(S(1)/(S(4)*c - x**S(2)), x), x, (b + S(2)*c*x)/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule13)


    def cons_f19(p):
        return LessEqual(S(3), Denominator(p), S(4))

    cons19 = CustomConstraint(cons_f19)
    pattern14 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons4, cons11, cons19)
    rule14 = ReplacementRule(pattern14, lambda p, c, x, b : (-c*(b*x + c*x**S(2))/b**S(2))**(-p)*(b*x + c*x**S(2))**p*Int((-c*x/b - c**S(2)*x**S(2)/b**S(2))**p, x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons7, cons11, )
    def With15(x, c, p, a, b):
        d = Denominator(p)
        if LessEqual(S(3), d, S(4)):
            return d*sqrt((b + S(2)*c*x)**S(2))*Subst(Int(x**(d*(p + S(1)) + S(-1))/sqrt(-S(4)*a*c + b**S(2) + S(4)*c*x**d), x), x, (a + b*x + c*x**S(2))**(S(1)/d))/(b + S(2)*c*x)
        print("Unable to Integrate")
    rule15 = ReplacementRule(pattern15, lambda x, c, p, a, b : With15(x, c, p, a, b))
    rubi.add(rule15)


    def cons_f20(p):
        return Not(IntegerQ(S(4)*p))

    cons20 = CustomConstraint(cons_f20)
    pattern16 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons6, cons7, cons20, )
    def With16(x, c, p, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -((-b - S(2)*c*x + q)/(S(2)*q))**(-p + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Hypergeometric2F1(-p, p + S(1), p + S(2), (b + S(2)*c*x + q)/(S(2)*q))/(q*(p + S(1)))
    rule16 = ReplacementRule(pattern16, lambda x, c, p, a, b : With16(x, c, p, a, b))
    rubi.add(rule16)


    def cons_f21(u, x):
        return LinearQ(u, x)

    cons21 = CustomConstraint(cons_f21)

    def cons_f22(u, x):
        return NonzeroQ(u - x)

    cons22 = CustomConstraint(cons_f22)
    pattern17 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons6, cons21, cons22)
    rule17 = ReplacementRule(pattern17, lambda x, c, u, p, a, b : Subst(Int((a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule17)


    def cons_f23(c, d, e, b):
        return ZeroQ(-b*e + S(2)*c*d)

    cons23 = CustomConstraint(cons_f23)

    def cons_f24(m):
        return IntegerQ(m/S(2) + S(1)/2)

    cons24 = CustomConstraint(cons_f24)

    def cons_f25(d, x):
        return FreeQ(d, x)

    cons25 = CustomConstraint(cons_f25)

    def cons_f26(e, x):
        return FreeQ(e, x)

    cons26 = CustomConstraint(cons_f26)

    def cons_f27(m, x):
        return FreeQ(m, x)

    cons27 = CustomConstraint(cons_f27)
    pattern18 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons23, cons24)
    rule18 = ReplacementRule(pattern18, lambda x, d, c, m, e, p, a, b : c**(-m/S(2) + S(-1)/2)*e**m*(a + b*x + c*x**S(2))**(m/S(2) + p + S(1)/2)/(m + S(2)*p + S(1)))
    rubi.add(rule18)


    def cons_f28(p, m):
        return ZeroQ(m + S(2)*p + S(1))

    cons28 = CustomConstraint(cons_f28)
    pattern19 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons23, cons28)
    rule19 = ReplacementRule(pattern19, lambda x, c, d, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*log(RemoveContent(d + e*x, x))/e)
    rubi.add(rule19)


    def cons_f29(p, m):
        return NonzeroQ(m + S(2)*p + S(1))

    cons29 = CustomConstraint(cons_f29)
    pattern20 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons23, cons29)
    rule20 = ReplacementRule(pattern20, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule20)


    def cons_f30(c, d, e, b):
        return NonzeroQ(-b*e + S(2)*c*d)

    cons30 = CustomConstraint(cons_f30)

    def cons_f31(p, m):
        return ZeroQ(m + S(2)*p + S(2))

    cons31 = CustomConstraint(cons_f31)

    def cons_f32(m):
        return NonzeroQ(m + S(1))

    cons32 = CustomConstraint(cons_f32)
    pattern21 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons30, cons31, cons32)
    rule21 = ReplacementRule(pattern21, lambda x, d, c, m, e, p, a, b : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0)))**S(2), x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30)
    rule22 = ReplacementRule(pattern22, lambda x, d, c, e, a, b : sqrt(a + b*x + c*x**S(2))*Int((b + S(2)*c*x)/(d + e*x)**S(2), x)/(b + S(2)*c*x))
    rubi.add(rule22)


    def cons_f33(m):
        return NonzeroQ(m + S(2))

    cons33 = CustomConstraint(cons_f33)
    pattern23 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons1, cons30, cons33)
    rule23 = ReplacementRule(pattern23, lambda x, d, c, m, e, a, b : (d + e*x)**(m + S(1))*sqrt(a + b*x + c*x**S(2))/(e*(m + S(2))) - (-b*e + S(2)*c*d)*sqrt(a + b*x + c*x**S(2))*Int((d + e*x)**m, x)/(e*(b + S(2)*c*x)*(m + S(2))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30)
    rule24 = ReplacementRule(pattern24, lambda x, d, c, e, a, b : -S(4)*c*e*sqrt(a + b*x + c*x**S(2))/((d + e*x)*(-b*e + S(2)*c*d)**S(2)) + S(2)*c*Int(S(1)/((d + e*x)*sqrt(a + b*x + c*x**S(2))), x)/(-b*e + S(2)*c*d))
    rubi.add(rule24)


    def cons_f34(p, m):
        return ZeroQ(m + S(2)*p + S(3))

    cons34 = CustomConstraint(cons_f34)
    pattern25 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons30, cons34, cons33)
    rule25 = ReplacementRule(pattern25, lambda x, d, c, m, e, p, a, b : -S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)**S(2)*(m*p + S(-1))) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(2))*(-b*e + S(2)*c*d)))
    rubi.add(rule25)


    def cons_f35(p):
        return NonzeroQ(p + S(3)/2)

    cons35 = CustomConstraint(cons_f35)
    pattern26 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons1, cons30, cons35)
    rule26 = ReplacementRule(pattern26, lambda x, d, c, e, p, a, b : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule26)


    def cons_f36(p, m):
        return RationalQ(m, p)

    cons36 = CustomConstraint(cons_f36)

    def cons_f37(p):
        return Greater(p, S(1))

    cons37 = CustomConstraint(cons_f37)

    def cons_f38(m):
        return Inequality(S(-2), LessEqual, m, Less, S(-1))

    cons38 = CustomConstraint(cons_f38)

    def cons_f39(p):
        return IntegerQ(S(2)*p)

    cons39 = CustomConstraint(cons_f39)
    pattern27 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30, cons36, cons37, cons38, cons39)
    rule27 = ReplacementRule(pattern27, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))))
    rubi.add(rule27)


    def cons_f40(m):
        return Less(m, S(-2))

    cons40 = CustomConstraint(cons_f40)

    def cons_f41(p, m):
        return Not(And(NegativeIntegerQ(m + S(2)*p + S(3)), Greater(m + S(3)*p + S(3), S(0))))

    cons41 = CustomConstraint(cons_f41)
    pattern28 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30, cons36, cons37, cons40, cons39, cons41)
    rule28 = ReplacementRule(pattern28, lambda x, d, c, m, e, p, a, b : S(2)*c*p*(S(2)*p + S(-1))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2))))
    rubi.add(rule28)


    def cons_f42(p, m):
        return NonzeroQ(m + S(2)*p)

    cons42 = CustomConstraint(cons_f42)

    def cons_f43(m):
        return Not(And(RationalQ(m), Less(m, S(-2))))

    cons43 = CustomConstraint(cons_f43)

    def cons_f44(p, m):
        return Not(And(IntegerQ(m), Less(S(0), m, S(2)*p)))

    cons44 = CustomConstraint(cons_f44)
    pattern29 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons1, cons30, cons11, cons12, cons42, cons29, cons41, cons43, cons44, cons39)
    rule29 = ReplacementRule(pattern29, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)**S(2)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))))
    rubi.add(rule29)


    def cons_f45(m):
        return Inequality(S(0), Less, m, LessEqual, S(1))

    cons45 = CustomConstraint(cons_f45)
    pattern30 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30, cons36, cons14, cons45, cons39)
    rule30 = ReplacementRule(pattern30, lambda x, d, c, m, e, p, a, b : e**S(2)*m*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule30)


    def cons_f46(m):
        return Greater(m, S(1))

    cons46 = CustomConstraint(cons_f46)
    pattern31 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons1, cons30, cons36, cons14, cons46, cons39)
    rule31 = ReplacementRule(pattern31, lambda x, d, c, m, e, p, a, b : e**S(2)*m*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) - e*m*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule31)


    def cons_f47(p, m):
        return NonzeroQ(m + p + S(1))

    cons47 = CustomConstraint(cons_f47)
    pattern32 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons1, cons30, cons36, cons14, cons47, cons39)
    rule32 = ReplacementRule(pattern32, lambda x, d, c, m, e, p, a, b : S(2)*c*e**S(2)*(m + S(2)*p + S(2))*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) - S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule32)


    def cons_f48(m):
        return RationalQ(m)

    cons48 = CustomConstraint(cons_f48)

    def cons_f49(m):
        return Greater(m, S(0))

    cons49 = CustomConstraint(cons_f49)

    def cons_f50(p, m):
        return Or(Not(RationalQ(p)), Inequality(S(-1), LessEqual, p, Less, S(0)), And(IntegerQ(m), Less(S(0), m, S(2)*p)), And(Equal(m, S(1)/2), Less(p, S(0))))

    cons50 = CustomConstraint(cons_f50)

    def cons_f51(p, m):
        return Or(IntegerQ(m), IntegerQ(S(2)*p))

    cons51 = CustomConstraint(cons_f51)
    pattern33 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons1, cons30, cons48, cons49, cons29, cons50, cons51)
    rule33 = ReplacementRule(pattern33, lambda x, d, c, m, e, p, a, b : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*(m + S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(m + S(2)*p + S(1))))
    rubi.add(rule33)


    def cons_f52(m):
        return Less(m, S(-1))

    cons52 = CustomConstraint(cons_f52)
    pattern34 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons1, cons30, cons48, cons52, cons39)
    rule34 = ReplacementRule(pattern34, lambda x, d, c, m, e, p, a, b : S(2)*c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-b*e + S(2)*c*d)) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule34)


    def cons_f53(p):
        return Not(IntegerQ(p))

    cons53 = CustomConstraint(cons_f53)
    pattern35 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons1, cons53, cons30)
    rule35 = ReplacementRule(pattern35, lambda x, d, c, m, e, p, a, b : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m, x))
    rubi.add(rule35)


    def cons_f54(c, d, e, a, b):
        return ZeroQ(a*e**S(2) - b*d*e + c*d**S(2))

    cons54 = CustomConstraint(cons_f54)

    def cons_f55(p):
        return IntegerQ(p)

    cons55 = CustomConstraint(cons_f55)
    pattern36 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons54, cons55)
    rule36 = ReplacementRule(pattern36, lambda x, c, d, m, e, p, a, b : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule36)


    def cons_f56(c, a, d, e):
        return ZeroQ(a*e**S(2) + c*d**S(2))

    cons56 = CustomConstraint(cons_f56)

    def cons_f57(p, d, a, m):
        return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), IntegerQ(m + p)))

    cons57 = CustomConstraint(cons_f57)
    pattern37 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons57)
    rule37 = ReplacementRule(pattern37, lambda x, c, d, m, e, p, a : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule37)


    def cons_f58(p, m):
        return ZeroQ(m + p)

    cons58 = CustomConstraint(cons_f58)
    pattern38 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons54, cons53, cons58)
    rule38 = ReplacementRule(pattern38, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons58)
    rule39 = ReplacementRule(pattern39, lambda x, c, d, m, e, p, a : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons54, cons53, cons31)
    rule40 = ReplacementRule(pattern40, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons31)
    rule41 = ReplacementRule(pattern41, lambda x, d, c, m, e, p, a : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(p + S(1))))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons54, cons53, cons11, cons14)
    rule42 = ReplacementRule(pattern42, lambda x, d, c, e, p, a, b : -e**S(2)*(p + S(2))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**S(2), x_), cons2, cons4, cons25, cons26, cons6, cons56, cons53, cons11, cons14)
    rule43 = ReplacementRule(pattern43, lambda x, c, d, e, p, a : -e**S(2)*(p + S(2))*Int((a + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)/(c*(p + S(1))))
    rubi.add(rule43)


    def cons_f59(m):
        return IntegerQ(m)

    cons59 = CustomConstraint(cons_f59)

    def cons_f60(p, m):
        return Or(Less(S(0), -m, p), Less(p, -m, S(0)))

    cons60 = CustomConstraint(cons_f60)

    def cons_f61(m):
        return Unequal(m, S(2))

    cons61 = CustomConstraint(cons_f61)

    def cons_f62(m):
        return Unequal(m, S(-1))

    cons62 = CustomConstraint(cons_f62)
    pattern44 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54, cons53, cons59, cons11, cons60, cons61, cons62)
    rule44 = ReplacementRule(pattern44, lambda x, c, d, m, e, p, a, b : Int((a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons59, cons11, cons60, cons61, cons62)
    rule45 = ReplacementRule(pattern45, lambda x, c, d, m, e, p, a : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m), x))
    rubi.add(rule45)


    def cons_f63(p, m):
        return PositiveIntegerQ(m + p)

    cons63 = CustomConstraint(cons_f63)
    pattern46 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons54, cons53, cons63)
    rule46 = ReplacementRule(pattern46, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons63)
    rule47 = ReplacementRule(pattern47, lambda x, d, c, m, e, p, a : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule47)


    def cons_f64(p, m):
        return NegativeIntegerQ(m + S(2)*p + S(2))

    cons64 = CustomConstraint(cons_f64)
    pattern48 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons54, cons53, cons64)
    rule48 = ReplacementRule(pattern48, lambda x, d, c, m, e, p, a, b : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons64)
    rule49 = ReplacementRule(pattern49, lambda x, c, d, m, e, p, a : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54)
    rule50 = ReplacementRule(pattern50, lambda x, d, c, e, a, b : S(2)*e*Subst(Int(S(1)/(-b*e + S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons56)
    rule51 = ReplacementRule(pattern51, lambda x, c, d, e, a : S(2)*e*Subst(Int(S(1)/(S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule51)


    def cons_f65(p, m):
        return Or(Less(m, S(-2)), ZeroQ(m + S(2)*p + S(1)))

    cons65 = CustomConstraint(cons_f65)
    pattern52 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54, cons36, cons12, cons65, cons47, cons39)
    rule52 = ReplacementRule(pattern52, lambda x, d, c, m, e, p, a, b : -c*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + p + S(1))))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons56, cons36, cons12, cons65, cons47, cons39)
    rule53 = ReplacementRule(pattern53, lambda x, c, d, m, e, p, a : -c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/(e**S(2)*(m + p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + p + S(1))))
    rubi.add(rule53)


    def cons_f66(p, m):
        return Or(Inequality(S(-2), LessEqual, m, Less, S(0)), Equal(m + p + S(1), S(0)))

    cons66 = CustomConstraint(cons_f66)
    pattern54 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54, cons36, cons12, cons66, cons29, cons39)
    rule54 = ReplacementRule(pattern54, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule54)

    pattern55 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons56, cons36, cons12, cons66, cons29, cons39)
    rule55 = ReplacementRule(pattern55, lambda x, c, d, m, e, p, a : -S(2)*c*d*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e**S(2)*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule55)

    pattern56 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54, cons36, cons14, cons45, cons39)
    rule56 = ReplacementRule(pattern56, lambda x, d, c, m, e, p, a, b : -(-b*e + S(2)*c*d)*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule56)

    pattern57 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons56, cons36, cons14, cons45, cons39)
    rule57 = ReplacementRule(pattern57, lambda x, d, c, m, e, p, a : d*(m + S(2)*p + S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - d*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*e*(p + S(1))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons54, cons36, cons14, cons46, cons39)
    rule58 = ReplacementRule(pattern58, lambda x, d, c, m, e, p, a, b : -e**S(2)*(m + p)*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons56, cons36, cons14, cons46, cons39)
    rule59 = ReplacementRule(pattern59, lambda x, c, d, m, e, p, a : -e**S(2)*(m + p)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule59)


    def cons_f67(m):
        return GreaterEqual(m, S(1))

    cons67 = CustomConstraint(cons_f67)
    pattern60 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons54, cons48, cons67, cons29, cons39)
    rule60 = ReplacementRule(pattern60, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule60)

    pattern61 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons6, cons56, cons48, cons67, cons29, cons39)
    rule61 = ReplacementRule(pattern61, lambda x, d, c, m, e, p, a : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule61)


    def cons_f68(m):
        return Less(m, S(0))

    cons68 = CustomConstraint(cons_f68)
    pattern62 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons54, cons48, cons68, cons47, cons39)
    rule62 = ReplacementRule(pattern62, lambda x, d, c, m, e, p, a, b : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons6, cons56, cons48, cons68, cons47, cons39)
    rule63 = ReplacementRule(pattern63, lambda x, c, d, m, e, p, a : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons4, cons26, cons27, cons53)
    rule64 = ReplacementRule(pattern64, lambda x, c, m, e, p, b : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p, x))
    rubi.add(rule64)


    def cons_f69(a):
        return PositiveQ(a)

    cons69 = CustomConstraint(cons_f69)

    def cons_f70(d):
        return PositiveQ(d)

    cons70 = CustomConstraint(cons_f70)
    pattern65 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons56, cons53, cons69, cons70)
    rule65 = ReplacementRule(pattern65, lambda x, d, c, m, e, p, a : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons54, cons53)
    rule66 = ReplacementRule(pattern66, lambda x, c, d, m, e, p, a, b : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons56, cons53)
    rule67 = ReplacementRule(pattern67, lambda x, d, c, m, e, p, a : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23)
    rule68 = ReplacementRule(pattern68, lambda x, c, d, e, a, b : b**S(2)*Int((d + e*x)/(a + b*x + c*x**S(2)), x)/(d**S(2)*(-S(4)*a*c + b**S(2))) - S(4)*b*c*Int(S(1)/(b + S(2)*c*x), x)/(d*(-S(4)*a*c + b**S(2))))
    rubi.add(rule68)


    def cons_f71(p):
        return NonzeroQ(p + S(1))

    cons71 = CustomConstraint(cons_f71)
    pattern69 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons23, cons34, cons71)
    rule69 = ReplacementRule(pattern69, lambda x, c, d, m, e, p, a, b : S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule69)


    def cons_f72(p, m):
        return Not(And(ZeroQ(m + S(-3)), Unequal(p, S(1))))

    cons72 = CustomConstraint(cons_f72)
    pattern70 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons23, cons8, cons72)
    rule70 = ReplacementRule(pattern70, lambda x, d, c, m, e, p, a, b : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule70)


    def cons_f73(p, m):
        return NonzeroQ(m + S(2)*p + S(3))

    cons73 = CustomConstraint(cons_f73)

    def cons_f74(p, m):
        return Not(And(EvenQ(m), Less(m + S(2)*p + S(3), S(0))))

    cons74 = CustomConstraint(cons_f74)
    pattern71 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23, cons73, cons36, cons12, cons52, cons74, cons39)
    rule71 = ReplacementRule(pattern71, lambda x, c, d, m, e, p, a, b : -b*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(d*e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule71)


    def cons_f75(m):
        return Not(And(RationalQ(m), Less(m, S(-1))))

    cons75 = CustomConstraint(cons_f75)

    def cons_f76(p, m):
        return Not(And(PositiveIntegerQ(m/S(2) + S(-1)/2), Or(Not(IntegerQ(p)), Less(m, S(2)*p))))

    cons76 = CustomConstraint(cons_f76)
    pattern72 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons23, cons73, cons11, cons12, cons75, cons76, cons48, cons39)
    rule72 = ReplacementRule(pattern72, lambda x, c, d, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - d*p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(b*e*(m + S(2)*p + S(1))))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23, cons73, cons36, cons14, cons46, cons39)
    rule73 = ReplacementRule(pattern73, lambda x, c, d, m, e, p, a, b : -d*e*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))) + d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))))
    rubi.add(rule73)


    def cons_f77(m):
        return Not(And(RationalQ(m), Greater(m, S(1))))

    cons77 = CustomConstraint(cons_f77)
    pattern74 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons23, cons73, cons11, cons14, cons77, cons48, cons39)
    rule74 = ReplacementRule(pattern74, lambda x, c, d, m, e, p, a, b : -S(2)*c*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule74)

    pattern75 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23)
    rule75 = ReplacementRule(pattern75, lambda x, c, d, e, a, b : S(4)*c*Subst(Int(S(1)/(-S(4)*a*c*e + b**S(2)*e + S(4)*c*e*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule75)


    def cons_f78(c, a, b):
        return NegativeQ(c/(-S(4)*a*c + b**S(2)))

    cons78 = CustomConstraint(cons_f78)
    pattern76 = Pattern(Integral(S(1)/(sqrt(d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23, cons78)
    rule76 = ReplacementRule(pattern76, lambda x, c, d, e, a, b : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(S(1)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule76)

    pattern77 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23, cons78)
    rule77 = ReplacementRule(pattern77, lambda x, c, d, e, a, b : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(x**S(2)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule77)


    def cons_f79(m):
        return EqQ(m**S(2), S(1)/4)

    cons79 = CustomConstraint(cons_f79)
    pattern78 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons23, cons79)
    rule78 = ReplacementRule(pattern78, lambda x, c, d, m, e, a, b : sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*Int((d + e*x)**m/sqrt(-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2))), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule78)


    def cons_f80(p, m):
        return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), OddQ(m))

    cons80 = CustomConstraint(cons_f80)
    pattern79 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons23, cons73, cons48, cons46, cons29, cons80)
    rule79 = ReplacementRule(pattern79, lambda x, c, d, m, e, p, a, b : S(2)*d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(2)*p + S(1))) + d**S(2)*(m + S(-1))*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p, x)/(b**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule79)


    def cons_f81(p, m):
        return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), IntegerQ(m/S(2) + p + S(3)/2))

    cons81 = CustomConstraint(cons_f81)
    pattern80 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons23, cons73, cons48, cons52, cons81)
    rule80 = ReplacementRule(pattern80, lambda x, c, d, m, e, p, a, b : b**S(2)*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/(d**S(2)*(m + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*b*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(d*(m + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons23)
    rule81 = ReplacementRule(pattern81, lambda x, c, d, m, e, p, a, b : Subst(Int(x**m*(a - b**S(2)/(S(4)*c) + c*x**S(2)/e**S(2))**p, x), x, d + e*x)/e)
    rubi.add(rule81)


    def cons_f82(c, d, e, a, b):
        return NonzeroQ(a*e**S(2) - b*d*e + c*d**S(2))

    cons82 = CustomConstraint(cons_f82)
    pattern82 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons8)
    rule82 = ReplacementRule(pattern82, lambda x, d, c, m, e, p, a, b : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule82)


    def cons_f83(c, a, d, e):
        return NonzeroQ(a*e**S(2) + c*d**S(2))

    cons83 = CustomConstraint(cons_f83)

    def cons_f84(p, m):
        return Not(And(ZeroQ(m + S(-1)), Greater(p, S(1))))

    cons84 = CustomConstraint(cons_f84)
    pattern83 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons83, cons8, cons84)
    rule83 = ReplacementRule(pattern83, lambda x, c, d, m, e, p, a : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m, x), x))
    rubi.add(rule83)


    def cons_f85(c, a, b):
        return NiceSqrtQ(-S(4)*a*c + b**S(2))

    cons85 = CustomConstraint(cons_f85)
    pattern84 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons85, )
    def With84(x, d, c, e, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (c*d - e*(b/S(2) - q/S(2)))*Int(S(1)/(b/S(2) + c*x - q/S(2)), x)/q - (c*d - e*(b/S(2) + q/S(2)))*Int(S(1)/(b/S(2) + c*x + q/S(2)), x)/q
    rule84 = ReplacementRule(pattern84, lambda x, d, c, e, a, b : With84(x, d, c, e, a, b))
    rubi.add(rule84)


    def cons_f86(c, a):
        return NiceSqrtQ(-a*c)

    cons86 = CustomConstraint(cons_f86)
    pattern85 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons86, )
    def With85(x, c, d, e, a):
        q = Rt(-a*c, S(2))
        return (-c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x + q), x) + (c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x - q), x)
    rule85 = ReplacementRule(pattern85, lambda x, c, d, e, a : With85(x, c, d, e, a))
    rubi.add(rule85)


    def cons_f87(c, a, b):
        return Not(NiceSqrtQ(-S(4)*a*c + b**S(2)))

    cons87 = CustomConstraint(cons_f87)
    pattern86 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons87)
    rule86 = ReplacementRule(pattern86, lambda x, d, c, e, a, b : e*Int((b + S(2)*c*x)/(a + b*x + c*x**S(2)), x)/(S(2)*c) + (-b*e + S(2)*c*d)*Int(S(1)/(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule86)


    def cons_f88(c, a):
        return Not(NiceSqrtQ(-a*c))

    cons88 = CustomConstraint(cons_f88)
    pattern87 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons88)
    rule87 = ReplacementRule(pattern87, lambda x, c, d, e, a : d*Int(S(1)/(a + c*x**S(2)), x) + e*Int(x/(a + c*x**S(2)), x))
    rubi.add(rule87)

    pattern88 = Pattern(Integral(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30)
    rule88 = ReplacementRule(pattern88, lambda x, d, c, e, a, b : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83)
    rule89 = ReplacementRule(pattern89, lambda x, c, d, e, a : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule89)


    def cons_f89(d, m):
        return Or(NonzeroQ(d), Greater(m, S(2)))

    cons89 = CustomConstraint(cons_f89)
    pattern90 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons59, cons46, cons89)
    rule90 = ReplacementRule(pattern90, lambda x, d, c, m, e, a, b : Int(PolynomialDivide((d + e*x)**m, a + b*x + c*x**S(2), x), x))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons59, cons46, cons89)
    rule91 = ReplacementRule(pattern91, lambda x, c, d, m, e, a : Int(PolynomialDivide((d + e*x)**m, a + c*x**S(2), x), x))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons48, cons46)
    rule92 = ReplacementRule(pattern92, lambda x, d, c, m, e, a, b : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + e*x*(-b*e + S(2)*c*d), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule92)

    pattern93 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons48, cons46)
    rule93 = ReplacementRule(pattern93, lambda x, c, d, m, e, a : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + S(2)*c*d*e*x, x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule93)

    pattern94 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30)
    rule94 = ReplacementRule(pattern94, lambda x, d, c, e, a, b : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((-b*e + c*d - c*e*x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons83)
    rule95 = ReplacementRule(pattern95, lambda x, c, d, e, a : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((c*d - c*e*x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule95)

    pattern96 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30)
    rule96 = ReplacementRule(pattern96, lambda x, d, c, e, a, b : S(2)*e*Subst(Int(S(1)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule96)

    pattern97 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons83)
    rule97 = ReplacementRule(pattern97, lambda x, c, d, e, a : S(2)*e*Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons48, cons52)
    rule98 = ReplacementRule(pattern98, lambda x, d, c, m, e, a, b : e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(-b*e + c*d - c*e*x, x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons27, cons83, cons48, cons52)
    rule99 = ReplacementRule(pattern99, lambda x, c, d, m, e, a : c*Int((d - e*x)*(d + e*x)**(m + S(1))/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)) + e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule99)


    def cons_f90(m):
        return Not(IntegerQ(m))

    cons90 = CustomConstraint(cons_f90)
    pattern100 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons90)
    rule100 = ReplacementRule(pattern100, lambda x, d, c, m, e, a, b : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons27, cons83, cons90)
    rule101 = ReplacementRule(pattern101, lambda x, c, d, m, e, a : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30)
    rule102 = ReplacementRule(pattern102, lambda x, d, c, e, a, b : -S(2)*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2), x_), cons2, cons4, cons25, cons26, cons83)
    rule103 = ReplacementRule(pattern103, lambda x, c, d, e, a : (-a*e + c*d*x)/(a*c*sqrt(a + c*x**S(2))))
    rubi.add(rule103)

    pattern104 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons11, cons14, cons15)
    rule104 = ReplacementRule(pattern104, lambda x, d, c, e, p, a, b : -(S(2)*p + S(3))*(-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons11, cons14, cons15)
    rule105 = ReplacementRule(pattern105, lambda x, c, d, e, p, a : d*(S(2)*p + S(3))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule105)


    def cons_f91(p):
        return Not(And(RationalQ(p), LessEqual(p, S(-1))))

    cons91 = CustomConstraint(cons_f91)
    pattern106 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons82, cons30, cons91)
    rule106 = ReplacementRule(pattern106, lambda x, d, c, e, p, a, b : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons6, cons83, cons91)
    rule107 = ReplacementRule(pattern107, lambda x, c, d, e, p, a : d*Int((a + c*x**S(2))**p, x) + e*(a + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule107)


    def cons_f92(d, a, e, b):
        return ZeroQ(a*e + b*d)

    cons92 = CustomConstraint(cons_f92)

    def cons_f93(c, d, e, b):
        return ZeroQ(b*e + c*d)

    cons93 = CustomConstraint(cons_f93)

    def cons_f94(p, m):
        return PositiveIntegerQ(m - p + S(1))

    cons94 = CustomConstraint(cons_f94)
    pattern108 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons92, cons93, cons94, cons53)
    rule108 = ReplacementRule(pattern108, lambda x, d, c, m, e, p, a, b : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m - p)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule108)


    def cons_f95(c, d, e, b):
        return NonzeroQ(-b*e + c*d)

    cons95 = CustomConstraint(cons_f95)

    def cons_f96(m):
        return Equal(m**S(2), S(1)/4)

    cons96 = CustomConstraint(cons_f96)

    def cons_f97(c):
        return NegativeQ(c)

    cons97 = CustomConstraint(cons_f97)

    def cons_f98(b):
        return RationalQ(b)

    cons98 = CustomConstraint(cons_f98)
    pattern109 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons4, cons25, cons26, cons95, cons30, cons48, cons96, cons97, cons98)
    rule109 = ReplacementRule(pattern109, lambda x, d, c, m, e, b : Int((d + e*x)**m/(sqrt(b*x)*sqrt(S(1) + c*x/b)), x))
    rubi.add(rule109)

    pattern110 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons3, cons4, cons25, cons26, cons95, cons30, cons48, cons96)
    rule110 = ReplacementRule(pattern110, lambda x, d, c, m, e, b : sqrt(x)*sqrt(b + c*x)*Int((d + e*x)**m/(sqrt(x)*sqrt(b + c*x)), x)/sqrt(b*x + c*x**S(2)))
    rubi.add(rule110)


    def cons_f99(m):
        return ZeroQ(m**S(2) + S(-1)/4)

    cons99 = CustomConstraint(cons_f99)
    pattern111 = Pattern(Integral(x_**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons7, cons99)
    rule111 = ReplacementRule(pattern111, lambda x, c, m, a, b : S(2)*Subst(Int(x**(S(2)*m + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((e_*x_)**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons26, cons7, cons99)
    rule112 = ReplacementRule(pattern112, lambda x, c, m, e, a, b : x**(-m)*(e*x)**m*Int(x**m/sqrt(a + b*x + c*x**S(2)), x))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons99)
    rule113 = ReplacementRule(pattern113, lambda x, d, c, m, e, a, b : S(2)*sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*(S(2)*c*(d + e*x)/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))))**(-m)*(d + e*x)**m*Rt(-S(4)*a*c + b**S(2), S(2))*Subst(Int((S(2)*e*x**S(2)*Rt(-S(4)*a*c + b**S(2), S(2))/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(S(2))*sqrt((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))/Rt(-S(4)*a*c + b**S(2), S(2)))/S(2))/(c*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons83, cons99)
    rule114 = ReplacementRule(pattern114, lambda x, c, d, m, e, a : S(2)*a*(c*(d + e*x)/(-a*e*Rt(-c/a, S(2)) + c*d))**(-m)*sqrt(S(1) + c*x**S(2)/a)*(d + e*x)**m*Rt(-c/a, S(2))*Subst(Int((S(2)*a*e*x**S(2)*Rt(-c/a, S(2))/(-a*e*Rt(-c/a, S(2)) + c*d) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(-x*Rt(-c/a, S(2))/S(2) + S(1)/2))/(c*sqrt(a + c*x**S(2))))
    rubi.add(rule114)


    def cons_f100(p, m):
        return Equal(m + S(2)*p + S(2), S(0))

    cons100 = CustomConstraint(cons_f100)
    pattern115 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons36, cons100, cons12)
    rule115 = ReplacementRule(pattern115, lambda x, d, c, m, e, p, a, b : p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons83, cons36, cons100, cons12)
    rule116 = ReplacementRule(pattern116, lambda x, c, d, m, e, p, a : -S(2)*a*c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*a*e + S(2)*c*d*x)/(S(2)*(m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons36, cons100, cons14)
    rule117 = ReplacementRule(pattern117, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*(S(2)*p + S(3))*(a*e**S(2) - b*d*e + c*d**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons83, cons36, cons100, cons14)
    rule118 = ReplacementRule(pattern118, lambda x, c, d, m, e, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) + (S(2)*p + S(3))*(a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule118)

    pattern119 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons30)
    rule119 = ReplacementRule(pattern119, lambda x, d, c, e, a, b : -S(2)*Subst(Int(S(1)/(S(4)*a*e**S(2) - S(4)*b*d*e + S(4)*c*d**S(2) - x**S(2)), x), x, (S(2)*a*e - b*d - x*(-b*e + S(2)*c*d))/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule119)


    def cons_f101(x, c, d, e, a):
        return FreeQ(List(a, c, d, e), x)

    cons101 = CustomConstraint(cons_f101)
    pattern120 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons101)
    rule120 = ReplacementRule(pattern120, lambda x, c, d, e, a : -Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - x**S(2)), x), x, (a*e - c*d*x)/sqrt(a + c*x**S(2))))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons53, cons31)
    rule121 = ReplacementRule(pattern121, lambda x, d, c, m, e, p, a, b : -((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))**(-p)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), -S(4)*c*(d + e*x)*Rt(-S(4)*a*c + b**S(2), S(2))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))/((m + S(1))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons53, cons31)
    rule122 = ReplacementRule(pattern122, lambda x, d, c, m, e, p, a : ((c*d + e*Rt(-a*c, S(2)))*(c*x + Rt(-a*c, S(2)))/((c*d - e*Rt(-a*c, S(2)))*(c*x - Rt(-a*c, S(2)))))**(-p)*(a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*x + Rt(-a*c, S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), S(2)*c*(d + e*x)*Rt(-a*c, S(2))/((c*d - e*Rt(-a*c, S(2)))*(-c*x + Rt(-a*c, S(2)))))/((m + S(1))*(c*d + e*Rt(-a*c, S(2)))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons34, cons11, cons14)
    rule123 = ReplacementRule(pattern123, lambda x, d, c, m, e, p, a, b : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons34, cons11, cons14)
    rule124 = ReplacementRule(pattern124, lambda x, c, d, m, e, p, a : -d*m*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons34)
    rule125 = ReplacementRule(pattern125, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + (-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons34)
    rule126 = ReplacementRule(pattern126, lambda x, c, d, m, e, p, a : c*d*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule126)


    def cons_f102(p, m):
        return Or(IntegerQ(p), And(RationalQ(m), Less(m, S(-1))))

    cons102 = CustomConstraint(cons_f102)

    def cons_f103(p, m):
        return Not(NegativeIntegerQ(m + S(2)*p + S(1)))

    cons103 = CustomConstraint(cons_f103)

    def cons_f104(x, c, d, m, e, p, a, b):
        return IntQuadraticQ(a, b, c, d, e, m, p, x)

    cons104 = CustomConstraint(cons_f104)
    pattern127 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons11, cons12, cons102, cons32, cons103, cons104)
    rule127 = ReplacementRule(pattern127, lambda x, d, c, m, e, p, a, b : -p*Int((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule127)


    def cons_f105(x, c, d, m, e, p, a):
        return IntQuadraticQ(a, S(0), c, d, e, m, p, x)

    cons105 = CustomConstraint(cons_f105)
    pattern128 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons83, cons11, cons12, cons102, cons32, cons103, cons105)
    rule128 = ReplacementRule(pattern128, lambda x, c, d, m, e, p, a : -S(2)*c*p*Int(x*(a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e*(m + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule128)


    def cons_f106(m):
        return Or(Not(RationalQ(m)), Less(m, S(1)))

    cons106 = CustomConstraint(cons_f106)

    def cons_f107(p, m):
        return Not(NegativeIntegerQ(m + S(2)*p))

    cons107 = CustomConstraint(cons_f107)
    pattern129 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons11, cons12, cons29, cons106, cons107, cons104)
    rule129 = ReplacementRule(pattern129, lambda x, d, c, m, e, p, a, b : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d), x), x)/(e*(m + S(2)*p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons83, cons11, cons12, cons29, cons106, cons107, cons105)
    rule130 = ReplacementRule(pattern130, lambda x, c, d, m, e, p, a : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*e - c*d*x, x), x)/(e*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule130)


    def cons_f108(p, m):
        return Or(Less(m, S(1)), And(NegativeIntegerQ(m + S(2)*p + S(3)), Unequal(m, S(2))))

    cons108 = CustomConstraint(cons_f108)
    pattern131 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons36, cons14, cons49, cons108, cons104)
    rule131 = ReplacementRule(pattern131, lambda x, d, c, m, e, p, a, b : (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(b*e*m + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(m + S(2)*p + S(3))), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule131)

    pattern132 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons83, cons36, cons14, cons49, cons108, cons105)
    rule132 = ReplacementRule(pattern132, lambda x, c, d, m, e, p, a : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(d*(S(2)*p + S(3)) + e*x*(m + S(2)*p + S(3))), x)/(S(2)*a*(p + S(1))))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons7, cons82, cons30, cons36, cons14, cons46, cons104)
    rule133 = ReplacementRule(pattern133, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*c*d**S(2)*(S(2)*p + S(3)) + e*x*(b*e - S(2)*c*d)*(m + S(2)*p + S(2)) + e*(S(2)*a*e*(m + S(-1)) + b*d*(-m + S(2)*p + S(4))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons83, cons36, cons14, cons46, cons105)
    rule134 = ReplacementRule(pattern134, lambda x, c, d, m, e, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(a*e**S(2)*(m + S(-1)) - c*d**S(2)*(S(2)*p + S(3)) - c*d*e*x*(m + S(2)*p + S(2)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons7, cons82, cons30, cons11, cons14, cons104)
    rule135 = ReplacementRule(pattern135, lambda x, d, c, m, e, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*e - b**S(2)*e + b*c*d + c*x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3)) - c*e*x*(-b*e + S(2)*c*d)*(m + S(2)*p + S(4)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule135)

    pattern136 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons83, cons11, cons14, cons105)
    rule136 = ReplacementRule(pattern136, lambda x, c, d, m, e, p, a : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(a*e + c*d*x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(a*e**S(2)*(m + S(2)*p + S(3)) + c*d**S(2)*(S(2)*p + S(3)) + c*d*e*x*(m + S(2)*p + S(4)), x), x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule136)


    def cons_f109(m):
        return If(RationalQ(m), Greater(m, S(1)), SumSimplerQ(m, S(-2)))

    cons109 = CustomConstraint(cons_f109)
    pattern137 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons109, cons29, cons104)
    rule137 = ReplacementRule(pattern137, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p*Simp(c*d**S(2)*(m + S(2)*p + S(1)) + e*x*(m + p)*(-b*e + S(2)*c*d) - e*(a*e*(m + S(-1)) + b*d*(p + S(1))), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule137)

    pattern138 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons109, cons29, cons105)
    rule138 = ReplacementRule(pattern138, lambda x, c, d, m, e, p, a : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-2))*Simp(-a*e**S(2)*(m + S(-1)) + c*d**S(2)*(m + S(2)*p + S(1)) + S(2)*c*d*e*x*(m + p), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule138)


    def cons_f110(x, d, c, m, e, p, a, b):
        return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, b, c, d, e, m, p, x)), And(SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

    cons110 = CustomConstraint(cons_f110)
    pattern139 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons110)
    rule139 = ReplacementRule(pattern139, lambda x, d, c, m, e, p, a, b : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*e*(m + p + S(2)) + c*d*(m + S(1)) - c*e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule139)


    def cons_f111(x, d, c, m, e, p, a):
        return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, S(0), c, d, e, m, p, x)), And(SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

    cons111 = CustomConstraint(cons_f111)
    pattern140 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons111)
    rule140 = ReplacementRule(pattern140, lambda x, c, d, m, e, p, a : c*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(d*(m + S(1)) - e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule140)


    def cons_f112(c, d, e, a, b):
        return ZeroQ(-S(3)*a*c*e**S(2) + b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

    cons112 = CustomConstraint(cons_f112)

    def cons_f113(c, d, e, b):
        return PosQ(c*e**S(2)*(-b*e + S(2)*c*d))

    cons113 = CustomConstraint(cons_f113)
    pattern141 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons4, cons25, cons26, cons30, cons112, cons113, )
    def With141(x, d, c, e, a, b):
        q = Rt(S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) + S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x - q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule141 = ReplacementRule(pattern141, lambda x, d, c, e, a, b : With141(x, d, c, e, a, b))
    rubi.add(rule141)


    def cons_f114(c, a, d, e):
        return ZeroQ(-S(3)*a*e**S(2) + c*d**S(2))

    cons114 = CustomConstraint(cons_f114)
    pattern142 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons114, )
    def With142(x, c, d, e, a):
        q = Rt(S(6)*c**S(2)*e**S(2)/d**S(2), S(3))
        return -sqrt(S(3))*c*e*ArcTan(S(2)*sqrt(S(3))*c*(d - e*x)/(S(3)*d*q*(a + c*x**S(2))**(S(1)/3)) + sqrt(S(3))/S(3))/(d**S(2)*q**S(2)) - S(3)*c*e*log(d + e*x)/(S(2)*d**S(2)*q**S(2)) + S(3)*c*e*log(c*d - c*e*x - d*q*(a + c*x**S(2))**(S(1)/3))/(S(2)*d**S(2)*q**S(2))
    rule142 = ReplacementRule(pattern142, lambda x, c, d, e, a : With142(x, c, d, e, a))
    rubi.add(rule142)


    def cons_f115(c, d, e, b):
        return NegQ(c*e**S(2)*(-b*e + S(2)*c*d))

    cons115 = CustomConstraint(cons_f115)
    pattern143 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons4, cons25, cons26, cons30, cons112, cons115, )
    def With143(x, d, c, e, a, b):
        q = Rt(-S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) - S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x + q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule143 = ReplacementRule(pattern143, lambda x, d, c, e, a, b : With143(x, d, c, e, a, b))
    rubi.add(rule143)


    def cons_f116(c, d, e, a, b):
        return ZeroQ(S(9)*a*c*e**S(2) - S(2)*b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

    cons116 = CustomConstraint(cons_f116)
    pattern144 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons2, cons3, cons4, cons25, cons26, cons7, cons116, )
    def With144(x, d, c, e, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)*Int(S(1)/((d + e*x)*(b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule144 = ReplacementRule(pattern144, lambda x, d, c, e, a, b : With144(x, d, c, e, a, b))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons83)
    rule145 = ReplacementRule(pattern145, lambda x, c, d, e, a : d*Int(S(1)/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/4)*(d_ + x_*WC('e', S(1)))), x_), cons2, cons4, cons25, cons26, cons83)
    rule146 = ReplacementRule(pattern146, lambda x, c, d, e, a : d*Int(S(1)/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons6, cons17, cons13)
    rule147 = ReplacementRule(pattern147, lambda x, d, c, e, p, a, b : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p/Simp(-b*e + S(2)*c*d + e*x, x), x), x, b + S(2)*c*x))
    rubi.add(rule147)


    def cons_f117(c, a, b):
        return Not(PositiveQ(S(4)*a - b**S(2)/c))

    cons117 = CustomConstraint(cons_f117)
    pattern148 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons6, cons117, cons13)
    rule148 = ReplacementRule(pattern148, lambda x, d, c, e, p, a, b : (-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))**(-p)*(a + b*x + c*x**S(2))**p*Int((-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2)))**p/(d + e*x), x))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons53, cons69, cons97)
    rule149 = ReplacementRule(pattern149, lambda x, d, c, m, e, p, a : Int((d + e*x)**m*(-x*Rt(-c, S(2)) + Rt(a, S(2)))**p*(x*Rt(-c, S(2)) + Rt(a, S(2)))**p, x))
    rubi.add(rule149)


    def cons_f118(m):
        return NegativeIntegerQ(m)

    cons118 = CustomConstraint(cons_f118)
    pattern150 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons82, cons30, cons53, cons118, )
    def With150(x, d, c, m, e, p, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -(e*(b + S(2)*c*x - q)/(S(2)*c*(d + e*x)))**(-p)*(e*(b + S(2)*c*x + q)/(S(2)*c*(d + e*x)))**(-p)*(a + b*x + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x*(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule150 = ReplacementRule(pattern150, lambda x, d, c, m, e, p, a, b : With150(x, d, c, m, e, p, a, b))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons6, cons83, cons53, cons118, )
    def With151(x, d, c, m, e, p, a):
        q = Rt(-a*c, S(2))
        return -(e*(c*x + q)/(c*(d + e*x)))**(-p)*(-e*(-c*x + q)/(c*(d + e*x)))**(-p)*(a + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*q/c) + S(1), x)**p*Simp(-x*(d + e*q/c) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule151 = ReplacementRule(pattern151, lambda x, d, c, m, e, p, a : With151(x, d, c, m, e, p, a))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons7, cons82, cons30, cons53, )
    def With152(x, d, c, m, e, p, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (-(d + e*x)/(d - e*(b - q)/(S(2)*c)) + S(1))**(-p)*(-(d + e*x)/(d - e*(b + q)/(S(2)*c)) + S(1))**(-p)*(a + b*x + c*x**S(2))**p*Subst(Int(x**m*Simp(-x/(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x/(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, d + e*x)/e
    rule152 = ReplacementRule(pattern152, lambda x, d, c, m, e, p, a, b : With152(x, d, c, m, e, p, a, b))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons83, cons53, )
    def With153(x, d, c, m, e, p, a):
        q = Rt(-a*c, S(2))
        return (a + c*x**S(2))**p*(-(d + e*x)/(d - e*q/c) + S(1))**(-p)*(-(d + e*x)/(d + e*q/c) + S(1))**(-p)*Subst(Int(x**m*Simp(-x/(d - e*q/c) + S(1), x)**p*Simp(-x/(d + e*q/c) + S(1), x)**p, x), x, d + e*x)/e
    rule153 = ReplacementRule(pattern153, lambda x, d, c, m, e, p, a : With153(x, d, c, m, e, p, a))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons27, cons6, cons21, cons22)
    rule154 = ReplacementRule(pattern154, lambda x, d, c, m, e, u, p, a, b : Subst(Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons27, cons6, cons21, cons22)
    rule155 = ReplacementRule(pattern155, lambda x, d, c, m, e, u, p, a : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule155)


    def cons_f119(n):
        return IntegerQ(n)

    cons119 = CustomConstraint(cons_f119)

    def cons_f120(p):
        return Not(IntegerQ(S(2)*p))

    cons120 = CustomConstraint(cons_f120)
    pattern156 = Pattern(Integral(x_**WC('n', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons6, cons119, cons120)
    rule156 = ReplacementRule(pattern156, lambda x, d, c, e, p, a, n : d*Int(x**n*(a + c*x**S(2))**p, x) + e*Int(x**(n + S(1))*(a + c*x**S(2))**p, x))
    rubi.add(rule156)


    def cons_f121(g, d, f, e):
        return NonzeroQ(-d*g + e*f)

    cons121 = CustomConstraint(cons_f121)

    def cons_f122(g, c, f, b):
        return ZeroQ(-b*g + S(2)*c*f)

    cons122 = CustomConstraint(cons_f122)

    def cons_f123(f, x):
        return FreeQ(f, x)

    cons123 = CustomConstraint(cons_f123)

    def cons_f124(g, x):
        return FreeQ(g, x)

    cons124 = CustomConstraint(cons_f124)
    pattern157 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons1, cons122)
    rule157 = ReplacementRule(pattern157, lambda g, x, d, c, m, e, f, a, b : (f + g*x)*Int((d + e*x)**m, x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons1, cons122, cons53, cons34)
    rule158 = ReplacementRule(pattern158, lambda g, x, d, c, m, e, f, p, a, b : -f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons1, cons122, cons53, cons36, cons14, cons49)
    rule159 = ReplacementRule(pattern159, lambda g, x, d, c, m, e, f, p, a, b : -e*g*m*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule159)


    def cons_f125(m):
        return Not(And(RationalQ(m), Greater(m, S(0))))

    cons125 = CustomConstraint(cons_f125)
    pattern160 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons1, cons122, cons53, cons11, cons14, cons125)
    rule160 = ReplacementRule(pattern160, lambda g, x, d, c, m, e, f, p, a, b : e*f*g*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))*(-d*g + e*f)) - f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule160)


    def cons_f126(p, m):
        return Or(Not(RationalQ(p)), And(Greater(p, S(0)), Or(Not(IntegerQ(m)), GreaterEqual(m, -S(2)*p + S(-2)), Less(m, -S(4)*p + S(-4)))))

    cons126 = CustomConstraint(cons_f126)
    pattern161 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons1, cons122, cons53, cons48, cons52, cons5, cons126)
    rule161 = ReplacementRule(pattern161, lambda g, x, d, c, m, e, f, p, a, b : -g*(S(2)*p + S(1))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule161)


    def cons_f127(p, m):
        return NonzeroQ(m + S(2)*p + S(2))

    cons127 = CustomConstraint(cons_f127)
    pattern162 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons1, cons122, cons53, cons48, cons52, cons127)
    rule162 = ReplacementRule(pattern162, lambda g, x, d, c, m, e, f, p, a, b : -g*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-d*g + e*f)) + S(2)*f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(1))*(-d*g + e*f)))
    rubi.add(rule162)


    def cons_f128(m):
        return PositiveIntegerQ(m)

    cons128 = CustomConstraint(cons_f128)

    def cons_f129(p, m):
        return Or(Not(RationalQ(p)), Less(m, S(2)*p + S(2)))

    cons129 = CustomConstraint(cons_f129)
    pattern163 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons1, cons122, cons53, cons128, cons127, cons129)
    rule163 = ReplacementRule(pattern163, lambda g, x, d, c, m, e, f, p, a, b : -b*m*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*f*(m + S(2)*p + S(2))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons1, cons122, cons53, cons127)
    rule164 = ReplacementRule(pattern164, lambda g, x, d, c, m, e, f, p, a, b : (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(2))) + (S(2)*p + S(1))*(-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(2)*p + S(2))))
    rubi.add(rule164)


    def cons_f130(g, c, f, b):
        return NonzeroQ(-b*g + S(2)*c*f)

    cons130 = CustomConstraint(cons_f130)

    def cons_f131(p):
        return Less(p, S(0))

    cons131 = CustomConstraint(cons_f131)
    pattern165 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons1, cons130, cons53, cons30, cons11, cons131)
    rule165 = ReplacementRule(pattern165, lambda g, x, d, c, f, e, p, a, b : (-b*g + S(2)*c*f)*Int((a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(d + e*x), x)/(-b*e + S(2)*c*d))
    rubi.add(rule165)


    def cons_f132(c, d, m, e, p, b):
        return Or(And(ZeroQ(m + S(2)*p + S(2)), NonzeroQ(m + S(1))), And(ZeroQ(-b*e + S(2)*c*d), NonzeroQ(m + S(-1))))

    cons132 = CustomConstraint(cons_f132)
    pattern166 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons1, cons130, cons53, cons132)
    rule166 = ReplacementRule(pattern166, lambda g, x, d, c, m, e, f, p, a, b : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule166)

    pattern167 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons1, cons130, cons53, cons30, cons34)
    rule167 = ReplacementRule(pattern167, lambda g, x, d, c, m, e, f, p, a, b : (-b*g + S(2)*c*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons1, cons130, cons53, cons30, cons127, cons73, cons48, cons52)
    rule168 = ReplacementRule(pattern168, lambda g, x, d, c, m, e, f, p, a, b : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))*(-b*e + S(2)*c*d)) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*d*(S(2)*p + S(1))))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule168)


    def cons_f133(g, x, d, f, m, e):
        return Not(And(ZeroQ(m + S(-1)), SimplerQ(f + g*x, d + e*x)))

    cons133 = CustomConstraint(cons_f133)
    pattern169 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons1, cons130, cons53, cons30, cons127, cons73, cons75, cons133)
    rule169 = ReplacementRule(pattern169, lambda g, x, d, c, m, e, f, p, a, b : g*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*c*e*(m + S(2)*p + S(2))) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*(S(2)*d*p + d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*e*(m + S(2)*p + S(2))))
    rubi.add(rule169)


    def cons_f134(n, x):
        return FreeQ(n, x)

    cons134 = CustomConstraint(cons_f134)
    pattern170 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons1, cons53)
    rule170 = ReplacementRule(pattern170, lambda g, x, d, c, f, e, m, p, a, b, n : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons7, cons54, cons55)
    rule171 = ReplacementRule(pattern171, lambda g, x, c, d, f, m, e, p, a, b, n : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule171)


    def cons_f135(p, d, a, m):
        return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), ZeroQ(m + p)))

    cons135 = CustomConstraint(cons_f135)
    pattern172 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons56, cons135)
    rule172 = ReplacementRule(pattern172, lambda g, x, c, d, f, m, e, p, a, n : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons134, cons6, cons54, cons118, cons120)
    rule173 = ReplacementRule(pattern173, lambda g, x, c, d, f, e, m, p, a, b, n : d**m*e**m*Int((f + g*x)**n*(a*e + c*d*x)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule173)

    pattern174 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons134, cons6, cons56, cons118, cons120)
    rule174 = ReplacementRule(pattern174, lambda g, x, c, d, f, e, m, p, a, n : d**m*e**m*Int((a + c*x**S(2))**(m + p)*(f + g*x)**n*(a*e + c*d*x)**(-m), x))
    rubi.add(rule174)


    def cons_f136(g, c, d, m, e, f, p, b):
        return ZeroQ(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))

    cons136 = CustomConstraint(cons_f136)
    pattern175 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons136)
    rule175 = ReplacementRule(pattern175, lambda g, x, d, c, m, e, f, p, a, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule175)


    def cons_f137(g, d, m, e, f, p):
        return ZeroQ(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))

    cons137 = CustomConstraint(cons_f137)
    pattern176 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons137)
    rule176 = ReplacementRule(pattern176, lambda g, x, c, d, f, m, e, p, a : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons36, cons14, cons49)
    rule177 = ReplacementRule(pattern177, lambda g, x, d, c, m, e, f, p, a, b : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons36, cons14, cons49)
    rule178 = ReplacementRule(pattern178, lambda g, x, d, c, m, e, f, p, a : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule178)


    def cons_f138(p):
        return SumSimplerQ(p, S(1))

    cons138 = CustomConstraint(cons_f138)

    def cons_f139(m):
        return SumSimplerQ(m, S(-1))

    cons139 = CustomConstraint(cons_f139)
    pattern179 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons138, cons139, cons71)
    rule179 = ReplacementRule(pattern179, lambda g, x, d, c, m, e, f, p, a, b : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons138, cons139, cons71)
    rule180 = ReplacementRule(pattern180, lambda g, x, c, d, f, m, e, p, a : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule180)


    def cons_f140(p, m):
        return Or(And(RationalQ(m), Less(m, S(-1)), Not(PositiveIntegerQ(m + p + S(1)))), And(RationalQ(m, p), Less(m, S(0)), Less(p, S(-1))), ZeroQ(m + S(2)*p + S(2)))

    cons140 = CustomConstraint(cons_f140)
    pattern181 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons140, cons47)
    rule181 = ReplacementRule(pattern181, lambda g, x, d, c, m, e, f, p, a, b : (d + e*x)**m*(d*g - e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons140, cons47)
    rule182 = ReplacementRule(pattern182, lambda g, x, c, d, f, m, e, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g - e*f)/(S(2)*c*d*(m + p + S(1))) + (S(2)*c*e*f*(p + S(1)) + m*(c*d*g + c*e*f))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*c*d*e*(m + p + S(1))))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons127)
    rule183 = ReplacementRule(pattern183, lambda g, x, d, c, m, e, f, p, a, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(c*e*(m + S(2)*p + S(2))))
    rubi.add(rule183)

    pattern184 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons127)
    rule184 = ReplacementRule(pattern184, lambda g, x, c, d, f, m, e, p, a : (S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/(e*(m + S(2)*p + S(2))) + g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule184)


    def cons_f141(f, a, g, c):
        return ZeroQ(a*g**S(2) + c*f**S(2))

    cons141 = CustomConstraint(cons_f141)

    def cons_f142(p):
        return Less(p, S(-2))

    cons142 = CustomConstraint(cons_f142)
    pattern185 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons4, cons123, cons124, cons141, cons11, cons142)
    rule185 = ReplacementRule(pattern185, lambda g, x, c, f, p, a : x**S(2)*(a + c*x**S(2))**(p + S(1))*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int(x*(a + c*x**S(2))**(p + S(1))*Simp(S(2)*a*g - c*f*x*(S(2)*p + S(5)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule185)

    pattern186 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons4, cons123, cons124, cons6, cons141)
    rule186 = ReplacementRule(pattern186, lambda g, x, c, f, p, a : -f**S(2)*Int((a + c*x**S(2))**(p + S(1))/(f - g*x), x)/c + Int((a + c*x**S(2))**(p + S(1))*(f + g*x), x)/c)
    rubi.add(rule186)


    def cons_f143(m, n):
        return IntegersQ(m, n)

    cons143 = CustomConstraint(cons_f143)

    def cons_f144(p, m):
        return Or(Less(S(0), -m, p + S(1)), Less(p, -m, S(0)))

    cons144 = CustomConstraint(cons_f144)
    pattern187 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons53, cons143, cons11, cons144)
    rule187 = ReplacementRule(pattern187, lambda g, x, c, d, f, e, m, p, a, b, n : Int((f + g*x)**n*(a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule187)

    pattern188 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons53, cons143, cons11, cons144)
    rule188 = ReplacementRule(pattern188, lambda g, x, c, d, f, e, m, p, a, n : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m)*(f + g*x)**n, x))
    rubi.add(rule188)


    def cons_f145(n):
        return PositiveIntegerQ(n)

    cons145 = CustomConstraint(cons_f145)

    def cons_f146(p, n):
        return NegativeIntegerQ(n + S(2)*p)

    cons146 = CustomConstraint(cons_f146)
    pattern189 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons53, cons145, cons146)
    rule189 = ReplacementRule(pattern189, lambda g, x, c, d, f, e, p, a, b, n : -(f + g*x)**n*(a*(-b*e + S(2)*c*d) + c*x*(-S(2)*a*e + b*d))*(a + b*x + c*x**S(2))**p/(d*e*p*(-S(4)*a*c + b**S(2))) - Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*Simp(-S(2)*a*c*(d*g*n - e*f*(S(2)*p + S(1))) + b*(a*e*g*n - c*d*f*(S(2)*p + S(1))) - c*g*x*(-S(2)*a*e + b*d)*(n + S(2)*p + S(1)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons53, cons145, cons146)
    rule190 = ReplacementRule(pattern190, lambda g, x, c, d, f, e, p, a, n : (a + c*x**S(2))**p*(d - e*x)*(f + g*x)**n/(S(2)*d*e*p) - Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*Simp(d*g*n - e*f*(S(2)*p + S(1)) - e*g*x*(n + S(2)*p + S(1)), x), x)/(S(2)*d*e*p))
    rubi.add(rule190)


    def cons_f147(n):
        return NegativeIntegerQ(n)

    cons147 = CustomConstraint(cons_f147)
    pattern191 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons53, cons147, cons146)
    rule191 = ReplacementRule(pattern191, lambda g, x, c, d, f, e, p, a, b, n : -(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p*(a*c*d*(-b*g + S(2)*c*f) - a*e*(S(2)*a*c*g - b**S(2)*g + b*c*f) + c*x*(-a*e*(-b*g + S(2)*c*f) + c*d*(-S(2)*a*g + b*f)))/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**p*Simp(S(2)*a*c*(a*e*g**S(2)*(n + S(2)*p + S(1)) + c*f*(-d*g*n + S(2)*e*f*p + e*f)) + b**S(2)*g*(-a*e*g*(n + p + S(1)) + c*d*f*p) + b*c*(a*g*(d*g*(n + S(1)) + e*f*(n - S(2)*p)) - c*d*f**S(2)*(S(2)*p + S(1))) + c*g*x*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f))*(n + S(2)*p + S(2)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons53, cons147, cons146)
    rule192 = ReplacementRule(pattern192, lambda g, x, c, d, f, e, p, a, n : (a + c*x**S(2))**p*(f + g*x)**(n + S(1))*(-a*e*g + c*d*f - c*x*(d*g + e*f))/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))) + Int((a + c*x**S(2))**p*(f + g*x)**n*Simp(a*e*g**S(2)*(n + S(2)*p + S(1)) - c*f*(d*g*n - e*(S(2)*f*p + f)) + c*g*x*(d*g + e*f)*(n + S(2)*p + S(2)), x), x)/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))))
    rubi.add(rule192)


    def cons_f148(g, d, c, f, e, b):
        return ZeroQ(-b*e*g + c*d*g + c*e*f)

    cons148 = CustomConstraint(cons_f148)

    def cons_f149(m, n):
        return NonzeroQ(m - n + S(-1))

    cons149 = CustomConstraint(cons_f149)
    pattern193 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons7, cons54, cons53, cons58, cons148, cons149)
    rule193 = ReplacementRule(pattern193, lambda g, x, c, d, f, e, m, p, a, b, n : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))))
    rubi.add(rule193)


    def cons_f150(g, d, f, e):
        return ZeroQ(d*g + e*f)

    cons150 = CustomConstraint(cons_f150)
    pattern194 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons56, cons53, cons58, cons150, cons149)
    rule194 = ReplacementRule(pattern194, lambda g, x, c, d, f, e, m, p, a, n : -e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule194)


    def cons_f151(m, n):
        return ZeroQ(m - n + S(-2))

    cons151 = CustomConstraint(cons_f151)
    pattern195 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons7, cons54, cons53, cons58, cons151)
    rule195 = ReplacementRule(pattern195, lambda g, x, c, d, f, e, m, p, a, b, n : -e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule195)

    pattern196 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons56, cons53, cons58, cons151)
    rule196 = ReplacementRule(pattern196, lambda g, x, c, d, f, e, m, p, a, n : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(n + S(1))*(d*g + e*f)))
    rubi.add(rule196)


    def cons_f152(p, n):
        return RationalQ(n, p)

    cons152 = CustomConstraint(cons_f152)

    def cons_f153(n):
        return Less(n, S(-1))

    cons153 = CustomConstraint(cons_f153)

    def cons_f154(p, n):
        return Not(And(IntegerQ(n + p), LessEqual(n + p + S(2), S(0))))

    cons154 = CustomConstraint(cons_f154)
    pattern197 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons53, cons58, cons152, cons12, cons153, cons154)
    rule197 = ReplacementRule(pattern197, lambda g, x, c, d, f, e, m, p, a, b, n : c*m*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*g*(n + S(1))) + (d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(n + S(1))))
    rubi.add(rule197)

    pattern198 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons53, cons58, cons152, cons12, cons153, cons154)
    rule198 = ReplacementRule(pattern198, lambda g, x, c, d, f, e, m, p, a, n : c*m*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(1)), x)/(e*g*(n + S(1))) + (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(n + S(1))))
    rubi.add(rule198)


    def cons_f155(n):
        return Not(PositiveIntegerQ(n))

    cons155 = CustomConstraint(cons_f155)

    def cons_f156(p, n):
        return Not(And(IntegerQ(n + p), Less(n + p + S(2), S(0))))

    cons156 = CustomConstraint(cons_f156)
    pattern199 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons134, cons121, cons7, cons54, cons53, cons58, cons152, cons12, cons149, cons155, cons156)
    rule199 = ReplacementRule(pattern199, lambda g, x, c, d, f, e, m, p, a, b, n : -(d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(m - n + S(-1))) - m*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**(m + S(1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*g*(m - n + S(-1))))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons134, cons121, cons56, cons53, cons58, cons152, cons12, cons149, cons155, cons156)
    rule200 = ReplacementRule(pattern200, lambda g, x, c, d, f, e, m, p, a, n : -c*m*(d*g + e*f)*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**n, x)/(e**S(2)*g*(m - n + S(-1))) - (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(m - n + S(-1))))
    rubi.add(rule200)


    def cons_f157(n):
        return Greater(n, S(0))

    cons157 = CustomConstraint(cons_f157)
    pattern201 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54, cons53, cons58, cons152, cons14, cons157)
    rule201 = ReplacementRule(pattern201, lambda g, x, c, d, f, e, m, p, a, b, n : -e*g*n*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56, cons53, cons58, cons152, cons14, cons157)
    rule202 = ReplacementRule(pattern202, lambda g, x, c, d, f, e, m, p, a, n : -e*g*n*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(p + S(1))))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons134, cons121, cons7, cons54, cons53, cons58, cons152, cons14)
    rule203 = ReplacementRule(pattern203, lambda g, x, c, d, f, e, m, p, a, b, n : e**S(2)*g*(m - n + S(-2))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-b*e*g + c*d*g + c*e*f)) + e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons134, cons121, cons56, cons53, cons58, cons152, cons14)
    rule204 = ReplacementRule(pattern204, lambda g, x, c, d, f, e, m, p, a, n : e**S(2)*g*(m - n + S(-2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(c*(p + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(p + S(1))*(d*g + e*f)))
    rubi.add(rule204)


    def cons_f158(n):
        return RationalQ(n)

    cons158 = CustomConstraint(cons_f158)

    def cons_f159(p, n):
        return Or(IntegerQ(S(2)*p), IntegerQ(n))

    cons159 = CustomConstraint(cons_f159)
    pattern205 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons53, cons58, cons158, cons157, cons149, cons159)
    rule205 = ReplacementRule(pattern205, lambda g, x, c, d, f, e, m, p, a, b, n : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))) - n*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*e*(m - n + S(-1))))
    rubi.add(rule205)

    pattern206 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons53, cons58, cons158, cons157, cons149, cons159)
    rule206 = ReplacementRule(pattern206, lambda g, x, c, d, f, e, m, p, a, n : -n*(d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/(e*(m - n + S(-1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule206)

    pattern207 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons53, cons58, cons158, cons153, cons39)
    rule207 = ReplacementRule(pattern207, lambda g, x, c, d, f, e, m, p, a, b, n : -c*e*(m - n + S(-2))*Int((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/((n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule207)

    pattern208 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons53, cons58, cons158, cons153, cons39)
    rule208 = ReplacementRule(pattern208, lambda g, x, c, d, f, e, m, p, a, n : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/((n + S(1))*(c*d*g + c*e*f)) - e*(m - n + S(-2))*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1)), x)/((n + S(1))*(d*g + e*f)))
    rubi.add(rule208)

    pattern209 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/((x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons54)
    rule209 = ReplacementRule(pattern209, lambda g, x, c, d, f, e, a, b : S(2)*e**S(2)*Subst(Int(S(1)/(-b*e*g + c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons56)
    rule210 = ReplacementRule(pattern210, lambda g, x, c, d, f, e, a : S(2)*e**S(2)*Subst(Int(S(1)/(c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule210)


    def cons_f160(p, m):
        return ZeroQ(m + p + S(-1))

    cons160 = CustomConstraint(cons_f160)

    def cons_f161(g, c, d, f, e, p, b, n):
        return ZeroQ(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))

    cons161 = CustomConstraint(cons_f161)

    def cons_f162(p, n):
        return NonzeroQ(n + p + S(2))

    cons162 = CustomConstraint(cons_f162)
    pattern211 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons7, cons54, cons53, cons160, cons161, cons162)
    rule211 = ReplacementRule(pattern211, lambda g, x, c, d, f, e, m, p, a, b, n : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule211)


    def cons_f163(g, d, f, e, p, n):
        return ZeroQ(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))

    cons163 = CustomConstraint(cons_f163)
    pattern212 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons56, cons53, cons160, cons163, cons162)
    rule212 = ReplacementRule(pattern212, lambda g, x, c, d, f, e, m, p, a, n : e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule212)

    pattern213 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons53, cons160, cons158, cons153, cons39)
    rule213 = ReplacementRule(pattern213, lambda g, x, c, d, f, e, m, p, a, b, n : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e*(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule213)

    pattern214 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons53, cons160, cons158, cons153, cons39)
    rule214 = ReplacementRule(pattern214, lambda g, x, c, d, f, e, m, p, a, n : -e*(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1)), x)/(g*(n + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)/(c*g*(n + S(1))*(d*g + e*f)))
    rubi.add(rule214)


    def cons_f164(n):
        return Not(And(RationalQ(n), Less(n, S(-1))))

    cons164 = CustomConstraint(cons_f164)
    pattern215 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons7, cons54, cons53, cons160, cons164, cons39)
    rule215 = ReplacementRule(pattern215, lambda g, x, c, d, f, e, m, p, a, b, n : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))) - (b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)/(c*g*(n + p + S(2))))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons121, cons56, cons53, cons160, cons164, cons39)
    rule216 = ReplacementRule(pattern216, lambda g, x, c, d, f, e, m, p, a, n : -(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(g*(n + p + S(2))) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule216)


    def cons_f165(m, n):
        return Or(PositiveIntegerQ(m), IntegersQ(m, n))

    cons165 = CustomConstraint(cons_f165)
    pattern217 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons134, cons6, cons121, cons7, cons54, cons53, cons165)
    rule217 = ReplacementRule(pattern217, lambda g, x, c, d, f, e, m, p, a, b, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule217)


    def cons_f166(p):
        return IntegerQ(p + S(-1)/2)

    cons166 = CustomConstraint(cons_f166)

    def cons_f167(p, m):
        return Not(And(Less(m, S(0)), Less(p, S(0))))

    cons167 = CustomConstraint(cons_f167)

    def cons_f168(p):
        return Unequal(p, S(1)/2)

    cons168 = CustomConstraint(cons_f168)
    pattern218 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons134, cons6, cons121, cons56, cons166, cons143, cons167, cons168)
    rule218 = ReplacementRule(pattern218, lambda g, x, c, d, f, e, m, p, a, n : Int(ExpandIntegrand(S(1)/sqrt(a + c*x**S(2)), (a + c*x**S(2))**(p + S(1)/2)*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule218)

    pattern219 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons134, cons6, cons121, cons56, cons53, cons165)
    rule219 = ReplacementRule(pattern219, lambda g, x, c, d, f, e, m, p, a, n : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral(x_**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons6, cons7, cons54)
    rule220 = ReplacementRule(pattern220, lambda x, c, d, e, p, a, b : d**S(2)*Int((a + b*x + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((d - e*x)*(a + b*x + c*x**S(2))**p, x)/e**S(2))
    rubi.add(rule220)

    pattern221 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_/(d_ + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons6, cons56)
    rule221 = ReplacementRule(pattern221, lambda x, c, d, e, p, a : d**S(2)*Int((a + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((a + c*x**S(2))**p*(d - e*x), x)/e**S(2))
    rubi.add(rule221)

    pattern222 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons54, cons53, cons73)
    rule222 = ReplacementRule(pattern222, lambda g, x, c, d, f, e, m, p, a, b : g*(d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(3))) - Int((d + e*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*e*g*(d*g + e*f*(m + p + S(1))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))) + e*g*x*(b*e*g*(m + p + S(2)) - c*(d*g*m + e*f*(m + S(2)*p + S(4)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule222)

    pattern223 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons56, cons53, cons73)
    rule223 = ReplacementRule(pattern223, lambda g, x, c, d, f, e, m, p, a : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(f + g*x)/(c*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(d + e*x)**m*Simp(-c*e*g*x*(d*g*m + e*f*(m + S(2)*p + S(4))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons3, cons4, cons26, cons123, cons124, cons27, cons134, cons53)
    rule224 = ReplacementRule(pattern224, lambda g, x, c, f, e, m, p, b, n : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p*(f + g*x)**n, x))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons56, cons53, cons69, cons70)
    rule225 = ReplacementRule(pattern225, lambda g, x, c, d, f, e, m, p, a, n : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule225)

    pattern226 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons7, cons54, cons53)
    rule226 = ReplacementRule(pattern226, lambda g, x, c, d, f, e, m, p, a, b, n : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons121, cons56, cons53)
    rule227 = ReplacementRule(pattern227, lambda g, x, c, d, f, e, m, p, a, n : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons7, cons82, cons8)
    rule228 = ReplacementRule(pattern228, lambda g, x, d, c, m, e, f, p, a, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule228)

    pattern229 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons83, cons8)
    rule229 = ReplacementRule(pattern229, lambda g, x, d, c, m, e, f, p, a : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x), x), x))
    rubi.add(rule229)

    pattern230 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82)
    rule230 = ReplacementRule(pattern230, lambda g, x, d, c, f, e, a, b : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int(Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83)
    rule231 = ReplacementRule(pattern231, lambda g, x, d, c, f, e, a : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int(Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule231)


    def cons_f169(g, d, c, f, e, a, b):
        return ZeroQ(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)

    cons169 = CustomConstraint(cons_f169)
    pattern232 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons82, cons34, cons169)
    rule232 = ReplacementRule(pattern232, lambda g, x, d, c, f, e, m, p, a, b : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule232)


    def cons_f170(g, c, d, f, e, a):
        return ZeroQ(a*e*g + c*d*f)

    cons170 = CustomConstraint(cons_f170)
    pattern233 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons83, cons34, cons170)
    rule233 = ReplacementRule(pattern233, lambda g, x, d, c, f, e, m, p, a : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule233)


    def cons_f171(d, c, m, e, b):
        return Not(And(Equal(m, S(1)), Or(ZeroQ(d), ZeroQ(-b*e + S(2)*c*d))))

    cons171 = CustomConstraint(cons_f171)
    pattern234 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons34, cons11, cons14, cons171)
    rule234 = ReplacementRule(pattern234, lambda g, x, d, c, m, e, f, p, a, b : -m*(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule234)


    def cons_f172(d, m):
        return Not(And(Equal(m, S(1)), ZeroQ(d)))

    cons172 = CustomConstraint(cons_f172)
    pattern235 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons34, cons11, cons14, cons172)
    rule235 = ReplacementRule(pattern235, lambda g, x, d, c, m, e, f, p, a : -m*(a*e*g + c*d*f)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*c*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons82, cons34)
    rule236 = ReplacementRule(pattern236, lambda g, x, d, c, f, e, m, p, a, b : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons83, cons34)
    rule237 = ReplacementRule(pattern237, lambda g, x, d, c, f, e, m, p, a : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e*g + c*d*f)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule237)


    def cons_f173(g, c, d, f, e, p, a, b):
        return ZeroQ(-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))

    cons173 = CustomConstraint(cons_f173)
    pattern238 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons7, cons82, cons173)
    rule238 = ReplacementRule(pattern238, lambda g, x, d, c, f, e, p, a, b : -(a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule238)


    def cons_f174(g, d, c, f, e, p, a):
        return ZeroQ(a*e*g - c*d*f*(S(2)*p + S(3)))

    cons174 = CustomConstraint(cons_f174)
    pattern239 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons83, cons174)
    rule239 = ReplacementRule(pattern239, lambda g, x, d, c, f, e, p, a : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons11, cons14)
    rule240 = ReplacementRule(pattern240, lambda g, x, d, c, f, e, p, a, b : -(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g - b*c*(d*g + e*f) + S(2)*c*(-a*e*g + c*d*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons11, cons14)
    rule241 = ReplacementRule(pattern241, lambda g, x, d, c, f, e, p, a : -(a + c*x**S(2))**(p + S(1))*(-a*(d*g + e*(f + g*x)) + c*d*f*x)/(S(2)*a*c*(p + S(1))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons7, cons82, cons91)
    rule242 = ReplacementRule(pattern242, lambda g, x, d, c, f, e, p, a, b : (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c**S(2)*(S(2)*p + S(3))) - (a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons83, cons91)
    rule243 = ReplacementRule(pattern243, lambda g, x, d, c, f, e, p, a : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**p, x)/(c*(S(2)*p + S(3))))
    rubi.add(rule243)


    def cons_f175(m):
        return Not(RationalQ(m))

    cons175 = CustomConstraint(cons_f175)

    def cons_f176(p):
        return Not(PositiveIntegerQ(p))

    cons176 = CustomConstraint(cons_f176)
    pattern244 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons2, cons4, cons26, cons123, cons124, cons6, cons175, cons176)
    rule244 = ReplacementRule(pattern244, lambda g, x, c, m, e, f, p, a : f*Int((e*x)**m*(a + c*x**S(2))**p, x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule244)


    def cons_f177(p, m):
        return ZeroQ(m - p)

    cons177 = CustomConstraint(cons_f177)
    pattern245 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons177, cons92, cons93)
    rule245 = ReplacementRule(pattern245, lambda g, x, d, c, f, e, m, p, a, b : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((f + g*x)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule245)


    def cons_f178(p, m):
        return Less(m + S(2)*p, S(0))

    cons178 = CustomConstraint(cons_f178)

    def cons_f179(p, m):
        return Not(NegativeIntegerQ(m + S(2)*p + S(3)))

    cons179 = CustomConstraint(cons_f179)
    pattern246 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons36, cons12, cons40, cons178, cons179)
    rule246 = ReplacementRule(pattern246, lambda g, x, d, c, f, e, m, p, a, b : -p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) + b**S(2)*e*(d*g*(p + S(1)) - e*f*(m + p + S(2))) + b*(a*e**S(2)*g*(m + S(1)) - c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))) - c*x*(S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2))) - e*(S(2)*a*e*g*(m + S(1)) - b*(d*g*(m - S(2)*p) + e*f*(m + S(2)*p + S(2))))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*p*(-b*e + S(2)*c*d)*(-d*g + e*f) - e*x*(g*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)) + p*(-b*e + S(2)*c*d)*(-d*g + e*f)) + (d*g - e*f*(m + S(2)))*(a*e**S(2) - b*d*e + c*d**S(2)))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule246)

    pattern247 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons36, cons12, cons40, cons178, cons179)
    rule247 = ReplacementRule(pattern247, lambda g, x, d, c, f, e, m, p, a : -p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) - c*x*(-S(2)*a*e**S(2)*g*(m + S(1)) + S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*c*d**S(2)*p*(-d*g + e*f) - e*x*(S(2)*c*d*p*(-d*g + e*f) + g*(m + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e**S(2) + c*d**S(2))*(d*g - e*f*(m + S(2))))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule247)


    def cons_f180(p, m):
        return Or(And(RationalQ(m), Less(m, S(-1))), Equal(p, S(1)), And(IntegerQ(p), Not(RationalQ(m))))

    cons180 = CustomConstraint(cons_f180)

    def cons_f181(p, m):
        return Or(IntegerQ(m), IntegerQ(p), IntegersQ(S(2)*m, S(2)*p))

    cons181 = CustomConstraint(cons_f181)
    pattern248 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons7, cons82, cons11, cons12, cons180, cons32, cons103, cons181)
    rule248 = ReplacementRule(pattern248, lambda g, x, d, c, f, e, m, p, a, b : p*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-b*e*f*(m + S(2)*p + S(2)) + g*(S(2)*a*e*m + S(2)*a*e + S(2)*b*d*p + b*d) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(b*e*m + b*e + S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule248)

    pattern249 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons83, cons11, cons12, cons180, cons32, cons103, cons181)
    rule249 = ReplacementRule(pattern249, lambda g, x, d, c, f, e, m, p, a : p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*Simp(g*(S(2)*a*e*m + S(2)*a*e) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule249)


    def cons_f182(p, m):
        return Or(IntegerQ(p), Not(RationalQ(m)), Inequality(S(-1), LessEqual, m, Less, S(0)))

    cons182 = CustomConstraint(cons_f182)
    pattern250 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons7, cons82, cons11, cons12, cons182, cons107, cons181)
    rule250 = ReplacementRule(pattern250, lambda g, x, d, c, f, e, m, p, a, b : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(c*e*f*(-S(2)*a*e + b*d)*(m + S(2)*p + S(2)) + g*(a*e*(b*e*m + b*e - S(2)*c*d*m) + b*d*(b*e*p - S(2)*c*d*p - c*d)) + x*(c*e*f*(-b*e + S(2)*c*d)*(m + S(2)*p + S(2)) + g*(b**S(2)*e**S(2)*(m + p + S(1)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(1)) - c*e*(S(2)*a*e*(m + S(2)*p + S(1)) + b*d*(m - S(2)*p)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)) - g*(-b*e*p + S(2)*c*d*p + c*d))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons83, cons11, cons12, cons182, cons107, cons181)
    rule251 = ReplacementRule(pattern251, lambda g, x, d, c, f, e, m, p, a : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*c*d*e*g*m + a*c*e**S(2)*f*(m + S(2)*p + S(2)) - x*(c**S(2)*d*e*f*(m + S(2)*p + S(2)) - g*(a*c*e**S(2)*(m + S(2)*p + S(1)) + c**S(2)*d**S(2)*(S(2)*p + S(1)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*d*g*(S(2)*p + S(1)) + c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule251)


    def cons_f183(g, c, d, m, e, f, p, a, b):
        return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, b, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

    cons183 = CustomConstraint(cons_f183)
    pattern252 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons36, cons14, cons46, cons183)
    rule252 = ReplacementRule(pattern252, lambda g, x, d, c, f, e, m, p, a, b : -(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g + S(2)*c**S(2)*d*f - c*(S(2)*a*e*g + b*d*g + b*e*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(b*e*g*(a*e*(m + S(-1)) + b*d*(p + S(2))) + S(2)*c**S(2)*d**S(2)*f*(S(2)*p + S(3)) - c*(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) + b*d*(d*g*(S(2)*p + S(3)) - e*f*(m - S(2)*p + S(-4)))) + e*x*(b**S(2)*e*g*(m + p + S(1)) + S(2)*c**S(2)*d*f*(m + S(2)*p + S(2)) - c*(S(2)*a*e*g*m + b*(d*g + e*f)*(m + S(2)*p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule252)


    def cons_f184(g, c, d, m, e, f, p, a):
        return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

    cons184 = CustomConstraint(cons_f184)
    pattern253 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons36, cons14, cons46, cons184)
    rule253 = ReplacementRule(pattern253, lambda g, x, d, c, f, e, m, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(S(2)*a*(d*g + e*f) - x*(-S(2)*a*e*g + S(2)*c*d*f))/(S(4)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) - S(2)*c*d**S(2)*f*(S(2)*p + S(3)) + e*x*(S(2)*a*e*g*m - S(2)*c*d*f*(m + S(2)*p + S(2))), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule253)


    def cons_f185(g, x, d, m, e, f):
        return Not(And(Equal(m, S(1)), SimplerQ(d + e*x, f + g*x)))

    cons185 = CustomConstraint(cons_f185)
    pattern254 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons36, cons14, cons49, cons185, cons181)
    rule254 = ReplacementRule(pattern254, lambda g, x, d, c, m, e, f, p, a, b : (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-e*x*(-b*g + S(2)*c*f)*(m + S(2)*p + S(3)) - f*(b*e*m + S(2)*c*d*(S(2)*p + S(3))) + g*(S(2)*a*e*m + b*d*(S(2)*p + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons36, cons14, cons49, cons185, cons181)
    rule255 = ReplacementRule(pattern255, lambda g, x, d, c, m, e, f, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*Simp(a*e*g*m - c*d*f*(S(2)*p + S(3)) - c*e*f*x*(m + S(2)*p + S(3)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons7, cons82, cons11, cons14, cons181)
    rule256 = ReplacementRule(pattern256, lambda g, x, d, c, f, e, m, p, a, b : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-a*g*(-b*e + S(2)*c*d) + c*x*(f*(-b*e + S(2)*c*d) - g*(-S(2)*a*e + b*d)) + f*(S(2)*a*c*e - b**S(2)*e + b*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*e*x*(-f*(-b*e + S(2)*c*d) + g*(-S(2)*a*e + b*d))*(m + S(2)*p + S(4)) + f*(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3))) - g*(a*e*(b*e*m + b*e - S(2)*c*d*m) - b*d*(-b*e*p - b*e + S(2)*c*d*p + S(3)*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule256)

    pattern257 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons11, cons14, cons181)
    rule257 = ReplacementRule(pattern257, lambda g, x, d, c, f, e, m, p, a : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-a*c*d*g + a*c*e*f + c*x*(a*e*g + c*d*f))/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(-a*c*d*e*g*m + c*e*x*(a*e*g + c*d*f)*(m + S(2)*p + S(4)) + f*(a*c*e**S(2)*(m + S(2)*p + S(3)) + c**S(2)*d**S(2)*(S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons59)
    rule258 = ReplacementRule(pattern258, lambda g, x, d, c, m, e, f, a, b : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons59)
    rule259 = ReplacementRule(pattern259, lambda g, x, d, c, m, e, f, a : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule259)


    def cons_f186(m):
        return FractionQ(m)

    cons186 = CustomConstraint(cons_f186)
    pattern260 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons186, cons49)
    rule260 = ReplacementRule(pattern260, lambda g, x, d, c, f, e, m, a, b : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule260)

    pattern261 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons186, cons49)
    rule261 = ReplacementRule(pattern261, lambda g, x, d, c, f, e, m, a : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule261)

    pattern262 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82)
    rule262 = ReplacementRule(pattern262, lambda g, x, d, c, f, e, a, b : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83)
    rule263 = ReplacementRule(pattern263, lambda g, x, d, c, f, e, a : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule263)

    pattern264 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons7, cons82, cons186, cons52)
    rule264 = ReplacementRule(pattern264, lambda g, x, d, c, f, e, m, a, b : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons121, cons83, cons186, cons52)
    rule265 = ReplacementRule(pattern265, lambda g, x, d, c, f, e, m, a : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons175)
    rule266 = ReplacementRule(pattern266, lambda g, x, d, c, f, e, m, a, b : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons175)
    rule267 = ReplacementRule(pattern267, lambda g, x, d, c, f, e, m, a : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule267)


    def cons_f187(g, x, d, f, m, e):
        return Not(And(Equal(m, S(1)), SimplerQ(f + g*x, d + e*x)))

    cons187 = CustomConstraint(cons_f187)
    pattern268 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons7, cons82, cons48, cons49, cons127, cons187, cons181)
    rule268 = ReplacementRule(pattern268, lambda g, x, d, c, m, e, f, p, a, b : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p*Simp(d*(p + S(1))*(-b*g + S(2)*c*f) + m*(-a*e*g + c*d*f) + x*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(-b*e*g + c*d*g + c*e*f)), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons83, cons48, cons49, cons127, cons187, cons181)
    rule269 = ReplacementRule(pattern269, lambda g, x, d, c, m, e, f, p, a : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*Simp(-a*e*g*m + c*d*f*(m + S(2)*p + S(2)) + c*x*(d*g*m + e*f*(m + S(2)*p + S(2))), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons7, cons82, cons48, cons52, cons181)
    rule270 = ReplacementRule(pattern270, lambda g, x, d, c, f, e, m, p, a, b : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons6, cons121, cons83, cons48, cons52, cons181)
    rule271 = ReplacementRule(pattern271, lambda g, x, d, c, f, e, m, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule271)


    def cons_f188(p, m):
        return NegativeIntegerQ(m + S(2)*p + S(3))

    cons188 = CustomConstraint(cons_f188)
    pattern272 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons82, cons188, cons32)
    rule272 = ReplacementRule(pattern272, lambda g, x, d, c, f, e, m, p, a, b : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons83, cons188, cons32)
    rule273 = ReplacementRule(pattern273, lambda g, x, d, c, f, e, m, p, a : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule273)


    def cons_f189(c, d, e, a, b):
        return ZeroQ(S(4)*c*(a - d) - (b - e)**S(2))

    cons189 = CustomConstraint(cons_f189)

    def cons_f190(g, d, f, e, a, b):
        return ZeroQ(e*f*(b - e) - S(2)*g*(-a*e + b*d))

    cons190 = CustomConstraint(cons_f190)

    def cons_f191(d, a, e, b):
        return NonzeroQ(-a*e + b*d)

    cons191 = CustomConstraint(cons_f191)
    pattern274 = Pattern(Integral((f_ + x_*WC('g', S(1)))/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons189, cons190, cons191)
    rule274 = ReplacementRule(pattern274, lambda g, x, d, c, f, e, a, b : S(4)*f*(a - d)*Subst(Int(S(1)/(S(4)*a - S(4)*d - x**S(2)), x), x, (S(2)*a - S(2)*d + x*(b - e))/sqrt(a + b*x + c*x**S(2)))/(-a*e + b*d))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons123, cons124, cons7)
    rule275 = ReplacementRule(pattern275, lambda g, x, c, f, a, b : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule275)


    def cons_f192(g, x, c, f, a):
        return FreeQ(List(a, c, f, g), x)

    cons192 = CustomConstraint(cons_f192)
    pattern276 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons4, cons123, cons124, cons192)
    rule276 = ReplacementRule(pattern276, lambda g, x, c, f, a : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons26, cons123, cons124, cons7)
    rule277 = ReplacementRule(pattern277, lambda g, x, c, f, e, a, b : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + b*x + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule277)


    def cons_f193(g, x, c, f, e, a):
        return FreeQ(List(a, c, e, f, g), x)

    cons193 = CustomConstraint(cons_f193)
    pattern278 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons2, cons4, cons26, cons123, cons124, cons193)
    rule278 = ReplacementRule(pattern278, lambda g, x, c, f, e, a : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons82)
    rule279 = ReplacementRule(pattern279, lambda g, x, d, c, f, e, m, p, a, b : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule279)

    pattern280 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons83)
    rule280 = ReplacementRule(pattern280, lambda g, x, d, c, f, e, m, p, a : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/e)
    rubi.add(rule280)


    def cons_f194(p, m, n):
        return IntegersQ(m, n, p)

    cons194 = CustomConstraint(cons_f194)
    pattern281 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons194)
    rule281 = ReplacementRule(pattern281, lambda g, x, d, c, f, e, m, p, a, b, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule281)

    pattern282 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons194)
    rule282 = ReplacementRule(pattern282, lambda g, x, d, c, f, e, m, p, a, n : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule282)


    def cons_f195(p):
        return FractionQ(p)

    cons195 = CustomConstraint(cons_f195)
    pattern283 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons195, cons12)
    rule283 = ReplacementRule(pattern283, lambda g, x, d, c, f, e, p, a, b : (a*e**S(2) - b*d*e + c*d**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + b*x + c*x**S(2))**(p + S(-1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons195, cons12)
    rule284 = ReplacementRule(pattern284, lambda g, x, d, c, f, e, p, a : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule284)


    def cons_f196(p, n):
        return IntegersQ(n, p)

    cons196 = CustomConstraint(cons_f196)
    pattern285 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons196, cons186, )
    def With285(g, x, d, c, f, e, m, p, a, b, n):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(c*x**(S(2)*q)/e**S(2) - x**q*(-b*e + S(2)*c*d)/e**S(2) + (a*e**S(2) - b*d*e + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule285 = ReplacementRule(pattern285, lambda g, x, d, c, f, e, m, p, a, b, n : With285(g, x, d, c, f, e, m, p, a, b, n))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons196, cons186, )
    def With286(g, x, d, c, f, e, m, p, a, n):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(-S(2)*c*d*x**q/e**S(2) + c*x**(S(2)*q)/e**S(2) + (a*e**S(2) + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule286 = ReplacementRule(pattern286, lambda g, x, d, c, f, e, m, p, a, n : With286(g, x, d, c, f, e, m, p, a, n))
    rubi.add(rule286)


    def cons_f197(m, n):
        return ZeroQ(m - n)

    cons197 = CustomConstraint(cons_f197)

    def cons_f198(d, m, f):
        return Or(IntegerQ(m), And(PositiveQ(d), PositiveQ(f)))

    cons198 = CustomConstraint(cons_f198)
    pattern287 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons197, cons150, cons198)
    rule287 = ReplacementRule(pattern287, lambda g, x, c, d, m, e, f, p, a, b, n : Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons197, cons150, cons198)
    rule288 = ReplacementRule(pattern288, lambda g, x, c, d, m, e, f, p, a, n : Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons197, cons150)
    rule289 = ReplacementRule(pattern289, lambda g, x, c, d, m, e, f, p, a, b, n : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons197, cons150)
    rule290 = ReplacementRule(pattern290, lambda g, x, d, c, m, e, f, p, a, n : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule290)


    def cons_f199(m, n):
        return RationalQ(m, n)

    cons199 = CustomConstraint(cons_f199)
    pattern291 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons199)
    rule291 = ReplacementRule(pattern291, lambda g, x, d, c, f, e, m, a, b, n : c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x) + Int((a + b*x)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons199)
    rule292 = ReplacementRule(pattern292, lambda g, x, d, c, f, e, m, a, n : a*Int((d + e*x)**m*(f + g*x)**n, x) + c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule292)


    def cons_f200(n):
        return Not(IntegerQ(n))

    cons200 = CustomConstraint(cons_f200)

    def cons_f201(n):
        return Greater(n, S(1))

    cons201 = CustomConstraint(cons_f201)
    pattern293 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons90, cons200, cons199, cons49, cons201)
    rule293 = ReplacementRule(pattern293, lambda g, x, d, c, f, e, m, a, b, n : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-b*e*g + c*d*g + S(2)*c*e*f + c*e*g*x, x), x)/c**S(2) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(a*b*e*g**S(2) - a*c*d*g**S(2) - S(2)*a*c*e*f*g + c**S(2)*d*f**S(2) + x*(-a*c*e*g**S(2) + b**S(2)*e*g**S(2) - b*c*d*g**S(2) - S(2)*b*c*e*f*g + S(2)*c**S(2)*d*f*g + c**S(2)*e*f**S(2)), x)/(a + b*x + c*x**S(2)), x)/c**S(2))
    rubi.add(rule293)

    pattern294 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons90, cons200, cons199, cons49, cons201)
    rule294 = ReplacementRule(pattern294, lambda g, x, d, c, f, e, m, a, n : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(d*g + S(2)*e*f + e*g*x, x), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-a*d*g**S(2) - S(2)*a*e*f*g + c*d*f**S(2) + x*(-a*e*g**S(2) + S(2)*c*d*f*g + c*e*f**S(2)), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule294)

    pattern295 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons90, cons200, cons199, cons49, cons157)
    rule295 = ReplacementRule(pattern295, lambda g, x, d, c, f, e, m, a, b, n : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule295)

    pattern296 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons90, cons200, cons199, cons49, cons157)
    rule296 = ReplacementRule(pattern296, lambda g, x, d, c, f, e, m, a, n : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule296)

    pattern297 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons90, cons200, cons199, cons49, cons153)
    rule297 = ReplacementRule(pattern297, lambda g, x, d, c, f, e, m, a, b, n : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) - b*f*g + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g - b*d*g + c*d*f + c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*g**S(2) - b*f*g + c*f**S(2)))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons90, cons200, cons199, cons49, cons153)
    rule298 = ReplacementRule(pattern298, lambda g, x, d, c, f, e, m, a, n : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g + c*d*f + c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*g**S(2) + c*f**S(2)))
    rubi.add(rule298)


    def cons_f202(m):
        return PositiveIntegerQ(m + S(1)/2)

    cons202 = CustomConstraint(cons_f202)
    pattern299 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons7, cons82, cons202)
    rule299 = ReplacementRule(pattern299, lambda g, x, d, c, f, e, m, a, b : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons83, cons202)
    rule300 = ReplacementRule(pattern300, lambda g, x, d, c, f, e, m, a : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + c*x**S(2)), x), x))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons7, cons82, cons90, cons200)
    rule301 = ReplacementRule(pattern301, lambda g, x, d, c, f, e, m, a, b, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule301)

    pattern302 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons83, cons90, cons200)
    rule302 = ReplacementRule(pattern302, lambda g, x, d, c, f, e, m, a, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule302)


    def cons_f203(p, m, n):
        return Or(IntegerQ(p), IntegersQ(m, n))

    cons203 = CustomConstraint(cons_f203)
    pattern303 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons203)
    rule303 = ReplacementRule(pattern303, lambda g, x, d, c, f, e, m, p, a, b, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons203)
    rule304 = ReplacementRule(pattern304, lambda g, x, d, c, f, e, m, p, a, n : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule304)

    pattern305 = Pattern(Integral((x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons26, cons124, cons27, cons134, cons6, cons177, cons92, cons93)
    rule305 = ReplacementRule(pattern305, lambda g, x, d, c, m, e, p, a, b, n : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((g*x)**n*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule305)

    pattern306 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons200, cons53, cons152, cons12, cons153)
    rule306 = ReplacementRule(pattern306, lambda g, x, d, c, f, e, p, a, b, n : (a*e**S(2) - b*d*e + c*d**S(2))*Int((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1))*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons200, cons53, cons152, cons12, cons153)
    rule307 = ReplacementRule(pattern307, lambda g, x, d, c, f, e, p, a, n : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**(n + S(1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**n*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons200, cons53, cons152, cons14, cons157)
    rule308 = ReplacementRule(pattern308, lambda g, x, d, c, f, e, p, a, b, n : e*(-d*g + e*f)*Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule308)

    pattern309 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons200, cons53, cons152, cons14, cons157)
    rule309 = ReplacementRule(pattern309, lambda g, x, d, c, f, e, p, a, n : e*(-d*g + e*f)*Int((a + c*x**S(2))**(p + S(1))*(f + g*x)**(n + S(-1))/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, )
    def With310(g, x, d, c, f, e, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -sqrt(S(2))*sqrt(-g*(b + S(2)*c*x - q)/(-b*g + S(2)*c*f + g*q))*sqrt(-g*(b + S(2)*c*x + q)/(-b*g + S(2)*c*f - g*q))*EllipticPi(e*(-b*g + S(2)*c*f + g*q)/(S(2)*c*(-d*g + e*f)), asin(sqrt(S(2))*sqrt(c/(-b*g + S(2)*c*f + g*q))*sqrt(f + g*x)), (-b*g + S(2)*c*f + g*q)/(-b*g + S(2)*c*f - g*q))/(sqrt(c/(-b*g + S(2)*c*f + g*q))*(-d*g + e*f)*sqrt(a + b*x + c*x**S(2)))
    rule310 = ReplacementRule(pattern310, lambda g, x, d, c, f, e, a, b : With310(g, x, d, c, f, e, a, b))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, )
    def With311(g, x, d, c, f, e, a):
        q = Rt(-a*c, S(2))
        return -S(2)*sqrt(g*(-c*x + q)/(c*f + g*q))*sqrt(-g*(c*x + q)/(c*f - g*q))*EllipticPi(e*(c*f + g*q)/(c*(-d*g + e*f)), asin(sqrt(c/(c*f + g*q))*sqrt(f + g*x)), (c*f + g*q)/(c*f - g*q))/(sqrt(c/(c*f + g*q))*sqrt(a + c*x**S(2))*(-d*g + e*f))
    rule311 = ReplacementRule(pattern311, lambda g, x, d, c, f, e, a : With311(g, x, d, c, f, e, a))
    rubi.add(rule311)


    def cons_f204(n):
        return IntegerQ(n + S(1)/2)

    cons204 = CustomConstraint(cons_f204)
    pattern312 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82, cons204)
    rule312 = ReplacementRule(pattern312, lambda g, x, d, c, f, e, a, b, n : Int(ExpandIntegrand(S(1)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**S(2))), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule312)

    pattern313 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83, cons204)
    rule313 = ReplacementRule(pattern313, lambda g, x, d, c, f, e, a, n : Int(ExpandIntegrand(S(1)/(sqrt(a + c*x**S(2))*sqrt(f + g*x)), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons121, cons7, cons82)
    rule314 = ReplacementRule(pattern314, lambda g, x, d, c, f, e, a, b : -S(2)*sqrt((-d*g + e*f)**S(2)*(a + b*x + c*x**S(2))/((d + e*x)**S(2)*(a*g**S(2) - b*f*g + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) - b*d*e + c*d**S(2))/(a*g**S(2) - b*f*g + c*f**S(2)) - x**S(2)*(S(2)*a*e*g - b*d*g - b*e*f + S(2)*c*d*f)/(a*g**S(2) - b*f*g + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/((-d*g + e*f)*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule314)

    pattern315 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons121, cons83)
    rule315 = ReplacementRule(pattern315, lambda g, x, d, c, f, e, a : -S(2)*sqrt((a + c*x**S(2))*(-d*g + e*f)**S(2)/((d + e*x)**S(2)*(a*g**S(2) + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) + c*d**S(2))/(a*g**S(2) + c*f**S(2)) - x**S(2)*(S(2)*a*e*g + S(2)*c*d*f)/(a*g**S(2) + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/(sqrt(a + c*x**S(2))*(-d*g + e*f)))
    rubi.add(rule315)


    def cons_f205(g, x, c, f, e, m, p, a):
        return FreeQ(List(a, c, e, f, g, m, p), x)

    cons205 = CustomConstraint(cons_f205)
    pattern316 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(2), x_), cons2, cons4, cons26, cons123, cons124, cons27, cons6, cons205)
    rule316 = ReplacementRule(pattern316, lambda g, x, c, m, e, f, p, a : Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + g**S(2)*x**S(2)), x) + S(2)*f*g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule316)

    pattern317 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(3), x_), cons2, cons4, cons26, cons123, cons124, cons27, cons6, cons205)
    rule317 = ReplacementRule(pattern317, lambda g, x, c, m, e, f, p, a : f*Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + S(3)*g**S(2)*x**S(2)), x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p*(S(3)*f**S(2) + g**S(2)*x**S(2)), x)/e)
    rubi.add(rule317)

    pattern318 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons7, cons82, cons145)
    rule318 = ReplacementRule(pattern318, lambda g, x, d, c, f, e, m, p, a, b, n : g*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule318)

    pattern319 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons6, cons121, cons83, cons145)
    rule319 = ReplacementRule(pattern319, lambda g, x, d, c, f, e, m, p, a, n : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/e)
    rubi.add(rule319)


    def cons_f206(g, x, c, d, f, e, m, p, a, b, n):
        return FreeQ(List(a, b, c, d, e, f, g, m, n, p), x)

    cons206 = CustomConstraint(cons_f206)
    pattern320 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons206)
    rule320 = ReplacementRule(pattern320, lambda g, x, d, c, m, e, f, p, a, b, n : Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule320)


    def cons_f207(g, x, c, d, f, e, m, p, a, n):
        return FreeQ(List(a, c, d, e, f, g, m, n, p), x)

    cons207 = CustomConstraint(cons_f207)
    pattern321 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons207)
    rule321 = ReplacementRule(pattern321, lambda g, x, d, c, m, e, f, p, a, n : Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons21, cons22)
    rule322 = ReplacementRule(pattern322, lambda g, x, d, c, m, e, f, u, p, a, b, n : Subst(Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons27, cons134, cons6, cons21, cons22)
    rule323 = ReplacementRule(pattern323, lambda g, x, d, c, m, e, f, u, p, a, n : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule323)


    def cons_f208(c, a, d, f):
        return ZeroQ(-a*f + c*d)

    cons208 = CustomConstraint(cons_f208)

    def cons_f209(d, a, e, b):
        return ZeroQ(-a*e + b*d)

    cons209 = CustomConstraint(cons_f209)

    def cons_f210(p, c, f):
        return Or(IntegerQ(p), PositiveQ(c/f))

    cons210 = CustomConstraint(cons_f210)

    def cons_f211(x, d, c, f, e, q, a, b):
        return Or(Not(IntegerQ(q)), LessEqual(LeafCount(d + e*x + f*x**S(2)), LeafCount(a + b*x + c*x**S(2))))

    cons211 = CustomConstraint(cons_f211)

    def cons_f212(q, x):
        return FreeQ(q, x)

    cons212 = CustomConstraint(cons_f212)
    pattern324 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons6, cons212, cons208, cons209, cons210, cons211)
    rule324 = ReplacementRule(pattern324, lambda x, c, d, f, e, q, p, a, b : (c/f)**p*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule324)


    def cons_f213(q):
        return Not(IntegerQ(q))

    cons213 = CustomConstraint(cons_f213)

    def cons_f214(c, f):
        return Not(PositiveQ(c/f))

    cons214 = CustomConstraint(cons_f214)
    pattern325 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons6, cons212, cons208, cons209, cons53, cons213, cons214)
    rule325 = ReplacementRule(pattern325, lambda x, c, d, f, e, q, p, a, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule325)

    pattern326 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons6, cons212, cons1, cons53)
    rule326 = ReplacementRule(pattern326, lambda x, c, d, f, e, q, p, a, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons6, cons212, cons1, cons53)
    rule327 = ReplacementRule(pattern327, lambda x, c, d, f, q, p, a, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule327)


    def cons_f215(c, d, f, e, q, a, b):
        return ZeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

    cons215 = CustomConstraint(cons_f215)

    def cons_f216(q):
        return NonzeroQ(q + S(1))

    cons216 = CustomConstraint(cons_f216)

    def cons_f217(q):
        return NonzeroQ(S(2)*q + S(3))

    cons217 = CustomConstraint(cons_f217)
    pattern328 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons212, cons215, cons216, cons217)
    rule328 = ReplacementRule(pattern328, lambda x, c, d, f, e, q, a, b : (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule328)


    def cons_f218(c, d, f, e, q, a):
        return ZeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

    cons218 = CustomConstraint(cons_f218)
    pattern329 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons212, cons218, cons216, cons217)
    rule329 = ReplacementRule(pattern329, lambda x, c, d, f, e, q, a : (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule329)


    def cons_f219(d, c, f, q, a):
        return ZeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

    cons219 = CustomConstraint(cons_f219)
    pattern330 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons123, cons212, cons216, cons219)
    rule330 = ReplacementRule(pattern330, lambda x, c, d, f, q, a, b : (d + f*x**S(2))**(q + S(1))*(S(2)*a*f*x*(q + S(1)) + b*d)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule330)


    def cons_f220(q):
        return PositiveIntegerQ(q + S(2))

    cons220 = CustomConstraint(cons_f220)
    pattern331 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons123, cons212, cons220)
    rule331 = ReplacementRule(pattern331, lambda x, c, d, f, q, a, b : b*Int(x*(d + f*x**S(2))**q, x) + Int((a + c*x**S(2))*(d + f*x**S(2))**q, x))
    rubi.add(rule331)


    def cons_f221(d, f, e):
        return NonzeroQ(-S(4)*d*f + e**S(2))

    cons221 = CustomConstraint(cons_f221)
    pattern332 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons221, cons220)
    rule332 = ReplacementRule(pattern332, lambda x, c, d, f, e, q, a, b : Int(ExpandIntegrand((a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule332)

    pattern333 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons221, cons220)
    rule333 = ReplacementRule(pattern333, lambda x, c, d, f, e, q, a : Int(ExpandIntegrand((a + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule333)


    def cons_f222(q):
        return RationalQ(q)

    cons222 = CustomConstraint(cons_f222)

    def cons_f223(q):
        return Less(q, S(-1))

    cons223 = CustomConstraint(cons_f223)

    def cons_f224(c, d, f, e, q, a, b):
        return NonzeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

    cons224 = CustomConstraint(cons_f224)
    pattern334 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons221, cons222, cons223, cons224)
    rule334 = ReplacementRule(pattern334, lambda x, c, d, f, e, q, a, b : -(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f - S(2)*b*d*f + c*d*e + x*(c*(-S(2)*d*f + e**S(2)) + f*(S(2)*a*f - b*e)))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule334)


    def cons_f225(c, d, f, e, q, a):
        return NonzeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

    cons225 = CustomConstraint(cons_f225)
    pattern335 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons221, cons222, cons223, cons225)
    rule335 = ReplacementRule(pattern335, lambda x, c, d, f, e, q, a : -(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f + c*d*e + x*(S(2)*a*f**S(2) + c*(-S(2)*d*f + e**S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule335)


    def cons_f226(d, c, f, q, a):
        return NonzeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

    cons226 = CustomConstraint(cons_f226)
    pattern336 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons2, cons3, cons4, cons25, cons123, cons222, cons223, cons226)
    rule336 = ReplacementRule(pattern336, lambda x, c, d, f, q, a, b : (d + f*x**S(2))**(q + S(1))*(b*d - x*(a*f - c*d))/(S(2)*d*f*(q + S(1))) + (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**(q + S(1)), x)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule336)


    def cons_f227(q):
        return Not(PositiveIntegerQ(q))

    cons227 = CustomConstraint(cons_f227)

    def cons_f228(q):
        return Not(And(RationalQ(q), LessEqual(q, S(-1))))

    cons228 = CustomConstraint(cons_f228)
    pattern337 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons25, cons26, cons123, cons2, cons3, cons4, cons212, cons221, cons227, cons228, cons224)
    rule337 = ReplacementRule(pattern337, lambda x, c, d, f, e, q, a, b : (c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons25, cons26, cons123, cons2, cons4, cons212, cons221, cons227, cons228, cons225)
    rule338 = ReplacementRule(pattern338, lambda x, c, d, f, e, q, a : (S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons25, cons123, cons2, cons3, cons4, cons212, cons227, cons228, cons226)
    rule339 = ReplacementRule(pattern339, lambda x, c, d, f, q, a, b : (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**q, x)/(f*(S(2)*q + S(3))) + (d + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule339)


    def cons_f229(p, q):
        return RationalQ(p, q)

    cons229 = CustomConstraint(cons_f229)

    def cons_f230(q):
        return Greater(q, S(0))

    cons230 = CustomConstraint(cons_f230)
    pattern340 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, cons229, cons14, cons230)
    rule340 = ReplacementRule(pattern340, lambda x, c, d, f, e, q, p, a, b : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(b*e*q + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(S(2)*b*f*q + S(2)*c*e*(S(2)*p + q + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule340)

    pattern341 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons7, cons229, cons14, cons230)
    rule341 = ReplacementRule(pattern341, lambda x, c, d, f, q, p, a, b : (b + S(2)*c*x)*(d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*b*f*q*x + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule341)

    pattern342 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons221, cons229, cons14, cons230)
    rule342 = ReplacementRule(pattern342, lambda x, c, d, f, e, q, p, a : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(S(2)*p + q + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule342)


    def cons_f231(c, d, f, e, a, b):
        return NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))

    cons231 = CustomConstraint(cons_f231)

    def cons_f232(p, q):
        return Not(And(Not(IntegerQ(p)), IntegerQ(q), Less(q, S(-1))))

    cons232 = CustomConstraint(cons_f232)
    pattern343 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons212, cons7, cons221, cons11, cons14, cons231, cons232)
    rule343 = ReplacementRule(pattern343, lambda x, c, d, f, e, q, p, a, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d) + c*x*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + S(2)*c*(p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) - e*(p + q + S(2))*(-S(2)*a*c**S(2)*e - b**S(3)*f + b**S(2)*c*e - b*c*(-S(3)*a*f + c*d)) + x*(S(2)*f*(p + q + S(2))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule343)


    def cons_f233(d, c, f, a, b):
        return NonzeroQ(b**S(2)*d*f + (-a*f + c*d)**S(2))

    cons233 = CustomConstraint(cons_f233)
    pattern344 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons212, cons7, cons11, cons14, cons233, cons232)
    rule344 = ReplacementRule(pattern344, lambda x, c, d, f, q, p, a, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(b**S(3)*f + b*c*(-S(3)*a*f + c*d) + c*x*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*c*(p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + x*(-b*f*(p + S(1))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*f*(b**S(3)*f + b*c*(-S(3)*a*f + c*d))*(p + q + S(2))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule344)


    def cons_f234(c, d, f, e, a):
        return NonzeroQ(a*c*e**S(2) + (-a*f + c*d)**S(2))

    cons234 = CustomConstraint(cons_f234)
    pattern345 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons212, cons221, cons11, cons14, cons234, cons232)
    rule345 = ReplacementRule(pattern345, lambda x, c, d, f, e, q, p, a : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c**S(2)*e + c*x*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(S(2)*a*c**S(2)*e**S(2)*(p + q + S(2)) + c*f*x**S(2)*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + S(2)*q + S(5)) + S(2)*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + x*(S(4)*a*c**S(2)*e*f*(p + q + S(2)) + c*e*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + q + S(4))) - (-S(2)*a*c*f + S(2)*c**S(2)*d)*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule345)


    def cons_f235(p, q):
        return NonzeroQ(p + q)

    cons235 = CustomConstraint(cons_f235)

    def cons_f236(p, q):
        return NonzeroQ(S(2)*p + S(2)*q + S(1))

    cons236 = CustomConstraint(cons_f236)
    pattern346 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons212, cons7, cons221, cons11, cons37, cons235, cons236)
    rule346 = ReplacementRule(pattern346, lambda x, c, d, f, e, q, p, a, b : (a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(3)*p + S(2)*q) - c*e*(S(2)*p + q) + S(2)*c*f*x*(p + q))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(x**S(2)*(c*(p + q)*(-c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))) + f*(-S(2)*a*f + b*e)*(S(4)*p + S(2)*q + S(-1))) + p*(-p + S(1))*(-b*f + c*e)**S(2)) + x*(S(2)*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)*(-b*f + c*e) - (p + q)*(b*(c*(S(2)*p + q)*(-S(4)*d*f + e**S(2)) + f*(S(2)*p + S(2)*q + S(1))*(S(2)*a*f - b*e + S(2)*c*d)) + e*f*(-p + S(1))*(-S(4)*a*c + b**S(2)))) + (-p + S(1))*(S(2)*p + q)*(-a*e + b*d)*(-b*f + c*e) - (p + q)*(-a*(c*(S(2)*d*f - e**S(2)*(S(2)*p + q)) + f*(-S(2)*a*f + b*e)*(S(2)*p + S(2)*q + S(1))) + b**S(2)*d*f*(-p + S(1))), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule346)

    pattern347 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons212, cons7, cons11, cons37, cons235, cons236)
    rule347 = ReplacementRule(pattern347, lambda x, c, d, f, q, p, a, b : (d + f*x**S(2))**(q + S(1))*(b*(S(3)*p + S(2)*q) + S(2)*c*x*(p + q))*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-2))*Simp(b**S(2)*d*(p + S(-1))*(S(2)*p + q) + x**S(2)*(b**S(2)*f*p*(-p + S(1)) + S(2)*c*(p + q)*(-a*f*(S(4)*p + S(2)*q + S(-1)) + c*d*(S(2)*p + S(-1)))) - x*(S(2)*b*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d) - S(2)*b*(p + q)*(S(2)*c*d*(S(2)*p + q) - (a*f + c*d)*(S(2)*p + S(2)*q + S(1)))) - (p + q)*(-S(2)*a*(-a*f*(S(2)*p + S(2)*q + S(1)) + c*d) + b**S(2)*d*(-p + S(1))), x), x)/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule347)

    pattern348 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons212, cons221, cons11, cons37, cons235, cons236)
    rule348 = ReplacementRule(pattern348, lambda x, c, d, f, e, q, p, a : -c*(a + c*x**S(2))**(p + S(-1))*(e*(S(2)*p + q) - S(2)*f*x*(p + q))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(-a*c*e**S(2)*(-p + S(1))*(S(2)*p + q) + a*(p + q)*(-S(2)*a*f**S(2)*(S(2)*p + S(2)*q + S(1)) + c*(S(2)*d*f - e**S(2)*(S(2)*p + q))) + x**S(2)*(c**S(2)*e**S(2)*p*(-p + S(1)) - c*(p + q)*(S(2)*a*f**S(2)*(S(4)*p + S(2)*q + S(-1)) + c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))))) + x*(S(4)*a*c*e*f*(-p + S(1))*(p + q) + S(2)*c*e*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule348)

    pattern349 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, )
    def With349(x, c, d, f, e, a, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-a*c*f + b**S(2)*f - b*c*e + c**S(2)*d - x*(-b*c*f + c**S(2)*e))/(a + b*x + c*x**S(2)), x)/q + Int((a*f**S(2) - b*e*f - c*d*f + c*e**S(2) + x*(-b*f**S(2) + c*e*f))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule349 = ReplacementRule(pattern349, lambda x, c, d, f, e, a, b : With349(x, c, d, f, e, a, b))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons7, )
    def With350(x, c, d, f, a, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return -Int((-a*f**S(2) + b*f**S(2)*x + c*d*f)/(d + f*x**S(2)), x)/q + Int((-a*c*f + b**S(2)*f + b*c*f*x + c**S(2)*d)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule350 = ReplacementRule(pattern350, lambda x, c, d, f, a, b : With350(x, c, d, f, a, b))
    rubi.add(rule350)


    def cons_f237(c, f, e, b):
        return ZeroQ(-b*f + c*e)

    cons237 = CustomConstraint(cons_f237)
    pattern351 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, cons237)
    rule351 = ReplacementRule(pattern351, lambda x, c, d, f, e, a, b : -S(2)*e*Subst(Int(S(1)/(e*(-S(4)*a*f + b*e) - x**S(2)*(-a*e + b*d)), x), x, (e + S(2)*f*x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule351)


    def cons_f238(c, f, e, b):
        return NonzeroQ(-b*f + c*e)

    cons238 = CustomConstraint(cons_f238)
    pattern352 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, cons238, cons16, )
    def With352(x, c, d, f, e, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - S(2)*c*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule352 = ReplacementRule(pattern352, lambda x, c, d, f, e, a, b : With352(x, c, d, f, e, a, b))
    rubi.add(rule352)


    def cons_f239(c, a):
        return PosQ(-a*c)

    cons239 = CustomConstraint(cons_f239)
    pattern353 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons221, cons239)
    rule353 = ReplacementRule(pattern353, lambda x, c, d, f, e, a : Int(S(1)/((a - x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2) + Int(S(1)/((a + x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons7, cons16, )
    def With354(x, c, d, f, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule354 = ReplacementRule(pattern354, lambda x, c, d, f, a, b : With354(x, c, d, f, a, b))
    rubi.add(rule354)


    def cons_f240(c, a, b):
        return NegQ(-S(4)*a*c + b**S(2))

    cons240 = CustomConstraint(cons_f240)
    pattern355 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, cons238, cons240, )
    def With355(x, c, d, f, e, a, b):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d - q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule355 = ReplacementRule(pattern355, lambda x, c, d, f, e, a, b : With355(x, c, d, f, e, a, b))
    rubi.add(rule355)


    def cons_f241(c, a):
        return NegQ(-a*c)

    cons241 = CustomConstraint(cons_f241)
    pattern356 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons221, cons241, )
    def With356(x, c, d, f, e, a):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d + c*e*x - q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + c*e*x + q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule356 = ReplacementRule(pattern356, lambda x, c, d, f, e, a : With356(x, c, d, f, e, a))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons123, cons7, cons240, )
    def With357(x, c, d, f, a, b):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f - b*f*x + c*d - q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) + Int((-a*f - b*f*x + c*d + q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule357 = ReplacementRule(pattern357, lambda x, c, d, f, a, b : With357(x, c, d, f, a, b))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221)
    rule358 = ReplacementRule(pattern358, lambda x, c, d, f, e, a, b : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f + c*d + x*(-b*f + c*e))/(sqrt(a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule358)

    pattern359 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1))), x_), cons2, cons3, cons4, cons25, cons123, cons7)
    rule359 = ReplacementRule(pattern359, lambda x, c, d, f, a, b : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f - b*f*x + c*d)/((d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), x)/f)
    rubi.add(rule359)

    pattern360 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons2, cons4, cons25, cons26, cons123, cons221)
    rule360 = ReplacementRule(pattern360, lambda x, c, d, f, e, a : c*Int(S(1)/sqrt(a + c*x**S(2)), x)/f - Int((-a*f + c*d + c*e*x)/(sqrt(a + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule360)

    pattern361 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, )
    def With361(x, c, d, f, e, a, b):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*sqrt(d + e*x + f*x**S(2))), x)/sqrt(a + b*x + c*x**S(2))
    rule361 = ReplacementRule(pattern361, lambda x, c, d, f, e, a, b : With361(x, c, d, f, e, a, b))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons7, )
    def With362(x, c, d, f, a, b):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(d + f*x**S(2))*sqrt(b + S(2)*c*x + r)), x)/sqrt(a + b*x + c*x**S(2))
    rule362 = ReplacementRule(pattern362, lambda x, c, d, f, a, b : With362(x, c, d, f, a, b))
    rubi.add(rule362)


    def cons_f242(x, c, d, f, e, q, p, a, b):
        return FreeQ(List(a, b, c, d, e, f, p, q), x)

    cons242 = CustomConstraint(cons_f242)
    pattern363 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons6, cons212, cons242)
    rule363 = ReplacementRule(pattern363, lambda x, c, d, f, e, q, p, a, b : Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule363)


    def cons_f243(x, c, d, f, e, q, p, a):
        return FreeQ(List(a, c, d, e, f, p, q), x)

    cons243 = CustomConstraint(cons_f243)
    pattern364 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons6, cons212, cons243)
    rule364 = ReplacementRule(pattern364, lambda x, c, d, f, e, q, p, a : Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule364)

    pattern365 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons6, cons212, cons21, cons22)
    rule365 = ReplacementRule(pattern365, lambda x, c, d, f, e, q, u, p, a, b : Subst(Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule365)

    pattern366 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons6, cons212, cons21, cons22)
    rule366 = ReplacementRule(pattern366, lambda x, c, d, f, e, q, u, p, a : Subst(Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule366)


    def cons_f244(h, x):
        return FreeQ(h, x)

    cons244 = CustomConstraint(cons_f244)
    pattern367 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons212, cons208, cons209, cons210, cons211)
    rule367 = ReplacementRule(pattern367, lambda g, x, c, d, m, e, f, q, p, h, a, b : (c/f)**p*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule367)

    pattern368 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons212, cons208, cons209, cons53, cons213, cons214)
    rule368 = ReplacementRule(pattern368, lambda g, x, c, d, m, e, f, q, p, h, a, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons212, cons1)
    rule369 = ReplacementRule(pattern369, lambda g, x, c, d, m, e, f, q, p, h, a, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons212, cons1)
    rule370 = ReplacementRule(pattern370, lambda g, x, c, d, m, f, q, p, h, a, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(g + h*x)**m, x))
    rubi.add(rule370)


    def cons_f245(g, c, h, a, b):
        return ZeroQ(a*h**S(2) - b*g*h + c*g**S(2))

    cons245 = CustomConstraint(cons_f245)

    def cons_f246(g, c, d, f, e, h, a):
        return ZeroQ(a**S(2)*f*h**S(2) - a*c*e*g*h + c**S(2)*d*g**S(2))

    cons246 = CustomConstraint(cons_f246)
    pattern371 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons245, cons246, cons59)
    rule371 = ReplacementRule(pattern371, lambda g, x, c, d, m, e, f, p, h, a, b : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule371)


    def cons_f247(c, a, g, h):
        return ZeroQ(a*h**S(2) + c*g**S(2))

    cons247 = CustomConstraint(cons_f247)
    pattern372 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons247, cons246, cons59)
    rule372 = ReplacementRule(pattern372, lambda g, x, c, d, m, e, f, p, h, a : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule372)


    def cons_f248(g, c, d, f, h, a):
        return ZeroQ(a**S(2)*f*h**S(2) + c**S(2)*d*g**S(2))

    cons248 = CustomConstraint(cons_f248)
    pattern373 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons6, cons245, cons248, cons59)
    rule373 = ReplacementRule(pattern373, lambda g, x, c, d, m, f, p, h, a, b : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons6, cons247, cons248, cons59)
    rule374 = ReplacementRule(pattern374, lambda g, x, c, d, m, f, p, h, a : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule374)


    def cons_f249(c, f, e, a, b):
        return ZeroQ(a*f**S(2) - b*e*f + c*e**S(2))

    cons249 = CustomConstraint(cons_f249)
    pattern375 = Pattern(Integral(x_**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons26, cons123, cons212, cons7, cons249, cons55)
    rule375 = ReplacementRule(pattern375, lambda x, c, f, e, q, p, a, b : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule375)


    def cons_f250(c, a, f, e):
        return ZeroQ(a*f**S(2) + c*e**S(2))

    cons250 = CustomConstraint(cons_f250)
    pattern376 = Pattern(Integral(x_**WC('p', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons26, cons123, cons212, cons250, cons55)
    rule376 = ReplacementRule(pattern376, lambda x, c, f, e, q, p, a : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule376)


    def cons_f251(g, c, m, f, e, p, h, b):
        return ZeroQ(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

    cons251 = CustomConstraint(cons_f251)

    def cons_f252(g, c, d, f, m, p, h, a, b):
        return ZeroQ(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

    cons252 = CustomConstraint(cons_f252)
    pattern377 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons251, cons252, cons73)
    rule377 = ReplacementRule(pattern377, lambda g, x, c, d, m, e, f, p, h, a, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule377)


    def cons_f253(g, c, f, e, m, p, h):
        return ZeroQ(c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

    cons253 = CustomConstraint(cons_f253)

    def cons_f254(c, d, m, f, p, h, a):
        return ZeroQ(h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

    cons254 = CustomConstraint(cons_f254)
    pattern378 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons253, cons254, cons73)
    rule378 = ReplacementRule(pattern378, lambda g, x, c, d, m, e, f, p, h, a : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule378)


    def cons_f255(g, c, m, f, p, h, b):
        return ZeroQ(b*f*h*(m + p + S(2)) + S(2)*c*f*g*(p + S(1)))

    cons255 = CustomConstraint(cons_f255)
    pattern379 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons255, cons252, cons73)
    rule379 = ReplacementRule(pattern379, lambda g, x, c, d, m, f, p, h, a, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule379)


    def cons_f256(p, m):
        return Or(IntegersQ(m, p), PositiveIntegerQ(p))

    cons256 = CustomConstraint(cons_f256)
    pattern380 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons7, cons221, cons256)
    rule380 = ReplacementRule(pattern380, lambda g, x, c, d, m, e, f, p, h, a, b : Int(ExpandIntegrand((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule380)

    pattern381 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons221, cons256)
    rule381 = ReplacementRule(pattern381, lambda g, x, c, d, m, e, f, p, h, a : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule381)

    pattern382 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons7, cons256)
    rule382 = ReplacementRule(pattern382, lambda g, x, c, d, m, f, p, h, a, b : Int(ExpandIntegrand((d + f*x**S(2))*(g + h*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule382)

    pattern383 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons27, cons256)
    rule383 = ReplacementRule(pattern383, lambda g, x, c, d, m, f, p, h, a : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + f*x**S(2))*(g + h*x)**m, x), x))
    rubi.add(rule383)


    def cons_f257(g, c, h, a, b):
        return NonzeroQ(a*h**S(2) - b*g*h + c*g**S(2))

    cons257 = CustomConstraint(cons_f257)
    pattern384 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons7, cons221, cons48, cons52, cons257)
    rule384 = ReplacementRule(pattern384, lambda g, x, c, d, f, e, m, p, h, a, b : (g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(f*g**S(2)*(p + S(1)) - h*(-d*h*(m + p + S(2)) + e*g*(p + S(1)))) + h*(m + S(1))*(a*e*h - a*f*g + c*d*g) - x*(c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule384)


    def cons_f258(c, a, g, h):
        return NonzeroQ(a*h**S(2) + c*g**S(2))

    cons258 = CustomConstraint(cons_f258)
    pattern385 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons221, cons48, cons52, cons258)
    rule385 = ReplacementRule(pattern385, lambda g, x, c, d, f, e, m, p, h, a : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(a*e*h - a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule385)

    pattern386 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons6, cons7, cons48, cons52, cons257)
    rule386 = ReplacementRule(pattern386, lambda g, x, c, d, f, m, p, h, a, b : (g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))*(a + b*x + c*x**S(2))**(p + S(1))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(d*h**S(2)*(m + p + S(2)) + f*g**S(2)*(p + S(1))) + h*(m + S(1))*(-a*f*g + c*d*g) - x*(c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule386)

    pattern387 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons6, cons48, cons52, cons258)
    rule387 = ReplacementRule(pattern387, lambda g, x, c, d, f, m, p, h, a : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(-a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule387)

    pattern388 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons257)
    rule388 = ReplacementRule(pattern388, lambda g, x, c, d, f, e, h, a, b : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((a*e*h - a*f*g - b*d*h + c*d*g + x*(a*f*h - b*f*g - c*d*h + c*e*g))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule388)

    pattern389 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons258)
    rule389 = ReplacementRule(pattern389, lambda g, x, d, c, f, e, h, a : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((a*e*h - a*f*g + c*d*g + x*(a*f*h - c*d*h + c*e*g))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule389)

    pattern390 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons257)
    rule390 = ReplacementRule(pattern390, lambda g, x, c, d, f, h, a, b : (d*h**S(2) + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((-a*f*g - b*d*h + c*d*g - x*(-a*f*h + b*f*g + c*d*h))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule390)

    pattern391 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(g_ + x_*WC('h', S(1)))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons258)
    rule391 = ReplacementRule(pattern391, lambda g, x, c, d, f, h, a : (d*h**S(2) + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((-a*f*g + c*d*g - x*(-a*f*h + c*d*h))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule391)

    pattern392 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons36, cons14, cons46)
    rule392 = ReplacementRule(pattern392, lambda g, x, c, d, f, e, m, p, h, a, b : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f - S(2)*a*c*e + b*c*d + x*(c*(-b*e + S(2)*c*d) + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(c*(S(2)*p + S(3))*(-b*e + S(2)*c*d) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f - S(2)*a*c*e + b*c*d) + h*x*(c*(-b*e + S(2)*c*d)*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule392)

    pattern393 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons36, cons14, cons46)
    rule393 = ReplacementRule(pattern393, lambda g, x, d, c, f, e, m, p, h, a : (a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(a*e - x*(-a*f + c*d))/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*(e*h*m + f*g) - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule393)

    pattern394 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons36, cons14, cons46)
    rule394 = ReplacementRule(pattern394, lambda g, x, c, d, f, m, p, h, a, b : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f + b*c*d + x*(S(2)*c**S(2)*d + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(S(2)*c**S(2)*d*(S(2)*p + S(3)) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f + b*c*d) + h*x*(S(2)*c**S(2)*d*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule394)

    pattern395 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons36, cons14, cons46)
    rule395 = ReplacementRule(pattern395, lambda g, x, c, d, f, m, p, h, a : -x*(a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(-a*f + c*d)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*f*g - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule395)


    def cons_f259(g, c, h, a, b):
        return NonzeroQ(c*g**S(2) - h*(-a*h + b*g))

    cons259 = CustomConstraint(cons_f259)
    pattern396 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons7, cons221, cons11, cons14, cons259)
    rule396 = ReplacementRule(pattern396, lambda g, x, c, d, m, e, f, p, h, a, b : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(-S(2)*a*e*h + S(2)*a*f*g + b*d*h + b*e*g)) - (-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) - h*(-(-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + (p + S(1))*(c*g**S(2) - h*(-a*h + b*g))*(S(2)*a*f - b*e + S(2)*c*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule396)

    pattern397 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons221, cons11, cons14, cons258)
    rule397 = ReplacementRule(pattern397, lambda g, x, c, d, m, e, f, p, h, a : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*c*e*g - a*h*(-a*f + c*d) - c*x*(a*e*h - a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(g*(p + S(2))*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g) + h*x*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g)*(m + S(2)*p + S(4)) - h*(a*c*e*g - a*h*(-a*f + c*d))*(m + p + S(2)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule397)

    pattern398 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons7, cons11, cons14, cons259)
    rule398 = ReplacementRule(pattern398, lambda g, x, c, d, m, f, p, h, a, b : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*d*(-b*h + S(2)*c*g) - x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(S(2)*a*f*g + b*d*h)) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) - h*(-b*d*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + S(2)*(p + S(1))*(a*f + c*d)*(c*g**S(2) - h*(-a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule398)

    pattern399 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons27, cons11, cons14, cons258)
    rule399 = ReplacementRule(pattern399, lambda g, x, d, c, m, f, p, h, a : -(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*h*(-a*f + c*d) + c*x*(-a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(a*h**S(2)*(-a*f + c*d)*(m + p + S(2)) + g*(p + S(2))*(-a*c*f*g + c**S(2)*d*g) + h*x*(-a*c*f*g + c**S(2)*d*g)*(m + S(2)*p + S(4)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule399)

    pattern400 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons7, cons221, cons34)
    rule400 = ReplacementRule(pattern400, lambda g, x, c, d, m, e, f, p, h, a, b : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons221, cons34)
    rule401 = ReplacementRule(pattern401, lambda g, x, c, d, m, e, f, p, h, a : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule401)

    pattern402 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons7, cons34)
    rule402 = ReplacementRule(pattern402, lambda g, x, c, d, m, f, p, h, a, b : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons34)
    rule403 = ReplacementRule(pattern403, lambda g, x, c, d, m, f, p, h, a : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule403)

    pattern404 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons7, cons221, cons73)
    rule404 = ReplacementRule(pattern404, lambda g, x, c, d, m, e, f, p, h, a, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))) + x*(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1)))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons221, cons73)
    rule405 = ReplacementRule(pattern405, lambda g, x, c, d, m, e, f, p, h, a : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(c*x*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule405)

    pattern406 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons7, cons73)
    rule406 = ReplacementRule(pattern406, lambda g, x, c, d, m, f, p, h, a, b : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + f*x*(b*h*(m + p + S(2)) + S(2)*c*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)))*(g_ + x_*WC('h', S(1)))**WC('m', S(1)), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons27, cons6, cons73)
    rule407 = ReplacementRule(pattern407, lambda g, x, c, d, m, f, p, h, a : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(S(2)*c*f*g*x*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule407)


    def cons_f260(p, q):
        return IntegersQ(p, q)

    cons260 = CustomConstraint(cons_f260)
    pattern408 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons260, cons12)
    rule408 = ReplacementRule(pattern408, lambda g, x, c, d, f, e, q, p, h, a, b : Int(ExpandIntegrand((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule408)


    def cons_f261(p, q):
        return Or(Greater(p, S(0)), Greater(q, S(0)))

    cons261 = CustomConstraint(cons_f261)
    pattern409 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons260, cons261)
    rule409 = ReplacementRule(pattern409, lambda g, x, c, d, f, e, q, p, h, a : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule409)

    pattern410 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons229, cons14, cons230)
    rule410 = ReplacementRule(pattern410, lambda g, x, c, d, f, e, q, p, h, a, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + e*q*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)) + x*(-e*(b*h - S(2)*c*g)*(S(2)*p + q + S(3)) + S(2)*f*q*(-S(2)*a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule410)

    pattern411 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons229, cons14, cons230)
    rule411 = ReplacementRule(pattern411, lambda g, x, c, d, f, e, q, p, h, a : (a + c*x**S(2))**(p + S(1))*(a*h - c*g*x)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-a*e*h*q + c*d*g*(S(2)*p + S(3)) + c*f*g*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(-S(2)*a*f*h*q + c*e*g*(S(2)*p + q + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule411)

    pattern412 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons229, cons14, cons230)
    rule412 = ReplacementRule(pattern412, lambda g, x, c, d, f, q, p, h, a, b : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + S(2)*f*q*x*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule412)

    pattern413 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons212, cons7, cons221, cons11, cons14, cons231, cons232)
    rule413 = ReplacementRule(pattern413, lambda g, x, c, d, f, e, q, p, h, a, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + c*x*(g*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - h*(a*b*f - S(2)*a*c*e + b*c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)) - e*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - x*(S(2)*f*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d))) + (p + S(1))*(b*h - S(2)*c*g)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule413)

    pattern414 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons212, cons221, cons11, cons14, cons234, cons232)
    rule414 = ReplacementRule(pattern414, lambda g, x, c, d, f, e, q, p, h, a : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d) + c*x*(S(2)*a*c*e*h + g*(-S(2)*a*c*f + S(2)*c**S(2)*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + S(2)*q + S(5)) - S(2)*c*g*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) - e*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2)) - x*(c*e*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + q + S(4)) + S(2)*f*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d)), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule414)

    pattern415 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons212, cons7, cons11, cons14, cons233, cons232)
    rule415 = ReplacementRule(pattern415, lambda g, x, c, d, f, q, p, h, a, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*c*g*(a*f + c*d) + c*x*(g*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - h*(a*b*f + b*c*d)) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) - x*(-b*f*(p + S(1))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) + S(2)*f*(-b*c*g*(a*f + c*d) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b*h - S(2)*c*g)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule415)


    def cons_f262(p, q):
        return NonzeroQ(p + q + S(1))

    cons262 = CustomConstraint(cons_f262)
    pattern416 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons212, cons7, cons221, cons11, cons12, cons262)
    rule416 = ReplacementRule(pattern416, lambda g, x, c, d, f, e, q, p, h, a, b : h*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-a*e + b*d) + x**S(2)*(c*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-b*f + c*e)) + x*(b*(e*h - S(2)*f*g)*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule416)

    pattern417 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons212, cons221, cons11, cons12, cons262)
    rule417 = ReplacementRule(pattern417, lambda g, x, c, d, f, e, q, p, h, a : h*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) + Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*e*h*p - a*(e*h - S(2)*f*g)*(p + q + S(1)) - S(2)*h*p*x*(-a*f + c*d) - x**S(2)*(c*e*h*p + c*(e*h - S(2)*f*g)*(p + q + S(1))), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule417)

    pattern418 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons212, cons7, cons11, cons12, cons262)
    rule418 = ReplacementRule(pattern418, lambda g, x, c, d, f, q, p, h, a, b : h*(d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*f*(p + q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*f*g*(p + q + S(1)) + b*d*h*p + x**S(2)*(-b*f*h*p - S(2)*c*f*g*(p + q + S(1))) + x*(-S(2)*b*f*g*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule418)

    pattern419 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, )
    def With419(g, x, c, d, f, e, h, a, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(-a*b*f*h + a*c*e*h - a*c*f*g + b**S(2)*f*g - b*c*e*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(a + b*x + c*x**S(2)), x)/q + Int(Simp(a*f**S(2)*g + b*d*f*h - b*e*f*g - c*d*e*h - c*d*f*g + c*e**S(2)*g - f*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule419 = ReplacementRule(pattern419, lambda g, x, c, d, f, e, h, a, b : With419(g, x, c, d, f, e, h, a, b))
    rubi.add(rule419)

    pattern420 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, )
    def With420(g, x, c, d, f, h, a, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(a*f**S(2)*g + b*d*f*h - c*d*f*g - f*x*(-a*f*h + b*f*g + c*d*h), x)/(d + f*x**S(2)), x)/q + Int(Simp(-a*b*f*h - a*c*f*g + b**S(2)*f*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h), x)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule420 = ReplacementRule(pattern420, lambda g, x, c, d, f, h, a, b : With420(g, x, c, d, f, h, a, b))
    rubi.add(rule420)


    def cons_f263(c, a):
        return PositiveQ(a*c)

    cons263 = CustomConstraint(cons_f263)
    pattern421 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons263)
    rule421 = ReplacementRule(pattern421, lambda g, x, c, d, f, h, a : g*Int(S(1)/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x) + h*Int(x/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x))
    rubi.add(rule421)


    def cons_f264(c, a):
        return Not(PositiveQ(a*c))

    cons264 = CustomConstraint(cons_f264)
    pattern422 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons264, )
    def With422(g, x, c, d, f, h, a):
        q = Rt(-a*c, S(2))
        return -(c*g - h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(c*x + q)), x)/(S(2)*q) - (c*g + h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(-c*x + q)), x)/(S(2)*q)
    rule422 = ReplacementRule(pattern422, lambda g, x, c, d, f, h, a : With422(g, x, c, d, f, h, a))
    rubi.add(rule422)


    def cons_f265(e, f, g, h):
        return ZeroQ(e*h - S(2)*f*g)

    cons265 = CustomConstraint(cons_f265)
    pattern423 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons237, cons265)
    rule423 = ReplacementRule(pattern423, lambda g, x, c, d, f, e, h, a, b : -S(2)*g*Subst(Int(S(1)/(-a*e + b*d - b*x**S(2)), x), x, sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule423)


    def cons_f266(e, f, g, h):
        return NonzeroQ(e*h - S(2)*f*g)

    cons266 = CustomConstraint(cons_f266)
    pattern424 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons237, cons266)
    rule424 = ReplacementRule(pattern424, lambda g, x, d, c, f, e, h, a, b : h*Int((e + S(2)*f*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f) - (e*h - S(2)*f*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f))
    rubi.add(rule424)

    pattern425 = Pattern(Integral(x_/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons7, cons221, cons209)
    rule425 = ReplacementRule(pattern425, lambda x, c, d, f, e, a, b : -S(2)*e*Subst(Int((-d*x**S(2) + S(1))/(-b*f + c*e + d**S(2)*x**S(4)*(-b*f + c*e) - e*x**S(2)*(S(2)*a*f - b*e + S(2)*c*d)), x), x, (S(1) + x*(e + sqrt(-S(4)*d*f + e**S(2)))/(S(2)*d))/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule425)


    def cons_f267(e, d, g, h):
        return ZeroQ(S(2)*d*h - e*g)

    cons267 = CustomConstraint(cons_f267)
    pattern426 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons209, cons267)
    rule426 = ReplacementRule(pattern426, lambda g, x, c, d, f, e, h, a, b : g*Subst(Int(S(1)/(a + x**S(2)*(-a*f + c*d)), x), x, x/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule426)


    def cons_f268(e, d, g, h):
        return NonzeroQ(S(2)*d*h - e*g)

    cons268 = CustomConstraint(cons_f268)
    pattern427 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons209, cons268)
    rule427 = ReplacementRule(pattern427, lambda g, x, c, d, f, e, h, a, b : h*Int((S(2)*d + e*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e - (S(2)*d*h - e*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e)
    rubi.add(rule427)


    def cons_f269(g, d, c, f, e, h, a, b):
        return ZeroQ(g**S(2)*(-b*f + c*e) - S(2)*g*h*(-a*f + c*d) + h**S(2)*(-a*e + b*d))

    cons269 = CustomConstraint(cons_f269)
    pattern428 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons191, cons269)
    rule428 = ReplacementRule(pattern428, lambda g, x, c, d, f, e, h, a, b : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g) - x**S(2)*(-a*e + b*d), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule428)


    def cons_f270(g, c, d, f, e, h, a):
        return ZeroQ(a*e*h**S(2) - c*e*g**S(2) + S(2)*g*h*(-a*f + c*d))

    cons270 = CustomConstraint(cons_f270)
    pattern429 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons270)
    rule429 = ReplacementRule(pattern429, lambda g, x, c, d, f, e, h, a : -S(2)*a*g*h*Subst(Int(S(1)/Simp(S(2)*a**S(2)*c*g*h + a*e*x**S(2), x), x), x, Simp(a*h - c*g*x, x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule429)


    def cons_f271(g, d, c, f, h, a, b):
        return ZeroQ(b*d*h**S(2) - b*f*g**S(2) - S(2)*g*h*(-a*f + c*d))

    cons271 = CustomConstraint(cons_f271)
    pattern430 = Pattern(Integral((g_ + x_*WC('h', S(1)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons271)
    rule430 = ReplacementRule(pattern430, lambda g, x, c, d, f, h, a, b : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(-b*d*x**S(2) + g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + f*x**S(2))))
    rubi.add(rule430)

    pattern431 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons16, )
    def With431(g, x, d, c, f, e, h, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule431 = ReplacementRule(pattern431, lambda g, x, d, c, f, e, h, a, b : With431(g, x, d, c, f, e, h, a, b))
    rubi.add(rule431)

    pattern432 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons239, )
    def With432(g, x, d, c, f, e, h, a):
        q = Rt(-a*c, S(2))
        return (-c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x + q)*sqrt(d + e*x + f*x**S(2))), x) + (c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x - q)*sqrt(d + e*x + f*x**S(2))), x)
    rule432 = ReplacementRule(pattern432, lambda g, x, d, c, f, e, h, a : With432(g, x, d, c, f, e, h, a))
    rubi.add(rule432)

    pattern433 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons16, )
    def With433(g, x, c, d, f, h, a, b):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule433 = ReplacementRule(pattern433, lambda g, x, c, d, f, h, a, b : With433(g, x, c, d, f, h, a, b))
    rubi.add(rule433)

    pattern434 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, cons191, cons240, )
    def With434(g, x, c, d, f, e, h, a, b):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-g*(-a*f + c*d - q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d + q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-g*(-a*f + c*d + q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d - q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule434 = ReplacementRule(pattern434, lambda g, x, c, d, f, e, h, a, b : With434(g, x, c, d, f, e, h, a, b))
    rubi.add(rule434)

    pattern435 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons221, cons241, )
    def With435(g, x, d, c, f, e, h, a):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-a*e*h - g*(-a*f + c*d - q) + x*(-c*e*g + h*(-a*f + c*d + q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-a*e*h - g*(-a*f + c*d + q) + x*(-c*e*g + h*(-a*f + c*d - q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule435 = ReplacementRule(pattern435, lambda g, x, d, c, f, e, h, a : With435(g, x, d, c, f, e, h, a))
    rubi.add(rule435)

    pattern436 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, cons240, )
    def With436(g, x, c, d, f, h, a, b):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(b*d*h - g*(-a*f + c*d - q) + x*(b*f*g + h*(-a*f + c*d + q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) - Int(Simp(b*d*h - g*(-a*f + c*d + q) + x*(b*f*g + h*(-a*f + c*d - q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule436 = ReplacementRule(pattern436, lambda g, x, c, d, f, h, a, b : With436(g, x, c, d, f, h, a, b))
    rubi.add(rule436)

    pattern437 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons7, cons221, )
    def With437(g, x, c, d, f, e, h, a, b):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f + e**S(2), S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)), x)/(sqrt(a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2)))
    rule437 = ReplacementRule(pattern437, lambda g, x, c, d, f, e, h, a, b : With437(g, x, c, d, f, e, h, a, b))
    rubi.add(rule437)

    pattern438 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons124, cons244, cons7, )
    def With438(g, x, c, d, f, h, a, b):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f, S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)), x)/(sqrt(d + f*x**S(2))*sqrt(a + b*x + c*x**S(2)))
    rule438 = ReplacementRule(pattern438, lambda g, x, c, d, f, h, a, b : With438(g, x, c, d, f, h, a, b))
    rubi.add(rule438)


    def cons_f272(c, d, f, a, b):
        return ZeroQ(c**S(2)*d - f*(-S(3)*a*c + b**S(2)))

    cons272 = CustomConstraint(cons_f272)

    def cons_f273(g, c, h, a, b):
        return ZeroQ(S(9)*a*c*h**S(2) - S(2)*b**S(2)*h**S(2) - b*c*g*h + c**S(2)*g**S(2))

    cons273 = CustomConstraint(cons_f273)

    def cons_f274(c, g, h, b):
        return PositiveQ(-S(9)*c*h**S(2)/(-b*h + S(2)*c*g)**S(2))

    cons274 = CustomConstraint(cons_f274)
    pattern439 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons237, cons272, cons273, cons274, )
    def With439(g, x, c, d, f, e, h, a, b):
        q = S(3)**(S(2)/3)*(-c*h**S(2)/(-b*h + S(2)*c*g)**S(2))**(S(1)/3)
        return sqrt(S(3))*h*q*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3)/(S(3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3)) + sqrt(S(3))/S(3))/f - S(3)*h*q*log((-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3) + S(2)**(S(1)/3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3))/(S(2)*f) + h*q*log(d + e*x + f*x**S(2))/(S(2)*f)
    rule439 = ReplacementRule(pattern439, lambda g, x, c, d, f, e, h, a, b : With439(g, x, c, d, f, e, h, a, b))
    rubi.add(rule439)


    def cons_f275(c, a, d, f):
        return ZeroQ(S(3)*a*f + c*d)

    cons275 = CustomConstraint(cons_f275)

    def cons_f276(c, a, g, h):
        return ZeroQ(S(9)*a*h**S(2) + c*g**S(2))

    cons276 = CustomConstraint(cons_f276)
    pattern440 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons275, cons276, cons69)
    rule440 = ReplacementRule(pattern440, lambda g, x, c, d, f, h, a : S(2)**(S(1)/3)*sqrt(S(3))*h*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(S(1) - S(3)*h*x/g)**(S(2)/3)/(S(3)*(S(1) + S(3)*h*x/g)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*a**(S(1)/3)*f) + S(2)**(S(1)/3)*h*log(d + f*x**S(2))/(S(4)*a**(S(1)/3)*f) - S(3)*S(2)**(S(1)/3)*h*log((S(1) - S(3)*h*x/g)**(S(2)/3) + S(2)**(S(1)/3)*(S(1) + S(3)*h*x/g)**(S(1)/3))/(S(4)*a**(S(1)/3)*f))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons237, cons272, cons273, cons117, )
    def With441(g, x, c, d, f, e, h, a, b):
        q = -c/(-S(4)*a*c + b**S(2))
        return (q*(a + b*x + c*x**S(2)))**(S(1)/3)*Int((g + h*x)/((d + e*x + f*x**S(2))*(a*q + b*q*x + c*q*x**S(2))**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule441 = ReplacementRule(pattern441, lambda g, x, c, d, f, e, h, a, b : With441(g, x, c, d, f, e, h, a, b))
    rubi.add(rule441)


    def cons_f277(a):
        return Not(PositiveQ(a))

    cons277 = CustomConstraint(cons_f277)
    pattern442 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons2, cons4, cons25, cons123, cons124, cons244, cons275, cons276, cons277)
    rule442 = ReplacementRule(pattern442, lambda g, x, c, d, f, h, a : (S(1) + c*x**S(2)/a)**(S(1)/3)*Int((g + h*x)/((S(1) + c*x**S(2)/a)**(S(1)/3)*(d + f*x**S(2))), x)/(a + c*x**S(2))**(S(1)/3))
    rubi.add(rule442)


    def cons_f278(g, x, c, d, f, e, q, p, h, a, b):
        return FreeQ(List(a, b, c, d, e, f, g, h, p, q), x)

    cons278 = CustomConstraint(cons_f278)
    pattern443 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons212, cons278)
    rule443 = ReplacementRule(pattern443, lambda g, x, c, d, f, e, q, p, h, a, b : Int((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule443)


    def cons_f279(g, x, c, d, f, e, q, p, h, a):
        return FreeQ(List(a, c, d, e, f, g, h, p, q), x)

    cons279 = CustomConstraint(cons_f279)
    pattern444 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons6, cons212, cons279)
    rule444 = ReplacementRule(pattern444, lambda g, x, d, c, f, e, q, p, h, a : Int((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule444)

    pattern445 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons212, cons21, cons22)
    rule445 = ReplacementRule(pattern445, lambda g, x, c, d, m, e, f, q, u, p, h, a, b : Subst(Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule445)

    pattern446 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons124, cons244, cons27, cons6, cons212, cons21, cons22)
    rule446 = ReplacementRule(pattern446, lambda g, x, c, d, m, e, f, q, u, p, h, a : Subst(Int((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule446)


    def cons_f280(z, x):
        return LinearQ(z, x)

    cons280 = CustomConstraint(cons_f280)

    def cons_f281(v, u, x):
        return QuadraticQ(List(u, v), x)

    cons281 = CustomConstraint(cons_f281)

    def cons_f282(z, v, u, x):
        return Not(And(LinearMatchQ(z, x), QuadraticMatchQ(List(u, v), x)))

    cons282 = CustomConstraint(cons_f282)
    pattern447 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*z_**WC('m', S(1)), x_), cons27, cons6, cons212, cons280, cons281, cons282)
    rule447 = ReplacementRule(pattern447, lambda z, x, m, q, u, p, v : Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**m, x))
    rubi.add(rule447)


    def cons_f283(i, x):
        return FreeQ(i, x)

    cons283 = CustomConstraint(cons_f283)
    pattern448 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons283, cons27, cons134, cons6, cons212, cons150, cons197, cons198)
    rule448 = ReplacementRule(pattern448, lambda g, x, c, d, m, e, q, f, p, h, a, i, b, n : Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule448)

    pattern449 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons283, cons27, cons134, cons6, cons212, cons8, cons118)
    rule449 = ReplacementRule(pattern449, lambda g, x, d, c, m, e, f, q, p, h, a, i, b, n : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(h + i*x)**q*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule449)

    pattern450 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons124, cons244, cons283, cons27, cons134, cons6, cons212, cons150, cons197)
    rule450 = ReplacementRule(pattern450, lambda g, x, c, d, m, e, q, f, p, h, a, i, b, n : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule450)


    def cons_f284(A, x):
        return FreeQ(A, x)

    cons284 = CustomConstraint(cons_f284)

    def cons_f285(B, x):
        return FreeQ(B, x)

    cons285 = CustomConstraint(cons_f285)

    def cons_f286(C, x):
        return FreeQ(C, x)

    cons286 = CustomConstraint(cons_f286)
    pattern451 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons208, cons209, cons210, cons211)
    rule451 = ReplacementRule(pattern451, lambda x, c, d, f, e, q, b, C, p, a, B, A : (c/f)**p*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule451)

    pattern452 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons208, cons209, cons210, cons211)
    rule452 = ReplacementRule(pattern452, lambda x, c, d, f, e, q, C, p, a, A, b : (c/f)**p*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule452)

    pattern453 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons208, cons209, cons53, cons213, cons214)
    rule453 = ReplacementRule(pattern453, lambda x, c, d, f, e, q, b, C, p, a, B, A : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule453)

    pattern454 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons208, cons209, cons53, cons213, cons214)
    rule454 = ReplacementRule(pattern454, lambda x, c, d, f, e, q, C, p, a, A, b : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule454)

    pattern455 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons1)
    rule455 = ReplacementRule(pattern455, lambda x, c, d, f, e, q, b, C, p, a, B, A : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule455)

    pattern456 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons1)
    rule456 = ReplacementRule(pattern456, lambda x, c, d, f, e, q, C, p, a, A, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule456)

    pattern457 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons6, cons212, cons1)
    rule457 = ReplacementRule(pattern457, lambda x, c, d, f, q, b, C, p, a, B, A : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(A + B*x + C*x**S(2)), x))
    rubi.add(rule457)

    pattern458 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons6, cons212, cons1)
    rule458 = ReplacementRule(pattern458, lambda x, c, d, f, q, C, p, a, A, b : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule458)

    pattern459 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons7, cons221, cons260, cons12)
    rule459 = ReplacementRule(pattern459, lambda x, c, d, f, e, q, b, C, p, a, B, A : Int(ExpandIntegrand((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule459)

    pattern460 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons7, cons221, cons260, cons12)
    rule460 = ReplacementRule(pattern460, lambda x, c, d, f, e, q, C, p, a, A, b : Int(ExpandIntegrand((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule460)

    pattern461 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons221, cons260, cons261)
    rule461 = ReplacementRule(pattern461, lambda x, c, d, f, e, q, C, p, a, B, A : Int(ExpandIntegrand((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule461)

    pattern462 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons221, cons260, cons261)
    rule462 = ReplacementRule(pattern462, lambda x, c, d, f, e, q, C, p, a, A : Int(ExpandIntegrand((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule462)

    pattern463 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons7, cons221, cons229, cons14, cons230)
    rule463 = ReplacementRule(pattern463, lambda x, c, d, f, e, q, b, C, p, a, B, A : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + e*q*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))) + x*(-e*(C*(S(2)*a*c*(q + S(1)) - b**S(2)*(p + q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + q + S(3))) + S(2)*f*q*(A*b*c - S(2)*B*a*c + C*a*b)), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule463)

    pattern464 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons7, cons221, cons229, cons14, cons230)
    rule464 = ReplacementRule(pattern464, lambda x, c, d, f, e, q, C, p, a, A, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*(b*e*q + S(2)*c*d*(S(2)*p + S(3))) - C*(-a*b*e*q + S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*c*(b*f*q + c*e*(S(2)*p + q + S(3))) + C*(S(2)*a*b*f*q - S(2)*a*c*e*(q + S(1)) + b**S(2)*e*(p + q + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule464)

    pattern465 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons221, cons229, cons14, cons230)
    rule465 = ReplacementRule(pattern465, lambda x, c, d, f, e, q, C, p, a, B, A : (a + c*x**S(2))**(p + S(1))*(B*a - x*(A*c - C*a))*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - a*(B*e*q + C*d) - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - a*(S(2)*B*f*q + C*e*(q + S(1)))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule465)

    pattern466 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons221, cons229, cons14, cons230)
    rule466 = ReplacementRule(pattern466, lambda x, c, d, f, e, q, C, p, a, A : -x*(a + c*x**S(2))**(p + S(1))*(A*c - C*a)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - C*a*d - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - C*a*e*(q + S(1))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule466)

    pattern467 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons7, cons229, cons14, cons230)
    rule467 = ReplacementRule(pattern467, lambda x, c, d, f, q, b, C, p, a, B, A : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + S(2)*f*q*x*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule467)

    pattern468 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons7, cons229, cons14, cons230)
    rule468 = ReplacementRule(pattern468, lambda x, c, d, f, q, C, p, a, A, b : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*A*c**S(2)*d*(S(2)*p + S(3)) - C*(S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*b*c*f*q + S(2)*C*a*b*f*q), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule468)

    pattern469 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons212, cons7, cons221, cons11, cons14, cons231, cons232)
    rule469 = ReplacementRule(pattern469, lambda x, c, d, f, e, q, b, C, p, a, B, A : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - B*(a*b*f - S(2)*a*c*e + b*c*d) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)) - e*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e))) + (p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule469)

    pattern470 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons212, cons7, cons221, cons11, cons14, cons231, cons232)
    rule470 = ReplacementRule(pattern470, lambda x, c, d, f, e, q, C, p, a, A, b : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)) - e*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule470)

    pattern471 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons212, cons221, cons11, cons14, cons234, cons232)
    rule471 = ReplacementRule(pattern471, lambda x, c, d, f, e, q, C, p, a, B, A : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*B*a*c*e - S(2)*C*a*(-a*f + c*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - e*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2)) - x*(c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + q + S(4)) + S(2)*f*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule471)

    pattern472 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons212, cons221, cons11, cons14, cons234, cons232)
    rule472 = ReplacementRule(pattern472, lambda x, c, d, f, e, q, C, p, a, A : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) - S(2)*C*a*(-a*f + c*d)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-S(2)*a*c*e**S(2)*(A*c - C*a)*(p + q + S(2)) - c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - x*(S(4)*a*c*e*f*(A*c - C*a)*(p + q + S(2)) + c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + q + S(4))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule472)

    pattern473 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons212, cons7, cons11, cons14, cons233, cons232)
    rule473 = ReplacementRule(pattern473, lambda x, c, d, f, q, b, C, p, a, B, A : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - B*(a*b*f + b*c*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) + S(2)*f*(-b*(A*c - C*a)*(a*f + c*d) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule473)

    pattern474 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons212, cons7, cons11, cons14, cons233, cons232)
    rule474 = ReplacementRule(pattern474, lambda x, c, d, f, q, C, p, a, A, b : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) + S(2)*f*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule474)


    def cons_f287(p, q):
        return NonzeroQ(S(2)*p + S(2)*q + S(3))

    cons287 = CustomConstraint(cons_f287)
    pattern475 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons212, cons7, cons221, cons11, cons12, cons262, cons287)
    rule475 = ReplacementRule(pattern475, lambda x, c, d, f, e, q, b, C, p, a, B, A : (a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) + S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(p*(-b*f + c*e)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule475)

    pattern476 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons212, cons7, cons221, cons11, cons12, cons262, cons287)
    rule476 = ReplacementRule(pattern476, lambda x, c, d, f, e, q, C, p, a, A, b : (S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + x**S(2)*(p*(-b*f + c*e)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2)))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule476)

    pattern477 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons212, cons221, cons11, cons12, cons262, cons287)
    rule477 = ReplacementRule(pattern477, lambda x, c, d, f, e, q, C, p, a, B, A : (a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule477)

    pattern478 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons212, cons221, cons11, cons12, cons262, cons287)
    rule478 = ReplacementRule(pattern478, lambda x, c, d, f, e, q, C, p, a, A : (a + c*x**S(2))**p*(-C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule478)

    pattern479 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons212, cons7, cons11, cons12, cons262, cons287)
    rule479 = ReplacementRule(pattern479, lambda x, c, d, f, q, b, C, p, a, B, A : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p*(B*c*f*(S(2)*p + S(2)*q + S(3)) + C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(b*d*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + x**S(2)*(-b*f*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1)))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule479)

    pattern480 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons212, cons7, cons11, cons12, cons262, cons287)
    rule480 = ReplacementRule(pattern480, lambda x, c, d, f, q, C, p, a, A, b : (d + f*x**S(2))**(q + S(1))*(C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))*(a + b*x + c*x**S(2))**p/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-C*b**S(2)*d*f*p*(q + S(1)) + x**S(2)*(C*b**S(2)*f**S(2)*p*(q + S(1)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(2)*C*b*f*p*(q + S(1))*(-a*f + c*d) - b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule480)

    pattern481 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons7, cons221, )
    def With481(x, c, d, f, e, C, b, a, B, A):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d - B*a*b*f + B*a*c*e + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) + B*b*d*f - B*c*d*e - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule481 = ReplacementRule(pattern481, lambda x, c, d, f, e, C, b, a, B, A : With481(x, c, d, f, e, C, b, a, B, A))
    rubi.add(rule481)

    pattern482 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons7, cons221, )
    def With482(x, c, d, f, e, C, a, A, b):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule482 = ReplacementRule(pattern482, lambda x, c, d, f, e, C, a, A, b : With482(x, c, d, f, e, C, a, A, b))
    rubi.add(rule482)

    pattern483 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons7, )
    def With483(x, c, d, f, C, b, a, B, A):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f + B*b*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d - B*a*b*f + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule483 = ReplacementRule(pattern483, lambda x, c, d, f, C, b, a, B, A : With483(x, c, d, f, C, b, a, B, A))
    rubi.add(rule483)

    pattern484 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons7, )
    def With484(x, c, d, f, C, a, A, b):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule484 = ReplacementRule(pattern484, lambda x, c, d, f, C, a, A, b : With484(x, c, d, f, C, a, A, b))
    rubi.add(rule484)

    pattern485 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons7, cons221)
    rule485 = ReplacementRule(pattern485, lambda x, c, d, f, e, C, b, a, B, A : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule485)

    pattern486 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons7, cons221)
    rule486 = ReplacementRule(pattern486, lambda x, c, d, f, e, C, a, A, b : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule486)

    pattern487 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons221)
    rule487 = ReplacementRule(pattern487, lambda x, d, c, f, e, C, a, B, A : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c + B*c*x - C*a)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule487)

    pattern488 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons221)
    rule488 = ReplacementRule(pattern488, lambda x, c, d, f, e, C, a, A : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + (A*c - C*a)*Int(S(1)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule488)

    pattern489 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons285, cons286, cons7)
    rule489 = ReplacementRule(pattern489, lambda x, c, d, f, C, b, a, B, A : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule489)

    pattern490 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons2, cons3, cons4, cons25, cons123, cons284, cons286, cons7)
    rule490 = ReplacementRule(pattern490, lambda x, c, d, f, C, a, A, b : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule490)


    def cons_f288(x, c, d, f, e, q, C, p, a, B, A, b):
        return FreeQ(List(a, b, c, d, e, f, A, B, C, p, q), x)

    cons288 = CustomConstraint(cons_f288)
    pattern491 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons288)
    rule491 = ReplacementRule(pattern491, lambda x, c, d, f, e, q, b, C, p, a, B, A : Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule491)


    def cons_f289(x, c, d, f, e, q, C, p, a, A, b):
        return FreeQ(List(a, b, c, d, e, f, A, C, p, q), x)

    cons289 = CustomConstraint(cons_f289)
    pattern492 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons289)
    rule492 = ReplacementRule(pattern492, lambda x, c, d, f, e, q, C, p, a, A, b : Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule492)


    def cons_f290(x, c, d, f, e, q, C, p, a, B, A):
        return FreeQ(List(a, c, d, e, f, A, B, C, p, q), x)

    cons290 = CustomConstraint(cons_f290)
    pattern493 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons290)
    rule493 = ReplacementRule(pattern493, lambda x, c, d, f, e, q, C, p, a, B, A : Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule493)


    def cons_f291(x, c, d, f, e, q, C, p, a, A):
        return FreeQ(List(a, c, d, e, f, A, C, p, q), x)

    cons291 = CustomConstraint(cons_f291)
    pattern494 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons291)
    rule494 = ReplacementRule(pattern494, lambda x, c, d, f, e, q, C, p, a, A : Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule494)

    pattern495 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons21, cons22)
    rule495 = ReplacementRule(pattern495, lambda x, c, d, f, e, q, C, u, p, a, B, A, b : Subst(Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule495)

    pattern496 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons21, cons22)
    rule496 = ReplacementRule(pattern496, lambda x, c, d, f, e, q, b, u, p, a, B, A : Subst(Int((A + B*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule496)

    pattern497 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons3, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons21, cons22)
    rule497 = ReplacementRule(pattern497, lambda x, c, d, f, e, q, C, u, p, a, A, b : Subst(Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule497)

    pattern498 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons21, cons22)
    rule498 = ReplacementRule(pattern498, lambda x, c, d, f, e, q, C, u, p, a, B, A : Subst(Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule498)

    pattern499 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons285, cons286, cons6, cons212, cons21, cons22)
    rule499 = ReplacementRule(pattern499, lambda x, c, d, f, e, q, u, p, a, B, A : Subst(Int((A + B*x)*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule499)

    pattern500 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons2, cons4, cons25, cons26, cons123, cons284, cons286, cons6, cons212, cons21, cons22)
    rule500 = ReplacementRule(pattern500, lambda x, c, d, f, e, q, C, u, p, a, A : Subst(Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule500)

    return rubi
