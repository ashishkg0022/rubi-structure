from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on

if matchpy:
    from matchpy import Pattern, ReplacementRule, CustomConstraint
    from sympy.integrals.rubi.utility_function import (
        sympy_op_factory, Int, Sum, Set, With, Module, Scan, MapAnd, FalseQ,
        ZeroQ, NegativeQ, NonzeroQ, FreeQ, NFreeQ, List, Log, PositiveQ,
        PositiveIntegerQ, NegativeIntegerQ, IntegerQ, IntegersQ,
        ComplexNumberQ, PureComplexNumberQ, RealNumericQ, PositiveOrZeroQ,
        NegativeOrZeroQ, FractionOrNegativeQ, NegQ, Equal, Unequal, IntPart,
        FracPart, RationalQ, ProductQ, SumQ, NonsumQ, Subst, First, Rest,
        SqrtNumberQ, SqrtNumberSumQ, LinearQ, Sqrt, ArcCosh, Coefficient,
        Denominator, Hypergeometric2F1, Not, Simplify, FractionalPart,
        IntegerPart, AppellF1, EllipticPi, EllipticE, EllipticF, ArcTan,
        ArcCot, ArcCoth, ArcTanh, ArcSin, ArcSinh, ArcCos, ArcCsc, ArcSec,
        ArcCsch, ArcSech, Sinh, Tanh, Cosh, Sech, Csch, Coth, LessEqual, Less,
        Greater, GreaterEqual, FractionQ, IntLinearcQ, Expand, IndependentQ,
        PowerQ, IntegerPowerQ, PositiveIntegerPowerQ, FractionalPowerQ, AtomQ,
        ExpQ, LogQ, Head, MemberQ, TrigQ, SinQ, CosQ, TanQ, CotQ, SecQ, CscQ,
        Sin, Cos, Tan, Cot, Sec, Csc, HyperbolicQ, SinhQ, CoshQ, TanhQ, CothQ,
        SechQ, CschQ, InverseTrigQ, SinCosQ, SinhCoshQ, LeafCount, Numerator,
        NumberQ, NumericQ, Length, ListQ, Im, Re, InverseHyperbolicQ,
        InverseFunctionQ, TrigHyperbolicFreeQ, InverseFunctionFreeQ, RealQ,
        EqQ, FractionalPowerFreeQ, ComplexFreeQ, PolynomialQ, FactorSquareFree,
        PowerOfLinearQ, Exponent, QuadraticQ, LinearPairQ, BinomialParts,
        TrinomialParts, PolyQ, EvenQ, OddQ, PerfectSquareQ, NiceSqrtAuxQ,
        NiceSqrtQ, Together, PosAux, PosQ, CoefficientList, ReplaceAll,
        ExpandLinearProduct, GCD, ContentFactor, NumericFactor,
        NonnumericFactors, MakeAssocList, GensymSubst, KernelSubst,
        ExpandExpression, Apart, SmartApart, MatchQ,
        PolynomialQuotientRemainder, FreeFactors, NonfreeFactors,
        RemoveContentAux, RemoveContent, FreeTerms, NonfreeTerms,
        ExpandAlgebraicFunction, CollectReciprocals, ExpandCleanup,
        AlgebraicFunctionQ, Coeff, LeadTerm, RemainingTerms, LeadFactor,
        RemainingFactors, LeadBase, LeadDegree, Numer, Denom, hypergeom, Expon,
        MergeMonomials, PolynomialDivide, BinomialQ, TrinomialQ,
        GeneralizedBinomialQ, GeneralizedTrinomialQ, FactorSquareFreeList,
        PerfectPowerTest, SquareFreeFactorTest, RationalFunctionQ,
        RationalFunctionFactors, NonrationalFunctionFactors, Reverse,
        RationalFunctionExponents, RationalFunctionExpand, ExpandIntegrand,
        SimplerQ, SimplerSqrtQ, SumSimplerQ, BinomialDegree, TrinomialDegree,
        CancelCommonFactors, SimplerIntegrandQ, GeneralizedBinomialDegree,
        GeneralizedBinomialParts, GeneralizedTrinomialDegree,
        GeneralizedTrinomialParts, MonomialQ, MonomialSumQ,
        MinimumMonomialExponent, MonomialExponent, LinearMatchQ,
        PowerOfLinearMatchQ, QuadraticMatchQ, CubicMatchQ, BinomialMatchQ,
        TrinomialMatchQ, GeneralizedBinomialMatchQ, GeneralizedTrinomialMatchQ,
        QuotientOfLinearsMatchQ, PolynomialTermQ, PolynomialTerms,
        NonpolynomialTerms, PseudoBinomialParts, NormalizePseudoBinomial,
        PseudoBinomialPairQ, PseudoBinomialQ, PolynomialGCD, PolyGCD,
        AlgebraicFunctionFactors, NonalgebraicFunctionFactors,
        QuotientOfLinearsP, QuotientOfLinearsParts, QuotientOfLinearsQ,
        Flatten, Sort, AbsurdNumberQ, AbsurdNumberFactors,
        NonabsurdNumberFactors, SumSimplerAuxQ, Prepend, Drop,
        CombineExponents, FactorInteger, FactorAbsurdNumber,
        SubstForInverseFunction, SubstForFractionalPower,
        SubstForFractionalPowerOfQuotientOfLinears,
        FractionalPowerOfQuotientOfLinears, SubstForFractionalPowerQ,
        SubstForFractionalPowerAuxQ, FractionalPowerOfSquareQ,
        FractionalPowerSubexpressionQ, Apply, FactorNumericGcd,
        MergeableFactorQ, MergeFactor, MergeFactors, TrigSimplifyQ,
        TrigSimplify, TrigSimplifyRecur, Order, FactorOrder, Smallest,
        OrderedQ, MinimumDegree, PositiveFactors, Sign, NonpositiveFactors,
        PolynomialInAuxQ, PolynomialInQ, ExponentInAux, ExponentIn,
        PolynomialInSubstAux, PolynomialInSubst, Distrib, DistributeDegree,
        FunctionOfPower, DivideDegreesOfFactors, MonomialFactor, FullSimplify,
        FunctionOfLinearSubst, FunctionOfLinear, NormalizeIntegrand,
        NormalizeIntegrandAux, NormalizeIntegrandFactor,
        NormalizeIntegrandFactorBase, NormalizeTogether,
        NormalizeLeadTermSigns, AbsorbMinusSign, NormalizeSumFactors,
        SignOfFactor, NormalizePowerOfLinear, SimplifyIntegrand, SimplifyTerm,
        TogetherSimplify, SmartSimplify, SubstForExpn, ExpandToSum, UnifySum,
        UnifyTerms, UnifyTerm, CalculusQ, FunctionOfInverseLinear,
        PureFunctionOfSinhQ, PureFunctionOfTanhQ, PureFunctionOfCoshQ,
        IntegerQuotientQ, OddQuotientQ, EvenQuotientQ, FindTrigFactor,
        FunctionOfSinhQ, FunctionOfCoshQ, OddHyperbolicPowerQ, FunctionOfTanhQ,
        FunctionOfTanhWeight, FunctionOfHyperbolicQ, SmartNumerator,
        SmartDenominator, SubstForAux, ActivateTrig, ExpandTrig, TrigExpand,
        SubstForTrig, SubstForHyperbolic, InertTrigFreeQ, LCM,
        SubstForFractionalPowerOfLinear, FractionalPowerOfLinear,
        InverseFunctionOfLinear, InertTrigQ, InertReciprocalQ, DeactivateTrig,
        FixInertTrigFunction, DeactivateTrigAux, PowerOfInertTrigSumQ,
        PiecewiseLinearQ, KnownTrigIntegrandQ, KnownSineIntegrandQ,
        KnownTangentIntegrandQ, KnownCotangentIntegrandQ,
        KnownSecantIntegrandQ, TryPureTanSubst, TryTanhSubst, TryPureTanhSubst,
        AbsurdNumberGCD, AbsurdNumberGCDList, ExpandTrigExpand,
        ExpandTrigReduce, ExpandTrigReduceAux, NormalizeTrig, TrigToExp,
        ExpandTrigToExp, TrigReduce, FunctionOfTrig, AlgebraicTrigFunctionQ,
        FunctionOfHyperbolic, FunctionOfQ, FunctionOfExpnQ, PureFunctionOfSinQ,
        PureFunctionOfCosQ, PureFunctionOfTanQ, PureFunctionOfCotQ,
        FunctionOfCosQ, FunctionOfSinQ, OddTrigPowerQ, FunctionOfTanQ,
        FunctionOfTanWeight, FunctionOfTrigQ, FunctionOfDensePolynomialsQ,
        FunctionOfLog, PowerVariableExpn, PowerVariableDegree,
        PowerVariableSubst, EulerIntegrandQ, FunctionOfSquareRootOfQuadratic,
        SquareRootOfQuadraticSubst, Divides, EasyDQ, ProductOfLinearPowersQ,
        Rt, NthRoot, AtomBaseQ, SumBaseQ, NegSumBaseQ, AllNegTermQ,
        SomeNegTermQ, TrigSquareQ, RtAux, TrigSquare, IntSum, IntTerm, Map2,
        ConstantFactor, SameQ, ReplacePart, CommonFactors,
        MostMainFactorPosition, FunctionOfExponentialQ, FunctionOfExponential,
        FunctionOfExponentialFunction, FunctionOfExponentialFunctionAux,
        FunctionOfExponentialTest, FunctionOfExponentialTestAux, stdev,
        rubi_test, If, IntQuadraticQ, IntBinomialQ, RectifyTangent,
        RectifyCotangent, Inequality, Condition, Simp, SimpHelp, SplitProduct,
        SplitSum, SubstFor, SubstForAux, FresnelS, FresnelC, Erfc, Erfi, Gamma,
        FunctionOfTrigOfLinearQ, ElementaryFunctionQ, Complex, UnsameQ,
        _SimpFixFactor, SimpFixFactor, _FixSimplify, FixSimplify,
        _SimplifyAntiderivativeSum, SimplifyAntiderivativeSum,
        _SimplifyAntiderivative, SimplifyAntiderivative, _TrigSimplifyAux,
        TrigSimplifyAux, Cancel, Part, PolyLog, D, Dist
    )
    from sympy import Integral, S, sqrt, And, Or

def cons_f1(m):
    return NonzeroQ(m + S(1))

cons1 = CustomConstraint(cons_f1)

def cons_f2(m, x):
    return FreeQ(m, x)

cons2 = CustomConstraint(cons_f2)

def cons_f3(x, a, b):
    return FreeQ(List(a, b), x)

cons3 = CustomConstraint(cons_f3)

def cons_f4(a, x):
    return FreeQ(a, x)

cons4 = CustomConstraint(cons_f4)

def cons_f5(b, x):
    return FreeQ(b, x)

cons5 = CustomConstraint(cons_f5)

def cons_f6(u, x):
    return LinearQ(u, x)

cons6 = CustomConstraint(cons_f6)

def cons_f7(u, x):
    return NonzeroQ(u - x)

cons7 = CustomConstraint(cons_f7)

def cons_f8(c, b, a, d):
    return ZeroQ(a*d + b*c)

cons8 = CustomConstraint(cons_f8)

def cons_f9(c, x):
    return FreeQ(c, x)

cons9 = CustomConstraint(cons_f9)

def cons_f10(d, x):
    return FreeQ(d, x)

cons10 = CustomConstraint(cons_f10)

def cons_f11(c, b, a, d):
    return NonzeroQ(-a*d + b*c)

cons11 = CustomConstraint(cons_f11)

def cons_f12(n, m):
    return ZeroQ(m + n + S(2))

cons12 = CustomConstraint(cons_f12)

def cons_f13(n, x):
    return FreeQ(n, x)

cons13 = CustomConstraint(cons_f13)

def cons_f14(m):
    return PositiveIntegerQ(m + S(1)/2)

cons14 = CustomConstraint(cons_f14)

def cons_f15(m):
    return NegativeIntegerQ(m + S(3)/2)

cons15 = CustomConstraint(cons_f15)

def cons_f16(a, m, c):
    return Or(IntegerQ(m), And(PositiveQ(a), PositiveQ(c)))

cons16 = CustomConstraint(cons_f16)

def cons_f17(a):
    return PositiveQ(a)

cons17 = CustomConstraint(cons_f17)

def cons_f18(a, c):
    return ZeroQ(a + c)

cons18 = CustomConstraint(cons_f18)

def cons_f19(m):
    return Not(IntegerQ(S(2)*m))

cons19 = CustomConstraint(cons_f19)

def cons_f20(d, b, a, c):
    return PosQ(b*d/(a*c))

cons20 = CustomConstraint(cons_f20)

def cons_f21(m):
    return IntegerQ(m + S(1)/2)

cons21 = CustomConstraint(cons_f21)

def cons_f22(n):
    return IntegerQ(n + S(1)/2)

cons22 = CustomConstraint(cons_f22)

def cons_f23(n, m):
    return Less(S(0), m, n)

cons23 = CustomConstraint(cons_f23)

def cons_f24(n, m):
    return Less(m, n, S(0))

cons24 = CustomConstraint(cons_f24)

def cons_f25(m):
    return PositiveIntegerQ(m)

cons25 = CustomConstraint(cons_f25)

def cons_f26(n, c, m):
    return Or(Not(IntegerQ(n)), And(ZeroQ(c), LessEqual(S(7)*m + S(4)*n, S(0))), Less(S(9)*m + S(5)*n + S(5), S(0)), Greater(m + n + S(2), S(0)))

cons26 = CustomConstraint(cons_f26)

def cons_f27(m):
    return NegativeIntegerQ(m)

cons27 = CustomConstraint(cons_f27)

def cons_f28(n):
    return IntegerQ(n)

cons28 = CustomConstraint(cons_f28)

def cons_f29(n, m):
    return Not(And(PositiveIntegerQ(n), Less(m + n + S(2), S(0))))

cons29 = CustomConstraint(cons_f29)

def cons_f30(n):
    return RationalQ(n)

cons30 = CustomConstraint(cons_f30)

def cons_f31(n):
    return Greater(n, S(0))

cons31 = CustomConstraint(cons_f31)

def cons_f32(n):
    return Less(n, S(-1))

cons32 = CustomConstraint(cons_f32)

def cons_f33(c, b, a, d):
    return PosQ((-a*d + b*c)/b)

cons33 = CustomConstraint(cons_f33)

def cons_f34(c, b, a, d):
    return NegQ((-a*d + b*c)/b)

cons34 = CustomConstraint(cons_f34)

def cons_f35(n):
    return Less(S(-1), n, S(0))

cons35 = CustomConstraint(cons_f35)

def cons_f36(n):
    return Not(IntegerQ(n))

cons36 = CustomConstraint(cons_f36)

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

def cons_f41(b, a, m, n, x, d, c):
    return IntLinearcQ(a, b, c, d, m, n, x)

cons41 = CustomConstraint(cons_f41)

def cons_f42(n, a, m, c):
    return Not(And(Less(n, S(-1)), Or(ZeroQ(a), And(NonzeroQ(c), Less(m, n), IntegerQ(n)))))

cons42 = CustomConstraint(cons_f42)

def cons_f43(n, m):
    return Unequal(m + n + S(1), S(0))

cons43 = CustomConstraint(cons_f43)

def cons_f44(n, m):
    return Not(And(PositiveIntegerQ(m), Or(Not(IntegerQ(n)), Less(S(0), m, n))))

cons44 = CustomConstraint(cons_f44)

def cons_f45(n, m):
    return Not(And(IntegerQ(m + n), Less(m + n + S(2), S(0))))

cons45 = CustomConstraint(cons_f45)

def cons_f46(d, b):
    return ZeroQ(b + d)

cons46 = CustomConstraint(cons_f46)

def cons_f47(a, c):
    return PositiveQ(a + c)

cons47 = CustomConstraint(cons_f47)

def cons_f48(c, b, a, d):
    return PositiveQ(-a*d + b*c)

cons48 = CustomConstraint(cons_f48)

def cons_f49(b):
    return PositiveQ(b)

cons49 = CustomConstraint(cons_f49)

def cons_f50(d, b):
    return ZeroQ(b - d)

cons50 = CustomConstraint(cons_f50)

def cons_f51(m):
    return RationalQ(m)

cons51 = CustomConstraint(cons_f51)

def cons_f52(m):
    return Less(S(-1), m, S(0))

cons52 = CustomConstraint(cons_f52)

def cons_f53(m):
    return LessEqual(S(3), Denominator(m), S(4))

cons53 = CustomConstraint(cons_f53)

def cons_f54(d, b):
    return PosQ(d/b)

cons54 = CustomConstraint(cons_f54)

def cons_f55(d, b):
    return NegQ(d/b)

cons55 = CustomConstraint(cons_f55)

def cons_f56(n, m):
    return Equal(m + n + S(1), S(0))

cons56 = CustomConstraint(cons_f56)

def cons_f57(n, m):
    return LessEqual(Denominator(n), Denominator(m))

cons57 = CustomConstraint(cons_f57)

def cons_f58(n, m):
    return NegativeIntegerQ(m + n + S(2))

cons58 = CustomConstraint(cons_f58)

def cons_f59(n, m):
    return Or(_SumSimplerQ(m, S(1)), Not(_SumSimplerQ(n, S(1))))

cons59 = CustomConstraint(cons_f59)

def cons_f60(m):
    return Not(IntegerQ(m))

cons60 = CustomConstraint(cons_f60)

def cons_f61(n, c, d, b):
    return Or(IntegerQ(n), And(PositiveQ(c), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2)), PositiveQ(-d/(b*c))))))

cons61 = CustomConstraint(cons_f61)

def cons_f62(d, m, b, c):
    return Or(IntegerQ(m), PositiveQ(-d/(b*c)))

cons62 = CustomConstraint(cons_f62)

def cons_f63(c):
    return Not(PositiveQ(c))

cons63 = CustomConstraint(cons_f63)

def cons_f64(d, b, c):
    return Not(PositiveQ(-d/(b*c)))

cons64 = CustomConstraint(cons_f64)

def cons_f65(n, c, m, d):
    return Or(And(RationalQ(m), Not(And(ZeroQ(n + S(1)/2), ZeroQ(c**S(2) - d**S(2))))), Not(RationalQ(n)))

cons65 = CustomConstraint(cons_f65)

def cons_f66(c, b, a, d):
    return PositiveQ(b/(-a*d + b*c))

cons66 = CustomConstraint(cons_f66)

def cons_f67(b, a, m, n, d, c):
    return Or(RationalQ(m), Not(And(RationalQ(n), PositiveQ(-d/(-a*d + b*c)))))

cons67 = CustomConstraint(cons_f67)

def cons_f68(n, m):
    return Or(RationalQ(m), Not(SimplerQ(n + S(1), m + S(1))))

cons68 = CustomConstraint(cons_f68)

def cons_f69(u, x):
    return NonzeroQ(Coefficient(u, x, S(0)))

cons69 = CustomConstraint(cons_f69)

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

def cons_f75(n, p):
    return NonzeroQ(n + p + S(2))

cons75 = CustomConstraint(cons_f75)

def cons_f76(b, f, p, e, a, n, d, c):
    return ZeroQ(a*d*f*(n + p + S(2)) - b*(c*f*(p + S(1)) + d*e*(n + S(1))))

cons76 = CustomConstraint(cons_f76)

def cons_f77(p):
    return PositiveIntegerQ(p)

cons77 = CustomConstraint(cons_f77)

def cons_f78(f, e, a, b):
    return ZeroQ(a*f + b*e)

cons78 = CustomConstraint(cons_f78)

def cons_f79(n, p):
    return Not(And(NegativeIntegerQ(n + p + S(2)), Greater(n + S(2)*p, S(0))))

cons79 = CustomConstraint(cons_f79)

def cons_f80(n, p):
    return Or(NonzeroQ(n + S(1)), Equal(p, S(1)))

cons80 = CustomConstraint(cons_f80)

def cons_f81(f, e, a, b):
    return NonzeroQ(a*f + b*e)

cons81 = CustomConstraint(cons_f81)

def cons_f82(b, p, f, e, a, n, d):
    return Or(Not(IntegerQ(n)), Less(S(5)*n + S(9)*p, S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, d, e, f)))

cons82 = CustomConstraint(cons_f82)

def cons_f83(b, d, p, f, e, a, n, c):
    return Or(NegativeIntegerQ(n, p), ZeroQ(p + S(-1)), And(PositiveIntegerQ(p), Or(Not(IntegerQ(n)), LessEqual(S(5)*n + S(9)*p + S(10), S(0)), GreaterEqual(n + p + S(1), S(0)), And(GreaterEqual(n + p + S(2), S(0)), RationalQ(a, b, c, d, e, f)))))

cons83 = CustomConstraint(cons_f83)

def cons_f84(n, p):
    return ZeroQ(n + p + S(2))

cons84 = CustomConstraint(cons_f84)

def cons_f85(p):
    return NonzeroQ(p + S(1))

cons85 = CustomConstraint(cons_f85)

def cons_f86(n, p):
    return Not(And(_SumSimplerQ(n, S(1)), Not(_SumSimplerQ(p, S(1)))))

cons86 = CustomConstraint(cons_f86)

def cons_f87(p):
    return RationalQ(p)

cons87 = CustomConstraint(cons_f87)

def cons_f88(p):
    return Less(p, S(-1))

cons88 = CustomConstraint(cons_f88)

def cons_f89(n, p, c, e):
    return Or(Not(And(RationalQ(n), Less(n, S(-1)))), IntegerQ(p), Not(Or(IntegerQ(n), Not(Or(ZeroQ(e), Not(Or(ZeroQ(c), Less(p, n))))))))

cons89 = CustomConstraint(cons_f89)

def cons_f90(p):
    return Not(RationalQ(p))

cons90 = CustomConstraint(cons_f90)

def cons_f91(p):
    return _SumSimplerQ(p, S(1))

cons91 = CustomConstraint(cons_f91)

def cons_f92(n, p):
    return NonzeroQ(n + p + S(3))

cons92 = CustomConstraint(cons_f92)

def cons_f93(b, f, p, e, a, n, d, c):
    return ZeroQ(-b*(c*f*(p + S(1)) + d*e*(n + S(1)))*(a*d*f*(n + p + S(4)) - b*(c*f*(p + S(2)) + d*e*(n + S(2)))) + d*f*(a**S(2)*d*f*(n + p + S(3)) - b*(a*(c*f*(p + S(1)) + d*e*(n + S(1))) + b*c*e))*(n + p + S(2)))

cons93 = CustomConstraint(cons_f93)

def cons_f94(n, m):
    return ZeroQ(m - n + S(-1))

cons94 = CustomConstraint(cons_f94)

def cons_f95(m):
    return Not(PositiveIntegerQ(m))

cons95 = CustomConstraint(cons_f95)

def cons_f96(n, p, m):
    return NonzeroQ(m + n + p + S(2))

cons96 = CustomConstraint(cons_f96)

def cons_f97(p):
    return IntegerQ(p)

cons97 = CustomConstraint(cons_f97)

def cons_f98(p):
    return Less(S(0), p, S(1))

cons98 = CustomConstraint(cons_f98)

def cons_f99(p):
    return Greater(p, S(1))

cons99 = CustomConstraint(cons_f99)

def cons_f100(p):
    return Not(IntegerQ(p))

cons100 = CustomConstraint(cons_f100)

def cons_f101(n):
    return PositiveIntegerQ(n)

cons101 = CustomConstraint(cons_f101)

def cons_f102(p):
    return FractionQ(p)

cons102 = CustomConstraint(cons_f102)

def cons_f103(n, m):
    return IntegersQ(m, n)

cons103 = CustomConstraint(cons_f103)

def cons_f104(n, p, m):
    return Or(IntegerQ(p), And(Greater(m, S(0)), GreaterEqual(n, S(-1))))

cons104 = CustomConstraint(cons_f104)

def cons_f105(n, p):
    return Or(And(RationalQ(n), Less(n, S(-1))), And(ZeroQ(n + p + S(3)), NonzeroQ(n + S(1)), Or(_SumSimplerQ(n, S(1)), Not(_SumSimplerQ(p, S(1))))))

cons105 = CustomConstraint(cons_f105)

def cons_f106(b, f, e, a, x, d, c):
    return FreeQ(List(a, b, c, d, e, f), x)

cons106 = CustomConstraint(cons_f106)

def cons_f107(b, f, e, a, d, c):
    return ZeroQ(S(2)*b*d*e - f*(a*d + b*c))

cons107 = CustomConstraint(cons_f107)

def cons_f108(n, m):
    return ZeroQ(m + n + S(1))

cons108 = CustomConstraint(cons_f108)

def cons_f109(b, a, x, d, c):
    return SimplerQ(a + b*x, c + d*x)

cons109 = CustomConstraint(cons_f109)

def cons_f110(n, p, m):
    return ZeroQ(m + n + p + S(2))

cons110 = CustomConstraint(cons_f110)

def cons_f111(p, m):
    return Not(And(_SumSimplerQ(p, S(1)), Not(_SumSimplerQ(m, S(1)))))

cons111 = CustomConstraint(cons_f111)

def cons_f112(n, p, m):
    return ZeroQ(m + n + p + S(3))

cons112 = CustomConstraint(cons_f112)

def cons_f113(b, f, p, e, a, m, n, d, c):
    return ZeroQ(a*d*f*(m + S(1)) + b*c*f*(n + S(1)) + b*d*e*(p + S(1)))

cons113 = CustomConstraint(cons_f113)

def cons_f114(m):
    return Or(And(RationalQ(m), Less(m, S(-1))), _SumSimplerQ(m, S(1)))

cons114 = CustomConstraint(cons_f114)

def cons_f115(n, p, m):
    return RationalQ(m, n, p)

cons115 = CustomConstraint(cons_f115)

def cons_f116(p):
    return Greater(p, S(0))

cons116 = CustomConstraint(cons_f116)

def cons_f117(n, p, m):
    return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), IntegersQ(m, n + p), IntegersQ(p, m + n))

cons117 = CustomConstraint(cons_f117)

def cons_f118(n):
    return Greater(n, S(1))

cons118 = CustomConstraint(cons_f118)

def cons_f119(m):
    return Greater(m, S(1))

cons119 = CustomConstraint(cons_f119)

def cons_f120(n, p, m):
    return NonzeroQ(m + n + p + S(1))

cons120 = CustomConstraint(cons_f120)

def cons_f121(m):
    return Greater(m, S(0))

cons121 = CustomConstraint(cons_f121)

def cons_f122(n, p, m):
    return Or(IntegersQ(S(2)*m, S(2)*n, S(2)*p), Or(IntegersQ(m, n + p), IntegersQ(p, m + n)))

cons122 = CustomConstraint(cons_f122)

def cons_f123(n, p, m):
    return IntegersQ(S(2)*m, S(2)*n, S(2)*p)

cons123 = CustomConstraint(cons_f123)

def cons_f124(n, p):
    return Or(IntegerQ(n), IntegersQ(S(2)*n, S(2)*p))

cons124 = CustomConstraint(cons_f124)

def cons_f125(n, m):
    return PositiveIntegerQ(m + n + S(1))

cons125 = CustomConstraint(cons_f125)

def cons_f126(n, m):
    return Or(And(RationalQ(m), Greater(m, S(0))), And(Not(RationalQ(m)), Or(_SumSimplerQ(m, S(-1)), Not(_SumSimplerQ(n, S(-1))))))

cons126 = CustomConstraint(cons_f126)

def cons_f127(f, d, e, c):
    return PositiveQ(-f/(-c*f + d*e))

cons127 = CustomConstraint(cons_f127)

def cons_f128(f, d, e, c):
    return Not(PositiveQ(-f/(-c*f + d*e)))

cons128 = CustomConstraint(cons_f128)

def cons_f129(f, e, d, c):
    return NonzeroQ(-c*f + d*e)

cons129 = CustomConstraint(cons_f129)

def cons_f130(c):
    return PositiveQ(c)

cons130 = CustomConstraint(cons_f130)

def cons_f131(e):
    return PositiveQ(e)

cons131 = CustomConstraint(cons_f131)

def cons_f132(d, b):
    return Not(NegativeQ(-b/d))

cons132 = CustomConstraint(cons_f132)

def cons_f133(d, b):
    return NegativeQ(-b/d)

cons133 = CustomConstraint(cons_f133)

def cons_f134(e, c):
    return Not(And(PositiveQ(c), PositiveQ(e)))

cons134 = CustomConstraint(cons_f134)

def cons_f135(f, e, a, b):
    return PositiveQ(b/(-a*f + b*e))

cons135 = CustomConstraint(cons_f135)

def cons_f136(c, b, a, d):
    return Not(NegativeQ(-(-a*d + b*c)/d))

cons136 = CustomConstraint(cons_f136)

def cons_f137(d, b, f, e, a, x, c):
    return Not(And(SimplerQ(c + d*x, a + b*x), PositiveQ(-d/(-a*d + b*c)), PositiveQ(d/(-c*f + d*e)), Not(NegativeQ((-a*d + b*c)/b))))

cons137 = CustomConstraint(cons_f137)

def cons_f138(b, d, f, e, a, c):
    return Not(And(PositiveQ(b/(-a*d + b*c)), PositiveQ(b/(-a*f + b*e))))

cons138 = CustomConstraint(cons_f138)

def cons_f139(f, d, b):
    return Or(PositiveQ(-b/d), NegativeQ(-b/f))

cons139 = CustomConstraint(cons_f139)

def cons_f140(f, d, b):
    return Or(PosQ(-b/d), NegQ(-b/f))

cons140 = CustomConstraint(cons_f140)

def cons_f141(b, f, e, a, x):
    return SimplerQ(a + b*x, e + f*x)

cons141 = CustomConstraint(cons_f141)

def cons_f142(b, d, f, e, a, c):
    return Or(PositiveQ(-(-a*d + b*c)/d), NegativeQ(-(-a*f + b*e)/f))

cons142 = CustomConstraint(cons_f142)

def cons_f143(b, d, f, e, a, c):
    return Or(PosQ(-(-a*d + b*c)/d), NegQ(-(-a*f + b*e)/f))

cons143 = CustomConstraint(cons_f143)

def cons_f144(b, f, e, a, d, c):
    return ZeroQ(-a*d*f - b*c*f + S(2)*b*d*e)

cons144 = CustomConstraint(cons_f144)

def cons_f145(n, m):
    return PositiveIntegerQ(m - n)

cons145 = CustomConstraint(cons_f145)

def cons_f146(n, m):
    return Or(PositiveIntegerQ(m), NegativeIntegerQ(m, n))

cons146 = CustomConstraint(cons_f146)

def cons_f147(n, p, m):
    return NegativeIntegerQ(m + n + p + S(2))

cons147 = CustomConstraint(cons_f147)

def cons_f148(n, p, m):
    return Or(_SumSimplerQ(m, S(1)), And(Not(And(NonzeroQ(n + S(1)), _SumSimplerQ(n, S(1)))), Not(And(NonzeroQ(p + S(1)), _SumSimplerQ(p, S(1))))))

cons148 = CustomConstraint(cons_f148)

def cons_f149(n):
    return NegativeIntegerQ(n)

cons149 = CustomConstraint(cons_f149)

def cons_f150(p, e):
    return Or(IntegerQ(p), PositiveQ(e))

cons150 = CustomConstraint(cons_f150)

def cons_f151(d, b, c):
    return PositiveQ(-d/(b*c))

cons151 = CustomConstraint(cons_f151)

def cons_f152(f, p, e, d, c):
    return Or(IntegerQ(p), PositiveQ(d/(-c*f + d*e)))

cons152 = CustomConstraint(cons_f152)

def cons_f153(b, a, x, d, c):
    return Not(And(PositiveQ(d/(a*d - b*c)), SimplerQ(c + d*x, a + b*x)))

cons153 = CustomConstraint(cons_f153)

def cons_f154(c, b, a, d):
    return Not(PositiveQ(b/(-a*d + b*c)))

cons154 = CustomConstraint(cons_f154)

def cons_f155(d, b, a, x, c):
    return Not(SimplerQ(c + d*x, a + b*x))

cons155 = CustomConstraint(cons_f155)

def cons_f156(b, f, e, a, x, d, c):
    return Not(And(PositiveQ(d/(a*d - b*c)), PositiveQ(d/(-c*f + d*e)), SimplerQ(c + d*x, a + b*x)))

cons156 = CustomConstraint(cons_f156)

def cons_f157(b, d, f, e, a, x, c):
    return Not(And(PositiveQ(f/(a*f - b*e)), PositiveQ(f/(c*f - d*e)), SimplerQ(e + f*x, a + b*x)))

cons157 = CustomConstraint(cons_f157)

def cons_f158(f, e, a, b):
    return Not(PositiveQ(b/(-a*f + b*e)))

cons158 = CustomConstraint(cons_f158)

def cons_f159(b, f, e, a, x):
    return Not(SimplerQ(e + f*x, a + b*x))

cons159 = CustomConstraint(cons_f159)

def cons_f160(n, m):
    return Or(PositiveIntegerQ(m), IntegersQ(m, n))

cons160 = CustomConstraint(cons_f160)

def cons_f161(g, x):
    return FreeQ(g, x)

cons161 = CustomConstraint(cons_f161)

def cons_f162(h, x):
    return FreeQ(h, x)

cons162 = CustomConstraint(cons_f162)

def cons_f163(n, m):
    return Not(And(_SumSimplerQ(n, S(1)), Not(_SumSimplerQ(m, S(1)))))

cons163 = CustomConstraint(cons_f163)

def cons_f164(n, m):
    return Or(And(RationalQ(m), Less(m, S(-2))), And(ZeroQ(m + n + S(3)), Not(And(RationalQ(n), Less(n, S(-2))))))

cons164 = CustomConstraint(cons_f164)

def cons_f165(m):
    return Or(And(RationalQ(m), Inequality(S(-2), LessEqual, m, Less, S(-1))), _SumSimplerQ(m, S(1)))

cons165 = CustomConstraint(cons_f165)

def cons_f166(n, m):
    return NonzeroQ(m + n + S(3))

cons166 = CustomConstraint(cons_f166)

def cons_f167(n, m):
    return NonzeroQ(m + n + S(2))

cons167 = CustomConstraint(cons_f167)

def cons_f168(n, p, m):
    return Or(IntegersQ(m, n, p), PositiveIntegerQ(n, p))

cons168 = CustomConstraint(cons_f168)

def cons_f169(g, b, h, f, e, a, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, h), x)

cons169 = CustomConstraint(cons_f169)

def cons_f170(g, b, h, f, p, e, a, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, h, n, p), x)

cons170 = CustomConstraint(cons_f170)

def cons_f171(d, f, e, x, c):
    return SimplerQ(c + d*x, e + f*x)

cons171 = CustomConstraint(cons_f171)

def cons_f172(n, p, m):
    return Or(_SumSimplerQ(m, S(1)), And(Not(_SumSimplerQ(n, S(1))), Not(_SumSimplerQ(p, S(1)))))

cons172 = CustomConstraint(cons_f172)

def cons_f173(q, x):
    return FreeQ(q, x)

cons173 = CustomConstraint(cons_f173)

def cons_f174(p, q):
    return IntegersQ(p, q)

cons174 = CustomConstraint(cons_f174)

def cons_f175(q):
    return PositiveIntegerQ(q)

cons175 = CustomConstraint(cons_f175)

def cons_f176(g, b, h, f, p, e, a, m, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, h, m, n, p, q), x)

cons176 = CustomConstraint(cons_f176)

def cons_f177(g, b, h, r, f, i, p, e, a, m, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, h, i, m, n, p, q, r), x)

cons177 = CustomConstraint(cons_f177)

def cons_f178(i, x):
    return FreeQ(i, x)

cons178 = CustomConstraint(cons_f178)

def cons_f179(r, x):
    return FreeQ(r, x)

cons179 = CustomConstraint(cons_f179)

def cons_f180(n, p, b, x):
    return FreeQ(List(b, n, p), x)

cons180 = CustomConstraint(cons_f180)

def cons_f181(n, p):
    return ZeroQ(p + S(1) + S(1)/n)

cons181 = CustomConstraint(cons_f181)

def cons_f182(n, p):
    return NegativeIntegerQ(p + S(1) + S(1)/n)

cons182 = CustomConstraint(cons_f182)

def cons_f183(n):
    return NonzeroQ(S(3)*n + S(1))

cons183 = CustomConstraint(cons_f183)

def cons_f184(n):
    return Less(n, S(0))

cons184 = CustomConstraint(cons_f184)

def cons_f185(n, p):
    return PositiveIntegerQ(n, p)

cons185 = CustomConstraint(cons_f185)

def cons_f186(n, p):
    return Or(IntegerQ(S(2)*p), And(Equal(n, S(2)), IntegerQ(S(4)*p)), And(Equal(n, S(2)), IntegerQ(S(3)*p)), Less(Denominator(p + S(1)/n), Denominator(p)))

cons186 = CustomConstraint(cons_f186)

def cons_f187(a, b):
    return PosQ(b/a)

cons187 = CustomConstraint(cons_f187)

def cons_f188(a):
    return Not(PositiveQ(a))

cons188 = CustomConstraint(cons_f188)

def cons_f189(n):
    return PositiveIntegerQ(n/S(2) + S(-3)/2)

cons189 = CustomConstraint(cons_f189)

def cons_f190(a, b):
    return PosQ(a/b)

cons190 = CustomConstraint(cons_f190)

def cons_f191(a, b):
    return NegQ(a/b)

cons191 = CustomConstraint(cons_f191)

def cons_f192(a, b):
    return Or(PositiveQ(a), PositiveQ(b))

cons192 = CustomConstraint(cons_f192)

def cons_f193(a, b):
    return Or(NegativeQ(a), NegativeQ(b))

cons193 = CustomConstraint(cons_f193)

def cons_f194(a, b):
    return Or(PositiveQ(a), NegativeQ(b))

cons194 = CustomConstraint(cons_f194)

def cons_f195(a, b):
    return Or(NegativeQ(a), PositiveQ(b))

cons195 = CustomConstraint(cons_f195)

def cons_f196(n):
    return PositiveIntegerQ(n/S(4) + S(-1)/2)

cons196 = CustomConstraint(cons_f196)

def cons_f197(a, b):
    return Or(PositiveQ(a/b), And(PosQ(a/b), AtomQ(SplitProduct(_SumBaseQ, a)), AtomQ(SplitProduct(_SumBaseQ, b))))

cons197 = CustomConstraint(cons_f197)

def cons_f198(a, b):
    return Not(PositiveQ(a/b))

cons198 = CustomConstraint(cons_f198)

def cons_f199(n):
    return PositiveIntegerQ(n/S(4) + S(-1))

cons199 = CustomConstraint(cons_f199)

def cons_f200(a, b):
    return PositiveQ(a/b)

cons200 = CustomConstraint(cons_f200)

def cons_f201(b):
    return PosQ(b)

cons201 = CustomConstraint(cons_f201)

def cons_f202(b):
    return NegQ(b)

cons202 = CustomConstraint(cons_f202)

def cons_f203(a):
    return PosQ(a)

cons203 = CustomConstraint(cons_f203)

def cons_f204(a):
    return NegQ(a)

cons204 = CustomConstraint(cons_f204)

def cons_f205(a, b):
    return NegQ(b/a)

cons205 = CustomConstraint(cons_f205)

def cons_f206(a):
    return NegativeQ(a)

cons206 = CustomConstraint(cons_f206)

def cons_f207(p):
    return Less(S(-1), p, S(0))

cons207 = CustomConstraint(cons_f207)

def cons_f208(p):
    return Unequal(p, S(-1)/2)

cons208 = CustomConstraint(cons_f208)

def cons_f209(n, p):
    return IntegerQ(p + S(1)/n)

cons209 = CustomConstraint(cons_f209)

def cons_f210(n, p):
    return Less(Denominator(p + S(1)/n), Denominator(p))

cons210 = CustomConstraint(cons_f210)

def cons_f211(n):
    return FractionQ(n)

cons211 = CustomConstraint(cons_f211)

def cons_f212(p):
    return Not(PositiveIntegerQ(p))

cons212 = CustomConstraint(cons_f212)

def cons_f213(n):
    return Not(IntegerQ(S(1)/n))

cons213 = CustomConstraint(cons_f213)

def cons_f214(n, p):
    return Not(NegativeIntegerQ(p + S(1)/n))

cons214 = CustomConstraint(cons_f214)

def cons_f215(p, a):
    return Or(IntegerQ(p), PositiveQ(a))

cons215 = CustomConstraint(cons_f215)

def cons_f216(p, a):
    return Not(Or(IntegerQ(p), PositiveQ(a)))

cons216 = CustomConstraint(cons_f216)

def cons_f217(a1, a2, b2, b1):
    return ZeroQ(a1*b2 + a2*b1)

cons217 = CustomConstraint(cons_f217)

def cons_f218(p, a2, a1):
    return Or(IntegerQ(p), And(PositiveQ(a1), PositiveQ(a2)))

cons218 = CustomConstraint(cons_f218)

def cons_f219(a1, x):
    return FreeQ(a1, x)

cons219 = CustomConstraint(cons_f219)

def cons_f220(b1, x):
    return FreeQ(b1, x)

cons220 = CustomConstraint(cons_f220)

def cons_f221(a2, x):
    return FreeQ(a2, x)

cons221 = CustomConstraint(cons_f221)

def cons_f222(b2, x):
    return FreeQ(b2, x)

cons222 = CustomConstraint(cons_f222)

def cons_f223(n):
    return PositiveIntegerQ(S(2)*n)

cons223 = CustomConstraint(cons_f223)

def cons_f224(n, p):
    return Or(IntegerQ(S(2)*p), Less(Denominator(p + S(1)/n), Denominator(p)))

cons224 = CustomConstraint(cons_f224)

def cons_f225(n):
    return NegativeIntegerQ(S(2)*n)

cons225 = CustomConstraint(cons_f225)

def cons_f226(n):
    return FractionQ(S(2)*n)

cons226 = CustomConstraint(cons_f226)

def cons_f227(c, m):
    return Or(IntegerQ(m), PositiveQ(c))

cons227 = CustomConstraint(cons_f227)

def cons_f228(n, m):
    return IntegerQ((m + S(1))/n)

cons228 = CustomConstraint(cons_f228)

def cons_f229(n, m):
    return Not(IntegerQ((m + S(1))/n))

cons229 = CustomConstraint(cons_f229)

def cons_f230(n):
    return NegQ(n)

cons230 = CustomConstraint(cons_f230)

def cons_f231(n, p, m):
    return ZeroQ(p + S(1) + (m + S(1))/n)

cons231 = CustomConstraint(cons_f231)

def cons_f232(n, p, m):
    return ZeroQ(p + S(1) + (m + S(1))/(S(2)*n))

cons232 = CustomConstraint(cons_f232)

def cons_f233(n, m):
    return IntegerQ((m + S(1))/(S(2)*n))

cons233 = CustomConstraint(cons_f233)

def cons_f234(n, p, m):
    return NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n)

cons234 = CustomConstraint(cons_f234)

def cons_f235(n, p, m):
    return NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n))

cons235 = CustomConstraint(cons_f235)

def cons_f236(p, m):
    return RationalQ(m, p)

cons236 = CustomConstraint(cons_f236)

def cons_f237(n, p, m):
    return Not(NegativeIntegerQ((m + n*p + n + S(1))/n))

cons237 = CustomConstraint(cons_f237)

def cons_f238(b, p, a, m, n, x, c):
    return IntBinomialQ(a, b, c, n, m, p, x)

cons238 = CustomConstraint(cons_f238)

def cons_f239(n, p, m):
    return NonzeroQ(m + S(2)*n*p + S(1))

cons239 = CustomConstraint(cons_f239)

def cons_f240(a1, p, m, b2, n, x, c, a2, b1):
    return IntBinomialQ(a1*a2, b1*b2, c, n, m, p, x)

cons240 = CustomConstraint(cons_f240)

def cons_f241(n, p, m):
    return NonzeroQ(m + n*p + S(1))

cons241 = CustomConstraint(cons_f241)

def cons_f242(m):
    return PositiveIntegerQ(m/S(4) + S(-1)/2)

cons242 = CustomConstraint(cons_f242)

def cons_f243(m):
    return NegativeIntegerQ(m/S(4) + S(-1)/2)

cons243 = CustomConstraint(cons_f243)

def cons_f244(m):
    return IntegerQ(S(2)*m)

cons244 = CustomConstraint(cons_f244)

def cons_f245(m):
    return Greater(m, S(3)/2)

cons245 = CustomConstraint(cons_f245)

def cons_f246(n, m):
    return Greater(m + S(1), n)

cons246 = CustomConstraint(cons_f246)

def cons_f247(n, p, m):
    return Not(NegativeIntegerQ((m + n*(p + S(1)) + S(1))/n))

cons247 = CustomConstraint(cons_f247)

def cons_f248(n, m):
    return Greater(m + S(1), S(2)*n)

cons248 = CustomConstraint(cons_f248)

def cons_f249(n, p, m):
    return Not(NegativeIntegerQ((m + S(2)*n*(p + S(1)) + S(1))/(S(2)*n)))

cons249 = CustomConstraint(cons_f249)

def cons_f250(n):
    return PositiveIntegerQ(n/S(2) + S(-1)/2)

cons250 = CustomConstraint(cons_f250)

def cons_f251(n, m):
    return Less(m, n + S(-1))

cons251 = CustomConstraint(cons_f251)

def cons_f252(n, m):
    return PositiveIntegerQ(m, n/S(2) + S(-1)/2)

cons252 = CustomConstraint(cons_f252)

def cons_f253(n, m):
    return PositiveIntegerQ(m, n/S(4) + S(-1)/2)

cons253 = CustomConstraint(cons_f253)

def cons_f254(n, m):
    return PositiveIntegerQ(m, n/S(4))

cons254 = CustomConstraint(cons_f254)

def cons_f255(n, m):
    return Less(m, n/S(2))

cons255 = CustomConstraint(cons_f255)

def cons_f256(n, m):
    return Inequality(n/S(2), LessEqual, m, Less, n)

cons256 = CustomConstraint(cons_f256)

def cons_f257(n, m):
    return PositiveIntegerQ(m, n)

cons257 = CustomConstraint(cons_f257)

def cons_f258(n, m):
    return Greater(m, S(2)*n + S(-1))

cons258 = CustomConstraint(cons_f258)

def cons_f259(n, m):
    return Greater(m, n + S(-1))

cons259 = CustomConstraint(cons_f259)

def cons_f260(n, m):
    return _SumSimplerQ(m, -n)

cons260 = CustomConstraint(cons_f260)

def cons_f261(n, p, m):
    return NegativeIntegerQ((m + n*p + S(1))/n)

cons261 = CustomConstraint(cons_f261)

def cons_f262(n, m):
    return _SumSimplerQ(m, -S(2)*n)

cons262 = CustomConstraint(cons_f262)

def cons_f263(n, p, m):
    return NegativeIntegerQ((m + S(2)*n*p + S(1))/(S(2)*n))

cons263 = CustomConstraint(cons_f263)

def cons_f264(n, m):
    return _SumSimplerQ(m, n)

cons264 = CustomConstraint(cons_f264)

def cons_f265(n, m):
    return _SumSimplerQ(m, S(2)*n)

cons265 = CustomConstraint(cons_f265)

def cons_f266(m):
    return FractionQ(m)

cons266 = CustomConstraint(cons_f266)

def cons_f267(n, p, m):
    return IntegersQ(m, p + (m + S(1))/n)

cons267 = CustomConstraint(cons_f267)

def cons_f268(n, p, m):
    return IntegersQ(m, p + (m + S(1))/(S(2)*n))

cons268 = CustomConstraint(cons_f268)

def cons_f269(n, p, m):
    return Less(Denominator(p + (m + S(1))/n), Denominator(p))

cons269 = CustomConstraint(cons_f269)

def cons_f270(n, p, m):
    return Less(Denominator(p + (m + S(1))/(S(2)*n)), Denominator(p))

cons270 = CustomConstraint(cons_f270)

def cons_f271(m):
    return Not(RationalQ(m))

cons271 = CustomConstraint(cons_f271)

def cons_f272(n, m):
    return IntegerQ(n/(m + S(1)))

cons272 = CustomConstraint(cons_f272)

def cons_f273(n, m):
    return IntegerQ(S(2)*n/(m + S(1)))

cons273 = CustomConstraint(cons_f273)

def cons_f274(n):
    return Not(IntegerQ(S(2)*n))

cons274 = CustomConstraint(cons_f274)

def cons_f275(n, p, m):
    return ZeroQ(p + (m + S(1))/n)

cons275 = CustomConstraint(cons_f275)

def cons_f276(n, p, m):
    return ZeroQ(p + (m + S(1))/(S(2)*n))

cons276 = CustomConstraint(cons_f276)

def cons_f277(n, p, m):
    return IntegerQ(p + (m + S(1))/n)

cons277 = CustomConstraint(cons_f277)

def cons_f278(n, p, m):
    return IntegerQ(p + (m + S(1))/(S(2)*n))

cons278 = CustomConstraint(cons_f278)

def cons_f279(n, m):
    return FractionQ((m + S(1))/n)

cons279 = CustomConstraint(cons_f279)

def cons_f280(n, m):
    return Or(_SumSimplerQ(m, n), _SumSimplerQ(m, -n))

cons280 = CustomConstraint(cons_f280)

def cons_f281(p, a):
    return Or(NegativeIntegerQ(p), PositiveQ(a))

cons281 = CustomConstraint(cons_f281)

def cons_f282(p, a):
    return Not(Or(NegativeIntegerQ(p), PositiveQ(a)))

cons282 = CustomConstraint(cons_f282)

def cons_f283(x, v):
    return LinearQ(v, x)

cons283 = CustomConstraint(cons_f283)

def cons_f284(x, v):
    return NonzeroQ(v - x)

cons284 = CustomConstraint(cons_f284)

def cons_f285(u, x, v):
    return LinearPairQ(u, v, x)

cons285 = CustomConstraint(cons_f285)

def cons_f286(p, q):
    return PositiveIntegerQ(p, q)

cons286 = CustomConstraint(cons_f286)

def cons_f287(n, p):
    return ZeroQ(n*p + S(1))

cons287 = CustomConstraint(cons_f287)

def cons_f288(n, p, q):
    return ZeroQ(n*(p + q + S(1)) + S(1))

cons288 = CustomConstraint(cons_f288)

def cons_f289(q):
    return RationalQ(q)

cons289 = CustomConstraint(cons_f289)

def cons_f290(q):
    return Greater(q, S(0))

cons290 = CustomConstraint(cons_f290)

def cons_f291(p):
    return NegativeIntegerQ(p)

cons291 = CustomConstraint(cons_f291)

def cons_f292(n, p, q):
    return ZeroQ(n*(p + q + S(2)) + S(1))

cons292 = CustomConstraint(cons_f292)

def cons_f293(b, p, a, q, d, c):
    return ZeroQ(a*d*(p + S(1)) + b*c*(q + S(1)))

cons293 = CustomConstraint(cons_f293)

def cons_f294(p, q):
    return Or(And(RationalQ(p), Less(p, S(-1))), Not(And(RationalQ(q), Less(q, S(-1)))))

cons294 = CustomConstraint(cons_f294)

def cons_f295(b, p, a, n, d, c):
    return ZeroQ(a*d - b*c*(n*(p + S(1)) + S(1)))

cons295 = CustomConstraint(cons_f295)

def cons_f296(n, p):
    return Or(And(RationalQ(p), Less(p, S(-1))), NegativeIntegerQ(p + S(1)/n))

cons296 = CustomConstraint(cons_f296)

def cons_f297(n, p):
    return NonzeroQ(n*(p + S(1)) + S(1))

cons297 = CustomConstraint(cons_f297)

def cons_f298(q):
    return NegativeIntegerQ(q)

cons298 = CustomConstraint(cons_f298)

def cons_f299(p, q):
    return GreaterEqual(p, -q)

cons299 = CustomConstraint(cons_f299)

def cons_f300(c, b, a, d):
    return ZeroQ(S(3)*a*d + b*c)

cons300 = CustomConstraint(cons_f300)

def cons_f301(p):
    return Or(Equal(p, S(1)/2), Equal(Denominator(p), S(4)))

cons301 = CustomConstraint(cons_f301)

def cons_f302(p):
    return Equal(Denominator(p), S(4))

cons302 = CustomConstraint(cons_f302)

def cons_f303(p):
    return Or(Equal(p, S(-5)/4), Equal(p, S(-7)/4))

cons303 = CustomConstraint(cons_f303)

def cons_f304(a, b):
    return PosQ(a*b)

cons304 = CustomConstraint(cons_f304)

def cons_f305(a, b):
    return NegQ(a*b)

cons305 = CustomConstraint(cons_f305)

def cons_f306(p):
    return Or(Equal(p, S(3)/4), Equal(p, S(5)/4))

cons306 = CustomConstraint(cons_f306)

def cons_f307(d, c):
    return PosQ(d/c)

cons307 = CustomConstraint(cons_f307)

def cons_f308(p, q):
    return RationalQ(p, q)

cons308 = CustomConstraint(cons_f308)

def cons_f309(q):
    return Less(S(0), q, S(1))

cons309 = CustomConstraint(cons_f309)

def cons_f310(b, p, a, q, n, x, d, c):
    return IntBinomialQ(a, b, c, d, n, p, q, x)

cons310 = CustomConstraint(cons_f310)

def cons_f311(q):
    return Greater(q, S(1))

cons311 = CustomConstraint(cons_f311)

def cons_f312(p, q):
    return Not(And(Not(IntegerQ(p)), IntegerQ(q), Less(q, S(-1))))

cons312 = CustomConstraint(cons_f312)

def cons_f313(p, q):
    return Greater(p + q, S(0))

cons313 = CustomConstraint(cons_f313)

def cons_f314(n, p, q):
    return NonzeroQ(n*(p + q) + S(1))

cons314 = CustomConstraint(cons_f314)

def cons_f315(p):
    return Not(And(IntegerQ(p), Greater(p, S(1))))

cons315 = CustomConstraint(cons_f315)

def cons_f316(c, a, b, d):
    return Not(SimplerSqrtQ(b/a, d/c))

cons316 = CustomConstraint(cons_f316)

def cons_f317(d, c):
    return NegQ(d/c)

cons317 = CustomConstraint(cons_f317)

def cons_f318(c, a, b, d):
    return Not(And(NegQ(b/a), SimplerSqrtQ(-b/a, -d/c)))

cons318 = CustomConstraint(cons_f318)

def cons_f319(a, b, d, c):
    return PositiveQ(a - b*c/d)

cons319 = CustomConstraint(cons_f319)

def cons_f320(n):
    return NonzeroQ(n + S(1))

cons320 = CustomConstraint(cons_f320)

def cons_f321(n, mn):
    return EqQ(mn, -n)

cons321 = CustomConstraint(cons_f321)

def cons_f322(q):
    return IntegerQ(q)

cons322 = CustomConstraint(cons_f322)

def cons_f323(n, p):
    return Or(PosQ(n), Not(IntegerQ(p)))

cons323 = CustomConstraint(cons_f323)

def cons_f324(q):
    return Not(IntegerQ(q))

cons324 = CustomConstraint(cons_f324)

def cons_f325(u, x, v):
    return PseudoBinomialPairQ(u, v, x)

cons325 = CustomConstraint(cons_f325)

def cons_f326(p, m):
    return IntegersQ(p, m/p)

cons326 = CustomConstraint(cons_f326)

def cons_f327(u, p, m, x, v):
    return PseudoBinomialPairQ(u*x**(m/p), v, x)

cons327 = CustomConstraint(cons_f327)

def cons_f328(e, m):
    return Or(IntegerQ(m), PositiveQ(e))

cons328 = CustomConstraint(cons_f328)

def cons_f329(n, m):
    return ZeroQ(m - n + S(1))

cons329 = CustomConstraint(cons_f329)

def cons_f330(d, b, p, a, m, n, c):
    return ZeroQ(a*d*(m + S(1)) - b*c*(m + n*(p + S(1)) + S(1)))

cons330 = CustomConstraint(cons_f330)

def cons_f331(n, non2):
    return ZeroQ(-n/S(2) + non2)

cons331 = CustomConstraint(cons_f331)

def cons_f332(a1, p, c, m, b2, n, d, a2, b1):
    return ZeroQ(a1*a2*d*(m + S(1)) - b1*b2*c*(m + n*(p + S(1)) + S(1)))

cons332 = CustomConstraint(cons_f332)

def cons_f333(n, p, m):
    return ZeroQ(m + n*(p + S(1)) + S(1))

cons333 = CustomConstraint(cons_f333)

def cons_f334(n, e):
    return Or(IntegerQ(n), PositiveQ(e))

cons334 = CustomConstraint(cons_f334)

def cons_f335(n, m):
    return Or(And(Greater(n, S(0)), Less(m, S(-1))), And(Less(n, S(0)), Greater(m + n, S(-1))))

cons335 = CustomConstraint(cons_f335)

def cons_f336(p):
    return Not(And(IntegerQ(p), Less(p, S(-1))))

cons336 = CustomConstraint(cons_f336)

def cons_f337(m):
    return PositiveIntegerQ(m/S(2))

cons337 = CustomConstraint(cons_f337)

def cons_f338(p, m):
    return Or(IntegerQ(p), Equal(m + S(2)*p + S(1), S(0)))

cons338 = CustomConstraint(cons_f338)

def cons_f339(m):
    return NegativeIntegerQ(m/S(2))

cons339 = CustomConstraint(cons_f339)

def cons_f340(n, p, m):
    return Or(IntegerQ(p), Not(RationalQ(m)), And(PositiveIntegerQ(n), NegativeIntegerQ(p + S(1)/2), LessEqual(S(-1), m, -n*(p + S(1)))))

cons340 = CustomConstraint(cons_f340)

def cons_f341(n, p, m):
    return NonzeroQ(m + n*(p + S(1)) + S(1))

cons341 = CustomConstraint(cons_f341)

def cons_f342(m):
    return Or(IntegerQ(m), PositiveIntegerQ(S(2)*m + S(2)), Not(RationalQ(m)))

cons342 = CustomConstraint(cons_f342)

def cons_f343(n, p, m):
    return NonzeroQ(m + n*(p + S(2)) + S(1))

cons343 = CustomConstraint(cons_f343)

def cons_f344(p, m, q):
    return RationalQ(m, p, q)

cons344 = CustomConstraint(cons_f344)

def cons_f345(n, m):
    return Greater(m - n + S(1), S(0))

cons345 = CustomConstraint(cons_f345)

def cons_f346(b, p, e, a, m, q, n, x, d, c):
    return IntBinomialQ(a, b, c, d, e, m, n, p, q, x)

cons346 = CustomConstraint(cons_f346)

def cons_f347(n, m):
    return Greater(m - n + S(1), n)

cons347 = CustomConstraint(cons_f347)

def cons_f348(n, m):
    return Inequality(n, GreaterEqual, m - n + S(1), Greater, S(0))

cons348 = CustomConstraint(cons_f348)

def cons_f349(m, q):
    return RationalQ(m, q)

cons349 = CustomConstraint(cons_f349)

def cons_f350(n, m):
    return LessEqual(n, m, S(2)*n + S(-1))

cons350 = CustomConstraint(cons_f350)

def cons_f351(n, m):
    return IntegersQ(m/S(2), n/S(2))

cons351 = CustomConstraint(cons_f351)

def cons_f352(n, m):
    return Less(S(0), m - n + S(1), n)

cons352 = CustomConstraint(cons_f352)

def cons_f353(n):
    return LessEqual(n, S(4))

cons353 = CustomConstraint(cons_f353)

def cons_f354(c, b, a, d):
    return ZeroQ(-a*d + S(4)*b*c)

cons354 = CustomConstraint(cons_f354)

def cons_f355(m):
    return PositiveIntegerQ(m/S(3) + S(-1)/3)

cons355 = CustomConstraint(cons_f355)

def cons_f356(m):
    return NegativeIntegerQ(m/S(3) + S(-1)/3)

cons356 = CustomConstraint(cons_f356)

def cons_f357(m):
    return IntegerQ(m/S(3) + S(-1)/3)

cons357 = CustomConstraint(cons_f357)

def cons_f358(n):
    return Or(EqQ(n, S(2)), EqQ(n, S(4)))

cons358 = CustomConstraint(cons_f358)

def cons_f359(b, d, a, n, c):
    return Not(And(EqQ(n, S(2)), SimplerSqrtQ(-b/a, -d/c)))

cons359 = CustomConstraint(cons_f359)

def cons_f360(n, p, m, q):
    return IntegersQ(p + (m + S(1))/n, q)

cons360 = CustomConstraint(cons_f360)

def cons_f361(n, m):
    return Or(ZeroQ(m - n), ZeroQ(m - S(2)*n + S(1)))

cons361 = CustomConstraint(cons_f361)

def cons_f362(p, m, q):
    return IntegersQ(m, p, q)

cons362 = CustomConstraint(cons_f362)

def cons_f363(p):
    return GreaterEqual(p, S(-2))

cons363 = CustomConstraint(cons_f363)

def cons_f364(m, q):
    return Or(GreaterEqual(q, S(-2)), And(Equal(q, S(-3)), IntegerQ(m/S(2) + S(-1)/2)))

cons364 = CustomConstraint(cons_f364)

def cons_f365(n, m):
    return NonzeroQ(m - n + S(1))

cons365 = CustomConstraint(cons_f365)

def cons_f366(n2, n):
    return ZeroQ(-S(2)*n + n2)

cons366 = CustomConstraint(cons_f366)

def cons_f367(p, q, r):
    return PositiveIntegerQ(p, q, r)

cons367 = CustomConstraint(cons_f367)

def cons_f368(b, f, e, a, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, n), x)

cons368 = CustomConstraint(cons_f368)

def cons_f369(b, d, a, n, c):
    return Not(And(ZeroQ(n + S(-2)), Or(And(PosQ(b/a), PosQ(d/c)), And(NegQ(b/a), Or(PosQ(d/c), And(PositiveQ(a), Or(Not(PositiveQ(c)), SimplerSqrtQ(-b/a, -d/c))))))))

cons369 = CustomConstraint(cons_f369)

def cons_f370(n, p, q):
    return NonzeroQ(n*(p + q + S(1)) + S(1))

cons370 = CustomConstraint(cons_f370)

def cons_f371(b, f, p, e, a, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, p, n), x)

cons371 = CustomConstraint(cons_f371)

def cons_f372(b, f, p, e, a, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, n, p, q), x)

cons372 = CustomConstraint(cons_f372)

def cons_f373(d, c):
    return PositiveQ(d/c)

cons373 = CustomConstraint(cons_f373)

def cons_f374(f, e):
    return PositiveQ(f/e)

cons374 = CustomConstraint(cons_f374)

def cons_f375(f, d, e, c):
    return Not(SimplerSqrtQ(d/c, f/e))

cons375 = CustomConstraint(cons_f375)

def cons_f376(f, d, e, c):
    return Not(SimplerSqrtQ(-f/e, -d/c))

cons376 = CustomConstraint(cons_f376)

def cons_f377(f, e):
    return PosQ(f/e)

cons377 = CustomConstraint(cons_f377)

def cons_f378(f, d, e, c):
    return Not(And(NegQ(f/e), SimplerSqrtQ(-f/e, -d/c)))

cons378 = CustomConstraint(cons_f378)

def cons_f379(q, r):
    return RationalQ(q, r)

cons379 = CustomConstraint(cons_f379)

def cons_f380(q):
    return Less(q, S(-1))

cons380 = CustomConstraint(cons_f380)

def cons_f381(r):
    return Greater(r, S(1))

cons381 = CustomConstraint(cons_f381)

def cons_f382(q):
    return LessEqual(q, S(-1))

cons382 = CustomConstraint(cons_f382)

def cons_f383(f, e, d, c):
    return PosQ((-c*f + d*e)/c)

cons383 = CustomConstraint(cons_f383)

def cons_f384(f, e, d, c):
    return NegQ((-c*f + d*e)/c)

cons384 = CustomConstraint(cons_f384)

def cons_f385(b, r, f, p, e, a, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, n, p, q, r), x)

cons385 = CustomConstraint(cons_f385)

def cons_f386(u, v):
    return ZeroQ(u - v)

cons386 = CustomConstraint(cons_f386)

def cons_f387(u, w):
    return ZeroQ(u - w)

cons387 = CustomConstraint(cons_f387)

def cons_f388(r):
    return IntegerQ(r)

cons388 = CustomConstraint(cons_f388)

def cons_f389(n2, n):
    return ZeroQ(-n/S(2) + n2)

cons389 = CustomConstraint(cons_f389)

def cons_f390(f1, e2, e1, f2):
    return ZeroQ(e1*f2 + e2*f1)

cons390 = CustomConstraint(cons_f390)

def cons_f391(e2, e1, r):
    return Or(IntegerQ(r), And(PositiveQ(e1), PositiveQ(e2)))

cons391 = CustomConstraint(cons_f391)

def cons_f392(e1, x):
    return FreeQ(e1, x)

cons392 = CustomConstraint(cons_f392)

def cons_f393(f1, x):
    return FreeQ(f1, x)

cons393 = CustomConstraint(cons_f393)

def cons_f394(e2, x):
    return FreeQ(e2, x)

cons394 = CustomConstraint(cons_f394)

def cons_f395(f2, x):
    return FreeQ(f2, x)

cons395 = CustomConstraint(cons_f395)

def cons_f396(g, m):
    return Or(IntegerQ(m), PositiveQ(g))

cons396 = CustomConstraint(cons_f396)

def cons_f397(p, q, r):
    return PositiveIntegerQ(p + S(2), q, r)

cons397 = CustomConstraint(cons_f397)

def cons_f398(p, q, r):
    return IntegersQ(p, q, r)

cons398 = CustomConstraint(cons_f398)

def cons_f399(b, d, f, e, a, q, c):
    return Not(And(Equal(q, S(1)), SimplerQ(-a*d + b*c, -a*f + b*e)))

cons399 = CustomConstraint(cons_f399)

def cons_f400(d, f, e, q, n, x, c):
    return Not(And(Equal(q, S(1)), SimplerQ(e + f*x**n, c + d*x**n)))

cons400 = CustomConstraint(cons_f400)

def cons_f401(r):
    return PositiveIntegerQ(r)

cons401 = CustomConstraint(cons_f401)

def cons_f402(g, b, f, p, e, a, m, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p), x)

cons402 = CustomConstraint(cons_f402)

def cons_f403(g, b, f, p, e, a, m, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q), x)

cons403 = CustomConstraint(cons_f403)

def cons_f404(g, b, r, f, p, e, a, m, q, n, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, m, n, p, q, r), x)

cons404 = CustomConstraint(cons_f404)

def cons_f405(a, b, c):
    return ZeroQ(-S(4)*a*c + b**S(2))

cons405 = CustomConstraint(cons_f405)

def cons_f406(p):
    return NonzeroQ(S(2)*p + S(1))

cons406 = CustomConstraint(cons_f406)

def cons_f407(a, b, c):
    return NonzeroQ(-S(4)*a*c + b**S(2))

cons407 = CustomConstraint(cons_f407)

def cons_f408(a, b, c):
    return PerfectSquareQ(-S(4)*a*c + b**S(2))

cons408 = CustomConstraint(cons_f408)

def cons_f409(a, b, c):
    return Not(PerfectSquareQ(-S(4)*a*c + b**S(2)))

cons409 = CustomConstraint(cons_f409)

def cons_f410(p):
    return IntegerQ(S(4)*p)

cons410 = CustomConstraint(cons_f410)

def cons_f411(p):
    return Unequal(p, S(-3)/2)

cons411 = CustomConstraint(cons_f411)

def cons_f412(a, b, c):
    return PosQ(-S(4)*a*c + b**S(2))

cons412 = CustomConstraint(cons_f412)

def cons_f413(a, b, c):
    return PositiveQ(S(4)*a - b**S(2)/c)

cons413 = CustomConstraint(cons_f413)

def cons_f414(x, c, b):
    return FreeQ(List(b, c), x)

cons414 = CustomConstraint(cons_f414)

def cons_f415(p):
    return LessEqual(S(3), Denominator(p), S(4))

cons415 = CustomConstraint(cons_f415)

def cons_f416(p):
    return Not(IntegerQ(S(4)*p))

cons416 = CustomConstraint(cons_f416)

def cons_f417(e, c, d, b):
    return ZeroQ(-b*e + S(2)*c*d)

cons417 = CustomConstraint(cons_f417)

def cons_f418(m):
    return IntegerQ(m/S(2) + S(1)/2)

cons418 = CustomConstraint(cons_f418)

def cons_f419(p, m):
    return ZeroQ(m + S(2)*p + S(1))

cons419 = CustomConstraint(cons_f419)

def cons_f420(p, m):
    return NonzeroQ(m + S(2)*p + S(1))

cons420 = CustomConstraint(cons_f420)

def cons_f421(e, c, d, b):
    return NonzeroQ(-b*e + S(2)*c*d)

cons421 = CustomConstraint(cons_f421)

def cons_f422(p, m):
    return ZeroQ(m + S(2)*p + S(2))

cons422 = CustomConstraint(cons_f422)

def cons_f423(m):
    return NonzeroQ(m + S(2))

cons423 = CustomConstraint(cons_f423)

def cons_f424(p, m):
    return ZeroQ(m + S(2)*p + S(3))

cons424 = CustomConstraint(cons_f424)

def cons_f425(p):
    return NonzeroQ(p + S(3)/2)

cons425 = CustomConstraint(cons_f425)

def cons_f426(m):
    return Inequality(S(-2), LessEqual, m, Less, S(-1))

cons426 = CustomConstraint(cons_f426)

def cons_f427(p):
    return IntegerQ(S(2)*p)

cons427 = CustomConstraint(cons_f427)

def cons_f428(m):
    return Less(m, S(-2))

cons428 = CustomConstraint(cons_f428)

def cons_f429(p, m):
    return Not(And(NegativeIntegerQ(m + S(2)*p + S(3)), Greater(m + S(3)*p + S(3), S(0))))

cons429 = CustomConstraint(cons_f429)

def cons_f430(p, m):
    return NonzeroQ(m + S(2)*p)

cons430 = CustomConstraint(cons_f430)

def cons_f431(m):
    return Not(And(RationalQ(m), Less(m, S(-2))))

cons431 = CustomConstraint(cons_f431)

def cons_f432(p, m):
    return Not(And(IntegerQ(m), Less(S(0), m, S(2)*p)))

cons432 = CustomConstraint(cons_f432)

def cons_f433(m):
    return Inequality(S(0), Less, m, LessEqual, S(1))

cons433 = CustomConstraint(cons_f433)

def cons_f434(p, m):
    return NonzeroQ(m + p + S(1))

cons434 = CustomConstraint(cons_f434)

def cons_f435(p, m):
    return Or(Not(RationalQ(p)), Inequality(S(-1), LessEqual, p, Less, S(0)), And(IntegerQ(m), Less(S(0), m, S(2)*p)), And(Equal(m, S(1)/2), Less(p, S(0))))

cons435 = CustomConstraint(cons_f435)

def cons_f436(p, m):
    return Or(IntegerQ(m), IntegerQ(S(2)*p))

cons436 = CustomConstraint(cons_f436)

def cons_f437(d, b, e, a, c):
    return ZeroQ(a*e**S(2) - b*d*e + c*d**S(2))

cons437 = CustomConstraint(cons_f437)

def cons_f438(e, c, d, a):
    return ZeroQ(a*e**S(2) + c*d**S(2))

cons438 = CustomConstraint(cons_f438)

def cons_f439(p, a, d, m):
    return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), IntegerQ(m + p)))

cons439 = CustomConstraint(cons_f439)

def cons_f440(p, m):
    return ZeroQ(m + p)

cons440 = CustomConstraint(cons_f440)

def cons_f441(p, m):
    return Or(Less(S(0), -m, p), Less(p, -m, S(0)))

cons441 = CustomConstraint(cons_f441)

def cons_f442(m):
    return Unequal(m, S(2))

cons442 = CustomConstraint(cons_f442)

def cons_f443(m):
    return Unequal(m, S(-1))

cons443 = CustomConstraint(cons_f443)

def cons_f444(p, m):
    return PositiveIntegerQ(m + p)

cons444 = CustomConstraint(cons_f444)

def cons_f445(p, m):
    return NegativeIntegerQ(m + S(2)*p + S(2))

cons445 = CustomConstraint(cons_f445)

def cons_f446(p, m):
    return Or(Less(m, S(-2)), ZeroQ(m + S(2)*p + S(1)))

cons446 = CustomConstraint(cons_f446)

def cons_f447(p, m):
    return Or(Inequality(S(-2), LessEqual, m, Less, S(0)), Equal(m + p + S(1), S(0)))

cons447 = CustomConstraint(cons_f447)

def cons_f448(m):
    return GreaterEqual(m, S(1))

cons448 = CustomConstraint(cons_f448)

def cons_f449(m):
    return Less(m, S(0))

cons449 = CustomConstraint(cons_f449)

def cons_f450(d):
    return PositiveQ(d)

cons450 = CustomConstraint(cons_f450)

def cons_f451(p, m):
    return Not(And(ZeroQ(m + S(-3)), Unequal(p, S(1))))

cons451 = CustomConstraint(cons_f451)

def cons_f452(p, m):
    return NonzeroQ(m + S(2)*p + S(3))

cons452 = CustomConstraint(cons_f452)

def cons_f453(p, m):
    return Not(And(EvenQ(m), Less(m + S(2)*p + S(3), S(0))))

cons453 = CustomConstraint(cons_f453)

def cons_f454(m):
    return Not(And(RationalQ(m), Less(m, S(-1))))

cons454 = CustomConstraint(cons_f454)

def cons_f455(p, m):
    return Not(And(PositiveIntegerQ(m/S(2) + S(-1)/2), Or(Not(IntegerQ(p)), Less(m, S(2)*p))))

cons455 = CustomConstraint(cons_f455)

def cons_f456(m):
    return Not(And(RationalQ(m), Greater(m, S(1))))

cons456 = CustomConstraint(cons_f456)

def cons_f457(c, b, a):
    return NegativeQ(c/(-S(4)*a*c + b**S(2)))

cons457 = CustomConstraint(cons_f457)

def cons_f458(m):
    return EqQ(m**S(2), S(1)/4)

cons458 = CustomConstraint(cons_f458)

def cons_f459(p, m):
    return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), OddQ(m))

cons459 = CustomConstraint(cons_f459)

def cons_f460(p, m):
    return Or(IntegerQ(S(2)*p), And(IntegerQ(m), RationalQ(p)), IntegerQ(m/S(2) + p + S(3)/2))

cons460 = CustomConstraint(cons_f460)

def cons_f461(d, b, e, a, c):
    return NonzeroQ(a*e**S(2) - b*d*e + c*d**S(2))

cons461 = CustomConstraint(cons_f461)

def cons_f462(e, c, d, a):
    return NonzeroQ(a*e**S(2) + c*d**S(2))

cons462 = CustomConstraint(cons_f462)

def cons_f463(p, m):
    return Not(And(ZeroQ(m + S(-1)), Greater(p, S(1))))

cons463 = CustomConstraint(cons_f463)

def cons_f464(a, b, c):
    return NiceSqrtQ(-S(4)*a*c + b**S(2))

cons464 = CustomConstraint(cons_f464)

def cons_f465(a, c):
    return NiceSqrtQ(-a*c)

cons465 = CustomConstraint(cons_f465)

def cons_f466(a, b, c):
    return Not(NiceSqrtQ(-S(4)*a*c + b**S(2)))

cons466 = CustomConstraint(cons_f466)

def cons_f467(a, c):
    return Not(NiceSqrtQ(-a*c))

cons467 = CustomConstraint(cons_f467)

def cons_f468(d, m):
    return Or(NonzeroQ(d), Greater(m, S(2)))

cons468 = CustomConstraint(cons_f468)

def cons_f469(p):
    return Not(And(RationalQ(p), LessEqual(p, S(-1))))

cons469 = CustomConstraint(cons_f469)

def cons_f470(e, d, b, a):
    return ZeroQ(a*e + b*d)

cons470 = CustomConstraint(cons_f470)

def cons_f471(e, c, d, b):
    return ZeroQ(b*e + c*d)

cons471 = CustomConstraint(cons_f471)

def cons_f472(p, m):
    return PositiveIntegerQ(m - p + S(1))

cons472 = CustomConstraint(cons_f472)

def cons_f473(e, c, d, b):
    return NonzeroQ(-b*e + c*d)

cons473 = CustomConstraint(cons_f473)

def cons_f474(m):
    return Equal(m**S(2), S(1)/4)

cons474 = CustomConstraint(cons_f474)

def cons_f475(c):
    return NegativeQ(c)

cons475 = CustomConstraint(cons_f475)

def cons_f476(b):
    return RationalQ(b)

cons476 = CustomConstraint(cons_f476)

def cons_f477(m):
    return ZeroQ(m**S(2) + S(-1)/4)

cons477 = CustomConstraint(cons_f477)

def cons_f478(p, m):
    return Equal(m + S(2)*p + S(2), S(0))

cons478 = CustomConstraint(cons_f478)

def cons_f479(e, a, x, d, c):
    return FreeQ(List(a, c, d, e), x)

cons479 = CustomConstraint(cons_f479)

def cons_f480(p, m):
    return Or(IntegerQ(p), And(RationalQ(m), Less(m, S(-1))))

cons480 = CustomConstraint(cons_f480)

def cons_f481(p, m):
    return Not(NegativeIntegerQ(m + S(2)*p + S(1)))

cons481 = CustomConstraint(cons_f481)

def cons_f482(b, p, e, a, m, x, d, c):
    return IntQuadraticQ(a, b, c, d, e, m, p, x)

cons482 = CustomConstraint(cons_f482)

def cons_f483(p, e, a, m, x, d, c):
    return IntQuadraticQ(a, S(0), c, d, e, m, p, x)

cons483 = CustomConstraint(cons_f483)

def cons_f484(m):
    return Or(Not(RationalQ(m)), Less(m, S(1)))

cons484 = CustomConstraint(cons_f484)

def cons_f485(p, m):
    return Not(NegativeIntegerQ(m + S(2)*p))

cons485 = CustomConstraint(cons_f485)

def cons_f486(p, m):
    return Or(Less(m, S(1)), And(NegativeIntegerQ(m + S(2)*p + S(3)), Unequal(m, S(2))))

cons486 = CustomConstraint(cons_f486)

def cons_f487(m):
    return If(RationalQ(m), Greater(m, S(1)), _SumSimplerQ(m, S(-2)))

cons487 = CustomConstraint(cons_f487)

def cons_f488(b, p, e, a, m, x, d, c):
    return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, b, c, d, e, m, p, x)), And(_SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

cons488 = CustomConstraint(cons_f488)

def cons_f489(p, e, a, m, x, d, c):
    return Or(And(RationalQ(m), Less(m, S(-1)), IntQuadraticQ(a, S(0), c, d, e, m, p, x)), And(_SumSimplerQ(m, S(1)), IntegerQ(p), NonzeroQ(m + S(1))), And(NegativeIntegerQ(m + S(2)*p + S(3)), NonzeroQ(m + S(1))))

cons489 = CustomConstraint(cons_f489)

def cons_f490(d, b, e, a, c):
    return ZeroQ(-S(3)*a*c*e**S(2) + b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

cons490 = CustomConstraint(cons_f490)

def cons_f491(e, c, d, b):
    return PosQ(c*e**S(2)*(-b*e + S(2)*c*d))

cons491 = CustomConstraint(cons_f491)

def cons_f492(e, c, d, a):
    return ZeroQ(-S(3)*a*e**S(2) + c*d**S(2))

cons492 = CustomConstraint(cons_f492)

def cons_f493(e, c, d, b):
    return NegQ(c*e**S(2)*(-b*e + S(2)*c*d))

cons493 = CustomConstraint(cons_f493)

def cons_f494(d, b, e, a, c):
    return ZeroQ(S(9)*a*c*e**S(2) - S(2)*b**S(2)*e**S(2) - b*c*d*e + c**S(2)*d**S(2))

cons494 = CustomConstraint(cons_f494)

def cons_f495(a, b, c):
    return Not(PositiveQ(S(4)*a - b**S(2)/c))

cons495 = CustomConstraint(cons_f495)

def cons_f496(p):
    return Not(IntegerQ(S(2)*p))

cons496 = CustomConstraint(cons_f496)

def cons_f497(f, e, d, g):
    return NonzeroQ(-d*g + e*f)

cons497 = CustomConstraint(cons_f497)

def cons_f498(g, f, c, b):
    return ZeroQ(-b*g + S(2)*c*f)

cons498 = CustomConstraint(cons_f498)

def cons_f499(m):
    return Not(And(RationalQ(m), Greater(m, S(0))))

cons499 = CustomConstraint(cons_f499)

def cons_f500(p, m):
    return Or(Not(RationalQ(p)), And(Greater(p, S(0)), Or(Not(IntegerQ(m)), GreaterEqual(m, -S(2)*p + S(-2)), Less(m, -S(4)*p + S(-4)))))

cons500 = CustomConstraint(cons_f500)

def cons_f501(p, m):
    return NonzeroQ(m + S(2)*p + S(2))

cons501 = CustomConstraint(cons_f501)

def cons_f502(p, m):
    return Or(Not(RationalQ(p)), Less(m, S(2)*p + S(2)))

cons502 = CustomConstraint(cons_f502)

def cons_f503(g, f, c, b):
    return NonzeroQ(-b*g + S(2)*c*f)

cons503 = CustomConstraint(cons_f503)

def cons_f504(p):
    return Less(p, S(0))

cons504 = CustomConstraint(cons_f504)

def cons_f505(b, d, p, e, m, c):
    return Or(And(ZeroQ(m + S(2)*p + S(2)), NonzeroQ(m + S(1))), And(ZeroQ(-b*e + S(2)*c*d), NonzeroQ(m + S(-1))))

cons505 = CustomConstraint(cons_f505)

def cons_f506(g, f, e, m, x, d):
    return Not(And(ZeroQ(m + S(-1)), SimplerQ(f + g*x, d + e*x)))

cons506 = CustomConstraint(cons_f506)

def cons_f507(p, a, d, m):
    return Or(IntegerQ(p), And(PositiveQ(a), PositiveQ(d), ZeroQ(m + p)))

cons507 = CustomConstraint(cons_f507)

def cons_f508(g, b, f, p, e, m, d, c):
    return ZeroQ(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))

cons508 = CustomConstraint(cons_f508)

def cons_f509(g, f, p, e, m, d):
    return ZeroQ(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))

cons509 = CustomConstraint(cons_f509)

def cons_f510(m):
    return _SumSimplerQ(m, S(-1))

cons510 = CustomConstraint(cons_f510)

def cons_f511(p, m):
    return Or(And(RationalQ(m), Less(m, S(-1)), Not(PositiveIntegerQ(m + p + S(1)))), And(RationalQ(m, p), Less(m, S(0)), Less(p, S(-1))), ZeroQ(m + S(2)*p + S(2)))

cons511 = CustomConstraint(cons_f511)

def cons_f512(g, f, a, c):
    return ZeroQ(a*g**S(2) + c*f**S(2))

cons512 = CustomConstraint(cons_f512)

def cons_f513(p):
    return Less(p, S(-2))

cons513 = CustomConstraint(cons_f513)

def cons_f514(p, m):
    return Or(Less(S(0), -m, p + S(1)), Less(p, -m, S(0)))

cons514 = CustomConstraint(cons_f514)

def cons_f515(n, p):
    return NegativeIntegerQ(n + S(2)*p)

cons515 = CustomConstraint(cons_f515)

def cons_f516(g, d, b, f, e, c):
    return ZeroQ(-b*e*g + c*d*g + c*e*f)

cons516 = CustomConstraint(cons_f516)

def cons_f517(n, m):
    return NonzeroQ(m - n + S(-1))

cons517 = CustomConstraint(cons_f517)

def cons_f518(f, e, d, g):
    return ZeroQ(d*g + e*f)

cons518 = CustomConstraint(cons_f518)

def cons_f519(n, m):
    return ZeroQ(m - n + S(-2))

cons519 = CustomConstraint(cons_f519)

def cons_f520(n, p):
    return RationalQ(n, p)

cons520 = CustomConstraint(cons_f520)

def cons_f521(n, p):
    return Not(And(IntegerQ(n + p), LessEqual(n + p + S(2), S(0))))

cons521 = CustomConstraint(cons_f521)

def cons_f522(n):
    return Not(PositiveIntegerQ(n))

cons522 = CustomConstraint(cons_f522)

def cons_f523(n, p):
    return Not(And(IntegerQ(n + p), Less(n + p + S(2), S(0))))

cons523 = CustomConstraint(cons_f523)

def cons_f524(n, p):
    return Or(IntegerQ(S(2)*p), IntegerQ(n))

cons524 = CustomConstraint(cons_f524)

def cons_f525(p, m):
    return ZeroQ(m + p + S(-1))

cons525 = CustomConstraint(cons_f525)

def cons_f526(g, b, d, f, p, e, n, c):
    return ZeroQ(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))

cons526 = CustomConstraint(cons_f526)

def cons_f527(g, f, p, e, n, d):
    return ZeroQ(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))

cons527 = CustomConstraint(cons_f527)

def cons_f528(n):
    return Not(And(RationalQ(n), Less(n, S(-1))))

cons528 = CustomConstraint(cons_f528)

def cons_f529(p):
    return IntegerQ(p + S(-1)/2)

cons529 = CustomConstraint(cons_f529)

def cons_f530(p, m):
    return Not(And(Less(m, S(0)), Less(p, S(0))))

cons530 = CustomConstraint(cons_f530)

def cons_f531(p):
    return Unequal(p, S(1)/2)

cons531 = CustomConstraint(cons_f531)

def cons_f532(g, b, f, e, a, d, c):
    return ZeroQ(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)

cons532 = CustomConstraint(cons_f532)

def cons_f533(g, d, f, e, a, c):
    return ZeroQ(a*e*g + c*d*f)

cons533 = CustomConstraint(cons_f533)

def cons_f534(b, e, m, d, c):
    return Not(And(Equal(m, S(1)), Or(ZeroQ(d), ZeroQ(-b*e + S(2)*c*d))))

cons534 = CustomConstraint(cons_f534)

def cons_f535(d, m):
    return Not(And(Equal(m, S(1)), ZeroQ(d)))

cons535 = CustomConstraint(cons_f535)

def cons_f536(g, b, d, p, f, e, a, c):
    return ZeroQ(-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))

cons536 = CustomConstraint(cons_f536)

def cons_f537(g, f, p, e, a, d, c):
    return ZeroQ(a*e*g - c*d*f*(S(2)*p + S(3)))

cons537 = CustomConstraint(cons_f537)

def cons_f538(p, m):
    return ZeroQ(m - p)

cons538 = CustomConstraint(cons_f538)

def cons_f539(p, m):
    return Less(m + S(2)*p, S(0))

cons539 = CustomConstraint(cons_f539)

def cons_f540(p, m):
    return Not(NegativeIntegerQ(m + S(2)*p + S(3)))

cons540 = CustomConstraint(cons_f540)

def cons_f541(p, m):
    return Or(And(RationalQ(m), Less(m, S(-1))), Equal(p, S(1)), And(IntegerQ(p), Not(RationalQ(m))))

cons541 = CustomConstraint(cons_f541)

def cons_f542(p, m):
    return Or(IntegerQ(m), IntegerQ(p), IntegersQ(S(2)*m, S(2)*p))

cons542 = CustomConstraint(cons_f542)

def cons_f543(p, m):
    return Or(IntegerQ(p), Not(RationalQ(m)), Inequality(S(-1), LessEqual, m, Less, S(0)))

cons543 = CustomConstraint(cons_f543)

def cons_f544(g, b, d, p, f, e, a, m, c):
    return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, b, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

cons544 = CustomConstraint(cons_f544)

def cons_f545(g, p, f, e, a, m, d, c):
    return Or(And(Equal(m, S(2)), Equal(p, S(-3)), RationalQ(a, c, d, e, f, g)), Not(NegativeIntegerQ(m + S(2)*p + S(3))))

cons545 = CustomConstraint(cons_f545)

def cons_f546(g, f, e, m, x, d):
    return Not(And(Equal(m, S(1)), SimplerQ(d + e*x, f + g*x)))

cons546 = CustomConstraint(cons_f546)

def cons_f547(g, f, e, m, x, d):
    return Not(And(Equal(m, S(1)), SimplerQ(f + g*x, d + e*x)))

cons547 = CustomConstraint(cons_f547)

def cons_f548(p, m):
    return NegativeIntegerQ(m + S(2)*p + S(3))

cons548 = CustomConstraint(cons_f548)

def cons_f549(d, b, e, a, c):
    return ZeroQ(S(4)*c*(a - d) - (b - e)**S(2))

cons549 = CustomConstraint(cons_f549)

def cons_f550(g, b, f, e, a, d):
    return ZeroQ(e*f*(b - e) - S(2)*g*(-a*e + b*d))

cons550 = CustomConstraint(cons_f550)

def cons_f551(e, d, b, a):
    return NonzeroQ(-a*e + b*d)

cons551 = CustomConstraint(cons_f551)

def cons_f552(g, f, a, x, c):
    return FreeQ(List(a, c, f, g), x)

cons552 = CustomConstraint(cons_f552)

def cons_f553(g, f, e, a, x, c):
    return FreeQ(List(a, c, e, f, g), x)

cons553 = CustomConstraint(cons_f553)

def cons_f554(n, p, m):
    return IntegersQ(m, n, p)

cons554 = CustomConstraint(cons_f554)

def cons_f555(n, p):
    return IntegersQ(n, p)

cons555 = CustomConstraint(cons_f555)

def cons_f556(f, d, m):
    return Or(IntegerQ(m), And(PositiveQ(d), PositiveQ(f)))

cons556 = CustomConstraint(cons_f556)

def cons_f557(n, p, m):
    return Or(IntegerQ(p), IntegersQ(m, n))

cons557 = CustomConstraint(cons_f557)

def cons_f558(g, f, p, e, a, m, x, c):
    return FreeQ(List(a, c, e, f, g, m, p), x)

cons558 = CustomConstraint(cons_f558)

def cons_f559(g, f, p, e, a, m, n, x, d, c):
    return FreeQ(List(a, c, d, e, f, g, m, n, p), x)

cons559 = CustomConstraint(cons_f559)

def cons_f560(f, c, d, a):
    return ZeroQ(-a*f + c*d)

cons560 = CustomConstraint(cons_f560)

def cons_f561(e, d, b, a):
    return ZeroQ(-a*e + b*d)

cons561 = CustomConstraint(cons_f561)

def cons_f562(f, p, c):
    return Or(IntegerQ(p), PositiveQ(c/f))

cons562 = CustomConstraint(cons_f562)

def cons_f563(b, f, e, a, q, x, d, c):
    return Or(Not(IntegerQ(q)), LessEqual(LeafCount(d + e*x + f*x**S(2)), LeafCount(a + b*x + c*x**S(2))))

cons563 = CustomConstraint(cons_f563)

def cons_f564(f, c):
    return Not(PositiveQ(c/f))

cons564 = CustomConstraint(cons_f564)

def cons_f565(b, d, f, e, a, q, c):
    return ZeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

cons565 = CustomConstraint(cons_f565)

def cons_f566(q):
    return NonzeroQ(q + S(1))

cons566 = CustomConstraint(cons_f566)

def cons_f567(q):
    return NonzeroQ(S(2)*q + S(3))

cons567 = CustomConstraint(cons_f567)

def cons_f568(d, f, e, a, q, c):
    return ZeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

cons568 = CustomConstraint(cons_f568)

def cons_f569(f, a, q, d, c):
    return ZeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

cons569 = CustomConstraint(cons_f569)

def cons_f570(q):
    return PositiveIntegerQ(q + S(2))

cons570 = CustomConstraint(cons_f570)

def cons_f571(f, e, d):
    return NonzeroQ(-S(4)*d*f + e**S(2))

cons571 = CustomConstraint(cons_f571)

def cons_f572(b, d, f, e, a, q, c):
    return NonzeroQ(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))

cons572 = CustomConstraint(cons_f572)

def cons_f573(d, f, e, a, q, c):
    return NonzeroQ(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))

cons573 = CustomConstraint(cons_f573)

def cons_f574(f, a, q, d, c):
    return NonzeroQ(S(2)*a*f*q + S(3)*a*f - c*d)

cons574 = CustomConstraint(cons_f574)

def cons_f575(q):
    return Not(PositiveIntegerQ(q))

cons575 = CustomConstraint(cons_f575)

def cons_f576(q):
    return Not(And(RationalQ(q), LessEqual(q, S(-1))))

cons576 = CustomConstraint(cons_f576)

def cons_f577(d, b, f, e, a, c):
    return NonzeroQ(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))

cons577 = CustomConstraint(cons_f577)

def cons_f578(b, f, a, d, c):
    return NonzeroQ(b**S(2)*d*f + (-a*f + c*d)**S(2))

cons578 = CustomConstraint(cons_f578)

def cons_f579(f, e, a, d, c):
    return NonzeroQ(a*c*e**S(2) + (-a*f + c*d)**S(2))

cons579 = CustomConstraint(cons_f579)

def cons_f580(p, q):
    return NonzeroQ(p + q)

cons580 = CustomConstraint(cons_f580)

def cons_f581(p, q):
    return NonzeroQ(S(2)*p + S(2)*q + S(1))

cons581 = CustomConstraint(cons_f581)

def cons_f582(f, e, c, b):
    return ZeroQ(-b*f + c*e)

cons582 = CustomConstraint(cons_f582)

def cons_f583(f, e, c, b):
    return NonzeroQ(-b*f + c*e)

cons583 = CustomConstraint(cons_f583)

def cons_f584(a, c):
    return PosQ(-a*c)

cons584 = CustomConstraint(cons_f584)

def cons_f585(a, b, c):
    return NegQ(-S(4)*a*c + b**S(2))

cons585 = CustomConstraint(cons_f585)

def cons_f586(a, c):
    return NegQ(-a*c)

cons586 = CustomConstraint(cons_f586)

def cons_f587(b, f, p, e, a, q, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, p, q), x)

cons587 = CustomConstraint(cons_f587)

def cons_f588(f, p, e, a, q, x, d, c):
    return FreeQ(List(a, c, d, e, f, p, q), x)

cons588 = CustomConstraint(cons_f588)

def cons_f589(g, b, h, a, c):
    return ZeroQ(a*h**S(2) - b*g*h + c*g**S(2))

cons589 = CustomConstraint(cons_f589)

def cons_f590(g, d, h, f, e, a, c):
    return ZeroQ(a**S(2)*f*h**S(2) - a*c*e*g*h + c**S(2)*d*g**S(2))

cons590 = CustomConstraint(cons_f590)

def cons_f591(g, c, h, a):
    return ZeroQ(a*h**S(2) + c*g**S(2))

cons591 = CustomConstraint(cons_f591)

def cons_f592(g, d, h, f, a, c):
    return ZeroQ(a**S(2)*f*h**S(2) + c**S(2)*d*g**S(2))

cons592 = CustomConstraint(cons_f592)

def cons_f593(b, f, e, a, c):
    return ZeroQ(a*f**S(2) - b*e*f + c*e**S(2))

cons593 = CustomConstraint(cons_f593)

def cons_f594(f, e, c, a):
    return ZeroQ(a*f**S(2) + c*e**S(2))

cons594 = CustomConstraint(cons_f594)

def cons_f595(g, b, h, f, p, e, m, c):
    return ZeroQ(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

cons595 = CustomConstraint(cons_f595)

def cons_f596(g, b, d, h, f, p, a, m, c):
    return ZeroQ(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

cons596 = CustomConstraint(cons_f596)

def cons_f597(g, h, f, p, e, m, c):
    return ZeroQ(c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))))

cons597 = CustomConstraint(cons_f597)

def cons_f598(d, h, f, p, a, m, c):
    return ZeroQ(h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))))

cons598 = CustomConstraint(cons_f598)

def cons_f599(g, b, h, f, p, m, c):
    return ZeroQ(b*f*h*(m + p + S(2)) + S(2)*c*f*g*(p + S(1)))

cons599 = CustomConstraint(cons_f599)

def cons_f600(p, m):
    return Or(IntegersQ(m, p), PositiveIntegerQ(p))

cons600 = CustomConstraint(cons_f600)

def cons_f601(g, b, h, a, c):
    return NonzeroQ(a*h**S(2) - b*g*h + c*g**S(2))

cons601 = CustomConstraint(cons_f601)

def cons_f602(g, c, h, a):
    return NonzeroQ(a*h**S(2) + c*g**S(2))

cons602 = CustomConstraint(cons_f602)

def cons_f603(g, b, h, a, c):
    return NonzeroQ(c*g**S(2) - h*(-a*h + b*g))

cons603 = CustomConstraint(cons_f603)

def cons_f604(p, q):
    return Or(Greater(p, S(0)), Greater(q, S(0)))

cons604 = CustomConstraint(cons_f604)

def cons_f605(p, q):
    return NonzeroQ(p + q + S(1))

cons605 = CustomConstraint(cons_f605)

def cons_f606(a, c):
    return PositiveQ(a*c)

cons606 = CustomConstraint(cons_f606)

def cons_f607(a, c):
    return Not(PositiveQ(a*c))

cons607 = CustomConstraint(cons_f607)

def cons_f608(f, e, g, h):
    return ZeroQ(e*h - S(2)*f*g)

cons608 = CustomConstraint(cons_f608)

def cons_f609(f, e, g, h):
    return NonzeroQ(e*h - S(2)*f*g)

cons609 = CustomConstraint(cons_f609)

def cons_f610(g, e, d, h):
    return ZeroQ(S(2)*d*h - e*g)

cons610 = CustomConstraint(cons_f610)

def cons_f611(g, e, d, h):
    return NonzeroQ(S(2)*d*h - e*g)

cons611 = CustomConstraint(cons_f611)

def cons_f612(g, b, h, f, e, a, d, c):
    return ZeroQ(g**S(2)*(-b*f + c*e) - S(2)*g*h*(-a*f + c*d) + h**S(2)*(-a*e + b*d))

cons612 = CustomConstraint(cons_f612)

def cons_f613(g, d, h, f, e, a, c):
    return ZeroQ(a*e*h**S(2) - c*e*g**S(2) + S(2)*g*h*(-a*f + c*d))

cons613 = CustomConstraint(cons_f613)

def cons_f614(g, b, h, f, a, d, c):
    return ZeroQ(b*d*h**S(2) - b*f*g**S(2) - S(2)*g*h*(-a*f + c*d))

cons614 = CustomConstraint(cons_f614)

def cons_f615(d, b, f, a, c):
    return ZeroQ(c**S(2)*d - f*(-S(3)*a*c + b**S(2)))

cons615 = CustomConstraint(cons_f615)

def cons_f616(g, b, h, a, c):
    return ZeroQ(S(9)*a*c*h**S(2) - S(2)*b**S(2)*h**S(2) - b*c*g*h + c**S(2)*g**S(2))

cons616 = CustomConstraint(cons_f616)

def cons_f617(g, c, b, h):
    return PositiveQ(-S(9)*c*h**S(2)/(-b*h + S(2)*c*g)**S(2))

cons617 = CustomConstraint(cons_f617)

def cons_f618(f, c, d, a):
    return ZeroQ(S(3)*a*f + c*d)

cons618 = CustomConstraint(cons_f618)

def cons_f619(g, c, h, a):
    return ZeroQ(S(9)*a*h**S(2) + c*g**S(2))

cons619 = CustomConstraint(cons_f619)

def cons_f620(g, b, h, f, p, e, a, q, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, g, h, p, q), x)

cons620 = CustomConstraint(cons_f620)

def cons_f621(g, h, f, p, e, a, q, x, d, c):
    return FreeQ(List(a, c, d, e, f, g, h, p, q), x)

cons621 = CustomConstraint(cons_f621)

def cons_f622(x, z):
    return LinearQ(z, x)

cons622 = CustomConstraint(cons_f622)

def cons_f623(u, x, v):
    return QuadraticQ(List(u, v), x)

cons623 = CustomConstraint(cons_f623)

def cons_f624(u, x, v, z):
    return Not(And(LinearMatchQ(z, x), QuadraticMatchQ(List(u, v), x)))

cons624 = CustomConstraint(cons_f624)

def cons_f625(A, x):
    return FreeQ(A, x)

cons625 = CustomConstraint(cons_f625)

def cons_f626(B, x):
    return FreeQ(B, x)

cons626 = CustomConstraint(cons_f626)

def cons_f627(C, x):
    return FreeQ(C, x)

cons627 = CustomConstraint(cons_f627)

def cons_f628(p, q):
    return NonzeroQ(S(2)*p + S(2)*q + S(3))

cons628 = CustomConstraint(cons_f628)

def cons_f629(b, f, A, C, p, B, e, a, q, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, A, B, C, p, q), x)

cons629 = CustomConstraint(cons_f629)

def cons_f630(b, f, A, C, p, e, a, q, x, d, c):
    return FreeQ(List(a, b, c, d, e, f, A, C, p, q), x)

cons630 = CustomConstraint(cons_f630)

def cons_f631(f, A, C, p, B, e, a, q, x, d, c):
    return FreeQ(List(a, c, d, e, f, A, B, C, p, q), x)

cons631 = CustomConstraint(cons_f631)

def cons_f632(f, A, C, p, e, a, q, x, d, c):
    return FreeQ(List(a, c, d, e, f, A, C, p, q), x)

cons632 = CustomConstraint(cons_f632)
