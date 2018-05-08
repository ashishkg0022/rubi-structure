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
    from sympy.integrals.rubi.symbol import WC
    from sympy.core.symbol import symbols, Symbol
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf, exp, log)
    from sympy.functions.elementary.hyperbolic import (acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh, tanh, coth, sech, csch)
    from sympy.functions.elementary.trigonometric import (atan, acsc, asin, acot, acos, asec)
    from sympy import pi as Pi


    A_, B_, C_, F_, G_, H_, a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, p_, q_, r_, t_, u_, v_, s_, w_, x_, y_, z_ = [WC(i) for i in 'ABCFGHabcdefghijklmnpqrtuvswxyz']
    a1_, a2_, b1_, b2_, c1_, c2_, d1_, d2_, n1_, n2_, e1_, e2_, f1_, f2_, g1_, g2_, n1_, n2_, n3_, Pq_, Pm_, Px_, Qm_, Qr_, Qx_, jn_, mn_, non2_, RFx_, RGx_ = [WC(i) for i in ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'd1', 'd2', 'n1', 'n2', 'e1', 'e2', 'f1', 'f2', 'g1', 'g2', 'n1', 'n2', 'n3', 'Pq', 'Pm', 'Px', 'Qm', 'Qr', 'Qx', 'jn', 'mn', 'non2', 'RFx', 'RGx']]

    _UseGamma = False

def quadratic_products(rubi):
    from sympy.integrals.rubi.constraints import cons405, cons4, cons5, cons9, cons406, cons74, cons407, cons77, cons408, cons409, cons87, cons116, cons410, cons88, cons411, cons412, cons413, cons414, cons415, cons416, cons6, cons7, cons417, cons418, cons10, cons72, cons2, cons419, cons420, cons421, cons422, cons1, cons423, cons424, cons425, cons236, cons99, cons426, cons427, cons428, cons429, cons430, cons431, cons432, cons433, cons119, cons434, cons51, cons121, cons435, cons436, cons38, cons100, cons437, cons97, cons438, cons439, cons440, cons71, cons441, cons442, cons443, cons444, cons445, cons446, cons447, cons448, cons449, cons17, cons450, cons85, cons451, cons452, cons453, cons454, cons455, cons456, cons457, cons458, cons459, cons460, cons461, cons462, cons463, cons464, cons465, cons466, cons467, cons468, cons60, cons469, cons470, cons471, cons472, cons473, cons474, cons475, cons476, cons477, cons478, cons479, cons480, cons481, cons482, cons483, cons484, cons485, cons486, cons487, cons488, cons489, cons490, cons491, cons492, cons493, cons494, cons495, cons27, cons28, cons496, cons497, cons498, cons73, cons161, cons499, cons500, cons501, cons25, cons502, cons503, cons504, cons505, cons506, cons13, cons507, cons508, cons509, cons91, cons510, cons511, cons512, cons513, cons103, cons514, cons101, cons515, cons149, cons516, cons517, cons518, cons519, cons520, cons32, cons521, cons522, cons523, cons31, cons30, cons524, cons525, cons526, cons75, cons527, cons528, cons160, cons529, cons530, cons531, cons532, cons533, cons534, cons535, cons536, cons537, cons271, cons212, cons538, cons539, cons540, cons541, cons542, cons543, cons544, cons545, cons546, cons266, cons547, cons548, cons549, cons550, cons551, cons552, cons553, cons554, cons102, cons555, cons70, cons556, cons37, cons36, cons118, cons14, cons557, cons22, cons558, cons402, cons559, cons560, cons561, cons562, cons563, cons173, cons324, cons564, cons565, cons566, cons567, cons568, cons569, cons570, cons571, cons289, cons380, cons572, cons573, cons574, cons575, cons576, cons308, cons290, cons577, cons312, cons578, cons579, cons580, cons581, cons582, cons583, cons584, cons585, cons586, cons587, cons588, cons162, cons589, cons590, cons591, cons592, cons593, cons594, cons595, cons596, cons597, cons598, cons599, cons600, cons601, cons602, cons603, cons174, cons604, cons605, cons606, cons607, cons608, cons609, cons610, cons611, cons612, cons613, cons614, cons615, cons616, cons617, cons618, cons619, cons188, cons620, cons621, cons622, cons623, cons624, cons178, cons625, cons626, cons627, cons628, cons629, cons630, cons631, cons632

    pattern1 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons405)
    rule1 = ReplacementRule(pattern1, lambda x, c, b, a : (b/S(2) + c*x)*Int(S(1)/(b/S(2) + c*x), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule1)

    pattern2 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons74, cons405, cons406)
    rule2 = ReplacementRule(pattern2, lambda b, p, a, x, c : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule2)

    pattern3 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons407, cons77, cons408, )
    def With3(b, p, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c**(-p)*Int(Simp(b/S(2) + c*x - q/S(2), x)**p*Simp(b/S(2) + c*x + q/S(2), x)**p, x)
    rule3 = ReplacementRule(pattern3, lambda b, p, a, x, c : With3(b, p, a, x, c))
    rubi.add(rule3)

    pattern4 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons407, cons77, cons409)
    rule4 = ReplacementRule(pattern4, lambda b, p, a, x, c : Int(ExpandIntegrand((a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule4)

    pattern5 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons407, cons87, cons116, cons410)
    rule5 = ReplacementRule(pattern5, lambda b, p, a, x, c : -p*(-S(4)*a*c + b**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*(S(2)*p + S(1))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule5)

    pattern6 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(-3)/2), x_), cons4, cons5, cons9, cons407)
    rule6 = ReplacementRule(pattern6, lambda x, a, b, c : -S(2)*(b + S(2)*c*x)/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule6)

    pattern7 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons407, cons87, cons88, cons411, cons410)
    rule7 = ReplacementRule(pattern7, lambda b, p, a, x, c : -S(2)*c*(S(2)*p + S(3))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule7)

    pattern8 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons407, cons412, cons408, )
    def With8(x, c, b, a):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return c*Int(S(1)/Simp(b/S(2) + c*x - q/S(2), x), x)/q - c*Int(S(1)/Simp(b/S(2) + c*x + q/S(2), x), x)/q
    rule8 = ReplacementRule(pattern8, lambda x, c, b, a : With8(x, c, b, a))
    rubi.add(rule8)

    pattern9 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons407, )
    def With9(x, c, b, a):
        q = -S(4)*a*c/b**S(2) + S(1)
        if And(RationalQ(q), Or(EqQ(q**S(2), S(1)), Not(RationalQ(-S(4)*a*c + b**S(2))))):
            return -S(2)*Subst(Int(S(1)/(q - x**S(2)), x), x, S(1) + S(2)*c*x/b)/b
        print("Unable to Integrate")
    rule9 = ReplacementRule(pattern9, lambda x, c, b, a : With9(x, c, b, a))
    rubi.add(rule9)

    pattern10 = Pattern(Integral(S(1)/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons407)
    rule10 = ReplacementRule(pattern10, lambda x, c, b, a : -S(2)*Subst(Int(S(1)/Simp(-S(4)*a*c + b**S(2) - x**S(2), x), x), x, b + S(2)*c*x))
    rubi.add(rule10)

    pattern11 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons74, cons413)
    rule11 = ReplacementRule(pattern11, lambda b, p, a, x, c : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p, x), x, b + S(2)*c*x)/(S(2)*c))
    rubi.add(rule11)

    pattern12 = Pattern(Integral(S(1)/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons5, cons9, cons414)
    rule12 = ReplacementRule(pattern12, lambda x, c, b : S(2)*Subst(Int(S(1)/(-c*x**S(2) + S(1)), x), x, x/sqrt(b*x + c*x**S(2))))
    rubi.add(rule12)

    pattern13 = Pattern(Integral(S(1)/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons407)
    rule13 = ReplacementRule(pattern13, lambda x, c, b, a : S(2)*Subst(Int(S(1)/(S(4)*c - x**S(2)), x), x, (b + S(2)*c*x)/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule13)

    pattern14 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons5, cons9, cons87, cons415)
    rule14 = ReplacementRule(pattern14, lambda x, c, b, p : (-c*(b*x + c*x**S(2))/b**S(2))**(-p)*(b*x + c*x**S(2))**p*Int((-c*x/b - c**S(2)*x**S(2)/b**S(2))**p, x))
    rubi.add(rule14)

    pattern15 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons407, cons87, )
    def With15(b, p, a, x, c):
        d = Denominator(p)
        if LessEqual(S(3), d, S(4)):
            return d*sqrt((b + S(2)*c*x)**S(2))*Subst(Int(x**(d*(p + S(1)) + S(-1))/sqrt(-S(4)*a*c + b**S(2) + S(4)*c*x**d), x), x, (a + b*x + c*x**S(2))**(S(1)/d))/(b + S(2)*c*x)
        print("Unable to Integrate")
    rule15 = ReplacementRule(pattern15, lambda b, p, a, x, c : With15(b, p, a, x, c))
    rubi.add(rule15)

    pattern16 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons74, cons407, cons416, )
    def With16(b, p, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -((-b - S(2)*c*x + q)/(S(2)*q))**(-p + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Hypergeometric2F1(-p, p + S(1), p + S(2), (b + S(2)*c*x + q)/(S(2)*q))/(q*(p + S(1)))
    rule16 = ReplacementRule(pattern16, lambda b, p, a, x, c : With16(b, p, a, x, c))
    rubi.add(rule16)

    pattern17 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons74, cons6, cons7)
    rule17 = ReplacementRule(pattern17, lambda b, u, p, a, x, c : Subst(Int((a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule17)

    pattern18 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons417, cons418)
    rule18 = ReplacementRule(pattern18, lambda b, d, p, e, a, m, x, c : c**(-m/S(2) + S(-1)/2)*e**m*(a + b*x + c*x**S(2))**(m/S(2) + p + S(1)/2)/(m + S(2)*p + S(1)))
    rubi.add(rule18)

    pattern19 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons417, cons419)
    rule19 = ReplacementRule(pattern19, lambda b, d, p, e, a, m, x, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*log(RemoveContent(d + e*x, x))/e)
    rubi.add(rule19)

    pattern20 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons417, cons420)
    rule20 = ReplacementRule(pattern20, lambda b, d, p, e, a, m, x, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule20)

    pattern21 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons421, cons422, cons1)
    rule21 = ReplacementRule(pattern21, lambda b, p, e, a, m, x, d, c : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule21)

    pattern22 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(x_*WC('e', S(1)) + WC('d', S(0)))**S(2), x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421)
    rule22 = ReplacementRule(pattern22, lambda b, e, a, x, d, c : sqrt(a + b*x + c*x**S(2))*Int((b + S(2)*c*x)/(d + e*x)**S(2), x)/(b + S(2)*c*x))
    rubi.add(rule22)

    pattern23 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons405, cons421, cons423)
    rule23 = ReplacementRule(pattern23, lambda b, e, a, m, x, d, c : (d + e*x)**(m + S(1))*sqrt(a + b*x + c*x**S(2))/(e*(m + S(2))) - (-b*e + S(2)*c*d)*sqrt(a + b*x + c*x**S(2))*Int((d + e*x)**m, x)/(e*(b + S(2)*c*x)*(m + S(2))))
    rubi.add(rule23)

    pattern24 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421)
    rule24 = ReplacementRule(pattern24, lambda b, e, a, x, d, c : -S(4)*c*e*sqrt(a + b*x + c*x**S(2))/((d + e*x)*(-b*e + S(2)*c*d)**S(2)) + S(2)*c*Int(S(1)/((d + e*x)*sqrt(a + b*x + c*x**S(2))), x)/(-b*e + S(2)*c*d))
    rubi.add(rule24)

    pattern25 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons421, cons424, cons423)
    rule25 = ReplacementRule(pattern25, lambda b, p, e, a, m, x, d, c : -S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)**S(2)*(m*p + S(-1))) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(2))*(-b*e + S(2)*c*d)))
    rubi.add(rule25)

    pattern26 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons405, cons421, cons425)
    rule26 = ReplacementRule(pattern26, lambda b, p, e, a, x, d, c : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule26)

    pattern27 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421, cons236, cons99, cons426, cons427)
    rule27 = ReplacementRule(pattern27, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(1))))
    rubi.add(rule27)

    pattern28 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421, cons236, cons99, cons428, cons427, cons429)
    rule28 = ReplacementRule(pattern28, lambda b, p, e, a, m, x, d, c : S(2)*c*p*(S(2)*p + S(-1))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(1))*(m + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))/(e**S(2)*(m + S(1))*(m + S(2))))
    rubi.add(rule28)

    pattern29 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons405, cons421, cons87, cons116, cons430, cons420, cons429, cons431, cons432, cons427)
    rule29 = ReplacementRule(pattern29, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))) + p*(S(2)*p + S(-1))*(-b*e + S(2)*c*d)**S(2)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*c*e**S(2)*(m + S(2)*p)*(m + S(2)*p + S(1))))
    rubi.add(rule29)

    pattern30 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421, cons236, cons88, cons433, cons427)
    rule30 = ReplacementRule(pattern30, lambda b, p, e, a, m, x, d, c : e**S(2)*m*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule30)

    pattern31 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons405, cons421, cons236, cons88, cons119, cons427)
    rule31 = ReplacementRule(pattern31, lambda b, p, e, a, m, x, d, c : e**S(2)*m*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) - e*m*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))*(S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(S(2)*p + S(1))))
    rubi.add(rule31)

    pattern32 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons405, cons421, cons236, cons88, cons434, cons427)
    rule32 = ReplacementRule(pattern32, lambda b, p, e, a, m, x, d, c : S(2)*c*e**S(2)*(m + S(2)*p + S(2))*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) - S(2)*c*e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(m + S(2)*p + S(2))/((p + S(1))*(S(2)*p + S(1))*(-b*e + S(2)*c*d)**S(2)) + (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((S(2)*p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule32)

    pattern33 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons405, cons421, cons51, cons121, cons420, cons435, cons436)
    rule33 = ReplacementRule(pattern33, lambda b, p, e, a, m, x, d, c : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*(m + S(2)*p + S(1))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p/(S(2)*c*(m + S(2)*p + S(1))))
    rubi.add(rule33)

    pattern34 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons405, cons421, cons51, cons38, cons427)
    rule34 = ReplacementRule(pattern34, lambda b, p, e, a, m, x, d, c : S(2)*c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-b*e + S(2)*c*d)) - (b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/((m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule34)

    pattern35 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons405, cons100, cons421)
    rule35 = ReplacementRule(pattern35, lambda b, p, e, a, m, x, d, c : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m, x))
    rubi.add(rule35)

    pattern36 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons437, cons97)
    rule36 = ReplacementRule(pattern36, lambda b, d, p, e, a, m, x, c : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule36)

    pattern37 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons439)
    rule37 = ReplacementRule(pattern37, lambda d, p, e, a, m, x, c : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule37)

    pattern38 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons437, cons100, cons440)
    rule38 = ReplacementRule(pattern38, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule38)

    pattern39 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons440)
    rule39 = ReplacementRule(pattern39, lambda d, p, e, a, m, x, c : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule39)

    pattern40 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons437, cons100, cons422)
    rule40 = ReplacementRule(pattern40, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule40)

    pattern41 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons422)
    rule41 = ReplacementRule(pattern41, lambda d, p, e, a, m, x, c : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(p + S(1))))
    rubi.add(rule41)

    pattern42 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons437, cons100, cons87, cons88)
    rule42 = ReplacementRule(pattern42, lambda b, p, e, a, x, d, c : -e**S(2)*(p + S(2))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule42)

    pattern43 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**S(2), x_), cons4, cons9, cons10, cons72, cons74, cons438, cons100, cons87, cons88)
    rule43 = ReplacementRule(pattern43, lambda d, p, e, a, x, c : -e**S(2)*(p + S(2))*Int((a + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)/(c*(p + S(1))))
    rubi.add(rule43)

    pattern44 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437, cons100, cons71, cons87, cons441, cons442, cons443)
    rule44 = ReplacementRule(pattern44, lambda b, p, e, a, m, x, d, c : Int((a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule44)

    pattern45 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons71, cons87, cons441, cons442, cons443)
    rule45 = ReplacementRule(pattern45, lambda d, p, e, a, m, x, c : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m), x))
    rubi.add(rule45)

    pattern46 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons437, cons100, cons444)
    rule46 = ReplacementRule(pattern46, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule46)

    pattern47 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons444)
    rule47 = ReplacementRule(pattern47, lambda d, p, e, a, m, x, c : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule47)

    pattern48 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons437, cons100, cons445)
    rule48 = ReplacementRule(pattern48, lambda b, p, e, a, m, x, d, c : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule48)

    pattern49 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons445)
    rule49 = ReplacementRule(pattern49, lambda d, p, e, a, m, x, c : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule49)

    pattern50 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437)
    rule50 = ReplacementRule(pattern50, lambda b, e, a, x, d, c : S(2)*e*Subst(Int(S(1)/(-b*e + S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule50)

    pattern51 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons438)
    rule51 = ReplacementRule(pattern51, lambda d, e, a, x, c : S(2)*e*Subst(Int(S(1)/(S(2)*c*d + e**S(2)*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule51)

    pattern52 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437, cons236, cons116, cons446, cons434, cons427)
    rule52 = ReplacementRule(pattern52, lambda b, p, e, a, m, x, d, c : -c*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + p + S(1))))
    rubi.add(rule52)

    pattern53 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons438, cons236, cons116, cons446, cons434, cons427)
    rule53 = ReplacementRule(pattern53, lambda d, p, e, a, m, x, c : -c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/(e**S(2)*(m + p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + p + S(1))))
    rubi.add(rule53)

    pattern54 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437, cons236, cons116, cons447, cons420, cons427)
    rule54 = ReplacementRule(pattern54, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - p*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule54)

    pattern55 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons438, cons236, cons116, cons447, cons420, cons427)
    rule55 = ReplacementRule(pattern55, lambda d, p, e, a, m, x, c : -S(2)*c*d*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e**S(2)*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule55)

    pattern56 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437, cons236, cons88, cons433, cons427)
    rule56 = ReplacementRule(pattern56, lambda b, p, e, a, m, x, d, c : -(-b*e + S(2)*c*d)*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(-b*e + S(2)*c*d)*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule56)

    pattern57 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons438, cons236, cons88, cons433, cons427)
    rule57 = ReplacementRule(pattern57, lambda d, p, e, a, m, x, c : d*(m + S(2)*p + S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - d*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*e*(p + S(1))))
    rubi.add(rule57)

    pattern58 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons437, cons236, cons88, cons119, cons427)
    rule58 = ReplacementRule(pattern58, lambda b, p, e, a, m, x, d, c : -e**S(2)*(m + p)*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule58)

    pattern59 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons438, cons236, cons88, cons119, cons427)
    rule59 = ReplacementRule(pattern59, lambda d, p, e, a, m, x, c : -e**S(2)*(m + p)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(p + S(1))))
    rubi.add(rule59)

    pattern60 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons437, cons51, cons448, cons420, cons427)
    rule60 = ReplacementRule(pattern60, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + (m + p)*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule60)

    pattern61 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons74, cons438, cons51, cons448, cons420, cons427)
    rule61 = ReplacementRule(pattern61, lambda d, p, e, a, m, x, c : S(2)*d*(m + p)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1)), x)/(m + S(2)*p + S(1)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))))
    rubi.add(rule61)

    pattern62 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons437, cons51, cons449, cons434, cons427)
    rule62 = ReplacementRule(pattern62, lambda b, p, e, a, m, x, d, c : c*(m + S(2)*p + S(2))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/((-b*e + S(2)*c*d)*(m + p + S(1))) - e*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule62)

    pattern63 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons74, cons438, cons51, cons449, cons434, cons427)
    rule63 = ReplacementRule(pattern63, lambda d, p, e, a, m, x, c : (m + S(2)*p + S(2))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*d*(m + p + S(1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*c*d*(m + p + S(1))))
    rubi.add(rule63)

    pattern64 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons5, cons9, cons72, cons2, cons100)
    rule64 = ReplacementRule(pattern64, lambda b, p, e, m, x, c : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p, x))
    rubi.add(rule64)

    pattern65 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons438, cons100, cons17, cons450)
    rule65 = ReplacementRule(pattern65, lambda d, p, e, a, m, x, c : Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule65)

    pattern66 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons437, cons100)
    rule66 = ReplacementRule(pattern66, lambda b, d, p, e, a, m, x, c : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule66)

    pattern67 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons438, cons100)
    rule67 = ReplacementRule(pattern67, lambda d, p, e, a, m, x, c : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(a/d + c*x/e)**p, x))
    rubi.add(rule67)

    pattern68 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417)
    rule68 = ReplacementRule(pattern68, lambda b, e, a, x, d, c : b**S(2)*Int((d + e*x)/(a + b*x + c*x**S(2)), x)/(d**S(2)*(-S(4)*a*c + b**S(2))) - S(4)*b*c*Int(S(1)/(b + S(2)*c*x), x)/(d*(-S(4)*a*c + b**S(2))))
    rubi.add(rule68)

    pattern69 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons417, cons424, cons85)
    rule69 = ReplacementRule(pattern69, lambda b, d, p, e, a, m, x, c : S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule69)

    pattern70 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons417, cons77, cons451)
    rule70 = ReplacementRule(pattern70, lambda b, p, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule70)

    pattern71 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417, cons452, cons236, cons116, cons38, cons453, cons427)
    rule71 = ReplacementRule(pattern71, lambda b, d, p, e, a, m, x, c : -b*p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(d*e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule71)

    pattern72 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons417, cons452, cons87, cons116, cons454, cons455, cons51, cons427)
    rule72 = ReplacementRule(pattern72, lambda b, d, p, e, a, m, x, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))) - d*p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(b*e*(m + S(2)*p + S(1))))
    rubi.add(rule72)

    pattern73 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417, cons452, cons236, cons88, cons119, cons427)
    rule73 = ReplacementRule(pattern73, lambda b, p, e, a, m, x, d, c : -d*e*(m + S(-1))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))) + d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))))
    rubi.add(rule73)

    pattern74 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons417, cons452, cons87, cons88, cons456, cons51, cons427)
    rule74 = ReplacementRule(pattern74, lambda b, d, p, e, a, m, x, c : -S(2)*c*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + S(2)*c*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(e*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule74)

    pattern75 = Pattern(Integral(S(1)/((d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417)
    rule75 = ReplacementRule(pattern75, lambda b, e, a, x, d, c : S(4)*c*Subst(Int(S(1)/(-S(4)*a*c*e + b**S(2)*e + S(4)*c*e*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule75)

    pattern76 = Pattern(Integral(S(1)/(sqrt(d_ + x_*WC('e', S(1)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417, cons457)
    rule76 = ReplacementRule(pattern76, lambda b, e, a, x, d, c : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(S(1)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule76)

    pattern77 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417, cons457)
    rule77 = ReplacementRule(pattern77, lambda b, e, a, x, d, c : S(4)*sqrt(-c/(-S(4)*a*c + b**S(2)))*Subst(Int(x**S(2)/sqrt(Simp(-b**S(2)*x**S(4)/(d**S(2)*(-S(4)*a*c + b**S(2))) + S(1), x)), x), x, sqrt(d + e*x))/e)
    rubi.add(rule77)

    pattern78 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons417, cons458)
    rule78 = ReplacementRule(pattern78, lambda b, e, a, m, x, d, c : sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*Int((d + e*x)**m/sqrt(-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2))), x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule78)

    pattern79 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons417, cons452, cons51, cons119, cons420, cons459)
    rule79 = ReplacementRule(pattern79, lambda b, d, p, e, a, m, x, c : S(2)*d*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(2)*p + S(1))) + d**S(2)*(m + S(-1))*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p, x)/(b**S(2)*(m + S(2)*p + S(1))))
    rubi.add(rule79)

    pattern80 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons417, cons452, cons51, cons38, cons460)
    rule80 = ReplacementRule(pattern80, lambda b, d, p, e, a, m, x, c : b**S(2)*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/(d**S(2)*(m + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*b*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(d*(m + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule80)

    pattern81 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons417)
    rule81 = ReplacementRule(pattern81, lambda b, d, p, e, a, m, x, c : Subst(Int(x**m*(a - b**S(2)/(S(4)*c) + c*x**S(2)/e**S(2))**p, x), x, d + e*x)/e)
    rubi.add(rule81)

    pattern82 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons77)
    rule82 = ReplacementRule(pattern82, lambda b, p, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule82)

    pattern83 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons462, cons77, cons463)
    rule83 = ReplacementRule(pattern83, lambda d, p, e, a, m, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m, x), x))
    rubi.add(rule83)

    pattern84 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons464, )
    def With84(b, e, a, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (c*d - e*(b/S(2) - q/S(2)))*Int(S(1)/(b/S(2) + c*x - q/S(2)), x)/q - (c*d - e*(b/S(2) + q/S(2)))*Int(S(1)/(b/S(2) + c*x + q/S(2)), x)/q
    rule84 = ReplacementRule(pattern84, lambda b, e, a, x, d, c : With84(b, e, a, x, d, c))
    rubi.add(rule84)

    pattern85 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons465, )
    def With85(d, e, a, x, c):
        q = Rt(-a*c, S(2))
        return (-c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x + q), x) + (c*d/(S(2)*q) + e/S(2))*Int(S(1)/(c*x - q), x)
    rule85 = ReplacementRule(pattern85, lambda d, e, a, x, c : With85(d, e, a, x, c))
    rubi.add(rule85)

    pattern86 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons466)
    rule86 = ReplacementRule(pattern86, lambda b, e, a, x, d, c : e*Int((b + S(2)*c*x)/(a + b*x + c*x**S(2)), x)/(S(2)*c) + (-b*e + S(2)*c*d)*Int(S(1)/(a + b*x + c*x**S(2)), x)/(S(2)*c))
    rubi.add(rule86)

    pattern87 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons467)
    rule87 = ReplacementRule(pattern87, lambda d, e, a, x, c : d*Int(S(1)/(a + c*x**S(2)), x) + e*Int(x/(a + c*x**S(2)), x))
    rubi.add(rule87)

    pattern88 = Pattern(Integral(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421)
    rule88 = ReplacementRule(pattern88, lambda b, e, a, x, d, c : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule88)

    pattern89 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462)
    rule89 = ReplacementRule(pattern89, lambda d, e, a, x, c : S(2)*e*Subst(Int(x**S(2)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule89)

    pattern90 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons71, cons119, cons468)
    rule90 = ReplacementRule(pattern90, lambda b, e, a, m, x, d, c : Int(PolynomialDivide((d + e*x)**m, a + b*x + c*x**S(2), x), x))
    rubi.add(rule90)

    pattern91 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons71, cons119, cons468)
    rule91 = ReplacementRule(pattern91, lambda d, e, a, m, x, c : Int(PolynomialDivide((d + e*x)**m, a + c*x**S(2), x), x))
    rubi.add(rule91)

    pattern92 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons51, cons119)
    rule92 = ReplacementRule(pattern92, lambda b, e, a, m, x, d, c : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + e*x*(-b*e + S(2)*c*d), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule92)

    pattern93 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons51, cons119)
    rule93 = ReplacementRule(pattern93, lambda d, e, a, m, x, c : e*(d + e*x)**(m + S(-1))/(c*(m + S(-1))) + Int((d + e*x)**(m + S(-2))*Simp(-a*e**S(2) + c*d**S(2) + S(2)*c*d*e*x, x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule93)

    pattern94 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421)
    rule94 = ReplacementRule(pattern94, lambda b, e, a, x, d, c : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((-b*e + c*d - c*e*x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule94)

    pattern95 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons462)
    rule95 = ReplacementRule(pattern95, lambda d, e, a, x, c : e**S(2)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((c*d - c*e*x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule95)

    pattern96 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421)
    rule96 = ReplacementRule(pattern96, lambda b, e, a, x, d, c : S(2)*e*Subst(Int(S(1)/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule96)

    pattern97 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons462)
    rule97 = ReplacementRule(pattern97, lambda d, e, a, x, c : S(2)*e*Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule97)

    pattern98 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons51, cons38)
    rule98 = ReplacementRule(pattern98, lambda b, e, a, m, x, d, c : e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(-b*e + c*d - c*e*x, x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule98)

    pattern99 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons2, cons462, cons51, cons38)
    rule99 = ReplacementRule(pattern99, lambda d, e, a, m, x, c : c*Int((d - e*x)*(d + e*x)**(m + S(1))/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)) + e*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule99)

    pattern100 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons60)
    rule100 = ReplacementRule(pattern100, lambda b, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule100)

    pattern101 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons2, cons462, cons60)
    rule101 = ReplacementRule(pattern101, lambda d, e, a, m, x, c : Int(ExpandIntegrand((d + e*x)**m, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule101)

    pattern102 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421)
    rule102 = ReplacementRule(pattern102, lambda b, e, a, x, d, c : -S(2)*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((-S(4)*a*c + b**S(2))*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule102)

    pattern103 = Pattern(Integral((d_ + x_*WC('e', S(1)))/(a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2), x_), cons4, cons9, cons10, cons72, cons462)
    rule103 = ReplacementRule(pattern103, lambda d, e, a, x, c : (-a*e + c*d*x)/(a*c*sqrt(a + c*x**S(2))))
    rubi.add(rule103)

    pattern104 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons87, cons88, cons411)
    rule104 = ReplacementRule(pattern104, lambda b, p, e, a, x, d, c : -(S(2)*p + S(3))*(-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule104)

    pattern105 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons87, cons88, cons411)
    rule105 = ReplacementRule(pattern105, lambda d, p, e, a, x, c : d*(S(2)*p + S(3))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule105)

    pattern106 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons461, cons421, cons469)
    rule106 = ReplacementRule(pattern106, lambda b, p, e, a, x, d, c : e*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))) + (-b*e + S(2)*c*d)*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c))
    rubi.add(rule106)

    pattern107 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons74, cons462, cons469)
    rule107 = ReplacementRule(pattern107, lambda d, p, e, a, x, c : d*Int((a + c*x**S(2))**p, x) + e*(a + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule107)

    pattern108 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons470, cons471, cons472, cons100)
    rule108 = ReplacementRule(pattern108, lambda b, p, e, a, m, x, d, c : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m - p)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule108)

    pattern109 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons5, cons9, cons10, cons72, cons473, cons421, cons51, cons474, cons475, cons476)
    rule109 = ReplacementRule(pattern109, lambda b, e, m, x, d, c : Int((d + e*x)**m/(sqrt(b*x)*sqrt(S(1) + c*x/b)), x))
    rubi.add(rule109)

    pattern110 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons5, cons9, cons10, cons72, cons473, cons421, cons51, cons474)
    rule110 = ReplacementRule(pattern110, lambda b, e, m, x, d, c : sqrt(x)*sqrt(b + c*x)*Int((d + e*x)**m/(sqrt(x)*sqrt(b + c*x)), x)/sqrt(b*x + c*x**S(2)))
    rubi.add(rule110)

    pattern111 = Pattern(Integral(x_**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons407, cons477)
    rule111 = ReplacementRule(pattern111, lambda b, a, m, x, c : S(2)*Subst(Int(x**(S(2)*m + S(1))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule111)

    pattern112 = Pattern(Integral((e_*x_)**m_/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons72, cons407, cons477)
    rule112 = ReplacementRule(pattern112, lambda b, e, a, m, x, c : x**(-m)*(e*x)**m*Int(x**m/sqrt(a + b*x + c*x**S(2)), x))
    rubi.add(rule112)

    pattern113 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons477)
    rule113 = ReplacementRule(pattern113, lambda b, e, a, m, x, d, c : S(2)*sqrt(-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))*(S(2)*c*(d + e*x)/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))))**(-m)*(d + e*x)**m*Rt(-S(4)*a*c + b**S(2), S(2))*Subst(Int((S(2)*e*x**S(2)*Rt(-S(4)*a*c + b**S(2), S(2))/(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2))) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(S(2))*sqrt((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))/Rt(-S(4)*a*c + b**S(2), S(2)))/S(2))/(c*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule113)

    pattern114 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_/sqrt(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons462, cons477)
    rule114 = ReplacementRule(pattern114, lambda d, e, a, m, x, c : S(2)*a*(c*(d + e*x)/(-a*e*Rt(-c/a, S(2)) + c*d))**(-m)*sqrt(S(1) + c*x**S(2)/a)*(d + e*x)**m*Rt(-c/a, S(2))*Subst(Int((S(2)*a*e*x**S(2)*Rt(-c/a, S(2))/(-a*e*Rt(-c/a, S(2)) + c*d) + S(1))**m/sqrt(-x**S(2) + S(1)), x), x, sqrt(-x*Rt(-c/a, S(2))/S(2) + S(1)/2))/(c*sqrt(a + c*x**S(2))))
    rubi.add(rule114)

    pattern115 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons236, cons478, cons116)
    rule115 = ReplacementRule(pattern115, lambda b, p, e, a, m, x, d, c : p*(-S(4)*a*c + b**S(2))*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/(S(2)*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule115)

    pattern116 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons462, cons236, cons478, cons116)
    rule116 = ReplacementRule(pattern116, lambda d, p, e, a, m, x, c : -S(2)*a*c*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2)), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*a*e + S(2)*c*d*x)/(S(2)*(m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule116)

    pattern117 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons236, cons478, cons88)
    rule117 = ReplacementRule(pattern117, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) - S(2)*(S(2)*p + S(3))*(a*e**S(2) - b*d*e + c*d**S(2))*Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule117)

    pattern118 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons462, cons236, cons478, cons88)
    rule118 = ReplacementRule(pattern118, lambda d, p, e, a, m, x, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) + (S(2)*p + S(3))*(a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule118)

    pattern119 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons421)
    rule119 = ReplacementRule(pattern119, lambda b, e, a, x, d, c : -S(2)*Subst(Int(S(1)/(S(4)*a*e**S(2) - S(4)*b*d*e + S(4)*c*d**S(2) - x**S(2)), x), x, (S(2)*a*e - b*d - x*(-b*e + S(2)*c*d))/sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule119)

    pattern120 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons479)
    rule120 = ReplacementRule(pattern120, lambda d, e, a, x, c : -Subst(Int(S(1)/(a*e**S(2) + c*d**S(2) - x**S(2)), x), x, (a*e - c*d*x)/sqrt(a + c*x**S(2))))
    rubi.add(rule120)

    pattern121 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons100, cons422)
    rule121 = ReplacementRule(pattern121, lambda b, p, e, a, m, x, d, c : -((b + S(2)*c*x + Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))**(-p)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), -S(4)*c*(d + e*x)*Rt(-S(4)*a*c + b**S(2), S(2))/((b + S(2)*c*x - Rt(-S(4)*a*c + b**S(2), S(2)))*(-b*e + S(2)*c*d - e*Rt(-S(4)*a*c + b**S(2), S(2)))))/((m + S(1))*(-b*e + S(2)*c*d + e*Rt(-S(4)*a*c + b**S(2), S(2)))))
    rubi.add(rule121)

    pattern122 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons100, cons422)
    rule122 = ReplacementRule(pattern122, lambda d, p, e, a, m, x, c : ((c*d + e*Rt(-a*c, S(2)))*(c*x + Rt(-a*c, S(2)))/((c*d - e*Rt(-a*c, S(2)))*(c*x - Rt(-a*c, S(2)))))**(-p)*(a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*x + Rt(-a*c, S(2)))*Hypergeometric2F1(m + S(1), -p, m + S(2), S(2)*c*(d + e*x)*Rt(-a*c, S(2))/((c*d - e*Rt(-a*c, S(2)))*(-c*x + Rt(-a*c, S(2)))))/((m + S(1))*(c*d + e*Rt(-a*c, S(2)))))
    rubi.add(rule122)

    pattern123 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons424, cons87, cons88)
    rule123 = ReplacementRule(pattern123, lambda b, p, e, a, m, x, d, c : m*(-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule123)

    pattern124 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons424, cons87, cons88)
    rule124 = ReplacementRule(pattern124, lambda d, p, e, a, m, x, c : -d*m*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*(p + S(1))) - x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))))
    rubi.add(rule124)

    pattern125 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons424)
    rule125 = ReplacementRule(pattern125, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + (-b*e + S(2)*c*d)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule125)

    pattern126 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons424)
    rule126 = ReplacementRule(pattern126, lambda d, p, e, a, m, x, c : c*d*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule126)

    pattern127 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons87, cons116, cons480, cons1, cons481, cons482)
    rule127 = ReplacementRule(pattern127, lambda b, p, e, a, m, x, d, c : -p*Int((b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule127)

    pattern128 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons462, cons87, cons116, cons480, cons1, cons481, cons483)
    rule128 = ReplacementRule(pattern128, lambda d, p, e, a, m, x, c : -S(2)*c*p*Int(x*(a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1)), x)/(e*(m + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(1))))
    rubi.add(rule128)

    pattern129 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons87, cons116, cons420, cons484, cons485, cons482)
    rule129 = ReplacementRule(pattern129, lambda b, p, e, a, m, x, d, c : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d), x), x)/(e*(m + S(2)*p + S(1))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(1))))
    rubi.add(rule129)

    pattern130 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons462, cons87, cons116, cons420, cons484, cons485, cons483)
    rule130 = ReplacementRule(pattern130, lambda d, p, e, a, m, x, c : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*e - c*d*x, x), x)/(e*(m + S(2)*p + S(1))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))/(e*(m + S(2)*p + S(1))))
    rubi.add(rule130)

    pattern131 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons236, cons88, cons121, cons486, cons482)
    rule131 = ReplacementRule(pattern131, lambda b, p, e, a, m, x, d, c : (b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(b*e*m + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(m + S(2)*p + S(3))), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule131)

    pattern132 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons462, cons236, cons88, cons121, cons486, cons483)
    rule132 = ReplacementRule(pattern132, lambda d, p, e, a, m, x, c : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(d*(S(2)*p + S(3)) + e*x*(m + S(2)*p + S(3))), x)/(S(2)*a*(p + S(1))))
    rubi.add(rule132)

    pattern133 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons407, cons461, cons421, cons236, cons88, cons119, cons482)
    rule133 = ReplacementRule(pattern133, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*e + b*d + x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*c*d**S(2)*(S(2)*p + S(3)) + e*x*(b*e - S(2)*c*d)*(m + S(2)*p + S(2)) + e*(S(2)*a*e*(m + S(-1)) + b*d*(-m + S(2)*p + S(4))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule133)

    pattern134 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons462, cons236, cons88, cons119, cons483)
    rule134 = ReplacementRule(pattern134, lambda d, p, e, a, m, x, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(a*e - c*d*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(a*e**S(2)*(m + S(-1)) - c*d**S(2)*(S(2)*p + S(3)) - c*d*e*x*(m + S(2)*p + S(2)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule134)

    pattern135 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons407, cons461, cons421, cons87, cons88, cons482)
    rule135 = ReplacementRule(pattern135, lambda b, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*e - b**S(2)*e + b*c*d + c*x*(-b*e + S(2)*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3)) - c*e*x*(-b*e + S(2)*c*d)*(m + S(2)*p + S(4)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule135)

    pattern136 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons462, cons87, cons88, cons483)
    rule136 = ReplacementRule(pattern136, lambda d, p, e, a, m, x, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(a*e + c*d*x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(a*e**S(2)*(m + S(2)*p + S(3)) + c*d**S(2)*(S(2)*p + S(3)) + c*d*e*x*(m + S(2)*p + S(4)), x), x)/(S(2)*a*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule136)

    pattern137 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons487, cons420, cons482)
    rule137 = ReplacementRule(pattern137, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(1))) + Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**p*Simp(c*d**S(2)*(m + S(2)*p + S(1)) + e*x*(m + p)*(-b*e + S(2)*c*d) - e*(a*e*(m + S(-1)) + b*d*(p + S(1))), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule137)

    pattern138 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons487, cons420, cons483)
    rule138 = ReplacementRule(pattern138, lambda d, p, e, a, m, x, c : e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))/(c*(m + S(2)*p + S(1))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-2))*Simp(-a*e**S(2)*(m + S(-1)) + c*d**S(2)*(m + S(2)*p + S(1)) + S(2)*c*d*e*x*(m + p), x), x)/(c*(m + S(2)*p + S(1))))
    rubi.add(rule138)

    pattern139 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons488)
    rule139 = ReplacementRule(pattern139, lambda b, p, e, a, m, x, d, c : e*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*e*(m + p + S(2)) + c*d*(m + S(1)) - c*e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule139)

    pattern140 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_, x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons489)
    rule140 = ReplacementRule(pattern140, lambda d, p, e, a, m, x, c : c*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(d*(m + S(1)) - e*x*(m + S(2)*p + S(3)), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule140)

    pattern141 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons72, cons421, cons490, cons491, )
    def With141(b, e, a, x, d, c):
        q = Rt(S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) + S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x - q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule141 = ReplacementRule(pattern141, lambda b, e, a, x, d, c : With141(b, e, a, x, d, c))
    rubi.add(rule141)

    pattern142 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons492, )
    def With142(d, e, a, x, c):
        q = Rt(S(6)*c**S(2)*e**S(2)/d**S(2), S(3))
        return -sqrt(S(3))*c*e*ArcTan(S(2)*sqrt(S(3))*c*(d - e*x)/(S(3)*d*q*(a + c*x**S(2))**(S(1)/3)) + sqrt(S(3))/S(3))/(d**S(2)*q**S(2)) - S(3)*c*e*log(d + e*x)/(S(2)*d**S(2)*q**S(2)) + S(3)*c*e*log(c*d - c*e*x - d*q*(a + c*x**S(2))**(S(1)/3))/(S(2)*d**S(2)*q**S(2))
    rule142 = ReplacementRule(pattern142, lambda d, e, a, x, c : With142(d, e, a, x, c))
    rubi.add(rule142)

    pattern143 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons72, cons421, cons490, cons493, )
    def With143(b, e, a, x, d, c):
        q = Rt(-S(3)*c*e**S(2)*(-b*e + S(2)*c*d), S(3))
        return -sqrt(S(3))*c*e*ArcTan(sqrt(S(3))/S(3) - S(2)*sqrt(S(3))*(-b*e + c*d - c*e*x)/(S(3)*q*(a + b*x + c*x**S(2))**(S(1)/3)))/q**S(2) - S(3)*c*e*log(d + e*x)/(S(2)*q**S(2)) + S(3)*c*e*log(-b*e + c*d - c*e*x + q*(a + b*x + c*x**S(2))**(S(1)/3))/(S(2)*q**S(2))
    rule143 = ReplacementRule(pattern143, lambda b, e, a, x, d, c : With143(b, e, a, x, d, c))
    rubi.add(rule143)

    pattern144 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**(S(1)/3)), x_), cons4, cons5, cons9, cons10, cons72, cons407, cons494, )
    def With144(b, e, a, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)*Int(S(1)/((d + e*x)*(b + S(2)*c*x - q)**(S(1)/3)*(b + S(2)*c*x + q)**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule144 = ReplacementRule(pattern144, lambda b, e, a, x, d, c : With144(b, e, a, x, d, c))
    rubi.add(rule144)

    pattern145 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/4)*(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons462)
    rule145 = ReplacementRule(pattern145, lambda d, e, a, x, c : d*Int(S(1)/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(1)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule145)

    pattern146 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/4)*(d_ + x_*WC('e', S(1)))), x_), cons4, cons9, cons10, cons72, cons462)
    rule146 = ReplacementRule(pattern146, lambda d, e, a, x, c : d*Int(S(1)/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x) - e*Int(x/((a + c*x**S(2))**(S(3)/4)*(d**S(2) - e**S(2)*x**S(2))), x))
    rubi.add(rule146)

    pattern147 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons74, cons413, cons410)
    rule147 = ReplacementRule(pattern147, lambda b, p, e, a, x, d, c : (-S(4)*c/(-S(4)*a*c + b**S(2)))**(-p)*Subst(Int(Simp(-x**S(2)/(-S(4)*a*c + b**S(2)) + S(1), x)**p/Simp(-b*e + S(2)*c*d + e*x, x), x), x, b + S(2)*c*x))
    rubi.add(rule147)

    pattern148 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons74, cons495, cons410)
    rule148 = ReplacementRule(pattern148, lambda b, p, e, a, x, d, c : (-c*(a + b*x + c*x**S(2))/(-S(4)*a*c + b**S(2)))**(-p)*(a + b*x + c*x**S(2))**p*Int((-a*c/(-S(4)*a*c + b**S(2)) - b*c*x/(-S(4)*a*c + b**S(2)) - c**S(2)*x**S(2)/(-S(4)*a*c + b**S(2)))**p/(d + e*x), x))
    rubi.add(rule148)

    pattern149 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons100, cons17, cons475)
    rule149 = ReplacementRule(pattern149, lambda d, p, e, a, m, x, c : Int((d + e*x)**m*(-x*Rt(-c, S(2)) + Rt(a, S(2)))**p*(x*Rt(-c, S(2)) + Rt(a, S(2)))**p, x))
    rubi.add(rule149)

    pattern150 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons461, cons421, cons100, cons27, )
    def With150(b, p, e, a, m, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -(e*(b + S(2)*c*x - q)/(S(2)*c*(d + e*x)))**(-p)*(e*(b + S(2)*c*x + q)/(S(2)*c*(d + e*x)))**(-p)*(a + b*x + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x*(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule150 = ReplacementRule(pattern150, lambda b, p, e, a, m, x, d, c : With150(b, p, e, a, m, x, d, c))
    rubi.add(rule150)

    pattern151 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons74, cons462, cons100, cons27, )
    def With151(d, p, e, a, m, x, c):
        q = Rt(-a*c, S(2))
        return -(e*(c*x + q)/(c*(d + e*x)))**(-p)*(-e*(-c*x + q)/(c*(d + e*x)))**(-p)*(a + c*x**S(2))**p*(S(1)/(d + e*x))**(S(2)*p)*Subst(Int(x**(-m - S(2)*p + S(-2))*Simp(-x*(d - e*q/c) + S(1), x)**p*Simp(-x*(d + e*q/c) + S(1), x)**p, x), x, S(1)/(d + e*x))/e
    rule151 = ReplacementRule(pattern151, lambda d, p, e, a, m, x, c : With151(d, p, e, a, m, x, c))
    rubi.add(rule151)

    pattern152 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons407, cons461, cons421, cons100, )
    def With152(b, p, e, a, m, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (-(d + e*x)/(d - e*(b - q)/(S(2)*c)) + S(1))**(-p)*(-(d + e*x)/(d - e*(b + q)/(S(2)*c)) + S(1))**(-p)*(a + b*x + c*x**S(2))**p*Subst(Int(x**m*Simp(-x/(d - e*(b - q)/(S(2)*c)) + S(1), x)**p*Simp(-x/(d - e*(b + q)/(S(2)*c)) + S(1), x)**p, x), x, d + e*x)/e
    rule152 = ReplacementRule(pattern152, lambda b, p, e, a, m, x, d, c : With152(b, p, e, a, m, x, d, c))
    rubi.add(rule152)

    pattern153 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons462, cons100, )
    def With153(d, p, e, a, m, x, c):
        q = Rt(-a*c, S(2))
        return (a + c*x**S(2))**p*(-(d + e*x)/(d - e*q/c) + S(1))**(-p)*(-(d + e*x)/(d + e*q/c) + S(1))**(-p)*Subst(Int(x**m*Simp(-x/(d - e*q/c) + S(1), x)**p*Simp(-x/(d + e*q/c) + S(1), x)**p, x), x, d + e*x)/e
    rule153 = ReplacementRule(pattern153, lambda d, p, e, a, m, x, c : With153(d, p, e, a, m, x, c))
    rubi.add(rule153)

    pattern154 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons2, cons74, cons6, cons7)
    rule154 = ReplacementRule(pattern154, lambda b, u, p, e, a, m, x, d, c : Subst(Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule154)

    pattern155 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons2, cons74, cons6, cons7)
    rule155 = ReplacementRule(pattern155, lambda u, p, e, a, m, x, d, c : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule155)

    pattern156 = Pattern(Integral(x_**WC('n', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons74, cons28, cons496)
    rule156 = ReplacementRule(pattern156, lambda p, e, a, n, x, d, c : d*Int(x**n*(a + c*x**S(2))**p, x) + e*Int(x**(n + S(1))*(a + c*x**S(2))**p, x))
    rubi.add(rule156)

    pattern157 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))/sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons405, cons498)
    rule157 = ReplacementRule(pattern157, lambda g, b, f, e, a, m, x, d, c : (f + g*x)*Int((d + e*x)**m, x)/sqrt(a + b*x + c*x**S(2)))
    rubi.add(rule157)

    pattern158 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons405, cons498, cons100, cons424)
    rule158 = ReplacementRule(pattern158, lambda g, b, p, f, e, a, m, x, d, c : -f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule158)

    pattern159 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons405, cons498, cons100, cons236, cons88, cons121)
    rule159 = ReplacementRule(pattern159, lambda g, b, f, p, e, a, m, x, d, c : -e*g*m*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(S(2)*c*(p + S(1))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*c*(p + S(1))))
    rubi.add(rule159)

    pattern160 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons405, cons498, cons100, cons87, cons88, cons499)
    rule160 = ReplacementRule(pattern160, lambda g, b, f, p, e, a, m, x, d, c : e*f*g*(m + S(2)*p + S(3))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1)), x)/(b*(p + S(1))*(-d*g + e*f)) - f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(p + S(1))*(-d*g + e*f)))
    rubi.add(rule160)

    pattern161 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons405, cons498, cons100, cons51, cons38, cons406, cons500)
    rule161 = ReplacementRule(pattern161, lambda g, b, f, p, e, a, m, x, d, c : -g*(S(2)*p + S(1))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))) + (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))))
    rubi.add(rule161)

    pattern162 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons405, cons498, cons100, cons51, cons38, cons501)
    rule162 = ReplacementRule(pattern162, lambda g, b, f, p, e, a, m, x, d, c : -g*(m + S(2)*p + S(3))*Int((d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/((m + S(1))*(-d*g + e*f)) + S(2)*f*g*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(b*(m + S(1))*(-d*g + e*f)))
    rubi.add(rule162)

    pattern163 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons405, cons498, cons100, cons25, cons501, cons502)
    rule163 = ReplacementRule(pattern163, lambda g, b, f, p, e, a, m, x, d, c : -b*m*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*f*(m + S(2)*p + S(2))) + g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule163)

    pattern164 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons405, cons498, cons100, cons501)
    rule164 = ReplacementRule(pattern164, lambda g, b, p, f, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(f + g*x)*(a + b*x + c*x**S(2))**p/(e*(m + S(2)*p + S(2))) + (S(2)*p + S(1))*(-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(2)*p + S(2))))
    rubi.add(rule164)

    pattern165 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons405, cons503, cons100, cons421, cons87, cons504)
    rule165 = ReplacementRule(pattern165, lambda g, b, f, p, e, a, x, d, c : (-b*g + S(2)*c*f)*Int((a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(a + b*x + c*x**S(2))**p/(d + e*x), x)/(-b*e + S(2)*c*d))
    rubi.add(rule165)

    pattern166 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons405, cons503, cons100, cons505)
    rule166 = ReplacementRule(pattern166, lambda g, b, f, p, e, a, m, x, d, c : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule166)

    pattern167 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons405, cons503, cons100, cons421, cons424)
    rule167 = ReplacementRule(pattern167, lambda g, b, f, p, e, a, m, x, d, c : (-b*g + S(2)*c*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d) - (-d*g + e*f)*Int((b + S(2)*c*x)*(d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(-b*e + S(2)*c*d))
    rubi.add(rule167)

    pattern168 = Pattern(Integral((f_ + x_*WC('g', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons405, cons503, cons100, cons421, cons501, cons452, cons51, cons38)
    rule168 = ReplacementRule(pattern168, lambda g, b, f, p, e, a, m, x, d, c : -(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**p/(e*(m + S(1))*(-b*e + S(2)*c*d)) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*d*(S(2)*p + S(1))))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(m + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule168)

    pattern169 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons405, cons503, cons100, cons421, cons501, cons452, cons454, cons506)
    rule169 = ReplacementRule(pattern169, lambda g, b, f, p, e, a, m, x, d, c : g*(b + S(2)*c*x)*(d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*c*e*(m + S(2)*p + S(2))) + (S(2)*c*e*f*(m + S(2)*p + S(2)) - g*(b*e*(m + S(1)) + S(2)*c*(S(2)*d*p + d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(S(2)*c*e*(m + S(2)*p + S(2))))
    rubi.add(rule169)

    pattern170 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons405, cons100)
    rule170 = ReplacementRule(pattern170, lambda g, b, f, p, e, a, m, n, x, d, c : c**(-IntPart(p))*(b/S(2) + c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b/S(2) + c*x)**(S(2)*p)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule170)

    pattern171 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons407, cons437, cons97)
    rule171 = ReplacementRule(pattern171, lambda g, b, d, f, p, e, a, m, n, x, c : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule171)

    pattern172 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons438, cons507)
    rule172 = ReplacementRule(pattern172, lambda g, d, f, p, e, a, m, n, x, c : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule172)

    pattern173 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons13, cons74, cons437, cons27, cons496)
    rule173 = ReplacementRule(pattern173, lambda g, b, d, f, p, e, a, m, n, x, c : d**m*e**m*Int((f + g*x)**n*(a*e + c*d*x)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule173)

    pattern174 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons13, cons74, cons438, cons27, cons496)
    rule174 = ReplacementRule(pattern174, lambda g, d, f, p, e, a, m, n, x, c : d**m*e**m*Int((a + c*x**S(2))**(m + p)*(f + g*x)**n*(a*e + c*d*x)**(-m), x))
    rubi.add(rule174)

    pattern175 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons508)
    rule175 = ReplacementRule(pattern175, lambda g, b, f, p, e, a, m, x, d, c : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))))
    rubi.add(rule175)

    pattern176 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons509)
    rule176 = ReplacementRule(pattern176, lambda g, d, f, p, e, a, m, x, c : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule176)

    pattern177 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons236, cons88, cons121)
    rule177 = ReplacementRule(pattern177, lambda g, b, f, p, e, a, m, x, d, c : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule177)

    pattern178 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons236, cons88, cons121)
    rule178 = ReplacementRule(pattern178, lambda g, f, p, e, a, m, x, d, c : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule178)

    pattern179 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons91, cons510, cons85)
    rule179 = ReplacementRule(pattern179, lambda g, b, f, p, e, a, m, x, d, c : -e*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-b*e + S(2)*c*d)) + (d + e*x)**m*(c*e*f + g*(-b*e + c*d))*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))*(-b*e + S(2)*c*d)))
    rubi.add(rule179)

    pattern180 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons91, cons510, cons85)
    rule180 = ReplacementRule(pattern180, lambda g, d, f, p, e, a, m, x, c : -e*(S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*c*d*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g + e*f)/(S(2)*c*d*(p + S(1))))
    rubi.add(rule180)

    pattern181 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons511, cons434)
    rule181 = ReplacementRule(pattern181, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**m*(d*g - e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((-b*e + S(2)*c*d)*(m + p + S(1))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(e*(-b*e + S(2)*c*d)*(m + p + S(1))))
    rubi.add(rule181)

    pattern182 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons511, cons434)
    rule182 = ReplacementRule(pattern182, lambda g, d, f, p, e, a, m, x, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(d*g - e*f)/(S(2)*c*d*(m + p + S(1))) + (S(2)*c*e*f*(p + S(1)) + m*(c*d*g + c*e*f))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(S(2)*c*d*e*(m + p + S(1))))
    rubi.add(rule182)

    pattern183 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons501)
    rule183 = ReplacementRule(pattern183, lambda g, b, f, p, e, a, m, x, d, c : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + (e*(p + S(1))*(-b*g + S(2)*c*f) + m*(c*e*f + g*(-b*e + c*d)))*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/(c*e*(m + S(2)*p + S(2))))
    rubi.add(rule183)

    pattern184 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons501)
    rule184 = ReplacementRule(pattern184, lambda g, d, f, p, e, a, m, x, c : (S(2)*e*f*(p + S(1)) + m*(d*g + e*f))*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/(e*(m + S(2)*p + S(2))) + g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))))
    rubi.add(rule184)

    pattern185 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons4, cons9, cons73, cons161, cons512, cons87, cons513)
    rule185 = ReplacementRule(pattern185, lambda g, f, p, a, x, c : x**S(2)*(a + c*x**S(2))**(p + S(1))*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int(x*(a + c*x**S(2))**(p + S(1))*Simp(S(2)*a*g - c*f*x*(S(2)*p + S(5)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule185)

    pattern186 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons4, cons9, cons73, cons161, cons74, cons512)
    rule186 = ReplacementRule(pattern186, lambda g, f, p, a, x, c : -f**S(2)*Int((a + c*x**S(2))**(p + S(1))/(f - g*x), x)/c + Int((a + c*x**S(2))**(p + S(1))*(f + g*x), x)/c)
    rubi.add(rule186)

    pattern187 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons100, cons103, cons87, cons514)
    rule187 = ReplacementRule(pattern187, lambda g, b, d, f, p, e, a, m, n, x, c : Int((f + g*x)**n*(a/d + c*x/e)**(-m)*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule187)

    pattern188 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons100, cons103, cons87, cons514)
    rule188 = ReplacementRule(pattern188, lambda g, d, f, p, e, a, m, n, x, c : a**(-m)*d**(S(2)*m)*Int((a + c*x**S(2))**(m + p)*(d - e*x)**(-m)*(f + g*x)**n, x))
    rubi.add(rule188)

    pattern189 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons100, cons101, cons515)
    rule189 = ReplacementRule(pattern189, lambda g, b, d, f, p, e, a, n, x, c : -(f + g*x)**n*(a*(-b*e + S(2)*c*d) + c*x*(-S(2)*a*e + b*d))*(a + b*x + c*x**S(2))**p/(d*e*p*(-S(4)*a*c + b**S(2))) - Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*Simp(-S(2)*a*c*(d*g*n - e*f*(S(2)*p + S(1))) + b*(a*e*g*n - c*d*f*(S(2)*p + S(1))) - c*g*x*(-S(2)*a*e + b*d)*(n + S(2)*p + S(1)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))))
    rubi.add(rule189)

    pattern190 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons100, cons101, cons515)
    rule190 = ReplacementRule(pattern190, lambda g, d, f, p, e, a, n, x, c : (a + c*x**S(2))**p*(d - e*x)*(f + g*x)**n/(S(2)*d*e*p) - Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*Simp(d*g*n - e*f*(S(2)*p + S(1)) - e*g*x*(n + S(2)*p + S(1)), x), x)/(S(2)*d*e*p))
    rubi.add(rule190)

    pattern191 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons100, cons149, cons515)
    rule191 = ReplacementRule(pattern191, lambda g, b, d, f, p, e, a, n, x, c : -(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p*(a*c*d*(-b*g + S(2)*c*f) - a*e*(S(2)*a*c*g - b**S(2)*g + b*c*f) + c*x*(-a*e*(-b*g + S(2)*c*f) + c*d*(-S(2)*a*g + b*f)))/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**p*Simp(S(2)*a*c*(a*e*g**S(2)*(n + S(2)*p + S(1)) + c*f*(-d*g*n + S(2)*e*f*p + e*f)) + b**S(2)*g*(-a*e*g*(n + p + S(1)) + c*d*f*p) + b*c*(a*g*(d*g*(n + S(1)) + e*f*(n - S(2)*p)) - c*d*f**S(2)*(S(2)*p + S(1))) + c*g*x*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f))*(n + S(2)*p + S(2)), x), x)/(d*e*p*(-S(4)*a*c + b**S(2))*(a*g**S(2) - b*f*g + c*f**S(2))))
    rubi.add(rule191)

    pattern192 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))/(d_ + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons100, cons149, cons515)
    rule192 = ReplacementRule(pattern192, lambda g, d, f, p, e, a, n, x, c : (a + c*x**S(2))**p*(f + g*x)**(n + S(1))*(-a*e*g + c*d*f - c*x*(d*g + e*f))/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))) + Int((a + c*x**S(2))**p*(f + g*x)**n*Simp(a*e*g**S(2)*(n + S(2)*p + S(1)) - c*f*(d*g*n - e*(S(2)*f*p + f)) + c*g*x*(d*g + e*f)*(n + S(2)*p + S(2)), x), x)/(S(2)*d*e*p*(a*g**S(2) + c*f**S(2))))
    rubi.add(rule192)

    pattern193 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons407, cons437, cons100, cons440, cons516, cons517)
    rule193 = ReplacementRule(pattern193, lambda g, b, d, f, p, e, a, m, n, x, c : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))))
    rubi.add(rule193)

    pattern194 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons438, cons100, cons440, cons518, cons517)
    rule194 = ReplacementRule(pattern194, lambda g, d, f, p, e, a, m, n, x, c : -e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule194)

    pattern195 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons407, cons437, cons100, cons440, cons519)
    rule195 = ReplacementRule(pattern195, lambda g, b, d, f, p, e, a, m, n, x, c : -e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule195)

    pattern196 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons438, cons100, cons440, cons519)
    rule196 = ReplacementRule(pattern196, lambda g, f, p, e, a, m, n, x, d, c : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(n + S(1))*(d*g + e*f)))
    rubi.add(rule196)

    pattern197 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons100, cons440, cons520, cons116, cons32, cons521)
    rule197 = ReplacementRule(pattern197, lambda g, b, d, f, p, e, a, m, n, x, c : c*m*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e*g*(n + S(1))) + (d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(n + S(1))))
    rubi.add(rule197)

    pattern198 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons100, cons440, cons520, cons116, cons32, cons521)
    rule198 = ReplacementRule(pattern198, lambda g, f, p, e, a, m, n, x, d, c : c*m*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(1)), x)/(e*g*(n + S(1))) + (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(n + S(1))))
    rubi.add(rule198)

    pattern199 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons13, cons497, cons407, cons437, cons100, cons440, cons520, cons116, cons517, cons522, cons523)
    rule199 = ReplacementRule(pattern199, lambda g, b, d, f, p, e, a, m, n, x, c : -(d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p/(g*(m - n + S(-1))) - m*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**(m + S(1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1)), x)/(e**S(2)*g*(m - n + S(-1))))
    rubi.add(rule199)

    pattern200 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons13, cons497, cons438, cons100, cons440, cons520, cons116, cons517, cons522, cons523)
    rule200 = ReplacementRule(pattern200, lambda g, d, f, p, e, a, m, n, x, c : -c*m*(d*g + e*f)*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*(f + g*x)**n, x)/(e**S(2)*g*(m - n + S(-1))) - (a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1))/(g*(m - n + S(-1))))
    rubi.add(rule200)

    pattern201 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437, cons100, cons440, cons520, cons88, cons31)
    rule201 = ReplacementRule(pattern201, lambda g, b, d, f, p, e, a, m, n, x, c : -e*g*n*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))) + e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(p + S(1))))
    rubi.add(rule201)

    pattern202 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438, cons100, cons440, cons520, cons88, cons31)
    rule202 = ReplacementRule(pattern202, lambda g, d, f, p, e, a, m, n, x, c : -e*g*n*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/(c*(p + S(1))) + e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(p + S(1))))
    rubi.add(rule202)

    pattern203 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons13, cons497, cons407, cons437, cons100, cons440, cons520, cons88)
    rule203 = ReplacementRule(pattern203, lambda g, b, d, f, p, e, a, m, n, x, c : e**S(2)*g*(m - n + S(-2))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-b*e*g + c*d*g + c*e*f)) + e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule203)

    pattern204 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons13, cons497, cons438, cons100, cons440, cons520, cons88)
    rule204 = ReplacementRule(pattern204, lambda g, d, f, p, e, a, m, n, x, c : e**S(2)*g*(m - n + S(-2))*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(c*(p + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/(c*(p + S(1))*(d*g + e*f)))
    rubi.add(rule204)

    pattern205 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons100, cons440, cons30, cons31, cons517, cons524)
    rule205 = ReplacementRule(pattern205, lambda g, b, d, f, p, e, a, m, n, x, c : -e*(d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m - n + S(-1))) - n*(-b*e*g + c*d*g + c*e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/(c*e*(m - n + S(-1))))
    rubi.add(rule205)

    pattern206 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons100, cons440, cons30, cons31, cons517, cons524)
    rule206 = ReplacementRule(pattern206, lambda g, d, f, p, e, a, m, n, x, c : -n*(d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/(e*(m - n + S(-1))) - e*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**n/(c*(m - n + S(-1))))
    rubi.add(rule206)

    pattern207 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons100, cons440, cons30, cons32, cons427)
    rule207 = ReplacementRule(pattern207, lambda g, b, d, f, p, e, a, m, n, x, c : -c*e*(m - n + S(-2))*Int((d + e*x)**m*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/((n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e**S(2)*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/((n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule207)

    pattern208 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons100, cons440, cons30, cons32, cons427)
    rule208 = ReplacementRule(pattern208, lambda g, f, p, e, a, m, n, x, d, c : -e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))/((n + S(1))*(c*d*g + c*e*f)) - e*(m - n + S(-2))*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(1)), x)/((n + S(1))*(d*g + e*f)))
    rubi.add(rule208)

    pattern209 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/((x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons437)
    rule209 = ReplacementRule(pattern209, lambda g, b, d, f, e, a, x, c : S(2)*e**S(2)*Subst(Int(S(1)/(-b*e*g + c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + b*x + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule209)

    pattern210 = Pattern(Integral(sqrt(d_ + x_*WC('e', S(1)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons438)
    rule210 = ReplacementRule(pattern210, lambda g, f, e, a, x, d, c : S(2)*e**S(2)*Subst(Int(S(1)/(c*(d*g + e*f) + e**S(2)*g*x**S(2)), x), x, sqrt(a + c*x**S(2))/sqrt(d + e*x)))
    rubi.add(rule210)

    pattern211 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons407, cons437, cons100, cons525, cons526, cons75)
    rule211 = ReplacementRule(pattern211, lambda g, b, d, f, p, e, a, m, n, x, c : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule211)

    pattern212 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons438, cons100, cons525, cons527, cons75)
    rule212 = ReplacementRule(pattern212, lambda g, d, f, p, e, a, m, n, x, c : e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule212)

    pattern213 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons100, cons525, cons30, cons32, cons427)
    rule213 = ReplacementRule(pattern213, lambda g, b, d, f, p, e, a, m, n, x, c : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)) - e*(b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**p, x)/(g*(n + S(1))*(-b*e*g + c*d*g + c*e*f)))
    rubi.add(rule213)

    pattern214 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons100, cons525, cons30, cons32, cons427)
    rule214 = ReplacementRule(pattern214, lambda g, f, p, e, a, m, n, x, d, c : -e*(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1)), x)/(g*(n + S(1))*(d*g + e*f)) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(-d*g + e*f)/(c*g*(n + S(1))*(d*g + e*f)))
    rubi.add(rule214)

    pattern215 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons407, cons437, cons100, cons525, cons528, cons427)
    rule215 = ReplacementRule(pattern215, lambda g, b, d, f, p, e, a, m, n, x, c : e**S(2)*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*g*(n + p + S(2))) - (b*e*g*(n + S(1)) - c*d*g*(S(2)*n + p + S(3)) + c*e*f*(p + S(1)))*Int((d + e*x)**(m + S(-1))*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x)/(c*g*(n + p + S(2))))
    rubi.add(rule215)

    pattern216 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons497, cons438, cons100, cons525, cons528, cons427)
    rule216 = ReplacementRule(pattern216, lambda g, d, f, p, e, a, m, n, x, c : -(-d*g*(S(2)*n + p + S(3)) + e*f*(p + S(1)))*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(g*(n + p + S(2))) + e**S(2)*(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*(f + g*x)**(n + S(1))/(c*g*(n + p + S(2))))
    rubi.add(rule216)

    pattern217 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons13, cons74, cons497, cons407, cons437, cons100, cons160)
    rule217 = ReplacementRule(pattern217, lambda g, b, d, f, p, e, a, m, n, x, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule217)

    pattern218 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons13, cons74, cons497, cons438, cons529, cons103, cons530, cons531)
    rule218 = ReplacementRule(pattern218, lambda g, d, f, p, e, a, m, n, x, c : Int(ExpandIntegrand(S(1)/sqrt(a + c*x**S(2)), (a + c*x**S(2))**(p + S(1)/2)*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule218)

    pattern219 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons13, cons74, cons497, cons438, cons100, cons160)
    rule219 = ReplacementRule(pattern219, lambda g, d, f, p, e, a, m, n, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule219)

    pattern220 = Pattern(Integral(x_**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(d_ + x_*WC('e', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons74, cons407, cons437)
    rule220 = ReplacementRule(pattern220, lambda b, d, p, e, a, x, c : d**S(2)*Int((a + b*x + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((d - e*x)*(a + b*x + c*x**S(2))**p, x)/e**S(2))
    rubi.add(rule220)

    pattern221 = Pattern(Integral(x_**S(2)*(a_ + x_**S(2)*WC('c', S(1)))**p_/(d_ + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons74, cons438)
    rule221 = ReplacementRule(pattern221, lambda d, p, e, a, x, c : d**S(2)*Int((a + c*x**S(2))**p/(d + e*x), x)/e**S(2) - Int((a + c*x**S(2))**p*(d - e*x), x)/e**S(2))
    rubi.add(rule221)

    pattern222 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2)*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons437, cons100, cons452)
    rule222 = ReplacementRule(pattern222, lambda g, b, d, f, p, e, a, m, x, c : g*(d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(3))) - Int((d + e*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*e*g*(d*g + e*f*(m + p + S(1))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))) + e*g*x*(b*e*g*(m + p + S(2)) - c*(d*g*m + e*f*(m + S(2)*p + S(4)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule222)

    pattern223 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**S(2), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons438, cons100, cons452)
    rule223 = ReplacementRule(pattern223, lambda g, f, p, e, a, m, x, d, c : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(f + g*x)/(c*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(d + e*x)**m*Simp(-c*e*g*x*(d*g*m + e*f*(m + S(2)*p + S(4))) - c*(d**S(2)*g**S(2) + d*e*f*g*m + e**S(2)*f**S(2)*(m + S(2)*p + S(3))), x), x)/(c*e**S(2)*(m + S(2)*p + S(3))))
    rubi.add(rule223)

    pattern224 = Pattern(Integral((x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons5, cons9, cons72, cons73, cons161, cons2, cons13, cons100)
    rule224 = ReplacementRule(pattern224, lambda g, b, f, p, e, m, n, x, c : x**(-m - p)*(e*x)**m*(b + c*x)**(-p)*(b*x + c*x**S(2))**p*Int(x**(m + p)*(b + c*x)**p*(f + g*x)**n, x))
    rubi.add(rule224)

    pattern225 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons438, cons100, cons17, cons450)
    rule225 = ReplacementRule(pattern225, lambda g, d, f, p, e, a, m, n, x, c : Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule225)

    pattern226 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons407, cons437, cons100)
    rule226 = ReplacementRule(pattern226, lambda g, b, d, f, p, e, a, m, n, x, c : (d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule226)

    pattern227 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(d_ + x_*WC('e', S(1)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons497, cons438, cons100)
    rule227 = ReplacementRule(pattern227, lambda g, d, f, p, e, a, m, n, x, c : (a + c*x**S(2))**FracPart(p)*(d + e*x)**(-FracPart(p))*(a/d + c*x/e)**(-FracPart(p))*Int((d + e*x)**(m + p)*(f + g*x)**n*(a/d + c*x/e)**p, x))
    rubi.add(rule227)

    pattern228 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons407, cons461, cons77)
    rule228 = ReplacementRule(pattern228, lambda g, b, f, p, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule228)

    pattern229 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons462, cons77)
    rule229 = ReplacementRule(pattern229, lambda g, f, p, e, a, m, x, d, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x), x), x))
    rubi.add(rule229)

    pattern230 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461)
    rule230 = ReplacementRule(pattern230, lambda g, b, f, e, a, x, d, c : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int(Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule230)

    pattern231 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462)
    rule231 = ReplacementRule(pattern231, lambda g, f, e, a, x, d, c : e*(-d*g + e*f)*Int(S(1)/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int(Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule231)

    pattern232 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons461, cons424, cons532)
    rule232 = ReplacementRule(pattern232, lambda g, b, f, p, e, a, m, x, d, c : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule232)

    pattern233 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons462, cons424, cons533)
    rule233 = ReplacementRule(pattern233, lambda g, f, p, e, a, m, x, d, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule233)

    pattern234 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons424, cons87, cons88, cons534)
    rule234 = ReplacementRule(pattern234, lambda g, b, f, p, e, a, m, x, d, c : -m*(-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1)), x)/((p + S(1))*(-S(4)*a*c + b**S(2))) + (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule234)

    pattern235 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons424, cons87, cons88, cons535)
    rule235 = ReplacementRule(pattern235, lambda g, f, p, e, a, m, x, d, c : -m*(a*e*g + c*d*f)*Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1)), x)/(S(2)*a*c*(p + S(1))) + (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule235)

    pattern236 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons461, cons424)
    rule236 = ReplacementRule(pattern236, lambda g, b, f, p, e, a, m, x, d, c : -(d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/(S(2)*(p + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) - (-S(2)*a*e*g + b*(d*g + e*f) - S(2)*c*d*f)*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/(S(2)*a*e**S(2) - S(2)*b*d*e + S(2)*c*d**S(2)))
    rubi.add(rule236)

    pattern237 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons462, cons424)
    rule237 = ReplacementRule(pattern237, lambda g, f, p, e, a, m, x, d, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/(S(2)*(p + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e*g + c*d*f)*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule237)

    pattern238 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons407, cons461, cons536)
    rule238 = ReplacementRule(pattern238, lambda g, b, f, p, e, a, x, d, c : -(a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule238)

    pattern239 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons462, cons537)
    rule239 = ReplacementRule(pattern239, lambda g, f, p, e, a, x, d, c : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule239)

    pattern240 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons87, cons88)
    rule240 = ReplacementRule(pattern240, lambda g, b, f, p, e, a, x, d, c : -(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g - b*c*(d*g + e*f) + S(2)*c*(-a*e*g + c*d*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**(p + S(1)), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule240)

    pattern241 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons87, cons88)
    rule241 = ReplacementRule(pattern241, lambda g, f, p, e, a, x, d, c : -(a + c*x**S(2))**(p + S(1))*(-a*(d*g + e*(f + g*x)) + c*d*f*x)/(S(2)*a*c*(p + S(1))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**(p + S(1)), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule241)

    pattern242 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons407, cons461, cons469)
    rule242 = ReplacementRule(pattern242, lambda g, b, f, p, e, a, x, d, c : (-S(2)*a*c*e*g + b**S(2)*e*g*(p + S(2)) + c*(S(2)*p + S(3))*(-b*(d*g + e*f) + S(2)*c*d*f))*Int((a + b*x + c*x**S(2))**p, x)/(S(2)*c**S(2)*(S(2)*p + S(3))) - (a + b*x + c*x**S(2))**(p + S(1))*(b*e*g*(p + S(2)) - S(2)*c*e*g*x*(p + S(1)) - c*(S(2)*p + S(3))*(d*g + e*f))/(S(2)*c**S(2)*(p + S(1))*(S(2)*p + S(3))))
    rubi.add(rule242)

    pattern243 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons462, cons469)
    rule243 = ReplacementRule(pattern243, lambda g, f, p, e, a, x, d, c : (a + c*x**S(2))**(p + S(1))*(S(2)*e*g*x*(p + S(1)) + (S(2)*p + S(3))*(d*g + e*f))/(S(2)*c*(p + S(1))*(S(2)*p + S(3))) - (a*e*g - c*d*f*(S(2)*p + S(3)))*Int((a + c*x**S(2))**p, x)/(c*(S(2)*p + S(3))))
    rubi.add(rule243)

    pattern244 = Pattern(Integral((x_*WC('e', S(1)))**m_*(a_ + x_**S(2)*WC('c', S(1)))**p_*(f_ + x_*WC('g', S(1))), x_), cons4, cons9, cons72, cons73, cons161, cons74, cons271, cons212)
    rule244 = ReplacementRule(pattern244, lambda g, f, p, e, a, m, x, c : f*Int((e*x)**m*(a + c*x**S(2))**p, x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule244)

    pattern245 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons538, cons470, cons471)
    rule245 = ReplacementRule(pattern245, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((f + g*x)*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule245)

    pattern246 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons236, cons116, cons428, cons539, cons540)
    rule246 = ReplacementRule(pattern246, lambda g, b, f, p, e, a, m, x, d, c : -p*Int((d + e*x)**(m + S(2))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) + b**S(2)*e*(d*g*(p + S(1)) - e*f*(m + p + S(2))) + b*(a*e**S(2)*g*(m + S(1)) - c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))) - c*x*(S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2))) - e*(S(2)*a*e*g*(m + S(1)) - b*(d*g*(m - S(2)*p) + e*f*(m + S(2)*p + S(2))))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) - (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*p*(-b*e + S(2)*c*d)*(-d*g + e*f) - e*x*(g*(m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2)) + p*(-b*e + S(2)*c*d)*(-d*g + e*f)) + (d*g - e*f*(m + S(2)))*(a*e**S(2) - b*d*e + c*d**S(2)))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule246)

    pattern247 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons236, cons116, cons428, cons539, cons540)
    rule247 = ReplacementRule(pattern247, lambda g, f, p, e, a, m, x, d, c : -p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(2))*Simp(S(2)*a*c*e*(m + S(2))*(-d*g + e*f) - c*x*(-S(2)*a*e**S(2)*g*(m + S(1)) + S(2)*c*d*(d*g*(S(2)*p + S(1)) - e*f*(m + S(2)*p + S(2)))), x), x)/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))) - (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-S(2)*c*d**S(2)*p*(-d*g + e*f) - e*x*(S(2)*c*d*p*(-d*g + e*f) + g*(m + S(1))*(a*e**S(2) + c*d**S(2))) + (a*e**S(2) + c*d**S(2))*(d*g - e*f*(m + S(2))))/(e**S(2)*(m + S(1))*(m + S(2))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule247)

    pattern248 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons407, cons461, cons87, cons116, cons541, cons1, cons481, cons542)
    rule248 = ReplacementRule(pattern248, lambda g, b, f, p, e, a, m, x, d, c : p*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-b*e*f*(m + S(2)*p + S(2)) + g*(S(2)*a*e*m + S(2)*a*e + S(2)*b*d*p + b*d) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(b*e*m + b*e + S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule248)

    pattern249 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons462, cons87, cons116, cons541, cons1, cons481, cons542)
    rule249 = ReplacementRule(pattern249, lambda g, f, p, e, a, m, x, d, c : p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**(m + S(1))*Simp(g*(S(2)*a*e*m + S(2)*a*e) + x*(-S(2)*c*e*f*(m + S(2)*p + S(2)) + g*(S(4)*c*d*p + S(2)*c*d)), x), x)/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-d*g*(S(2)*p + S(1)) + e*f*(m + S(2)*p + S(2)) + e*g*x*(m + S(1)))/(e**S(2)*(m + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule249)

    pattern250 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons407, cons461, cons87, cons116, cons543, cons485, cons542)
    rule250 = ReplacementRule(pattern250, lambda g, b, f, p, e, a, m, x, d, c : -p*Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(c*e*f*(-S(2)*a*e + b*d)*(m + S(2)*p + S(2)) + g*(a*e*(b*e*m + b*e - S(2)*c*d*m) + b*d*(b*e*p - S(2)*c*d*p - c*d)) + x*(c*e*f*(-b*e + S(2)*c*d)*(m + S(2)*p + S(2)) + g*(b**S(2)*e**S(2)*(m + p + S(1)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(1)) - c*e*(S(2)*a*e*(m + S(2)*p + S(1)) + b*d*(m - S(2)*p)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*(c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)) - g*(-b*e*p + S(2)*c*d*p + c*d))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule250)

    pattern251 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons462, cons87, cons116, cons543, cons485, cons542)
    rule251 = ReplacementRule(pattern251, lambda g, f, p, e, a, m, x, d, c : S(2)*p*Int((a + c*x**S(2))**(p + S(-1))*(d + e*x)**m*Simp(a*c*d*e*g*m + a*c*e**S(2)*f*(m + S(2)*p + S(2)) - x*(c**S(2)*d*e*f*(m + S(2)*p + S(2)) - g*(a*c*e**S(2)*(m + S(2)*p + S(1)) + c**S(2)*d**S(2)*(S(2)*p + S(1)))), x), x)/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))) + (a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(-c*d*g*(S(2)*p + S(1)) + c*e*f*(m + S(2)*p + S(2)) + c*e*g*x*(m + S(2)*p + S(1)))/(c*e**S(2)*(m + S(2)*p + S(1))*(m + S(2)*p + S(2))))
    rubi.add(rule251)

    pattern252 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons236, cons88, cons119, cons544)
    rule252 = ReplacementRule(pattern252, lambda g, b, f, p, e, a, m, x, d, c : -(d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*(S(2)*a*c*(d*g + e*f) - b*(a*e*g + c*d*f) - x*(b**S(2)*e*g + S(2)*c**S(2)*d*f - c*(S(2)*a*e*g + b*d*g + b*e*f)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + e*x)**(m + S(-2))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(b*e*g*(a*e*(m + S(-1)) + b*d*(p + S(2))) + S(2)*c**S(2)*d**S(2)*f*(S(2)*p + S(3)) - c*(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) + b*d*(d*g*(S(2)*p + S(3)) - e*f*(m - S(2)*p + S(-4)))) + e*x*(b**S(2)*e*g*(m + p + S(1)) + S(2)*c**S(2)*d*f*(m + S(2)*p + S(2)) - c*(S(2)*a*e*g*m + b*(d*g + e*f)*(m + S(2)*p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule252)

    pattern253 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons236, cons88, cons119, cons545)
    rule253 = ReplacementRule(pattern253, lambda g, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*(S(2)*a*(d*g + e*f) - x*(-S(2)*a*e*g + S(2)*c*d*f))/(S(4)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-2))*Simp(S(2)*a*e*(d*g*m + e*f*(m + S(-1))) - S(2)*c*d**S(2)*f*(S(2)*p + S(3)) + e*x*(S(2)*a*e*g*m - S(2)*c*d*f*(m + S(2)*p + S(2))), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule253)

    pattern254 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons236, cons88, cons121, cons546, cons542)
    rule254 = ReplacementRule(pattern254, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*g + b*f + x*(-b*g + S(2)*c*f))/((p + S(1))*(-S(4)*a*c + b**S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-e*x*(-b*g + S(2)*c*f)*(m + S(2)*p + S(3)) - f*(b*e*m + S(2)*c*d*(S(2)*p + S(3))) + g*(S(2)*a*e*m + b*d*(S(2)*p + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule254)

    pattern255 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons236, cons88, cons121, cons546, cons542)
    rule255 = ReplacementRule(pattern255, lambda g, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**m*(a*g - c*f*x)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(-1))*Simp(a*e*g*m - c*d*f*(S(2)*p + S(3)) - c*e*f*x*(m + S(2)*p + S(3)), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule255)

    pattern256 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons407, cons461, cons87, cons88, cons542)
    rule256 = ReplacementRule(pattern256, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-a*g*(-b*e + S(2)*c*d) + c*x*(f*(-b*e + S(2)*c*d) - g*(-S(2)*a*e + b*d)) + f*(S(2)*a*c*e - b**S(2)*e + b*c*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*e*x*(-f*(-b*e + S(2)*c*d) + g*(-S(2)*a*e + b*d))*(m + S(2)*p + S(4)) + f*(-S(2)*a*c*e**S(2)*(m + S(2)*p + S(3)) + b**S(2)*e**S(2)*(m + p + S(2)) + b*c*d*e*(-m + S(2)*p + S(2)) - S(2)*c**S(2)*d**S(2)*(S(2)*p + S(3))) - g*(a*e*(b*e*m + b*e - S(2)*c*d*m) - b*d*(-b*e*p - b*e + S(2)*c*d*p + S(3)*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule256)

    pattern257 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons87, cons88, cons542)
    rule257 = ReplacementRule(pattern257, lambda g, f, p, e, a, m, x, d, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-a*c*d*g + a*c*e*f + c*x*(a*e*g + c*d*f))/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x)**m*Simp(-a*c*d*e*g*m + c*e*x*(a*e*g + c*d*f)*(m + S(2)*p + S(4)) + f*(a*c*e**S(2)*(m + S(2)*p + S(3)) + c**S(2)*d**S(2)*(S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule257)

    pattern258 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons71)
    rule258 = ReplacementRule(pattern258, lambda g, b, f, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule258)

    pattern259 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons71)
    rule259 = ReplacementRule(pattern259, lambda g, f, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule259)

    pattern260 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons266, cons121)
    rule260 = ReplacementRule(pattern260, lambda g, b, f, e, a, m, x, d, c : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule260)

    pattern261 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons266, cons121)
    rule261 = ReplacementRule(pattern261, lambda g, f, e, a, m, x, d, c : g*(d + e*x)**m/(c*m) + Int((d + e*x)**(m + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule261)

    pattern262 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461)
    rule262 = ReplacementRule(pattern262, lambda g, b, f, e, a, x, d, c : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) - b*d*e + c*d**S(2) + c*x**S(4) - x**S(2)*(-b*e + S(2)*c*d)), x), x, sqrt(d + e*x)))
    rubi.add(rule262)

    pattern263 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462)
    rule263 = ReplacementRule(pattern263, lambda g, f, e, a, x, d, c : S(2)*Subst(Int((-d*g + e*f + g*x**S(2))/(a*e**S(2) + c*d**S(2) - S(2)*c*d*x**S(2) + c*x**S(4)), x), x, sqrt(d + e*x)))
    rubi.add(rule263)

    pattern264 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons407, cons461, cons266, cons38)
    rule264 = ReplacementRule(pattern264, lambda g, b, f, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule264)

    pattern265 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons497, cons462, cons266, cons38)
    rule265 = ReplacementRule(pattern265, lambda g, f, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((d + e*x)**(m + S(1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule265)

    pattern266 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons271)
    rule266 = ReplacementRule(pattern266, lambda g, b, f, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule266)

    pattern267 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons271)
    rule267 = ReplacementRule(pattern267, lambda g, f, e, a, m, x, d, c : Int(ExpandIntegrand((d + e*x)**m, (f + g*x)/(a + c*x**S(2)), x), x))
    rubi.add(rule267)

    pattern268 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons407, cons461, cons51, cons121, cons501, cons547, cons542)
    rule268 = ReplacementRule(pattern268, lambda g, b, f, p, e, a, m, x, d, c : g*(d + e*x)**m*(a + b*x + c*x**S(2))**(p + S(1))/(c*(m + S(2)*p + S(2))) + Int((d + e*x)**(m + S(-1))*(a + b*x + c*x**S(2))**p*Simp(d*(p + S(1))*(-b*g + S(2)*c*f) + m*(-a*e*g + c*d*f) + x*(e*(p + S(1))*(-b*g + S(2)*c*f) + m*(-b*e*g + c*d*g + c*e*f)), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule268)

    pattern269 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons462, cons51, cons121, cons501, cons547, cons542)
    rule269 = ReplacementRule(pattern269, lambda g, f, p, e, a, m, x, d, c : g*(a + c*x**S(2))**(p + S(1))*(d + e*x)**m/(c*(m + S(2)*p + S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(-1))*Simp(-a*e*g*m + c*d*f*(m + S(2)*p + S(2)) + c*x*(d*g*m + e*f*(m + S(2)*p + S(2))), x), x)/(c*(m + S(2)*p + S(2))))
    rubi.add(rule269)

    pattern270 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons407, cons461, cons51, cons38, cons542)
    rule270 = ReplacementRule(pattern270, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule270)

    pattern271 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons74, cons497, cons462, cons51, cons38, cons542)
    rule271 = ReplacementRule(pattern271, lambda g, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule271)

    pattern272 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons461, cons548, cons1)
    rule272 = ReplacementRule(pattern272, lambda g, b, f, p, e, a, m, x, d, c : (d + e*x)**(m + S(1))*(-d*g + e*f)*(a + b*x + c*x**S(2))**(p + S(1))/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))) + Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(b*(p + S(1))*(d*g - e*f) - c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g - b*e*f + c*d*f), x), x)/((m + S(1))*(a*e**S(2) - b*d*e + c*d**S(2))))
    rubi.add(rule272)

    pattern273 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons462, cons548, cons1)
    rule273 = ReplacementRule(pattern273, lambda g, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(d + e*x)**(m + S(1))*(-d*g + e*f)/((m + S(1))*(a*e**S(2) + c*d**S(2))) + Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*Simp(-c*x*(-d*g + e*f)*(m + S(2)*p + S(3)) + (m + S(1))*(a*e*g + c*d*f), x), x)/((m + S(1))*(a*e**S(2) + c*d**S(2))))
    rubi.add(rule273)

    pattern274 = Pattern(Integral((f_ + x_*WC('g', S(1)))/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons549, cons550, cons551)
    rule274 = ReplacementRule(pattern274, lambda g, b, f, e, a, x, d, c : S(4)*f*(a - d)*Subst(Int(S(1)/(S(4)*a - S(4)*d - x**S(2)), x), x, (S(2)*a - S(2)*d + x*(b - e))/sqrt(a + b*x + c*x**S(2)))/(-a*e + b*d))
    rubi.add(rule274)

    pattern275 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons73, cons161, cons407)
    rule275 = ReplacementRule(pattern275, lambda g, b, f, a, x, c : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + b*x**S(2) + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule275)

    pattern276 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons4, cons9, cons73, cons161, cons552)
    rule276 = ReplacementRule(pattern276, lambda g, f, a, x, c : S(2)*Subst(Int((f + g*x**S(2))/sqrt(a + c*x**S(4)), x), x, sqrt(x)))
    rubi.add(rule276)

    pattern277 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons72, cons73, cons161, cons407)
    rule277 = ReplacementRule(pattern277, lambda g, b, f, e, a, x, c : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + b*x + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule277)

    pattern278 = Pattern(Integral((f_ + x_*WC('g', S(1)))/(sqrt(e_*x_)*sqrt(a_ + x_**S(2)*WC('c', S(1)))), x_), cons4, cons9, cons72, cons73, cons161, cons553)
    rule278 = ReplacementRule(pattern278, lambda g, f, e, a, x, c : sqrt(x)*Int((f + g*x)/(sqrt(x)*sqrt(a + c*x**S(2))), x)/sqrt(e*x))
    rubi.add(rule278)

    pattern279 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons461)
    rule279 = ReplacementRule(pattern279, lambda g, b, f, p, e, a, m, x, d, c : g*Int((d + e*x)**(m + S(1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule279)

    pattern280 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons462)
    rule280 = ReplacementRule(pattern280, lambda g, f, p, e, a, m, x, d, c : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m, x)/e)
    rubi.add(rule280)

    pattern281 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons554)
    rule281 = ReplacementRule(pattern281, lambda g, b, f, p, e, a, m, n, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule281)

    pattern282 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons554)
    rule282 = ReplacementRule(pattern282, lambda g, f, p, e, a, m, n, x, d, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule282)

    pattern283 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons102, cons116)
    rule283 = ReplacementRule(pattern283, lambda g, b, f, p, e, a, x, d, c : (a*e**S(2) - b*d*e + c*d**S(2))*Int((a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + b*x + c*x**S(2))**(p + S(-1))*Simp(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule283)

    pattern284 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_/((x_*WC('e', S(1)) + WC('d', S(0)))*(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons102, cons116)
    rule284 = ReplacementRule(pattern284, lambda g, f, p, e, a, x, d, c : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*Simp(a*e*g + c*d*f - c*x*(-d*g + e*f), x)/(f + g*x), x)/(e*(-d*g + e*f)))
    rubi.add(rule284)

    pattern285 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons555, cons266, )
    def With285(g, b, f, p, e, a, m, n, x, d, c):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(c*x**(S(2)*q)/e**S(2) - x**q*(-b*e + S(2)*c*d)/e**S(2) + (a*e**S(2) - b*d*e + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule285 = ReplacementRule(pattern285, lambda g, b, f, p, e, a, m, n, x, d, c : With285(g, b, f, p, e, a, m, n, x, d, c))
    rubi.add(rule285)

    pattern286 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons555, cons266, )
    def With286(g, f, p, e, a, m, n, x, d, c):
        q = Denominator(m)
        return q*Subst(Int(x**(q*(m + S(1)) + S(-1))*(g*x**q/e + (-d*g + e*f)/e)**n*(-S(2)*c*d*x**q/e**S(2) + c*x**(S(2)*q)/e**S(2) + (a*e**S(2) + c*d**S(2))/e**S(2))**p, x), x, (d + e*x)**(S(1)/q))/e
    rule286 = ReplacementRule(pattern286, lambda g, f, p, e, a, m, n, x, d, c : With286(g, f, p, e, a, m, n, x, d, c))
    rubi.add(rule286)

    pattern287 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons70, cons518, cons556)
    rule287 = ReplacementRule(pattern287, lambda g, b, d, p, f, e, a, m, n, x, c : Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule287)

    pattern288 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons70, cons518, cons556)
    rule288 = ReplacementRule(pattern288, lambda g, d, p, f, e, a, m, n, x, c : Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule288)

    pattern289 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons70, cons518)
    rule289 = ReplacementRule(pattern289, lambda g, b, d, p, f, e, a, m, n, x, c : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule289)

    pattern290 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons70, cons518)
    rule290 = ReplacementRule(pattern290, lambda g, f, p, e, a, m, n, x, d, c : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((a + c*x**S(2))**p*(d*f + e*g*x**S(2))**m, x))
    rubi.add(rule290)

    pattern291 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons37)
    rule291 = ReplacementRule(pattern291, lambda g, b, f, e, a, m, n, x, d, c : c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x) + Int((a + b*x)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule291)

    pattern292 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons37)
    rule292 = ReplacementRule(pattern292, lambda g, f, e, a, m, n, x, d, c : a*Int((d + e*x)**m*(f + g*x)**n, x) + c*Int(x**S(2)*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule292)

    pattern293 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons60, cons36, cons37, cons121, cons118)
    rule293 = ReplacementRule(pattern293, lambda g, b, f, e, a, m, n, x, d, c : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-b*e*g + c*d*g + S(2)*c*e*f + c*e*g*x, x), x)/c**S(2) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(a*b*e*g**S(2) - a*c*d*g**S(2) - S(2)*a*c*e*f*g + c**S(2)*d*f**S(2) + x*(-a*c*e*g**S(2) + b**S(2)*e*g**S(2) - b*c*d*g**S(2) - S(2)*b*c*e*f*g + S(2)*c**S(2)*d*f*g + c**S(2)*e*f**S(2)), x)/(a + b*x + c*x**S(2)), x)/c**S(2))
    rubi.add(rule293)

    pattern294 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons60, cons36, cons37, cons121, cons118)
    rule294 = ReplacementRule(pattern294, lambda g, f, e, a, m, n, x, d, c : g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(d*g + S(2)*e*f + e*g*x, x), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-2))*Simp(-a*d*g**S(2) - S(2)*a*e*f*g + c*d*f**S(2) + x*(-a*e*g**S(2) + S(2)*c*d*f*g + c*e*f**S(2)), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule294)

    pattern295 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons60, cons36, cons37, cons121, cons31)
    rule295 = ReplacementRule(pattern295, lambda g, b, f, e, a, m, n, x, d, c : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(-b*e*g + c*d*g + c*e*f), x)/(a + b*x + c*x**S(2)), x)/c)
    rubi.add(rule295)

    pattern296 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons60, cons36, cons37, cons121, cons31)
    rule296 = ReplacementRule(pattern296, lambda g, f, e, a, m, n, x, d, c : e*g*Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1)), x)/c + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(-1))*Simp(-a*e*g + c*d*f + x*(c*d*g + c*e*f), x)/(a + c*x**S(2)), x)/c)
    rubi.add(rule296)

    pattern297 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons60, cons36, cons37, cons121, cons32)
    rule297 = ReplacementRule(pattern297, lambda g, b, f, e, a, m, n, x, d, c : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) - b*f*g + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g - b*d*g + c*d*f + c*x*(-d*g + e*f), x)/(a + b*x + c*x**S(2)), x)/(a*g**S(2) - b*f*g + c*f**S(2)))
    rubi.add(rule297)

    pattern298 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons60, cons36, cons37, cons121, cons32)
    rule298 = ReplacementRule(pattern298, lambda g, f, e, a, m, n, x, d, c : -g*(-d*g + e*f)*Int((d + e*x)**(m + S(-1))*(f + g*x)**n, x)/(a*g**S(2) + c*f**S(2)) + Int((d + e*x)**(m + S(-1))*(f + g*x)**(n + S(1))*Simp(a*e*g + c*d*f + c*x*(-d*g + e*f), x)/(a + c*x**S(2)), x)/(a*g**S(2) + c*f**S(2)))
    rubi.add(rule298)

    pattern299 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons407, cons461, cons14)
    rule299 = ReplacementRule(pattern299, lambda g, b, f, e, a, m, x, d, c : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule299)

    pattern300 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_/(sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons462, cons14)
    rule300 = ReplacementRule(pattern300, lambda g, f, e, a, m, x, d, c : Int(ExpandIntegrand(S(1)/(sqrt(d + e*x)*sqrt(f + g*x)), (d + e*x)**(m + S(1)/2)/(a + c*x**S(2)), x), x))
    rubi.add(rule300)

    pattern301 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons407, cons461, cons60, cons36)
    rule301 = ReplacementRule(pattern301, lambda g, b, f, e, a, m, n, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + b*x + c*x**S(2)), x), x))
    rubi.add(rule301)

    pattern302 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(a_ + x_**S(2)*WC('c', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons462, cons60, cons36)
    rule302 = ReplacementRule(pattern302, lambda g, f, e, a, m, n, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n, S(1)/(a + c*x**S(2)), x), x))
    rubi.add(rule302)

    pattern303 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons557)
    rule303 = ReplacementRule(pattern303, lambda g, b, f, p, e, a, m, n, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule303)

    pattern304 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons557)
    rule304 = ReplacementRule(pattern304, lambda g, f, p, e, a, m, n, x, d, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x))
    rubi.add(rule304)

    pattern305 = Pattern(Integral((x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons72, cons161, cons2, cons13, cons74, cons538, cons470, cons471)
    rule305 = ReplacementRule(pattern305, lambda g, b, p, e, a, m, n, x, d, c : (d + e*x)**FracPart(p)*(a*d + c*e*x**S(3))**(-FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((g*x)**n*(a*d + c*e*x**S(3))**p, x))
    rubi.add(rule305)

    pattern306 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons36, cons100, cons520, cons116, cons32)
    rule306 = ReplacementRule(pattern306, lambda g, b, f, p, e, a, n, x, d, c : (a*e**S(2) - b*d*e + c*d**S(2))*Int((f + g*x)**(n + S(1))*(a + b*x + c*x**S(2))**(p + S(-1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((f + g*x)**n*(a + b*x + c*x**S(2))**(p + S(-1))*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule306)

    pattern307 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons36, cons100, cons520, cons116, cons32)
    rule307 = ReplacementRule(pattern307, lambda g, f, p, e, a, n, x, d, c : (a*e**S(2) + c*d**S(2))*Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**(n + S(1))/(d + e*x), x)/(e*(-d*g + e*f)) - Int((a + c*x**S(2))**(p + S(-1))*(f + g*x)**n*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(e*(-d*g + e*f)))
    rubi.add(rule307)

    pattern308 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons36, cons100, cons520, cons88, cons31)
    rule308 = ReplacementRule(pattern308, lambda g, b, f, p, e, a, n, x, d, c : e*(-d*g + e*f)*Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))/(d + e*x), x)/(a*e**S(2) - b*d*e + c*d**S(2)) + Int((f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p*(a*e*g - b*e*f + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) - b*d*e + c*d**S(2)))
    rubi.add(rule308)

    pattern309 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_/(x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons36, cons100, cons520, cons88, cons31)
    rule309 = ReplacementRule(pattern309, lambda g, f, p, e, a, n, x, d, c : e*(-d*g + e*f)*Int((a + c*x**S(2))**(p + S(1))*(f + g*x)**(n + S(-1))/(d + e*x), x)/(a*e**S(2) + c*d**S(2)) + Int((a + c*x**S(2))**p*(f + g*x)**(n + S(-1))*(a*e*g + c*d*f - c*x*(-d*g + e*f)), x)/(a*e**S(2) + c*d**S(2)))
    rubi.add(rule309)

    pattern310 = Pattern(Integral(S(1)/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, )
    def With310(g, b, f, e, a, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return -sqrt(S(2))*sqrt(-g*(b + S(2)*c*x - q)/(-b*g + S(2)*c*f + g*q))*sqrt(-g*(b + S(2)*c*x + q)/(-b*g + S(2)*c*f - g*q))*EllipticPi(e*(-b*g + S(2)*c*f + g*q)/(S(2)*c*(-d*g + e*f)), asin(sqrt(S(2))*sqrt(c/(-b*g + S(2)*c*f + g*q))*sqrt(f + g*x)), (-b*g + S(2)*c*f + g*q)/(-b*g + S(2)*c*f - g*q))/(sqrt(c/(-b*g + S(2)*c*f + g*q))*(-d*g + e*f)*sqrt(a + b*x + c*x**S(2)))
    rule310 = ReplacementRule(pattern310, lambda g, b, f, e, a, x, d, c : With310(g, b, f, e, a, x, d, c))
    rubi.add(rule310)

    pattern311 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, )
    def With311(g, f, e, a, x, d, c):
        q = Rt(-a*c, S(2))
        return -S(2)*sqrt(g*(-c*x + q)/(c*f + g*q))*sqrt(-g*(c*x + q)/(c*f - g*q))*EllipticPi(e*(c*f + g*q)/(c*(-d*g + e*f)), asin(sqrt(c/(c*f + g*q))*sqrt(f + g*x)), (c*f + g*q)/(c*f - g*q))/(sqrt(c/(c*f + g*q))*sqrt(a + c*x**S(2))*(-d*g + e*f))
    rule311 = ReplacementRule(pattern311, lambda g, f, e, a, x, d, c : With311(g, f, e, a, x, d, c))
    rubi.add(rule311)

    pattern312 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/((x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461, cons22)
    rule312 = ReplacementRule(pattern312, lambda g, b, f, e, a, n, x, d, c : Int(ExpandIntegrand(S(1)/(sqrt(f + g*x)*sqrt(a + b*x + c*x**S(2))), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule312)

    pattern313 = Pattern(Integral((x_*WC('g', S(1)) + WC('f', S(0)))**n_/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*(x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462, cons22)
    rule313 = ReplacementRule(pattern313, lambda g, f, e, a, n, x, d, c : Int(ExpandIntegrand(S(1)/(sqrt(a + c*x**S(2))*sqrt(f + g*x)), (f + g*x)**(n + S(1)/2)/(d + e*x), x), x))
    rubi.add(rule313)

    pattern314 = Pattern(Integral(S(1)/(sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))*sqrt(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons497, cons407, cons461)
    rule314 = ReplacementRule(pattern314, lambda g, b, f, e, a, x, d, c : -S(2)*sqrt((-d*g + e*f)**S(2)*(a + b*x + c*x**S(2))/((d + e*x)**S(2)*(a*g**S(2) - b*f*g + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) - b*d*e + c*d**S(2))/(a*g**S(2) - b*f*g + c*f**S(2)) - x**S(2)*(S(2)*a*e*g - b*d*g - b*e*f + S(2)*c*d*f)/(a*g**S(2) - b*f*g + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/((-d*g + e*f)*sqrt(a + b*x + c*x**S(2))))
    rubi.add(rule314)

    pattern315 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_*WC('e', S(1)) + WC('d', S(0)))*sqrt(x_*WC('g', S(1)) + WC('f', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons497, cons462)
    rule315 = ReplacementRule(pattern315, lambda g, f, e, a, x, d, c : -S(2)*sqrt((a + c*x**S(2))*(-d*g + e*f)**S(2)/((d + e*x)**S(2)*(a*g**S(2) + c*f**S(2))))*(d + e*x)*Subst(Int(S(1)/sqrt(x**S(4)*(a*e**S(2) + c*d**S(2))/(a*g**S(2) + c*f**S(2)) - x**S(2)*(S(2)*a*e*g + S(2)*c*d*f)/(a*g**S(2) + c*f**S(2)) + S(1)), x), x, sqrt(f + g*x)/sqrt(d + e*x))/(sqrt(a + c*x**S(2))*(-d*g + e*f)))
    rubi.add(rule315)

    pattern316 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(2), x_), cons4, cons9, cons72, cons73, cons161, cons2, cons74, cons558)
    rule316 = ReplacementRule(pattern316, lambda g, p, f, e, a, m, x, c : Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + g**S(2)*x**S(2)), x) + S(2)*f*g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p, x)/e)
    rubi.add(rule316)

    pattern317 = Pattern(Integral((x_*WC('e', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(f_ + x_*WC('g', S(1)))**S(3), x_), cons4, cons9, cons72, cons73, cons161, cons2, cons74, cons558)
    rule317 = ReplacementRule(pattern317, lambda g, p, f, e, a, m, x, c : f*Int((e*x)**m*(a + c*x**S(2))**p*(f**S(2) + S(3)*g**S(2)*x**S(2)), x) + g*Int((e*x)**(m + S(1))*(a + c*x**S(2))**p*(S(3)*f**S(2) + g**S(2)*x**S(2)), x)/e)
    rubi.add(rule317)

    pattern318 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons407, cons461, cons101)
    rule318 = ReplacementRule(pattern318, lambda g, b, f, p, e, a, m, n, x, d, c : g*Int((d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e + (-d*g + e*f)*Int((d + e*x)**m*(f + g*x)**(n + S(-1))*(a + b*x + c*x**S(2))**p, x)/e)
    rubi.add(rule318)

    pattern319 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**m_*(x_*WC('g', S(1)) + WC('f', S(0)))**n_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons74, cons497, cons462, cons101)
    rule319 = ReplacementRule(pattern319, lambda g, f, p, e, a, m, n, x, d, c : g*Int((a + c*x**S(2))**p*(d + e*x)**(m + S(1))*(f + g*x)**(n + S(-1)), x)/e + (-d*g + e*f)*Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**(n + S(-1)), x)/e)
    rubi.add(rule319)

    pattern320 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons402)
    rule320 = ReplacementRule(pattern320, lambda g, b, f, p, e, a, m, n, x, d, c : Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule320)

    pattern321 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons559)
    rule321 = ReplacementRule(pattern321, lambda g, f, p, e, a, m, n, x, d, c : Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x))
    rubi.add(rule321)

    pattern322 = Pattern(Integral((u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(a_ + u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons6, cons7)
    rule322 = ReplacementRule(pattern322, lambda g, b, u, f, p, e, a, m, n, x, d, c : Subst(Int((d + e*x)**m*(f + g*x)**n*(a + b*x + c*x**S(2))**p, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule322)

    pattern323 = Pattern(Integral((a_ + u_**S(2)*WC('c', S(1)))**WC('p', S(1))*(u_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(u_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons2, cons13, cons74, cons6, cons7)
    rule323 = ReplacementRule(pattern323, lambda g, u, f, p, e, a, m, n, x, d, c : Subst(Int((a + c*x**S(2))**p*(d + e*x)**m*(f + g*x)**n, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule323)

    pattern324 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons173, cons560, cons561, cons562, cons563)
    rule324 = ReplacementRule(pattern324, lambda b, d, f, p, e, a, q, x, c : (c/f)**p*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule324)

    pattern325 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons173, cons560, cons561, cons100, cons324, cons564)
    rule325 = ReplacementRule(pattern325, lambda b, d, f, p, e, a, q, x, c : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule325)

    pattern326 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons173, cons405, cons100)
    rule326 = ReplacementRule(pattern326, lambda b, d, f, p, e, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule326)

    pattern327 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons74, cons173, cons405, cons100)
    rule327 = ReplacementRule(pattern327, lambda b, d, f, p, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule327)

    pattern328 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons173, cons565, cons566, cons567)
    rule328 = ReplacementRule(pattern328, lambda b, d, f, e, a, q, x, c : (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule328)

    pattern329 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons173, cons568, cons566, cons567)
    rule329 = ReplacementRule(pattern329, lambda d, f, e, a, q, x, c : (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule329)

    pattern330 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons73, cons173, cons566, cons569)
    rule330 = ReplacementRule(pattern330, lambda b, d, f, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(S(2)*a*f*x*(q + S(1)) + b*d)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule330)

    pattern331 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons73, cons173, cons570)
    rule331 = ReplacementRule(pattern331, lambda b, d, f, a, q, x, c : b*Int(x*(d + f*x**S(2))**q, x) + Int((a + c*x**S(2))*(d + f*x**S(2))**q, x))
    rubi.add(rule331)

    pattern332 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons571, cons570)
    rule332 = ReplacementRule(pattern332, lambda b, d, f, e, a, q, x, c : Int(ExpandIntegrand((a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule332)

    pattern333 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons571, cons570)
    rule333 = ReplacementRule(pattern333, lambda d, f, e, a, q, x, c : Int(ExpandIntegrand((a + c*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule333)

    pattern334 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons571, cons289, cons380, cons572)
    rule334 = ReplacementRule(pattern334, lambda b, d, f, e, a, q, x, c : -(c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f - S(2)*b*d*f + c*d*e + x*(c*(-S(2)*d*f + e**S(2)) + f*(S(2)*a*f - b*e)))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule334)

    pattern335 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons571, cons289, cons380, cons573)
    rule335 = ReplacementRule(pattern335, lambda d, f, e, a, q, x, c : -(S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**(q + S(1)), x)/(f*(q + S(1))*(-S(4)*d*f + e**S(2))) + (d + e*x + f*x**S(2))**(q + S(1))*(a*e*f + c*d*e + x*(S(2)*a*f**S(2) + c*(-S(2)*d*f + e**S(2))))/(f*(q + S(1))*(-S(4)*d*f + e**S(2))))
    rubi.add(rule335)

    pattern336 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons4, cons5, cons9, cons10, cons73, cons289, cons380, cons574)
    rule336 = ReplacementRule(pattern336, lambda b, d, f, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(b*d - x*(a*f - c*d))/(S(2)*d*f*(q + S(1))) + (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**(q + S(1)), x)/(S(2)*d*f*(q + S(1))))
    rubi.add(rule336)

    pattern337 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons10, cons72, cons73, cons4, cons5, cons9, cons173, cons571, cons575, cons576, cons572)
    rule337 = ReplacementRule(pattern337, lambda b, d, f, e, a, q, x, c : (c*(-S(2)*d*f + e**S(2)*(q + S(2))) + f*(S(2)*q + S(3))*(S(2)*a*f - b*e))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) - c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule337)

    pattern338 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons10, cons72, cons73, cons4, cons9, cons173, cons571, cons575, cons576, cons573)
    rule338 = ReplacementRule(pattern338, lambda d, f, e, a, q, x, c : (S(2)*a*f**S(2)*(S(2)*q + S(3)) + c*(-S(2)*d*f + e**S(2)*(q + S(2))))*Int((d + e*x + f*x**S(2))**q, x)/(S(2)*f**S(2)*(S(2)*q + S(3))) + (-c*e*(q + S(2)) + S(2)*c*f*x*(q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule338)

    pattern339 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**q_*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1))), x_), cons10, cons73, cons4, cons5, cons9, cons173, cons575, cons576, cons574)
    rule339 = ReplacementRule(pattern339, lambda b, d, f, a, q, x, c : (S(2)*a*f*q + S(3)*a*f - c*d)*Int((d + f*x**S(2))**q, x)/(f*(S(2)*q + S(3))) + (d + f*x**S(2))**(q + S(1))*(b*f*(S(2)*q + S(3)) + S(2)*c*f*x*(q + S(1)))/(S(2)*f**S(2)*(q + S(1))*(S(2)*q + S(3))))
    rubi.add(rule339)

    pattern340 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, cons308, cons88, cons290)
    rule340 = ReplacementRule(pattern340, lambda b, f, p, e, a, q, x, d, c : (b + S(2)*c*x)*(a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(b*e*q + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(S(2)*b*f*q + S(2)*c*e*(S(2)*p + q + S(3))), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule340)

    pattern341 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons407, cons308, cons88, cons290)
    rule341 = ReplacementRule(pattern341, lambda b, f, p, a, q, x, d, c : (b + S(2)*c*x)*(d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*b*f*q*x + S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule341)

    pattern342 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons571, cons308, cons88, cons290)
    rule342 = ReplacementRule(pattern342, lambda f, p, e, a, q, x, d, c : -x*(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q/(S(2)*a*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(S(2)*c*d*(S(2)*p + S(3)) + S(2)*c*e*x*(S(2)*p + q + S(3)) + S(2)*c*f*x**S(2)*(S(2)*p + S(2)*q + S(3)), x), x)/(S(4)*a*c*(p + S(1))))
    rubi.add(rule342)

    pattern343 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons173, cons407, cons571, cons87, cons88, cons577, cons312)
    rule343 = ReplacementRule(pattern343, lambda b, f, p, e, a, q, x, d, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d) + c*x*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + S(2)*c*(p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) - e*(p + q + S(2))*(-S(2)*a*c**S(2)*e - b**S(3)*f + b**S(2)*c*e - b*c*(-S(3)*a*f + c*d)) + x*(S(2)*f*(p + q + S(2))*(S(2)*a*c**S(2)*e + b**S(3)*f - b**S(2)*c*e + b*c*(-S(3)*a*f + c*d)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule343)

    pattern344 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons173, cons407, cons87, cons88, cons578, cons312)
    rule344 = ReplacementRule(pattern344, lambda b, f, p, a, q, x, d, c : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(b**S(3)*f + b*c*(-S(3)*a*f + c*d) + c*x*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*c*(p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + x*(-b*f*(p + S(1))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + S(2)*f*(b**S(3)*f + b*c*(-S(3)*a*f + c*d))*(p + q + S(2))) - (a*f*(p + S(1)) - c*d*(p + S(2)))*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule344)

    pattern345 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons173, cons571, cons87, cons88, cons579, cons312)
    rule345 = ReplacementRule(pattern345, lambda f, p, e, a, q, x, d, c : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c**S(2)*e + c*x*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(S(2)*a*c**S(2)*e**S(2)*(p + q + S(2)) + c*f*x**S(2)*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + S(2)*q + S(5)) + S(2)*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + x*(S(4)*a*c**S(2)*e*f*(p + q + S(2)) + c*e*(-S(2)*a*c*f + S(2)*c**S(2)*d)*(S(2)*p + q + S(4))) - (-S(2)*a*c*f + S(2)*c**S(2)*d)*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule345)

    pattern346 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons173, cons407, cons571, cons87, cons99, cons580, cons581)
    rule346 = ReplacementRule(pattern346, lambda b, f, p, e, a, q, x, d, c : (a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**(q + S(1))*(b*f*(S(3)*p + S(2)*q) - c*e*(S(2)*p + q) + S(2)*c*f*x*(p + q))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(x**S(2)*(c*(p + q)*(-c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))) + f*(-S(2)*a*f + b*e)*(S(4)*p + S(2)*q + S(-1))) + p*(-p + S(1))*(-b*f + c*e)**S(2)) + x*(S(2)*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)*(-b*f + c*e) - (p + q)*(b*(c*(S(2)*p + q)*(-S(4)*d*f + e**S(2)) + f*(S(2)*p + S(2)*q + S(1))*(S(2)*a*f - b*e + S(2)*c*d)) + e*f*(-p + S(1))*(-S(4)*a*c + b**S(2)))) + (-p + S(1))*(S(2)*p + q)*(-a*e + b*d)*(-b*f + c*e) - (p + q)*(-a*(c*(S(2)*d*f - e**S(2)*(S(2)*p + q)) + f*(-S(2)*a*f + b*e)*(S(2)*p + S(2)*q + S(1))) + b**S(2)*d*f*(-p + S(1))), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule346)

    pattern347 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons173, cons407, cons87, cons99, cons580, cons581)
    rule347 = ReplacementRule(pattern347, lambda b, f, p, a, q, x, d, c : (d + f*x**S(2))**(q + S(1))*(b*(S(3)*p + S(2)*q) + S(2)*c*x*(p + q))*(a + b*x + c*x**S(2))**(p + S(-1))/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-2))*Simp(b**S(2)*d*(p + S(-1))*(S(2)*p + q) + x**S(2)*(b**S(2)*f*p*(-p + S(1)) + S(2)*c*(p + q)*(-a*f*(S(4)*p + S(2)*q + S(-1)) + c*d*(S(2)*p + S(-1)))) - x*(S(2)*b*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d) - S(2)*b*(p + q)*(S(2)*c*d*(S(2)*p + q) - (a*f + c*d)*(S(2)*p + S(2)*q + S(1)))) - (p + q)*(-S(2)*a*(-a*f*(S(2)*p + S(2)*q + S(1)) + c*d) + b**S(2)*d*(-p + S(1))), x), x)/(S(2)*f*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule347)

    pattern348 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons173, cons571, cons87, cons99, cons580, cons581)
    rule348 = ReplacementRule(pattern348, lambda f, p, e, a, q, x, d, c : -c*(a + c*x**S(2))**(p + S(-1))*(e*(S(2)*p + q) - S(2)*f*x*(p + q))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))) - Int((a + c*x**S(2))**(p + S(-2))*(d + e*x + f*x**S(2))**q*Simp(-a*c*e**S(2)*(-p + S(1))*(S(2)*p + q) + a*(p + q)*(-S(2)*a*f**S(2)*(S(2)*p + S(2)*q + S(1)) + c*(S(2)*d*f - e**S(2)*(S(2)*p + q))) + x**S(2)*(c**S(2)*e**S(2)*p*(-p + S(1)) - c*(p + q)*(S(2)*a*f**S(2)*(S(4)*p + S(2)*q + S(-1)) + c*(S(2)*d*f*(-S(2)*p + S(1)) + e**S(2)*(S(3)*p + q + S(-1))))) + x*(S(4)*a*c*e*f*(-p + S(1))*(p + q) + S(2)*c*e*(-p + S(1))*(S(2)*p + q)*(-a*f + c*d)), x), x)/(S(2)*f**S(2)*(p + q)*(S(2)*p + S(2)*q + S(1))))
    rubi.add(rule348)

    pattern349 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, )
    def With349(b, d, f, e, a, x, c):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-a*c*f + b**S(2)*f - b*c*e + c**S(2)*d - x*(-b*c*f + c**S(2)*e))/(a + b*x + c*x**S(2)), x)/q + Int((a*f**S(2) - b*e*f - c*d*f + c*e**S(2) + x*(-b*f**S(2) + c*e*f))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule349 = ReplacementRule(pattern349, lambda b, d, f, e, a, x, c : With349(b, d, f, e, a, x, c))
    rubi.add(rule349)

    pattern350 = Pattern(Integral(S(1)/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons407, )
    def With350(b, d, f, a, x, c):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return -Int((-a*f**S(2) + b*f**S(2)*x + c*d*f)/(d + f*x**S(2)), x)/q + Int((-a*c*f + b**S(2)*f + b*c*f*x + c**S(2)*d)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule350 = ReplacementRule(pattern350, lambda b, d, f, a, x, c : With350(b, d, f, a, x, c))
    rubi.add(rule350)

    pattern351 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, cons582)
    rule351 = ReplacementRule(pattern351, lambda b, d, f, e, a, x, c : -S(2)*e*Subst(Int(S(1)/(e*(-S(4)*a*f + b*e) - x**S(2)*(-a*e + b*d)), x), x, (e + S(2)*f*x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule351)

    pattern352 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, cons583, cons412, )
    def With352(b, d, f, e, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - S(2)*c*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule352 = ReplacementRule(pattern352, lambda b, d, f, e, a, x, c : With352(b, d, f, e, a, x, c))
    rubi.add(rule352)

    pattern353 = Pattern(Integral(S(1)/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons571, cons584)
    rule353 = ReplacementRule(pattern353, lambda d, f, e, a, x, c : Int(S(1)/((a - x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2) + Int(S(1)/((a + x*Rt(-a*c, S(2)))*sqrt(d + e*x + f*x**S(2))), x)/S(2))
    rubi.add(rule353)

    pattern354 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons407, cons412, )
    def With354(b, d, f, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - S(2)*c*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule354 = ReplacementRule(pattern354, lambda b, d, f, a, x, c : With354(b, d, f, a, x, c))
    rubi.add(rule354)

    pattern355 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, cons583, cons585, )
    def With355(b, f, e, a, x, d, c):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d - q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + q + x*(-b*f + c*e))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule355 = ReplacementRule(pattern355, lambda b, f, e, a, x, d, c : With355(b, f, e, a, x, d, c))
    rubi.add(rule355)

    pattern356 = Pattern(Integral(S(1)/((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons571, cons586, )
    def With356(f, e, a, x, d, c):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f + c*d + c*e*x - q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) + Int((-a*f + c*d + c*e*x + q)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule356 = ReplacementRule(pattern356, lambda f, e, a, x, d, c : With356(f, e, a, x, d, c))
    rubi.add(rule356)

    pattern357 = Pattern(Integral(S(1)/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons73, cons407, cons585, )
    def With357(b, f, a, x, d, c):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return -Int((-a*f - b*f*x + c*d - q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) + Int((-a*f - b*f*x + c*d + q)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule357 = ReplacementRule(pattern357, lambda b, f, a, x, d, c : With357(b, f, a, x, d, c))
    rubi.add(rule357)

    pattern358 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571)
    rule358 = ReplacementRule(pattern358, lambda b, d, f, e, a, x, c : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f + c*d + x*(-b*f + c*e))/(sqrt(a + b*x + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule358)

    pattern359 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))/(d_ + x_**S(2)*WC('f', S(1))), x_), cons4, cons5, cons9, cons10, cons73, cons407)
    rule359 = ReplacementRule(pattern359, lambda b, d, f, a, x, c : c*Int(S(1)/sqrt(a + b*x + c*x**S(2)), x)/f - Int((-a*f - b*f*x + c*d)/((d + f*x**S(2))*sqrt(a + b*x + c*x**S(2))), x)/f)
    rubi.add(rule359)

    pattern360 = Pattern(Integral(sqrt(a_ + x_**S(2)*WC('c', S(1)))/(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1))), x_), cons4, cons9, cons10, cons72, cons73, cons571)
    rule360 = ReplacementRule(pattern360, lambda d, f, e, a, x, c : c*Int(S(1)/sqrt(a + c*x**S(2)), x)/f - Int((-a*f + c*d + c*e*x)/(sqrt(a + c*x**S(2))*(d + e*x + f*x**S(2))), x)/f)
    rubi.add(rule360)

    pattern361 = Pattern(Integral(S(1)/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, )
    def With361(b, d, f, e, a, x, c):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*sqrt(d + e*x + f*x**S(2))), x)/sqrt(a + b*x + c*x**S(2))
    rule361 = ReplacementRule(pattern361, lambda b, d, f, e, a, x, c : With361(b, d, f, e, a, x, c))
    rubi.add(rule361)

    pattern362 = Pattern(Integral(S(1)/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons407, )
    def With362(b, d, f, a, x, c):
        r = Rt(-S(4)*a*c + b**S(2), S(2))
        return sqrt(S(2)*a + x*(b + r))*sqrt(b + S(2)*c*x + r)*Int(S(1)/(sqrt(S(2)*a + x*(b + r))*sqrt(d + f*x**S(2))*sqrt(b + S(2)*c*x + r)), x)/sqrt(a + b*x + c*x**S(2))
    rule362 = ReplacementRule(pattern362, lambda b, d, f, a, x, c : With362(b, d, f, a, x, c))
    rubi.add(rule362)

    pattern363 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons173, cons587)
    rule363 = ReplacementRule(pattern363, lambda b, f, p, e, a, q, x, d, c : Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule363)

    pattern364 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons74, cons173, cons588)
    rule364 = ReplacementRule(pattern364, lambda d, f, p, e, a, q, x, c : Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule364)

    pattern365 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons74, cons173, cons6, cons7)
    rule365 = ReplacementRule(pattern365, lambda b, u, p, f, e, a, q, x, d, c : Subst(Int((a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule365)

    pattern366 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons74, cons173, cons6, cons7)
    rule366 = ReplacementRule(pattern366, lambda u, p, f, e, a, q, x, d, c : Subst(Int((a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule366)

    pattern367 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons173, cons560, cons561, cons562, cons563)
    rule367 = ReplacementRule(pattern367, lambda g, b, d, h, p, f, e, a, m, q, x, c : (c/f)**p*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule367)

    pattern368 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons173, cons560, cons561, cons100, cons324, cons564)
    rule368 = ReplacementRule(pattern368, lambda g, b, d, h, p, f, e, a, m, q, x, c : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((g + h*x)**m*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule368)

    pattern369 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons173, cons405)
    rule369 = ReplacementRule(pattern369, lambda g, b, d, h, f, p, e, a, m, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule369)

    pattern370 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons173, cons405)
    rule370 = ReplacementRule(pattern370, lambda g, b, d, h, f, p, a, m, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(g + h*x)**m, x))
    rubi.add(rule370)

    pattern371 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons589, cons590, cons71)
    rule371 = ReplacementRule(pattern371, lambda g, b, d, h, p, f, e, a, m, x, c : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule371)

    pattern372 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons591, cons590, cons71)
    rule372 = ReplacementRule(pattern372, lambda g, d, h, p, f, e, a, m, x, c : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule372)

    pattern373 = Pattern(Integral((g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons74, cons589, cons592, cons71)
    rule373 = ReplacementRule(pattern373, lambda g, b, d, h, p, f, a, m, x, c : Int((f*h*x/c + d*g/a)**m*(a + b*x + c*x**S(2))**(m + p), x))
    rubi.add(rule373)

    pattern374 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons74, cons591, cons592, cons71)
    rule374 = ReplacementRule(pattern374, lambda g, d, h, p, f, a, m, x, c : Int((a + c*x**S(2))**(m + p)*(f*h*x/c + d*g/a)**m, x))
    rubi.add(rule374)

    pattern375 = Pattern(Integral(x_**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons72, cons73, cons173, cons407, cons593, cons97)
    rule375 = ReplacementRule(pattern375, lambda b, p, f, e, a, q, x, c : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule375)

    pattern376 = Pattern(Integral(x_**WC('p', S(1))*(a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons72, cons73, cons173, cons594, cons97)
    rule376 = ReplacementRule(pattern376, lambda f, p, e, a, q, x, c : Int((a/e + c*x/f)**p*(e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule376)

    pattern377 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons595, cons596, cons452)
    rule377 = ReplacementRule(pattern377, lambda g, b, d, h, p, f, e, a, m, x, c : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule377)

    pattern378 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons597, cons598, cons452)
    rule378 = ReplacementRule(pattern378, lambda g, d, h, p, f, e, a, m, x, c : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule378)

    pattern379 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons599, cons596, cons452)
    rule379 = ReplacementRule(pattern379, lambda g, b, d, h, p, f, a, m, x, c : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule379)

    pattern380 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons407, cons571, cons600)
    rule380 = ReplacementRule(pattern380, lambda g, b, d, h, p, f, e, a, m, x, c : Int(ExpandIntegrand((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule380)

    pattern381 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons571, cons600)
    rule381 = ReplacementRule(pattern381, lambda g, d, h, p, f, e, a, m, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2)), x), x))
    rubi.add(rule381)

    pattern382 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons407, cons600)
    rule382 = ReplacementRule(pattern382, lambda g, b, d, h, p, f, a, m, x, c : Int(ExpandIntegrand((d + f*x**S(2))*(g + h*x)**m*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule382)

    pattern383 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons2, cons600)
    rule383 = ReplacementRule(pattern383, lambda g, d, h, p, f, a, m, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(d + f*x**S(2))*(g + h*x)**m, x), x))
    rubi.add(rule383)

    pattern384 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons407, cons571, cons51, cons38, cons601)
    rule384 = ReplacementRule(pattern384, lambda g, b, d, h, p, f, e, a, m, x, c : (g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(f*g**S(2)*(p + S(1)) - h*(-d*h*(m + p + S(2)) + e*g*(p + S(1)))) + h*(m + S(1))*(a*e*h - a*f*g + c*d*g) - x*(c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule384)

    pattern385 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons571, cons51, cons38, cons602)
    rule385 = ReplacementRule(pattern385, lambda g, h, p, f, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) - e*g*h + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(a*e*h - a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(S(2)*f*g**S(2)*(p + S(1)) - h*(-d*h + e*g)*(m + S(2)*p + S(3)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule385)

    pattern386 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons74, cons407, cons51, cons38, cons601)
    rule386 = ReplacementRule(pattern386, lambda g, b, d, h, p, f, a, m, x, c : (g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))*(a + b*x + c*x**S(2))**(p + S(1))/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))) + Int((g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**p*Simp(-b*(d*h**S(2)*(m + p + S(2)) + f*g**S(2)*(p + S(1))) + h*(m + S(1))*(-a*f*g + c*d*g) - x*(c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1))) + f*h*(m + S(1))*(-a*h + b*g)), x), x)/(h*(m + S(1))*(a*h**S(2) - b*g*h + c*g**S(2))))
    rubi.add(rule386)

    pattern387 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons74, cons51, cons38, cons602)
    rule387 = ReplacementRule(pattern387, lambda g, d, h, f, p, a, m, x, c : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(d*h**S(2) + f*g**S(2))/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**p*(g + h*x)**(m + S(1))*Simp(h*(m + S(1))*(-a*f*g + c*d*g) + x*(a*f*h**S(2)*(m + S(1)) - c*(d*h**S(2)*(m + S(2)*p + S(3)) + S(2)*f*g**S(2)*(p + S(1)))), x), x)/(h*(m + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule387)

    pattern388 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons601)
    rule388 = ReplacementRule(pattern388, lambda g, b, d, h, f, e, a, x, c : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((a*e*h - a*f*g - b*d*h + c*d*g + x*(a*f*h - b*f*g - c*d*h + c*e*g))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule388)

    pattern389 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(x_*WC('h', S(1)) + WC('g', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons602)
    rule389 = ReplacementRule(pattern389, lambda g, h, f, e, a, x, d, c : (d*h**S(2) - e*g*h + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((a*e*h - a*f*g + c*d*g + x*(a*f*h - c*d*h + c*e*g))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule389)

    pattern390 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(3)/2)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons601)
    rule390 = ReplacementRule(pattern390, lambda g, b, d, h, f, a, x, c : (d*h**S(2) + f*g**S(2))*Int(S(1)/((g + h*x)*sqrt(a + b*x + c*x**S(2))), x)/(a*h**S(2) - b*g*h + c*g**S(2)) + Int((-a*f*g - b*d*h + c*d*g - x*(-a*f*h + b*f*g + c*d*h))/(a + b*x + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) - b*g*h + c*g**S(2)))
    rubi.add(rule390)

    pattern391 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(3)/2)*(g_ + x_*WC('h', S(1)))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons602)
    rule391 = ReplacementRule(pattern391, lambda g, d, h, f, a, x, c : (d*h**S(2) + f*g**S(2))*Int(S(1)/(sqrt(a + c*x**S(2))*(g + h*x)), x)/(a*h**S(2) + c*g**S(2)) + Int((-a*f*g + c*d*g - x*(-a*f*h + c*d*h))/(a + c*x**S(2))**(S(3)/2), x)/(a*h**S(2) + c*g**S(2)))
    rubi.add(rule391)

    pattern392 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons236, cons88, cons119)
    rule392 = ReplacementRule(pattern392, lambda g, b, d, h, f, p, e, a, m, x, c : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f - S(2)*a*c*e + b*c*d + x*(c*(-b*e + S(2)*c*d) + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(c*(S(2)*p + S(3))*(-b*e + S(2)*c*d) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f - S(2)*a*c*e + b*c*d) + h*x*(c*(-b*e + S(2)*c*d)*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule392)

    pattern393 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons236, cons88, cons119)
    rule393 = ReplacementRule(pattern393, lambda g, h, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(a*e - x*(-a*f + c*d))/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*(e*h*m + f*g) - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule393)

    pattern394 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons236, cons88, cons119)
    rule394 = ReplacementRule(pattern394, lambda g, b, d, h, f, p, a, m, x, c : (g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*(a*b*f + b*c*d + x*(S(2)*c**S(2)*d + f*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((g + h*x)**(m + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(S(2)*c**S(2)*d*(S(2)*p + S(3)) - f*(S(2)*a*c - b**S(2)*(p + S(2)))) + h*m*(a*b*f + b*c*d) + h*x*(S(2)*c**S(2)*d*(m + S(2)*p + S(3)) - f*(S(2)*a*c*(m + S(1)) - b**S(2)*(m + p + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule394)

    pattern395 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**m_*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons236, cons88, cons119)
    rule395 = ReplacementRule(pattern395, lambda g, d, h, f, p, a, m, x, c : -x*(a + c*x**S(2))**(p + S(1))*(g + h*x)**m*(-a*f + c*d)/(S(2)*a*c*(p + S(1))) - Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(-1))*Simp(a*f*g - c*d*g*(S(2)*p + S(3)) + h*x*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule395)

    pattern396 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons407, cons571, cons87, cons88, cons603)
    rule396 = ReplacementRule(pattern396, lambda g, b, d, h, f, p, e, a, m, x, c : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(-S(2)*a*e*h + S(2)*a*f*g + b*d*h + b*e*g)) - (-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*(-c*e*h + c*f*g) + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c**S(2)*d*g) - h*(-(-a*e + b*d)*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + (p + S(1))*(c*g**S(2) - h*(-a*h + b*g))*(S(2)*a*f - b*e + S(2)*c*d), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule396)

    pattern397 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons571, cons87, cons88, cons602)
    rule397 = ReplacementRule(pattern397, lambda g, h, f, p, e, a, m, x, d, c : (a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*c*e*g - a*h*(-a*f + c*d) - c*x*(a*e*h - a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(g*(p + S(2))*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g) + h*x*(-a*(-c*e*h + c*f*g) + c**S(2)*d*g)*(m + S(2)*p + S(4)) - h*(a*c*e*g - a*h*(-a*f + c*d))*(m + p + S(2)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule397)

    pattern398 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_, x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons407, cons87, cons88, cons603)
    rule398 = ReplacementRule(pattern398, lambda g, b, d, h, f, p, a, m, x, c : -(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*d*(-b*h + S(2)*c*g) - x*(b*f*(-a*h + b*g) + S(2)*c**S(2)*d*g - c*(S(2)*a*f*g + b*d*h)) + (-a*f + c*d)*(-S(2)*a*h + b*g))/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**(p + S(1))*Simp(g*(p + S(2))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) + h*x*(m + S(2)*p + S(4))*(-S(2)*a*c*f*g + b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c**S(2)*d*g) - h*(-b*d*(-b*h + S(2)*c*g) + (-a*f + c*d)*(-S(2)*a*h + b*g))*(m + p + S(2)) + S(2)*(p + S(1))*(a*f + c*d)*(c*g**S(2) - h*(-a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(c*g**S(2) - h*(-a*h + b*g))))
    rubi.add(rule398)

    pattern399 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**p_*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons2, cons87, cons88, cons602)
    rule399 = ReplacementRule(pattern399, lambda g, d, h, f, p, a, m, x, c : -(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))*(a*h*(-a*f + c*d) + c*x*(-a*f*g + c*d*g))/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))) + Int((a + c*x**S(2))**(p + S(1))*(g + h*x)**m*Simp(a*h**S(2)*(-a*f + c*d)*(m + p + S(2)) + g*(p + S(2))*(-a*c*f*g + c**S(2)*d*g) + h*x*(-a*c*f*g + c**S(2)*d*g)*(m + S(2)*p + S(4)) + (p + S(1))*(a*f + c*d)*(a*h**S(2) + c*g**S(2)), x), x)/(S(2)*a*c*(p + S(1))*(a*h**S(2) + c*g**S(2))))
    rubi.add(rule399)

    pattern400 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons407, cons571, cons424)
    rule400 = ReplacementRule(pattern400, lambda g, b, d, h, p, f, e, a, m, x, c : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule400)

    pattern401 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons571, cons424)
    rule401 = ReplacementRule(pattern401, lambda g, d, h, p, f, e, a, m, x, c : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + h*x*(-e*h + S(2)*f*g)), x)/h**S(2))
    rubi.add(rule401)

    pattern402 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons407, cons424)
    rule402 = ReplacementRule(pattern402, lambda g, b, d, h, p, f, a, m, x, c : f*Int((g + h*x)**(m + S(2))*(a + b*x + c*x**S(2))**p, x)/h**S(2) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule402)

    pattern403 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(g_ + x_*WC('h', S(1)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons424)
    rule403 = ReplacementRule(pattern403, lambda g, d, h, p, f, a, m, x, c : f*Int((a + c*x**S(2))**p*(g + h*x)**(m + S(2)), x)/h**S(2) - Int((a + c*x**S(2))**p*(g + h*x)**m*(-d*h**S(2) + f*g**S(2) + S(2)*f*g*h*x), x)/h**S(2))
    rubi.add(rule403)

    pattern404 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons407, cons571, cons452)
    rule404 = ReplacementRule(pattern404, lambda g, b, d, h, p, f, e, a, m, x, c : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))) + x*(b*f*h*(m + p + S(2)) + c*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1)))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule404)

    pattern405 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons571, cons452)
    rule405 = ReplacementRule(pattern405, lambda g, d, h, p, f, e, a, m, x, c : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(c*x*(-e*h*(m + S(2)*p + S(3)) + S(2)*f*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule405)

    pattern406 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))*(x_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons407, cons452)
    rule406 = ReplacementRule(pattern406, lambda g, b, d, h, p, f, a, m, x, c : f*(g + h*x)**(m + S(1))*(a + b*x + c*x**S(2))**(p + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*Simp(b*f*g*(p + S(1)) + f*x*(b*h*(m + p + S(2)) + S(2)*c*g*(p + S(1))) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule406)

    pattern407 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)))*(g_ + x_*WC('h', S(1)))**WC('m', S(1)), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons2, cons74, cons452)
    rule407 = ReplacementRule(pattern407, lambda g, d, h, f, p, a, m, x, c : f*(a + c*x**S(2))**(p + S(1))*(g + h*x)**(m + S(1))/(c*h*(m + S(2)*p + S(3))) - Int((a + c*x**S(2))**p*(g + h*x)**m*Simp(S(2)*c*f*g*x*(p + S(1)) + h*(a*f*(m + S(1)) - c*d*(m + S(2)*p + S(3))), x), x)/(c*h*(m + S(2)*p + S(3))))
    rubi.add(rule407)

    pattern408 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons174, cons116)
    rule408 = ReplacementRule(pattern408, lambda g, b, d, h, p, f, e, a, q, x, c : Int(ExpandIntegrand((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule408)

    pattern409 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons174, cons604)
    rule409 = ReplacementRule(pattern409, lambda g, d, h, p, f, e, a, q, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule409)

    pattern410 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons308, cons88, cons290)
    rule410 = ReplacementRule(pattern410, lambda g, b, d, h, p, f, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + e*q*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)) + x*(-e*(b*h - S(2)*c*g)*(S(2)*p + q + S(3)) + S(2)*f*q*(-S(2)*a*h + b*g)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule410)

    pattern411 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons308, cons88, cons290)
    rule411 = ReplacementRule(pattern411, lambda g, d, h, p, f, e, a, q, x, c : (a + c*x**S(2))**(p + S(1))*(a*h - c*g*x)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-a*e*h*q + c*d*g*(S(2)*p + S(3)) + c*f*g*x**S(2)*(S(2)*p + S(2)*q + S(3)) + x*(-S(2)*a*f*h*q + c*e*g*(S(2)*p + q + S(3))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule411)

    pattern412 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons308, cons88, cons290)
    rule412 = ReplacementRule(pattern412, lambda g, b, d, h, p, f, a, q, x, c : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g))/((p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(S(2)*p + S(3))*(b*h - S(2)*c*g) + S(2)*f*q*x*(-S(2)*a*h + b*g) - f*x**S(2)*(b*h - S(2)*c*g)*(S(2)*p + S(2)*q + S(3)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule412)

    pattern413 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons173, cons407, cons571, cons87, cons88, cons577, cons312)
    rule413 = ReplacementRule(pattern413, lambda g, b, d, h, p, f, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + c*x*(g*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - h*(a*b*f - S(2)*a*c*e + b*c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)) - e*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - x*(S(2)*f*(c*g*(S(2)*a*c*e - b*(a*f + c*d)) + (-a*h + b*g)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d))) + (p + S(1))*(b*h - S(2)*c*g)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + b**S(2)*f*g - b*(a*f*h + c*d*h + c*e*g) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule413)

    pattern414 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons173, cons571, cons87, cons88, cons579, cons312)
    rule414 = ReplacementRule(pattern414, lambda g, d, h, p, f, e, a, q, x, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d) + c*x*(S(2)*a*c*e*h + g*(-S(2)*a*c*f + S(2)*c**S(2)*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + S(2)*q + S(5)) - S(2)*c*g*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2)) - e*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2)) - x*(c*e*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d))*(S(2)*p + q + S(4)) + S(2)*f*(S(2)*a*c**S(2)*e*g - a*h*(-S(2)*a*c*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*a*c*e*h + S(2)*c*g*(-a*f + c*d)), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule414)

    pattern415 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons173, cons407, cons87, cons88, cons578, cons312)
    rule415 = ReplacementRule(pattern415, lambda g, b, d, h, p, f, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*c*g*(a*f + c*d) + c*x*(g*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - h*(a*b*f + b*c*d)) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) - x*(-b*f*(p + S(1))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)) + S(2)*f*(-b*c*g*(a*f + c*d) + (-a*h + b*g)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b*h - S(2)*c*g)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(b**S(2)*f*g - b*(a*f*h + c*d*h) + S(2)*c*g*(-a*f + c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule415)

    pattern416 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons173, cons407, cons571, cons87, cons116, cons605)
    rule416 = ReplacementRule(pattern416, lambda g, b, d, h, p, f, e, a, q, x, c : h*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-a*e + b*d) + x**S(2)*(c*(e*h - S(2)*f*g)*(p + q + S(1)) + h*p*(-b*f + c*e)) + x*(b*(e*h - S(2)*f*g)*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule416)

    pattern417 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons173, cons571, cons87, cons116, cons605)
    rule417 = ReplacementRule(pattern417, lambda g, d, h, p, f, e, a, q, x, c : h*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*f*(p + q + S(1))) + Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*e*h*p - a*(e*h - S(2)*f*g)*(p + q + S(1)) - S(2)*h*p*x*(-a*f + c*d) - x**S(2)*(c*e*h*p + c*(e*h - S(2)*f*g)*(p + q + S(1))), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule417)

    pattern418 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_*WC('h', S(1)) + WC('g', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons173, cons407, cons87, cons116, cons605)
    rule418 = ReplacementRule(pattern418, lambda g, b, d, h, p, f, a, q, x, c : h*(d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p/(S(2)*f*(p + q + S(1))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-S(2)*a*f*g*(p + q + S(1)) + b*d*h*p + x**S(2)*(-b*f*h*p - S(2)*c*f*g*(p + q + S(1))) + x*(-S(2)*b*f*g*(p + q + S(1)) + S(2)*h*p*(-a*f + c*d)), x), x)/(S(2)*f*(p + q + S(1))))
    rubi.add(rule418)

    pattern419 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, )
    def With419(g, b, d, h, f, e, a, x, c):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(-a*b*f*h + a*c*e*h - a*c*f*g + b**S(2)*f*g - b*c*e*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(a + b*x + c*x**S(2)), x)/q + Int(Simp(a*f**S(2)*g + b*d*f*h - b*e*f*g - c*d*e*h - c*d*f*g + c*e**S(2)*g - f*x*(-a*f*h + b*f*g + c*d*h - c*e*g), x)/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule419 = ReplacementRule(pattern419, lambda g, b, d, h, f, e, a, x, c : With419(g, b, d, h, f, e, a, x, c))
    rubi.add(rule419)

    pattern420 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, )
    def With420(g, b, d, h, f, a, x, c):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int(Simp(a*f**S(2)*g + b*d*f*h - c*d*f*g - f*x*(-a*f*h + b*f*g + c*d*h), x)/(d + f*x**S(2)), x)/q + Int(Simp(-a*b*f*h - a*c*f*g + b**S(2)*f*g + c**S(2)*d*g + c*x*(-a*f*h + b*f*g + c*d*h), x)/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule420 = ReplacementRule(pattern420, lambda g, b, d, h, f, a, x, c : With420(g, b, d, h, f, a, x, c))
    rubi.add(rule420)

    pattern421 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons606)
    rule421 = ReplacementRule(pattern421, lambda g, d, h, f, a, x, c : g*Int(S(1)/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x) + h*Int(x/((a + c*x**S(2))*sqrt(d + f*x**S(2))), x))
    rubi.add(rule421)

    pattern422 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons607, )
    def With422(g, d, h, f, a, x, c):
        q = Rt(-a*c, S(2))
        return -(c*g - h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(c*x + q)), x)/(S(2)*q) - (c*g + h*q)*Int(S(1)/(sqrt(d + f*x**S(2))*(-c*x + q)), x)/(S(2)*q)
    rule422 = ReplacementRule(pattern422, lambda g, d, h, f, a, x, c : With422(g, d, h, f, a, x, c))
    rubi.add(rule422)

    pattern423 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons582, cons608)
    rule423 = ReplacementRule(pattern423, lambda g, b, d, h, f, e, a, x, c : -S(2)*g*Subst(Int(S(1)/(-a*e + b*d - b*x**S(2)), x), x, sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule423)

    pattern424 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons582, cons609)
    rule424 = ReplacementRule(pattern424, lambda g, b, h, f, e, a, x, d, c : h*Int((e + S(2)*f*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f) - (e*h - S(2)*f*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*f))
    rubi.add(rule424)

    pattern425 = Pattern(Integral(x_/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons407, cons571, cons561)
    rule425 = ReplacementRule(pattern425, lambda b, d, f, e, a, x, c : -S(2)*e*Subst(Int((-d*x**S(2) + S(1))/(-b*f + c*e + d**S(2)*x**S(4)*(-b*f + c*e) - e*x**S(2)*(S(2)*a*f - b*e + S(2)*c*d)), x), x, (S(1) + x*(e + sqrt(-S(4)*d*f + e**S(2)))/(S(2)*d))/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule425)

    pattern426 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons561, cons610)
    rule426 = ReplacementRule(pattern426, lambda g, b, d, h, f, e, a, x, c : g*Subst(Int(S(1)/(a + x**S(2)*(-a*f + c*d)), x), x, x/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule426)

    pattern427 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons561, cons611)
    rule427 = ReplacementRule(pattern427, lambda g, b, d, h, f, e, a, x, c : h*Int((S(2)*d + e*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e - (S(2)*d*h - e*g)*Int(S(1)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/e)
    rubi.add(rule427)

    pattern428 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons551, cons612)
    rule428 = ReplacementRule(pattern428, lambda g, b, d, h, f, e, a, x, c : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g) - x**S(2)*(-a*e + b*d), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule428)

    pattern429 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons613)
    rule429 = ReplacementRule(pattern429, lambda g, d, h, f, e, a, x, c : -S(2)*a*g*h*Subst(Int(S(1)/Simp(S(2)*a**S(2)*c*g*h + a*e*x**S(2), x), x), x, Simp(a*h - c*g*x, x)/sqrt(d + e*x + f*x**S(2))))
    rubi.add(rule429)

    pattern430 = Pattern(Integral((g_ + x_*WC('h', S(1)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons614)
    rule430 = ReplacementRule(pattern430, lambda g, b, d, h, f, a, x, c : -S(2)*g*(-S(2)*a*h + b*g)*Subst(Int(S(1)/Simp(-b*d*x**S(2) + g*(-S(4)*a*c + b**S(2))*(-S(2)*a*h + b*g), x), x), x, Simp(-S(2)*a*h + b*g - x*(b*h - S(2)*c*g), x)/sqrt(d + f*x**S(2))))
    rubi.add(rule430)

    pattern431 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons412, )
    def With431(g, b, h, f, e, a, x, d, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/((b + S(2)*c*x - q)*sqrt(d + e*x + f*x**S(2))), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/((b + S(2)*c*x + q)*sqrt(d + e*x + f*x**S(2))), x)/q
    rule431 = ReplacementRule(pattern431, lambda g, b, h, f, e, a, x, d, c : With431(g, b, h, f, e, a, x, d, c))
    rubi.add(rule431)

    pattern432 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons584, )
    def With432(g, h, f, e, a, x, d, c):
        q = Rt(-a*c, S(2))
        return (-c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x + q)*sqrt(d + e*x + f*x**S(2))), x) + (c*g/(S(2)*q) + h/S(2))*Int(S(1)/((c*x - q)*sqrt(d + e*x + f*x**S(2))), x)
    rule432 = ReplacementRule(pattern432, lambda g, h, f, e, a, x, d, c : With432(g, h, f, e, a, x, d, c))
    rubi.add(rule432)

    pattern433 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons412, )
    def With433(g, b, d, h, f, a, x, c):
        q = Rt(-S(4)*a*c + b**S(2), S(2))
        return (S(2)*c*g - h*(b - q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x - q)), x)/q - (S(2)*c*g - h*(b + q))*Int(S(1)/(sqrt(d + f*x**S(2))*(b + S(2)*c*x + q)), x)/q
    rule433 = ReplacementRule(pattern433, lambda g, b, d, h, f, a, x, c : With433(g, b, d, h, f, a, x, c))
    rubi.add(rule433)

    pattern434 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, cons551, cons585, )
    def With434(g, b, d, h, f, e, a, x, c):
        q = Rt(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-g*(-a*f + c*d - q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d + q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-g*(-a*f + c*d + q) + h*(-a*e + b*d) - x*(g*(-b*f + c*e) - h*(-a*f + c*d - q)), x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule434 = ReplacementRule(pattern434, lambda g, b, d, h, f, e, a, x, c : With434(g, b, d, h, f, e, a, x, c))
    rubi.add(rule434)

    pattern435 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons571, cons586, )
    def With435(g, h, f, e, a, x, d, c):
        q = Rt(a*c*e**S(2) + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(-a*e*h - g*(-a*f + c*d - q) + x*(-c*e*g + h*(-a*f + c*d + q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q) - Int(Simp(-a*e*h - g*(-a*f + c*d + q) + x*(-c*e*g + h*(-a*f + c*d - q)), x)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/(S(2)*q)
    rule435 = ReplacementRule(pattern435, lambda g, h, f, e, a, x, d, c : With435(g, h, f, e, a, x, d, c))
    rubi.add(rule435)

    pattern436 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, cons585, )
    def With436(g, b, d, h, f, a, x, c):
        q = Rt(b**S(2)*d*f + (-a*f + c*d)**S(2), S(2))
        return Int(Simp(b*d*h - g*(-a*f + c*d - q) + x*(b*f*g + h*(-a*f + c*d + q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q) - Int(Simp(b*d*h - g*(-a*f + c*d + q) + x*(b*f*g + h*(-a*f + c*d - q)), x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/(S(2)*q)
    rule436 = ReplacementRule(pattern436, lambda g, b, d, h, f, a, x, c : With436(g, b, d, h, f, a, x, c))
    rubi.add(rule436)

    pattern437 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons407, cons571, )
    def With437(g, b, d, h, f, e, a, x, c):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f + e**S(2), S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + x*(e + t))*sqrt(b + S(2)*c*x + s)*sqrt(e + S(2)*f*x + t)), x)/(sqrt(a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2)))
    rule437 = ReplacementRule(pattern437, lambda g, b, d, h, f, e, a, x, c : With437(g, b, d, h, f, e, a, x, c))
    rubi.add(rule437)

    pattern438 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/(sqrt(d_ + x_**S(2)*WC('f', S(1)))*sqrt(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons161, cons162, cons407, )
    def With438(g, b, d, h, f, a, x, c):
        s = Rt(-S(4)*a*c + b**S(2), S(2))
        t = Rt(-S(4)*d*f, S(2))
        return sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)*Int((g + h*x)/(sqrt(S(2)*a + x*(b + s))*sqrt(S(2)*d + t*x)*sqrt(S(2)*f*x + t)*sqrt(b + S(2)*c*x + s)), x)/(sqrt(d + f*x**S(2))*sqrt(a + b*x + c*x**S(2)))
    rule438 = ReplacementRule(pattern438, lambda g, b, d, h, f, a, x, c : With438(g, b, d, h, f, a, x, c))
    rubi.add(rule438)

    pattern439 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons582, cons615, cons616, cons617, )
    def With439(g, b, d, h, f, e, a, x, c):
        q = S(3)**(S(2)/3)*(-c*h**S(2)/(-b*h + S(2)*c*g)**S(2))**(S(1)/3)
        return sqrt(S(3))*h*q*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3)/(S(3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3)) + sqrt(S(3))/S(3))/f - S(3)*h*q*log((-S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(2)/3) + S(2)**(S(1)/3)*(S(3)*h*(b + S(2)*c*x)/(-b*h + S(2)*c*g) + S(1))**(S(1)/3))/(S(2)*f) + h*q*log(d + e*x + f*x**S(2))/(S(2)*f)
    rule439 = ReplacementRule(pattern439, lambda g, b, d, h, f, e, a, x, c : With439(g, b, d, h, f, e, a, x, c))
    rubi.add(rule439)

    pattern440 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons618, cons619, cons17)
    rule440 = ReplacementRule(pattern440, lambda g, d, h, f, a, x, c : S(2)**(S(1)/3)*sqrt(S(3))*h*ArcTan(-S(2)**(S(2)/3)*sqrt(S(3))*(S(1) - S(3)*h*x/g)**(S(2)/3)/(S(3)*(S(1) + S(3)*h*x/g)**(S(1)/3)) + sqrt(S(3))/S(3))/(S(2)*a**(S(1)/3)*f) + S(2)**(S(1)/3)*h*log(d + f*x**S(2))/(S(4)*a**(S(1)/3)*f) - S(3)*S(2)**(S(1)/3)*h*log((S(1) - S(3)*h*x/g)**(S(2)/3) + S(2)**(S(1)/3)*(S(1) + S(3)*h*x/g)**(S(1)/3))/(S(4)*a**(S(1)/3)*f))
    rubi.add(rule440)

    pattern441 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))/((x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**(S(1)/3)*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons582, cons615, cons616, cons495, )
    def With441(g, b, d, h, f, e, a, x, c):
        q = -c/(-S(4)*a*c + b**S(2))
        return (q*(a + b*x + c*x**S(2)))**(S(1)/3)*Int((g + h*x)/((d + e*x + f*x**S(2))*(a*q + b*q*x + c*q*x**S(2))**(S(1)/3)), x)/(a + b*x + c*x**S(2))**(S(1)/3)
    rule441 = ReplacementRule(pattern441, lambda g, b, d, h, f, e, a, x, c : With441(g, b, d, h, f, e, a, x, c))
    rubi.add(rule441)

    pattern442 = Pattern(Integral((g_ + x_*WC('h', S(1)))/((a_ + x_**S(2)*WC('c', S(1)))**(S(1)/3)*(d_ + x_**S(2)*WC('f', S(1)))), x_), cons4, cons9, cons10, cons73, cons161, cons162, cons618, cons619, cons188)
    rule442 = ReplacementRule(pattern442, lambda g, d, h, f, a, x, c : (S(1) + c*x**S(2)/a)**(S(1)/3)*Int((g + h*x)/((S(1) + c*x**S(2)/a)**(S(1)/3)*(d + f*x**S(2))), x)/(a + c*x**S(2))**(S(1)/3))
    rubi.add(rule442)

    pattern443 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons173, cons620)
    rule443 = ReplacementRule(pattern443, lambda g, b, d, h, f, p, e, a, q, x, c : Int((g + h*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule443)

    pattern444 = Pattern(Integral((x_*WC('h', S(1)) + WC('g', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons74, cons173, cons621)
    rule444 = ReplacementRule(pattern444, lambda g, h, f, p, e, a, q, x, d, c : Int((a + c*x**S(2))**p*(g + h*x)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule444)

    pattern445 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons173, cons6, cons7)
    rule445 = ReplacementRule(pattern445, lambda g, b, d, h, u, p, f, e, a, m, q, x, c : Subst(Int((g + h*x)**m*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule445)

    pattern446 = Pattern(Integral((u_*WC('h', S(1)) + WC('g', S(0)))**WC('m', S(1))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons161, cons162, cons2, cons74, cons173, cons6, cons7)
    rule446 = ReplacementRule(pattern446, lambda g, d, h, u, p, f, e, a, m, q, x, c : Subst(Int((a + c*x**S(2))**p*(g + h*x)**m*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule446)

    pattern447 = Pattern(Integral(u_**WC('p', S(1))*v_**WC('q', S(1))*z_**WC('m', S(1)), x_), cons2, cons74, cons173, cons622, cons623, cons624)
    rule447 = ReplacementRule(pattern447, lambda u, p, z, m, q, x, v : Int(ExpandToSum(u, x)**p*ExpandToSum(v, x)**q*ExpandToSum(z, x)**m, x))
    rubi.add(rule447)

    pattern448 = Pattern(Integral((d_ + x_*WC('e', S(1)))**WC('m', S(1))*(f_ + x_*WC('g', S(1)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons178, cons2, cons13, cons74, cons173, cons518, cons70, cons556)
    rule448 = ReplacementRule(pattern448, lambda g, b, d, h, p, i, f, e, a, m, q, n, x, c : Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule448)

    pattern449 = Pattern(Integral((x_*WC('e', S(1)) + WC('d', S(0)))**WC('m', S(1))*(x_*WC('g', S(1)) + WC('f', S(0)))**WC('n', S(1))*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons178, cons2, cons13, cons74, cons173, cons77, cons27)
    rule449 = ReplacementRule(pattern449, lambda g, b, h, f, i, p, e, a, m, q, n, x, d, c : Int(ExpandIntegrand((d + e*x)**m*(f + g*x)**n*(h + i*x)**q*(a + b*x + c*x**S(2))**p, x), x))
    rubi.add(rule449)

    pattern450 = Pattern(Integral((d_ + x_*WC('e', S(1)))**m_*(f_ + x_*WC('g', S(1)))**n_*(x_*WC('i', S(1)) + WC('h', S(0)))**WC('q', S(1))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons161, cons162, cons178, cons2, cons13, cons74, cons173, cons518, cons70)
    rule450 = ReplacementRule(pattern450, lambda g, b, d, h, p, i, f, e, a, m, q, n, x, c : (d + e*x)**FracPart(m)*(f + g*x)**FracPart(m)*(d*f + e*g*x**S(2))**(-FracPart(m))*Int((h + i*x)**q*(d*f + e*g*x**S(2))**m*(a + b*x + c*x**S(2))**p, x))
    rubi.add(rule450)

    pattern451 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons560, cons561, cons562, cons563)
    rule451 = ReplacementRule(pattern451, lambda b, d, A, p, C, f, B, e, a, q, x, c : (c/f)**p*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule451)

    pattern452 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons560, cons561, cons562, cons563)
    rule452 = ReplacementRule(pattern452, lambda b, d, A, p, C, f, e, a, q, x, c : (c/f)**p*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule452)

    pattern453 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons560, cons561, cons100, cons324, cons564)
    rule453 = ReplacementRule(pattern453, lambda b, d, A, p, C, f, B, e, a, q, x, c : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule453)

    pattern454 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons560, cons561, cons100, cons324, cons564)
    rule454 = ReplacementRule(pattern454, lambda b, d, A, p, C, f, e, a, q, x, c : a**IntPart(p)*d**(-IntPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*(d + e*x + f*x**S(2))**(-FracPart(p))*Int((A + C*x**S(2))*(d + e*x + f*x**S(2))**(p + q), x))
    rubi.add(rule454)

    pattern455 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons405)
    rule455 = ReplacementRule(pattern455, lambda b, d, A, p, C, f, B, e, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule455)

    pattern456 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons405)
    rule456 = ReplacementRule(pattern456, lambda b, d, A, p, C, f, e, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule456)

    pattern457 = Pattern(Integral((x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons74, cons173, cons405)
    rule457 = ReplacementRule(pattern457, lambda b, d, A, p, C, f, B, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q*(A + B*x + C*x**S(2)), x))
    rubi.add(rule457)

    pattern458 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons74, cons173, cons405)
    rule458 = ReplacementRule(pattern458, lambda b, d, A, p, C, f, a, q, x, c : (S(4)*c)**(-IntPart(p))*(b + S(2)*c*x)**(-S(2)*FracPart(p))*(a + b*x + c*x**S(2))**FracPart(p)*Int((A + C*x**S(2))*(b + S(2)*c*x)**(S(2)*p)*(d + f*x**S(2))**q, x))
    rubi.add(rule458)

    pattern459 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons407, cons571, cons174, cons116)
    rule459 = ReplacementRule(pattern459, lambda b, d, A, p, C, f, B, e, a, q, x, c : Int(ExpandIntegrand((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule459)

    pattern460 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons407, cons571, cons174, cons116)
    rule460 = ReplacementRule(pattern460, lambda b, d, A, p, C, f, e, a, q, x, c : Int(ExpandIntegrand((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule460)

    pattern461 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons571, cons174, cons604)
    rule461 = ReplacementRule(pattern461, lambda d, A, p, C, f, B, e, a, q, x, c : Int(ExpandIntegrand((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule461)

    pattern462 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons571, cons174, cons604)
    rule462 = ReplacementRule(pattern462, lambda d, A, p, C, f, e, a, q, x, c : Int(ExpandIntegrand((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x))
    rubi.add(rule462)

    pattern463 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons407, cons571, cons308, cons88, cons290)
    rule463 = ReplacementRule(pattern463, lambda b, d, A, p, C, f, B, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + e*q*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))) + x*(-e*(C*(S(2)*a*c*(q + S(1)) - b**S(2)*(p + q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + q + S(3))) + S(2)*f*q*(A*b*c - S(2)*B*a*c + C*a*b)), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule463)

    pattern464 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons407, cons571, cons308, cons88, cons290)
    rule464 = ReplacementRule(pattern464, lambda b, d, A, p, C, f, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*(b*e*q + S(2)*c*d*(S(2)*p + S(3))) - C*(-a*b*e*q + S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*c*(b*f*q + c*e*(S(2)*p + q + S(3))) + C*(S(2)*a*b*f*q - S(2)*a*c*e*(q + S(1)) + b**S(2)*e*(p + q + S(2)))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule464)

    pattern465 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons571, cons308, cons88, cons290)
    rule465 = ReplacementRule(pattern465, lambda d, A, p, C, f, B, e, a, q, x, c : (a + c*x**S(2))**(p + S(1))*(B*a - x*(A*c - C*a))*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - a*(B*e*q + C*d) - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - a*(S(2)*B*f*q + C*e*(q + S(1)))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule465)

    pattern466 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons571, cons308, cons88, cons290)
    rule466 = ReplacementRule(pattern466, lambda d, A, p, C, f, e, a, q, x, c : -x*(a + c*x**S(2))**(p + S(1))*(A*c - C*a)*(d + e*x + f*x**S(2))**q/(S(2)*a*c*(p + S(1))) + Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(-1))*Simp(A*c*d*(S(2)*p + S(3)) - C*a*d - f*x**S(2)*(-A*c*(S(2)*p + S(2)*q + S(3)) + C*a*(S(2)*q + S(1))) + x*(A*c*e*(S(2)*p + q + S(3)) - C*a*e*(q + S(1))), x), x)/(S(2)*a*c*(p + S(1))))
    rubi.add(rule466)

    pattern467 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons407, cons308, cons88, cons290)
    rule467 = ReplacementRule(pattern467, lambda b, d, A, p, C, f, B, a, q, x, c : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c - S(2)*B*a*c + C*a*b - x*(-C*(-S(2)*a*c + b**S(2)) + c*(-S(2)*A*c + B*b)))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-d*(C*(S(2)*a*c - b**S(2)*(p + S(2))) + c*(S(2)*p + S(3))*(-S(2)*A*c + B*b)) + S(2)*f*q*x*(A*b*c - S(2)*B*a*c + C*a*b) - f*x**S(2)*(C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2))) + c*(-S(2)*A*c + B*b)*(S(2)*p + S(2)*q + S(3))), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule467)

    pattern468 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons407, cons308, cons88, cons290)
    rule468 = ReplacementRule(pattern468, lambda b, d, A, p, C, f, a, q, x, c : (d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*c + C*a*b + x*(S(2)*A*c**S(2) + C*(-S(2)*a*c + b**S(2))))/(c*(p + S(1))*(-S(4)*a*c + b**S(2))) - Int((d + f*x**S(2))**(q + S(-1))*(a + b*x + c*x**S(2))**(p + S(1))*Simp(S(2)*A*c**S(2)*d*(S(2)*p + S(3)) - C*(S(2)*a*c*d - b**S(2)*d*(p + S(2))) - f*x**S(2)*(-S(2)*A*c**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*a*c*(S(2)*q + S(1)) - b**S(2)*(p + S(2)*q + S(2)))) + x*(S(2)*A*b*c*f*q + S(2)*C*a*b*f*q), x), x)/(c*(p + S(1))*(-S(4)*a*c + b**S(2))))
    rubi.add(rule468)

    pattern469 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons173, cons407, cons571, cons87, cons88, cons577, cons312)
    rule469 = ReplacementRule(pattern469, lambda b, d, A, p, C, f, B, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) - B*(a*b*f - S(2)*a*c*e + b*c*d) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)) - e*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*((A*b - B*a)*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e))) + (p + S(1))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + B*a*f + B*c*d + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule469)

    pattern470 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons173, cons407, cons571, cons87, cons88, cons577, cons312)
    rule470 = ReplacementRule(pattern470, lambda b, d, A, p, C, f, e, a, q, x, c : (a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + c*x*(A*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + C*(-a*b*e - S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))) + Int((a + b*x + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)) - e*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - x*(S(2)*f*(A*b*(b**S(2)*f + S(2)*c**S(2)*d - c*(S(2)*a*f + b*e)) + (A*c - C*a)*(S(2)*a*c*e - b*(a*f + c*d)))*(p + q + S(2)) - (b*f*(p + S(1)) - c*e*(S(2)*p + q + S(4)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(A*c*e + C*a*e)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(-(-a*e + b*d)*(-b*f + c*e) + (-a*f + c*d)**S(2))))
    rubi.add(rule470)

    pattern471 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons173, cons571, cons87, cons88, cons579, cons312)
    rule471 = ReplacementRule(pattern471, lambda d, A, p, C, f, B, e, a, q, x, c : -(a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**(q + S(1))*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*B*a*c*e - S(2)*C*a*(-a*f + c*d)))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - e*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2)) - x*(c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(S(2)*p + q + S(4)) + S(2)*f*(-B*a*(-S(2)*a*c*f + S(2)*c**S(2)*d) + S(2)*a*c*e*(A*c - C*a))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-B*c*e - C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule471)

    pattern472 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons173, cons571, cons87, cons88, cons579, cons312)
    rule472 = ReplacementRule(pattern472, lambda d, A, p, C, f, e, a, q, x, c : -(a + c*x**S(2))**(p + S(1))*(S(2)*a*c*e*(A*c - C*a) + c*x*(A*(-S(2)*a*c*f + S(2)*c**S(2)*d) - S(2)*C*a*(-a*f + c*d)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))) - Int((a + c*x**S(2))**(p + S(1))*(d + e*x + f*x**S(2))**q*Simp(-S(2)*a*c*e**S(2)*(A*c - C*a)*(p + q + S(2)) - c*f*x**S(2)*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + S(2)*q + S(5)) - x*(S(4)*a*c*e*f*(A*c - C*a)*(p + q + S(2)) + c*e*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(S(2)*p + q + S(4))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(a*c*e**S(2) + (-a*f + c*d)**S(2)) + (S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d))*(a*f*(p + S(1)) - c*d*(p + S(2))), x), x)/(S(4)*a*c*(p + S(1))*(a*c*e**S(2) + (-a*f + c*d)**S(2))))
    rubi.add(rule472)

    pattern473 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons173, cons407, cons87, cons88, cons578, cons312)
    rule473 = ReplacementRule(pattern473, lambda b, d, A, p, C, f, B, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(-b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - B*(a*b*f + b*c*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)) + S(2)*f*(-b*(A*c - C*a)*(a*f + c*d) + (A*b - B*a)*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d))*(p + q + S(2))) + (p + S(1))*(b**S(2)*d*f + (-a*f + c*d)**S(2))*(-S(2)*A*c + B*b - S(2)*C*a) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d) - b*(B*a*f + B*c*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule473)

    pattern474 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons173, cons407, cons87, cons88, cons578, cons312)
    rule474 = ReplacementRule(pattern474, lambda b, d, A, p, C, f, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**(p + S(1))*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d) + c*x*(A*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) + C*(-S(2)*a*(-a*f + c*d) + b**S(2)*d)))/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))) + Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(1))*Simp(-c*f*x**S(2)*(S(2)*p + S(2)*q + S(5))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) - x*(-b*f*(p + S(1))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)) + S(2)*f*(A*b*(-S(2)*a*c*f + b**S(2)*f + S(2)*c**S(2)*d) - b*(A*c - C*a)*(a*f + c*d))*(p + q + S(2))) + (p + S(1))*(-S(2)*A*c - S(2)*C*a)*(b**S(2)*d*f + (-a*f + c*d)**S(2)) + (a*f*(p + S(1)) - c*d*(p + S(2)))*(S(2)*A*c*(-a*f + c*d) - S(2)*a*(-C*a*f + C*c*d) + b**S(2)*(A*f + C*d)), x), x)/((p + S(1))*(-S(4)*a*c + b**S(2))*(b**S(2)*d*f + (-a*f + c*d)**S(2))))
    rubi.add(rule474)

    pattern475 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons173, cons407, cons571, cons87, cons116, cons605, cons628)
    rule475 = ReplacementRule(pattern475, lambda b, d, A, p, C, f, B, e, a, q, x, c : (a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) + S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(p*(-b*f + c*e)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(C*(q + S(1))*(-b*f + c*e) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule475)

    pattern476 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons173, cons407, cons571, cons87, cons116, cons605, cons628)
    rule476 = ReplacementRule(pattern476, lambda b, d, A, p, C, f, e, a, q, x, c : (S(2)*C*c*f*x*(p + q + S(1)) + C*(b*f*p - c*e*(S(2)*p + q + S(2))))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + b*x + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(p*(-a*e + b*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + x**S(2)*(p*(-b*f + c*e)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(S(2)*p*(-a*f + c*d)*(-C*c*e*(S(2)*p + S(2)*q + S(3)) + C*(q + S(1))*(-b*f + c*e)) + (C*e*f*p*(-S(4)*a*c + b**S(2)) - b*c*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2)))))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule476)

    pattern477 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons173, cons571, cons87, cons116, cons605, cons628)
    rule477 = ReplacementRule(pattern477, lambda d, A, p, C, f, B, e, a, q, x, c : (a + c*x**S(2))**p*(d + e*x + f*x**S(2))**(q + S(1))*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))) + f*(-S(2)*A*f + B*e)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*p + S(2)*q + S(3))*(S(2)*A*f - B*e + S(2)*C*d)))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - c*(-B*f + C*e)*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule477)

    pattern478 = Pattern(Integral((a_ + x_**S(2)*WC('c', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons173, cons571, cons87, cons116, cons605, cons628)
    rule478 = ReplacementRule(pattern478, lambda d, A, p, C, f, e, a, q, x, c : (a + c*x**S(2))**p*(-C*c*e*(S(2)*p + q + S(2)) + S(2)*C*c*f*x*(p + q + S(1)))*(d + e*x + f*x**S(2))**(q + S(1))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((a + c*x**S(2))**(p + S(-1))*(d + e*x + f*x**S(2))**q*Simp(a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + C*(S(2)*d*f - e**S(2)*(S(2)*p + q + S(2))))*(p + q + S(1)) - a*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + x**S(2)*(c*e*p*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3))) + (-S(4)*C*a*c*f**S(2)*p - c**S(2)*(C*(-S(4)*d*f + e**S(2))*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(4)*C*a*c*e*f*p*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(C*c*e*(q + S(1)) - C*c*e*(S(2)*p + S(2)*q + S(3)))), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule478)

    pattern479 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1))*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons173, cons407, cons87, cons116, cons605, cons628)
    rule479 = ReplacementRule(pattern479, lambda b, d, A, p, C, f, B, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(a + b*x + c*x**S(2))**p*(B*c*f*(S(2)*p + S(2)*q + S(3)) + C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(b*d*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + x**S(2)*(-b*f*p*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1))) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1)) + S(2)*p*(-a*f + c*d)*(B*c*f*(S(2)*p + S(2)*q + S(3)) - C*b*f*(q + S(1)))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule479)

    pattern480 = Pattern(Integral((d_ + x_**S(2)*WC('f', S(1)))**WC('q', S(1))*(x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))**WC('p', S(1)), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons173, cons407, cons87, cons116, cons605, cons628)
    rule480 = ReplacementRule(pattern480, lambda b, d, A, p, C, f, a, q, x, c : (d + f*x**S(2))**(q + S(1))*(C*b*f*p + S(2)*C*c*f*x*(p + q + S(1)))*(a + b*x + c*x**S(2))**p/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))) - Int((d + f*x**S(2))**q*(a + b*x + c*x**S(2))**(p + S(-1))*Simp(-C*b**S(2)*d*f*p*(q + S(1)) + x**S(2)*(C*b**S(2)*f**S(2)*p*(q + S(1)) + (C*f**S(2)*p*(-S(4)*a*c + b**S(2)) - c**S(2)*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3))))*(p + q + S(1))) + x*(-S(2)*C*b*f*p*(q + S(1))*(-a*f + c*d) - b*c*(-S(4)*C*d*f*(S(2)*p + q + S(2)) + f*(S(2)*A*f + S(2)*C*d)*(S(2)*p + S(2)*q + S(3)))*(p + q + S(1))) + (C*b**S(2)*d*f*p + a*c*(-S(2)*A*f**S(2)*(S(2)*p + S(2)*q + S(3)) + S(2)*C*d*f))*(p + q + S(1)), x), x)/(S(2)*c*f**S(2)*(p + q + S(1))*(S(2)*p + S(2)*q + S(3))))
    rubi.add(rule480)

    pattern481 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons407, cons571, )
    def With481(b, d, A, C, f, B, e, a, x, c):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d - B*a*b*f + B*a*c*e + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) + B*b*d*f - B*c*d*e - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e - B*a*f + B*c*d + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule481 = ReplacementRule(pattern481, lambda b, d, A, C, f, B, e, a, x, c : With481(b, d, A, C, f, B, e, a, x, c))
    rubi.add(rule481)

    pattern482 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*(d_ + x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons407, cons571, )
    def With482(b, d, A, C, f, e, a, x, c):
        q = a**S(2)*f**S(2) - a*b*e*f - S(2)*a*c*d*f + a*c*e**S(2) + b**S(2)*d*f - b*c*d*e + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((-A*a*c*f + A*b**S(2)*f - A*b*c*e + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(a + b*x + c*x**S(2)), x)/q + Int((A*a*f**S(2) - A*b*e*f - A*c*d*f + A*c*e**S(2) - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - A*c*e + C*a*e - C*b*d))/(d + e*x + f*x**S(2)), x)/q
        print("Unable to Integrate")
    rule482 = ReplacementRule(pattern482, lambda b, d, A, C, f, e, a, x, c : With482(b, d, A, C, f, e, a, x, c))
    rubi.add(rule482)

    pattern483 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons407, )
    def With483(b, d, A, C, f, B, a, x, c):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f + B*b*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d - B*a*b*f + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - B*a*f + B*c*d - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule483 = ReplacementRule(pattern483, lambda b, d, A, C, f, B, a, x, c : With483(b, d, A, C, f, B, a, x, c))
    rubi.add(rule483)

    pattern484 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((d_ + x_**S(2)*WC('f', S(1)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons407, )
    def With484(b, d, f, A, C, a, x, c):
        q = a**S(2)*f**S(2) - S(2)*a*c*d*f + b**S(2)*d*f + c**S(2)*d**S(2)
        if NonzeroQ(q):
            return Int((A*a*f**S(2) - A*c*d*f - C*a*d*f + C*c*d**S(2) - f*x*(A*b*f - C*b*d))/(d + f*x**S(2)), x)/q + Int((-A*a*c*f + A*b**S(2)*f + A*c**S(2)*d + C*a**S(2)*f - C*a*c*d + c*x*(A*b*f - C*b*d))/(a + b*x + c*x**S(2)), x)/q
        print("Unable to Integrate")
    rule484 = ReplacementRule(pattern484, lambda b, d, f, A, C, a, x, c : With484(b, d, f, A, C, a, x, c))
    rubi.add(rule484)

    pattern485 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons407, cons571)
    rule485 = ReplacementRule(pattern485, lambda b, d, A, C, f, B, e, a, x, c : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule485)

    pattern486 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons407, cons571)
    rule486 = ReplacementRule(pattern486, lambda b, d, A, C, f, e, a, x, c : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/((a + b*x + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule486)

    pattern487 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons571)
    rule487 = ReplacementRule(pattern487, lambda d, A, C, f, B, e, a, x, c : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + Int((A*c + B*c*x - C*a)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule487)

    pattern488 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/((a_ + x_**S(2)*WC('c', S(1)))*sqrt(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons571)
    rule488 = ReplacementRule(pattern488, lambda d, A, C, f, e, a, x, c : C*Int(S(1)/sqrt(d + e*x + f*x**S(2)), x)/c + (A*c - C*a)*Int(S(1)/((a + c*x**S(2))*sqrt(d + e*x + f*x**S(2))), x)/c)
    rubi.add(rule488)

    pattern489 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons626, cons627, cons407)
    rule489 = ReplacementRule(pattern489, lambda b, d, A, C, f, B, a, x, c : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a + x*(B*c - C*b))/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule489)

    pattern490 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))/(sqrt(x_**S(2)*WC('f', S(1)) + WC('d', S(0)))*(a_ + x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)))), x_), cons4, cons5, cons9, cons10, cons73, cons625, cons627, cons407)
    rule490 = ReplacementRule(pattern490, lambda b, d, A, C, f, a, x, c : C*Int(S(1)/sqrt(d + f*x**S(2)), x)/c + Int((A*c - C*a - C*b*x)/(sqrt(d + f*x**S(2))*(a + b*x + c*x**S(2))), x)/c)
    rubi.add(rule490)

    pattern491 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons629)
    rule491 = ReplacementRule(pattern491, lambda b, d, A, C, f, p, B, e, a, q, x, c : Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule491)

    pattern492 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + x_*WC('b', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons630)
    rule492 = ReplacementRule(pattern492, lambda b, d, A, C, f, p, e, a, q, x, c : Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule492)

    pattern493 = Pattern(Integral((x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('C', S(1)) + x_*WC('B', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons631)
    rule493 = ReplacementRule(pattern493, lambda d, A, C, f, p, B, e, a, q, x, c : Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule493)

    pattern494 = Pattern(Integral((x_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(x_**S(2)*WC('c', S(1)) + WC('a', S(0)))**p_*(x_**S(2)*WC('f', S(1)) + x_*WC('e', S(1)) + WC('d', S(0)))**q_, x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons632)
    rule494 = ReplacementRule(pattern494, lambda A, C, f, p, e, a, q, x, d, c : Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x))
    rubi.add(rule494)

    pattern495 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons6, cons7)
    rule495 = ReplacementRule(pattern495, lambda b, d, u, A, p, C, f, B, e, a, q, x, c : Subst(Int((A + B*x + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule495)

    pattern496 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons6, cons7)
    rule496 = ReplacementRule(pattern496, lambda b, d, u, A, p, f, B, e, a, q, x, c : Subst(Int((A + B*x)*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule496)

    pattern497 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + u_*WC('b', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons5, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons6, cons7)
    rule497 = ReplacementRule(pattern497, lambda b, d, u, A, p, C, f, e, a, q, x, c : Subst(Int((A + C*x**S(2))*(a + b*x + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule497)

    pattern498 = Pattern(Integral((u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('C', S(1)) + u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons6, cons7)
    rule498 = ReplacementRule(pattern498, lambda d, u, A, p, C, f, B, e, a, q, x, c : Subst(Int((a + c*x**S(2))**p*(A + B*x + C*x**S(2))*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule498)

    pattern499 = Pattern(Integral((u_*WC('B', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons626, cons627, cons74, cons173, cons6, cons7)
    rule499 = ReplacementRule(pattern499, lambda d, u, A, p, f, B, e, a, q, x, c : Subst(Int((A + B*x)*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule499)

    pattern500 = Pattern(Integral((u_**S(2)*WC('C', S(1)) + WC('A', S(0)))*(u_**S(2)*WC('c', S(1)) + WC('a', S(0)))**WC('p', S(1))*(u_**S(2)*WC('f', S(1)) + u_*WC('e', S(1)) + WC('d', S(0)))**WC('q', S(1)), x_), cons4, cons9, cons10, cons72, cons73, cons625, cons627, cons74, cons173, cons6, cons7)
    rule500 = ReplacementRule(pattern500, lambda d, u, A, p, C, f, e, a, q, x, c : Subst(Int((A + C*x**S(2))*(a + c*x**S(2))**p*(d + e*x + f*x**S(2))**q, x), x, u)/Coefficient(u, x, S(1)))
    rubi.add(rule500)

    return rubi
