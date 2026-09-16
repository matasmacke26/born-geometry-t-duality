(* ::Package:: *)

(* Corrected Born T-duality implementation. Mathematica 13+.
   Source audit performed; full execution requires a Wolfram kernel.
   No forced tensor symmetries or assumed successful checks. *)
BeginPackage["BornTDual`"];
RunBornExample::usage = "RunBornExample - Born T-duality notebook API.";
HorizonReport::usage = "HorizonReport - Born T-duality notebook API.";
ShowValidation::usage = "ShowValidation - Born T-duality notebook API.";
ShowNonZeroConnection::usage = "ShowNonZeroConnection - Born T-duality notebook API.";
ShowNonZeroBornTorsion::usage = "ShowNonZeroBornTorsion - Born T-duality notebook API.";
ShowNonZeroRicci::usage = "ShowNonZeroRicci - Born T-duality notebook API.";
ShowNonZeroGeneralisedRicci::usage = "ShowNonZeroGeneralisedRicci - Born T-duality notebook API.";
ShowNonZeroIndependentRiemann::usage = "ShowNonZeroIndependentRiemann - Born T-duality notebook API.";
PrintGeometrySummary::usage = "PrintGeometrySummary - Born T-duality notebook API.";
r::usage = "r - Born T-duality notebook API.";
mBH::usage = "mBH - Born T-duality notebook API.";
qCharge::usage = "qCharge - Born T-duality notebook API.";
th::usage = "th - Born T-duality notebook API.";
ph::usage = "ph - Born T-duality notebook API.";
v::usage = "v - Born T-duality notebook API.";
f::usage = "f - Born T-duality notebook API.";
Begin["`Private`"];
ClearAll[CurvatureFromFrame];
CurvatureFromFrame[name_String,x_List,lcOmega_,etaDD_,geDD_,anchor_] := Module[
    {d=Length[x],etaUU,geUU,lcOmegaLow,anchorD,firstBlock,lastTerm,
     rmOrderedLocal,rmFun,rmArray,rawChecks,ricFull,rcGen,rScalar,
     pairsLocal,pairMetric,rmPair,kScalar},
    etaUU=qMap[Inverse[etaDD],2]; geUU=qMap[Inverse[geDD],2];
    lcOmegaLow=Table[qSimp[Total[Flatten[Table[lcOmega[[A,B,E]] etaDD[[E,C]],{E,2 d}]]]],
        {A,2 d},{B,2 d},{C,2 d}];
    anchorD[A_Integer,expr_] := Total[Flatten[Table[anchor[[A,mu]] D[expr,x[[mu]]],{mu,d}]]];
    Clear[firstBlock, lastTerm, rmOrderedLocal, rmFun];

    firstBlock[A_Integer, B_Integer, C_Integer, Dd_Integer] :=
        anchorD[A, lcOmegaLow[[B, C, Dd]]] -
        Total[Flatten[Table[
            lcOmega[[A, Dd, E]] etaDD[[E, F]]
                lcOmega[[B, C, F]],
            {E, 2 d}, {F, 2 d}
        ]]] -
        anchorD[B, lcOmegaLow[[A, C, Dd]]] +
        Total[Flatten[Table[
            lcOmega[[B, Dd, E]] etaDD[[E, F]]
                lcOmega[[A, C, F]],
            {E, 2 d}, {F, 2 d}
        ]]] -
        Total[Flatten[Table[
            lcOmegaLow[[E, C, Dd]] lcOmega[[A, B, E]],
            {E, 2 d}
        ]]] +
        Total[Flatten[Table[
            lcOmegaLow[[E, C, Dd]] lcOmega[[B, A, E]],
            {E, 2 d}
        ]]];

    lastTerm[A_Integer, B_Integer, C_Integer, Dd_Integer] :=
        Total[Flatten[Table[
            lcOmegaLow[[E, A, B]] etaUU[[E, F]]
                lcOmegaLow[[F, C, Dd]],
            {E, 2 d}, {F, 2 d}
        ]]];

    rmOrderedLocal[A_Integer, B_Integer, C_Integer, Dd_Integer] :=
        rmOrderedLocal[A, B, C, Dd] = qSimp[
            1/2 (
                firstBlock[A, B, C, Dd] +
                firstBlock[C, Dd, A, B] -
                lastTerm[A, B, C, Dd]
            )
        ];

    Print["[", name, "] Evaluating all raw curvature components (no imposed symmetries)..."];
    rmArray = Table[rmOrderedLocal[A,B,C,Dd],
        {A,2 d},{B,2 d},{C,2 d},{Dd,2 d}];
    rawChecks = <|
        "FirstPairAntisymmetry" -> ZeroTensorQ[Table[rmArray[[A,B,C,Dd]]+rmArray[[B,A,C,Dd]],{A,2 d},{B,2 d},{C,2 d},{Dd,2 d}]],
        "SecondPairAntisymmetry" -> ZeroTensorQ[Table[rmArray[[A,B,C,Dd]]+rmArray[[A,B,Dd,C]],{A,2 d},{B,2 d},{C,2 d},{Dd,2 d}]],
        "PairExchange" -> ZeroTensorQ[Table[rmArray[[A,B,C,Dd]]-rmArray[[C,Dd,A,B]],{A,2 d},{B,2 d},{C,2 d},{Dd,2 d}]],
        "AlgebraicBianchi" -> ZeroTensorQ[Table[
            rmArray[[A,B,C,Dd]] + rmArray[[B,C,A,Dd]] + rmArray[[C,A,B,Dd]],
            {A,2 d},{B,2 d},{C,2 d},{Dd,2 d}]]
    |>;
    If[!And@@Values[rawChecks],
        Return[Failure["CurvatureSymmetriesNotProved",<|"Frame"->name,"Checks"->rawChecks|>]]];
    rmFun[A_Integer,B_Integer,C_Integer,Dd_Integer] := rmArray[[A,B,C,Dd]];

    (* ------------------------------------------------------------ *)
    (* Ricci, scalar and Kretschmann                                *)
    (* ------------------------------------------------------------ *)

    ricFull = Table[
        qSimp[
            Total[Flatten[Table[
                etaUU[[A, E]] rmFun[E, B, A, Dd],
                {A, 2 d}, {E, 2 d}
            ]]]
        ],
        {B, 2 d}, {Dd, 2 d}
    ];

    rcGen = Table[
        If[
            sectorSign[B, d] != sectorSign[Dd, d],
            ricFull[[B, Dd]],
            0
        ],
        {B, 2 d}, {Dd, 2 d}
    ];

    rScalar = qSimp[
        1/2 Total[Flatten[Table[
            geUU[[A, B]] ricFull[[A, B]],
            {A, 2 d}, {B, 2 d}
        ]]]
    ];

    pairsLocal = Subsets[Range[2 d], {2}];

    pairMetric = Table[
        With[
            {
                a = pairsLocal[[i, 1]],
                b = pairsLocal[[i, 2]],
                c = pairsLocal[[j, 1]],
                dd = pairsLocal[[j, 2]]
            },
            qSimp[
                geUU[[a, c]] geUU[[b, dd]] -
                geUU[[a, dd]] geUU[[b, c]]
            ]
        ],
        {i, Length[pairsLocal]},
        {j, Length[pairsLocal]}
    ];

    Print[
        "[", name, "] Building the ",
        Length[pairsLocal], " x ", Length[pairsLocal],
        " independent-pair Riemann matrix..."
    ];

    rmPair = Table[
        With[
            {p = pairsLocal[[i]], q = pairsLocal[[j]]},
            rmFun[p[[1]], p[[2]], q[[1]], q[[2]]]
        ],
        {i, Length[pairsLocal]},
        {j, Length[pairsLocal]}
    ];

    kScalar = qSimp[
        4 Tr[
            pairMetric .
            rmPair .
            pairMetric .
            Transpose[rmPair]
        ]
    ];


    <|"RmArray"->rmArray,"RawCurvatureChecks"->rawChecks,
      "RicFull"->ricFull,"RcGen"->rcGen,"RScalar"->rScalar,
      "Pairs"->pairsLocal,"PairMetricUU"->pairMetric,"RmPair"->rmPair,"KGen"->kScalar|>
];

ClearAll[RawFrameBracketLow,TorsionFromLow];
RawFrameBracketLow[frame_,x_List,h_] := Module[{d=Length[x],sz=Length[frame],br,vec,one},
    br=Table[
        vec=Table[qSimp[Total[Flatten[Table[frame[[i,A]] D[frame[[j,B]],x[[i]]] -
            frame[[i,B]] D[frame[[j,A]],x[[i]]],{i,d}]]]],{j,d}];
        one=Table[qSimp[Total[Flatten[Table[
            frame[[i,A]] D[frame[[d+j,B]],x[[i]]] +
            frame[[d+i,B]] D[frame[[i,A]],x[[j]]] -
            frame[[i,B]] (D[frame[[d+j,A]],x[[i]]] - D[frame[[d+i,A]],x[[j]]]),{i,d}]]]+
            Total[Flatten[Table[frame[[i,A]] frame[[k,B]] h[[i,k,j]],{i,d},{k,d}]]]],{j,d}];
        Join[vec,one],{A,sz},{B,sz}];
    Table[qSimp[1/2 Total[Flatten[Table[br[[A,B,i]] frame[[d+i,C]]+
        br[[A,B,d+i]] frame[[i,C]],{i,d}]]]],{A,sz},{B,sz},{C,sz}]
];
TorsionFromLow[omega_,bracket_] := Module[{sz=Length[omega]},
    Table[qSimp[omega[[A,B,C]]-omega[[B,A,C]]-bracket[[A,B,C]]+omega[[C,A,B]]],
        {A,sz},{B,sz},{C,sz}]];

(* Finite means proved finite and real under the given parameter assumptions.
   Unknown is not silently accepted. *)
ClearAll[LimitStatus];
LimitStatus[expr_,assumptions_:True] := Which[
    !FreeQ[expr,Indeterminate|ComplexInfinity|DirectedInfinity[_]],"NonFinite",
    !FreeQ[expr,_Limit|_ConditionalExpression|_Piecewise|_Missing|_Failure|_Inactive|Inactive[_][___]|$Aborted],"Inconclusive",
    TrueQ[FullSimplify[Element[expr,Reals],Assumptions->assumptions]],"Finite",
    True,"Inconclusive"
];

ClearAll[HorizonReport];
HorizonReport[result_Association] := Module[
    {rh=result["Horizon"],ass=result["ParameterAssumptions"],ge=result["SpecialisedDualRegular"],
     keys,rows,vals,limitsAbove,limitsBelow,statusA,statusB,agreement,smooth,
     frame,det,detAt,checks},
    keys={"RawFrame","GEDD","EtaDD","Anchor","BornOmega","LCOmega",
          "BornTorsionLow","RicFull","RcGen","RmArray","RScalar","KGen"};
    rows=Table[
        vals=DeleteDuplicates[Flatten[{ge[key]}]];
        limitsAbove=HorizonLimit[#,rh,"FromAbove",ass]& /@ vals;
        limitsBelow=HorizonLimit[#,rh,"FromBelow",ass]& /@ vals;
        statusA=LimitStatus[#,ass]& /@ limitsAbove;
        statusB=LimitStatus[#,ass]& /@ limitsBelow;
        agreement=And@@MapThread[TrueQ[FullSimplify[#1==#2,Assumptions->ass]]&,
            {limitsAbove,limitsBelow}];
        (* For these examples every component is rational in r with smooth
           angular coefficients. A nonzero denominator at the horizon proves
           smoothness locally, not just boundedness along a radial limit. *)
        smooth=And@@(Function[e,TrueQ[FullSimplify[
            (Denominator[Together[e]]/.r->rh)!=0,Assumptions->ass]] &&
            LimitStatus[Numerator[Together[e]]/.r->rh,ass]=="Finite"] /@ vals);
        <|"Object"->key,"DistinctExpressions"->Length[vals],
          "Above"->If[MemberQ[statusA,"NonFinite"],"NonFinite",If[AllTrue[statusA,#=="Finite"&],"Finite","Inconclusive"]],
          "Below"->If[MemberQ[statusB,"NonFinite"],"NonFinite",If[AllTrue[statusB,#=="Finite"&],"Finite","Inconclusive"]],
          "LimitsAgree"->agreement,"SmoothRationalExtensionProved"->smooth|>,{key,keys}];
    frame=ge["RawFrame"];det=FullSimplify[Det[frame],Assumptions->ass];
    detAt=FullSimplify[det/.r->rh,Assumptions->ass];
    checks=<|"HorizonIsRoot"->TrueQ[FullSimplify[(result["fExpression"]/.r->rh)==0,Assumptions->ass]],
        "HorizonIsSimple"->TrueQ[FullSimplify[(D[result["fExpression"],r]/.r->rh)!=0,Assumptions->ass]],
        "FrameNonDegenerateAtHorizon"->TrueQ[FullSimplify[detAt!=0,Assumptions->ass]],
        "GeneralisedMetricNonDegenerateAtHorizon"->TrueQ[FullSimplify[(Det[ge["GEDD"]]/.r->rh)!=0,Assumptions->ass]]|>;
    <|"Example"->result["Example"],"Horizon"->rh,"Checks"->checks,
      "RawFrameDeterminantAtHorizon"->detAt,"Objects"->rows,
      "AllRegularityChecksProved"->(And@@Values[checks] &&
         AllTrue[rows,(#["Above"]=="Finite" && #["Below"]=="Finite" && TrueQ[#["LimitsAgree"]] && TrueQ[#["SmoothRationalExtensionProved"]])&])|>
];

ClearAll[ShowValidation];
ShowValidation[result_Association] := Column[{
    Style[result["Example"]<>" \[LongDash] verification",Bold],
    Grid[Prepend[KeyValueMap[{#1,If[TrueQ[#2],"PASS","NOT PROVED"]}&,result["Checks"]],
        {"Check","Result"}],Frame->All,Alignment->Left],
    Row[{"All checks proved: ",result["AllChecksProved"]}],
    Row[{"Generalised scalar: ",result["SpecialisedDualRegular"]["RScalar"]}],
    Row[{"Generalised Kretschmann: ",result["SpecialisedDualRegular"]["KGen"]}]
}];


(* Independent ordinary curvature computation for the factor-four check.
   The square is unchanged by the overall ordinary Riemann sign convention. *)
ClearAll[OrdinaryKretschmann];
OrdinaryKretschmann[geom_Association] := Module[
    {x=geom["Coordinates"],d=geom["Dimension"],g=geom["gDD"],gi=geom["gUU"],
     ga=geom["GammaLC"],ru,rd,t1,t2,t3,t4},
    ru=Table[qSimp[D[ga[[a,mu,b]],x[[nu]]]-D[ga[[a,nu,b]],x[[mu]]]+
        Total[Flatten[Table[ga[[a,nu,c]] ga[[c,mu,b]]-ga[[a,mu,c]] ga[[c,nu,b]],{c,d}]]]],
        {a,d},{b,d},{mu,d},{nu,d}];
    rd=Table[qSimp[Total[Flatten[Table[g[[a,c]] ru[[c,b,mu,nu]],{c,d}]]]],{a,d},{b,d},{mu,d},{nu,d}];
    t1=Table[qSimp[Total[Flatten[Table[gi[[a,c]] rd[[c,b,mu,nu]],{c,d}]]]],{a,d},{b,d},{mu,d},{nu,d}];
    t2=Table[qSimp[Total[Flatten[Table[gi[[b,c]] t1[[a,c,mu,nu]],{c,d}]]]],{a,d},{b,d},{mu,d},{nu,d}];
    t3=Table[qSimp[Total[Flatten[Table[gi[[mu,c]] t2[[a,b,c,nu]],{c,d}]]]],{a,d},{b,d},{mu,d},{nu,d}];
    t4=Table[qSimp[Total[Flatten[Table[gi[[nu,c]] t3[[a,b,mu,c]],{c,d}]]]],{a,d},{b,d},{mu,d},{nu,d}];
    qSimp[Total[Flatten[Table[rd[[a,b,mu,nu]] t4[[a,b,mu,nu]],{a,d},{b,d},{mu,d},{nu,d}]]]]
];

RunBornExample[choice_String] := Module[
    {exampleF,exampleHorizon,parameterAssumptions,checks,metricDivOriginal,
     metricDivDual,dPhiDual,kTransverse,ordinaryK,result},
    If[!MemberQ[{"Schwarzschild","PlanarRN"},choice],
        Return[Failure["UnknownExample",<|"Allowed"->{"Schwarzschild","PlanarRN"}|>]]];



(* ::Title:: *)
(* Born-Selected Generalised Geometry and its T-Dual *)


(*
   PURPOSE
   -------
   This source implements the purely Born analysis discussed in Chapter 4.

   The geometric input is the ORIGINAL Born geometry:
       - a local isotropic splitting of TM + T*M,
       - an S^1-invariant metric g,
       - a two-form b in that splitting,
       - the raw Courant three-form twist H_twist,
       - a generalised Born structure represented on TM by the
         g-isometry Psi.

   The program then:
       1. constructs the original unique Generalised Born connection D^B;
       2. computes its intrinsic generalised torsion T_{D^B};
       3. constructs the Born-selected Generalised Levi-Civita connection
              D^B_LC = D^B - (1/3) T_{D^B};
       4. performs a local factorised T-duality along the chosen S^1 direction;
       5. extracts the T-dual metric and b-field from the transported
          generalised-metric eigenbundles;
       6. transports the Born structure using
              PsiHat = TMinus . Psi . Inverse[TPlus];
       7. constructs the unique T-dual Generalised Born connection directly
          from (gHat, HhatPreferred, PsiHat);
       8. removes its intrinsic torsion to obtain the T-dual Born-selected
          Generalised Levi-Civita connection;
       9. computes, independently on BOTH sides:
              - full generalised Riemann tensor,
              - full generalised Ricci curvature,
              - mixed generalised Ricci tensor,
              - generalised Ricci scalar,
              - generalised Kretschmann scalar;
      10. provides readable tables of all non-zero connection, torsion,
          Ricci and independent Riemann components.

   IMPORTANT
   ---------
   There is NO "naive dual canonical connection" in this program.
   Both sides are selected purely by their Generalised Born structures.

   SPLITTING
   ---------
   Mathematica needs a concrete local representation.  We therefore work in
   a local isotropic coordinate splitting.  The T-duality map is implemented
   by the standard local factorised O(d,d) transformation that swaps the
   chosen circle vector with the corresponding one-form.  The dual metric and
   b-field are then EXTRACTED from the transported eigenbundles.

   The Born/curvature calculations themselves are carried out in the preferred
   splitting determined by the corresponding generalised metric.  Thus the
   three-form entering the Bismut connections is

       HPreferred = HTwist + db.

   In the default horizon example HTwist = 0 on both sides, b = 0 originally,
   and the T-dual b-field is locally pure gauge.

   COMPONENT OUTPUT
   ----------------
   All display routines now:
       - suppress zero components completely;
       - remove exact duplicate rows;
       - use only independent components where a tensor symmetry is known;
       - print "all components are zero" instead of an empty Grid.

   CURRENT SCOPE
   -------------
   The automatic T-duality layer is designed for the local d=1 Buscher setting
   used in the horizon analysis.  It assumes the raw Courant twist is zero on
   the local patch.  Non-zero raw twists require a separate bracket-intertwining duality construction
   and are rejected by this driver. They cannot be inserted independently.
*)

$HistoryLength = 0; (* Do not clear the user's Global context. *)
$HistoryLength = 0;

(* ================================================================ *)
(* 0. User controls                                                 *)
(* ================================================================ *)

coords = {v, r, th, ph};
n = Length[coords];
genDim = 2 n;

(* The S^1 / T-duality direction.  In the EF examples this is v. *)
dualIndex = 1;

(* Choose the transverse geometry used by the default example. *)
transverseModel = If[choice == "Schwarzschild", "Spherical", "Flat2"];

transverseBlock = Switch[
    transverseModel,
    "Spherical", DiagonalMatrix[{r^2, r^2 Sin[th]^2}],
    "Flat2",     DiagonalMatrix[{r^2, r^2}],
    _, (
        Print["Unknown transverseModel; using Spherical."];
        DiagonalMatrix[{r^2, r^2 Sin[th]^2}]
    )
];

(* The symbolic engine works on the open patch where the Buscher
   parametrisation is invertible.  Horizon limits are taken only later. *)
bornAssumptions = r > 0 && 0 < th < Pi && f[r] != 0;

(* This is deliberately the simplifier used in the current thesis notebook. *)
ClearAll[qSimp];
qSimp[expr_] := FullSimplify[
    TrigReduce[expr],
    Assumptions -> bornAssumptions
];

ClearAll[qMap];
qMap[array_, level_Integer] := Map[qSimp, array, {level}];

(* ================================================================ *)
(* 1. Original Born-geometry input                                  *)
(* ================================================================ *)

(*
   DEFAULT: the generic Eddington-Finkelstein horizon geometry

       ds^2 = -f(r) dv^2 + 2 dv dr + r^2 dOmega_2^2.

   The original b-field and raw Courant twist vanish.

   The tangent representative PsiOriginal of the Born structure satisfies

       PsiOriginal^T . g . PsiOriginal = g.

   For the anchor-induced Born structure, PsiOriginal = IdentityMatrix[n].

   Replace the objects in this section to study another invariant Born
   geometry.  Everything below is generated from them.
*)

gOriginalDD = ArrayFlatten[{
    {{{-f[r], 1}, {1, 0}}, ConstantArray[0, {2, 2}]},
    {ConstantArray[0, {2, 2}], transverseBlock}
}];

bOriginalDD = ConstantArray[0, {n, n}];

(* Raw twist in the chosen isotropic splitting. *)
HTwistOriginalDDD = ConstantArray[0, {n, n, n}];

(* Tangent-bundle representative of the original Generalised Born structure. *)
PsiOriginal = IdentityMatrix[n];

Print["Original metric g ="];
Print[MatrixForm[gOriginalDD]];
Print["Original b-field ="];
Print[MatrixForm[bOriginalDD]];
Print["Original Born tangent representative Psi ="];
Print[MatrixForm[PsiOriginal]];


(* ================================================================ *)
(* 2. Basic tensor utilities                                        *)
(* ================================================================ *)

ClearAll[ExteriorDerivative2Form];
ExteriorDerivative2Form[b_, x_List] := Module[{d = Length[x]},
    Table[
        qSimp[
            D[b[[nu, rho]], x[[mu]]] +
            D[b[[rho, mu]], x[[nu]]] +
            D[b[[mu, nu]], x[[rho]]]
        ],
        {mu, d}, {nu, d}, {rho, d}
    ]
];

ClearAll[ChristoffelFromMetric];
ChristoffelFromMetric[gdd_, guu_, x_List] := Module[{d = Length[x]},
    Table[
        qSimp[
            1/2 Total[Flatten[Table[
                guu[[rho, sig]] (
                    D[gdd[[sig, nu]], x[[mu]]] +
                    D[gdd[[sig, mu]], x[[nu]]] -
                    D[gdd[[mu, nu]], x[[sig]]]
                ),
                {sig, d}
            ]]]
        ],
        {rho, d}, {mu, d}, {nu, d}
    ]
];

ClearAll[RaiseLast3FormIndex];
RaiseLast3FormIndex[hddd_, guu_] := Module[{d = Length[guu]},
    Table[
        qSimp[
            Total[Flatten[Table[
                guu[[rho, sig]] hddd[[mu, nu, sig]],
                {sig, d}
            ]]]
        ],
        {rho, d}, {mu, d}, {nu, d}
    ]
];

ClearAll[sectorSign, spacetimeIndex];
sectorSign[A_Integer, d_Integer] := If[A <= d, 1, -1];
spacetimeIndex[A_Integer, d_Integer] := 1 + Mod[A - 1, d];

ClearAll[ZeroTensorQ];
ZeroTensorQ[array_] := And @@ (
    TrueQ[qSimp[#] == 0] & /@ Flatten[array]
);

ClearAll[nonZeroQ];
nonZeroQ[expr_] := ! TrueQ[qSimp[expr] == 0];


(* ================================================================ *)
(* 3. Local factorised T-duality                                    *)
(* ================================================================ *)

(*
   We use the local coordinate splitting TM + T*M.

   The raw generalised-metric eigenbundle lifts are

       j_+(X) = X + (b+g)(X),
       j_-(X) = X + (b-g)(X).

   TGen is the factorised O(d,d) transformation which swaps
   partial_y <-> dy along the chosen duality coordinate y.

   If
       TGen . j_+ = jHat_+ . TPlus,
       TGen . j_- = jHat_- . TMinus,

   then the vector blocks of the transported lifts are precisely TPlus and
   TMinus.  The dual metric and b-field are reconstructed from the one-form
   blocks, while the transported Born representative is

       PsiHat = TMinus . Psi . Inverse[TPlus].

   This is the matrix implementation of the Chapter-4 tangent-bundle
   transport law.
*)

ClearAll[BuildFactorisedTDual];
BuildFactorisedTDual[gdd_, bdd_, psi_, k_Integer] := Module[
    {
        d, tgen, jPlus, jMinus, jHatPlus, jHatMinus,
        tPlus, tMinus, qPlus, qMinus, tPlusInv, tMinusInv,
        mPlusT, mMinusT, mPlus, mMinus, gHat, bHat, psiHat,
        tFrame
    },

    d = Length[gdd];

    tgen = IdentityMatrix[2 d];
    tgen[[{k, d + k}]] = tgen[[{d + k, k}]];

    (* Columns are the lifts of the coordinate tangent basis. *)
    jPlus  = Join[IdentityMatrix[d], Transpose[bdd + gdd]];
    jMinus = Join[IdentityMatrix[d], Transpose[bdd - gdd]];

    jHatPlus  = qMap[tgen . jPlus, 2];
    jHatMinus = qMap[tgen . jMinus, 2];

    tPlus  = jHatPlus[[1 ;; d, All]];
    tMinus = jHatMinus[[1 ;; d, All]];
    qPlus  = jHatPlus[[d + 1 ;; 2 d, All]];
    qMinus = jHatMinus[[d + 1 ;; 2 d, All]];

    tPlusInv  = qMap[Inverse[tPlus], 2];
    tMinusInv = qMap[Inverse[tMinus], 2];

    (* qPlus = Transpose[bHat+gHat] . TPlus, and similarly for -. *)
    mPlusT  = qMap[qPlus . tPlusInv, 2];
    mMinusT = qMap[qMinus . tMinusInv, 2];

    mPlus  = Transpose[mPlusT];
    mMinus = Transpose[mMinusT];

    gHat = qMap[(mPlus - mMinus)/2, 2];
    bHat = qMap[(mPlus + mMinus)/2, 2];

    psiHat = qMap[tMinus . psi . tPlusInv, 2];

    tFrame = ArrayFlatten[{
        {tPlus, ConstantArray[0, {d, d}]},
        {ConstantArray[0, {d, d}], tMinus}
    }];

    <|
        "TGen" -> tgen,
        "JPlusOriginalRaw" -> jPlus,
        "JMinusOriginalRaw" -> jMinus,
        "JPlusTDualRaw" -> jHatPlus,
        "JMinusTDualRaw" -> jHatMinus,
        "TPlus" -> tPlus,
        "TMinus" -> tMinus,
        "TFrame" -> tFrame,
        "gHatDD" -> gHat,
        "bHatDD" -> bHat,
        "PsiHat" -> psiHat
    |>
];

If[!ZeroTensorQ[HTwistOriginalDDD] ||
   !ZeroTensorQ[D[{gOriginalDD,bOriginalDD,PsiOriginal},coords[[dualIndex]]]] ||
   !ZeroTensorQ[Transpose[PsiOriginal] . gOriginalDD . PsiOriginal-gOriginalDD],
   Return[Failure["UnsupportedInput",<|"Reason"->"Requires invariant data, zero raw twist and a metric-isometric Psi."|>]]];

TDualData = BuildFactorisedTDual[
    gOriginalDD,
    bOriginalDD,
    PsiOriginal,
    dualIndex
];

gTDualDD = TDualData["gHatDD"];
bTDualDD = TDualData["bHatDD"];
PsiTDual = TDualData["PsiHat"];
TPlus = TDualData["TPlus"];
TMinus = TDualData["TMinus"];
TFrame = TDualData["TFrame"];

(* In the local Buscher examples used in the horizon analysis the raw
   Courant twist is zero on both patches.  The preferred flux on the dual
   side is therefore dbHat. Non-zero raw twists are outside this driver. *)
HTwistTDualDDD = ConstantArray[0, {n, n, n}];

Print["T-dual metric gHat ="];
Print[MatrixForm[gTDualDD]];
Print["T-dual b-field bHat ="];
Print[MatrixForm[bTDualDD]];
Print["TPlus ="];
Print[MatrixForm[TPlus]];
Print["TMinus ="];
Print[MatrixForm[TMinus]];
Print["Transported Born representative PsiHat ="];
Print[MatrixForm[PsiTDual]];

Print[
    "Original Born-isometry check Psi^T g Psi - g = 0: ",
    ZeroTensorQ[
        qMap[Transpose[PsiOriginal] . gOriginalDD . PsiOriginal - gOriginalDD, 2]
    ]
];

Print[
    "T-dual Born-isometry check PsiHat^T gHat PsiHat - gHat = 0: ",
    ZeroTensorQ[
        qMap[Transpose[PsiTDual] . gTDualDD . PsiTDual - gTDualDD, 2]
    ]
];


(* ================================================================ *)
(* 4. Optional KK diagnostics                                       *)
(* ================================================================ *)

(*
   These quantities are not needed by the core engine, but they make contact
   directly with the Chapter-4 Kaluza-Klein notation

       g = g0 kappa^2 + h,
       A = g1/g0,
       kappa = theta + A.

   Here theta is the coordinate one-form along dualIndex.
*)

ClearAll[BuildKKData];
BuildKKData[gdd_, k_Integer] := Module[
    {d, g0, g1, g2, a, h, theta},
    d = Length[gdd];
    theta = UnitVector[d, k];
    g0 = qSimp[gdd[[k, k]]];
    g1 = Table[If[i == k, 0, qSimp[gdd[[k, i]]]], {i, d}];
    g2 = Table[
        If[i == k || j == k, 0, qSimp[gdd[[i, j]]]],
        {i, d}, {j, d}
    ];
    a = Map[qSimp, g1/g0];
    h = qMap[g2 - Outer[Times, g1, g1]/g0, 2];
    <|
        "g0" -> g0,
        "g1D" -> g1,
        "g2DD" -> g2,
        "AD" -> a,
        "thetaD" -> theta,
        "kappaD" -> Map[qSimp, theta + a],
        "hDD" -> h
    |>
];

KKOriginal = BuildKKData[gOriginalDD, dualIndex];

Print["KK g0 = ", KKOriginal["g0"]];
Print["KK g1 = ", KKOriginal["g1D"]];
Print["KK A = g1/g0 = ", KKOriginal["AD"]];
Print["KK h ="];
Print[MatrixForm[KKOriginal["hDD"]]];


(* ================================================================ *)
(* 5. Generalised Born connection and Born-selected LC engine       *)
(* ================================================================ *)

(*
   Mathematics implemented in this section.

   Let Psi be the tangent representative of the Generalised Born structure and

       nabla^+_X Y = nabla_X Y + (1/2) g^{-1} H(X,Y,.),
       nabla^-_X Y = nabla_X Y - (1/2) g^{-1} H(X,Y,.).

   Proposition 4.13 gives the pure-type Born components

       pi_+(D^B_{u+} v+)
          = Psi^{-1} nabla^-_X(Psi Y),

       pi_-(D^B_{u-} v-)
          = Psi nabla^+_X(Psi^{-1} Y),

   while the mixed components are

       pi_+(D^B_{u-} v+) = nabla^+_X Y,
       pi_-(D^B_{u+} v-) = nabla^-_X Y.

   The intrinsic Generalised Born torsion is computed from its DEFINITION:

       T_D(u,v,w)
          = <D_u v,w> - <D_v u,w>
            - <[u,v]_H,w> + <v,D_w u>.

   We then construct

       D^B_LC = D^B - (1/3) T_{D^B}.

   No canonical Chapter-3 dual connection is introduced anywhere.
*)

ClearAll[BuildBornGeometry];
BuildBornGeometry[
    name_String,
    x_List,
    gRaw_,
    bRaw_,
    hTwistRaw_,
    psiRaw_
] := Module[
    {
        d, gd, bd, gu, db, hPref, hUp, gamma,
        psi, psiInv, gammaPlus, gammaMinus,
        etaDD, etaUU, geDD, geUU, anchor,
        bornOmega, bornOmegaLow,
        bracketLow, torsionBornLow, torsionBornVec,
        lcOmega, lcOmegaLow, torsionLCLow,
        courantResidualBorn, metricResidualBorn,
        courantResidualLC, metricResidualLC,
        anchorD, firstBlock, lastTerm, rmOrderedLocal, rmFun,
        ricFull, rcGen, rScalar,
        pairsLocal, pairMetric, rmPair, kScalar,
        rmArray, bornTorsionPureQ, lcTorsionFreeQ,
        bornCourantQ, bornMetricQ, lcCourantQ, lcMetricQ,
        jPlusRaw, jMinusRaw, curvResult, geometryChecks, bornK, bornKResidual
    },

    d = Length[x];
    gd = qMap[gRaw, 2];
    bd = qMap[bRaw, 2];
    gu = qMap[Inverse[gd], 2];

    db = ExteriorDerivative2Form[bd, x];
    hPref = qMap[hTwistRaw + db, 3];
    hUp = RaiseLast3FormIndex[hPref, gu];

    gamma = ChristoffelFromMetric[gd, gu, x];

    psi = qMap[psiRaw, 2];
    psiInv = qMap[Inverse[psi], 2];

    gammaPlus = Table[
        qSimp[gamma[[rho, mu, nu]] + 1/2 hUp[[rho, mu, nu]]],
        {rho, d}, {mu, d}, {nu, d}
    ];

    gammaMinus = Table[
        qSimp[gamma[[rho, mu, nu]] - 1/2 hUp[[rho, mu, nu]]],
        {rho, d}, {mu, d}, {nu, d}
    ];

    etaDD = ArrayFlatten[{
        {gd, ConstantArray[0, {d, d}]},
        {ConstantArray[0, {d, d}], -gd}
    }];

    etaUU = ArrayFlatten[{
        {gu, ConstantArray[0, {d, d}]},
        {ConstantArray[0, {d, d}], -gu}
    }];

    geDD = ArrayFlatten[{
        {gd, ConstantArray[0, {d, d}]},
        {ConstantArray[0, {d, d}], gd}
    }];

    geUU = ArrayFlatten[{
        {gu, ConstantArray[0, {d, d}]},
        {ConstantArray[0, {d, d}], gu}
    }];

    anchor = Join[IdentityMatrix[d], IdentityMatrix[d]];

    (* Raw eigenbundle lifts retained for inspection. *)
    jPlusRaw = Join[IdentityMatrix[d], Transpose[bd + gd]];
    jMinusRaw = Join[IdentityMatrix[d], Transpose[bd - gd]];

    (* ------------------------------------------------------------ *)
    (* Unique Generalised Born connection                           *)
    (* ------------------------------------------------------------ *)

    bornOmega = Table[
        Module[{sA, sB, sC, mu, nu, rho},
            sA = sectorSign[A, d];
            sB = sectorSign[B, d];
            sC = sectorSign[C, d];
            mu = spacetimeIndex[A, d];
            nu = spacetimeIndex[B, d];
            rho = spacetimeIndex[C, d];

            If[
                sC != sB,
                0,

                Which[
                    (* Pure + : Psi^{-1} nabla^-_mu(Psi e_nu). *)
                    sA == 1 && sB == 1,
                    qSimp[
                        Total[Flatten[Table[
                            psiInv[[rho, sig]] (
                                D[psi[[sig, nu]], x[[mu]]] +
                                Total[Flatten[Table[
                                    gammaMinus[[sig, mu, lam]] psi[[lam, nu]],
                                    {lam, d}
                                ]]]
                            ),
                            {sig, d}
                        ]]]
                    ],

                    (* Pure - : Psi nabla^+_mu(Psi^{-1} e_nu). *)
                    sA == -1 && sB == -1,
                    qSimp[
                        Total[Flatten[Table[
                            psi[[rho, sig]] (
                                D[psiInv[[sig, nu]], x[[mu]]] +
                                Total[Flatten[Table[
                                    gammaPlus[[sig, mu, lam]] psiInv[[lam, nu]],
                                    {lam, d}
                                ]]]
                            ),
                            {sig, d}
                        ]]]
                    ],

                    (* Mixed: output + and first argument -. *)
                    sA == -1 && sB == 1,
                    gammaPlus[[rho, mu, nu]],

                    (* Mixed: output - and first argument +. *)
                    sA == 1 && sB == -1,
                    gammaMinus[[rho, mu, nu]],

                    True,
                    0
                ]
            ]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    bornOmegaLow = Table[
        qSimp[
            Total[Flatten[Table[
                bornOmega[[A, B, E]] etaDD[[E, C]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    (* ------------------------------------------------------------ *)
    (* Dorfman bracket in the preferred splitting                   *)
    (* ------------------------------------------------------------ *)

    (*
       The preferred eigenframe is

           e_{mu,+} = partial_mu + g_{mu nu} dx^nu,
           e_{mu,-} = partial_mu - g_{mu nu} dx^nu.

       For coordinate vector fields the vector part of their Dorfman bracket
       vanishes.  Its one-form part is evaluated directly from

           [X+xi,Y+eta]_H
             = [X,Y] + L_X eta - i_Y dxi + i_Y i_X H.

       Pairing that pure one-form with e_{rho,+/-} gives one half of its
       rho-component because the Courant pairing carries the factor 1/2.
    *)

    bracketLow = Table[
        Module[{sA, sB, mu, nu, rho, qOne},
            sA = sectorSign[A, d];
            sB = sectorSign[B, d];
            mu = spacetimeIndex[A, d];
            nu = spacetimeIndex[B, d];
            rho = spacetimeIndex[C, d];

            qOne = Table[
                qSimp[
                    sB D[gd[[nu, j]], x[[mu]]] -
                    sA (
                        D[gd[[mu, j]], x[[nu]]] -
                        D[gd[[mu, nu]], x[[j]]]
                    ) +
                    hPref[[mu, nu, j]]
                ],
                {j, d}
            ];

            qSimp[1/2 qOne[[rho]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    (* ------------------------------------------------------------ *)
    (* Intrinsic Born torsion and Born-selected LC connection       *)
    (* ------------------------------------------------------------ *)

    torsionBornLow = Table[
        qSimp[
            bornOmegaLow[[A, B, C]] -
            bornOmegaLow[[B, A, C]] -
            bracketLow[[A, B, C]] +
            bornOmegaLow[[C, A, B]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    torsionBornVec = Table[
        qSimp[
            Total[Flatten[Table[
                etaUU[[C, Dd]] torsionBornLow[[A, B, Dd]],
                {Dd, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    lcOmega = Table[
        qSimp[
            bornOmega[[A, B, C]] -
            1/3 torsionBornVec[[A, B, C]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    lcOmegaLow = Table[
        qSimp[
            Total[Flatten[Table[
                lcOmega[[A, B, E]] etaDD[[E, C]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    torsionLCLow = Table[
        qSimp[
            lcOmegaLow[[A, B, C]] -
            lcOmegaLow[[B, A, C]] -
            bracketLow[[A, B, C]] +
            lcOmegaLow[[C, A, B]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    bornTorsionPureQ = And @@ Flatten[
        Table[
            If[
                sectorSign[A, d] == sectorSign[B, d] ==
                    sectorSign[C, d],
                True,
                TrueQ[qSimp[torsionBornLow[[A, B, C]]] == 0]
            ],
            {A, 2 d}, {B, 2 d}, {C, 2 d}
        ]
    ];

    lcTorsionFreeQ = And @@ Flatten[
        Map[TrueQ[qSimp[#] == 0] &, torsionLCLow, {3}]
    ];

    (* ------------------------------------------------------------ *)
    (* Compatibility checks                                        *)
    (* ------------------------------------------------------------ *)

    Clear[anchorD];
    anchorD[A_Integer, expr_] := Total[Flatten[Table[
        anchor[[A, mu]] D[expr, x[[mu]]],
        {mu, d}
    ]]];

    courantResidualBorn = Table[
        qSimp[
            anchorD[A, etaDD[[B, C]]] -
            Total[Flatten[Table[
                bornOmega[[A, B, E]] etaDD[[E, C]] +
                bornOmega[[A, C, E]] etaDD[[B, E]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    metricResidualBorn = Table[
        qSimp[
            anchorD[A, geDD[[B, C]]] -
            Total[Flatten[Table[
                bornOmega[[A, B, E]] geDD[[E, C]] +
                bornOmega[[A, C, E]] geDD[[B, E]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    courantResidualLC = Table[
        qSimp[
            anchorD[A, etaDD[[B, C]]] -
            Total[Flatten[Table[
                lcOmega[[A, B, E]] etaDD[[E, C]] +
                lcOmega[[A, C, E]] etaDD[[B, E]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    metricResidualLC = Table[
        qSimp[
            anchorD[A, geDD[[B, C]]] -
            Total[Flatten[Table[
                lcOmega[[A, B, E]] geDD[[E, C]] +
                lcOmega[[A, C, E]] geDD[[B, E]],
                {E, 2 d}
            ]]]
        ],
        {A, 2 d}, {B, 2 d}, {C, 2 d}
    ];

    bornCourantQ = And @@ Flatten[
        Map[TrueQ[qSimp[#] == 0] &, courantResidualBorn, {3}]
    ];
    bornMetricQ = And @@ Flatten[
        Map[TrueQ[qSimp[#] == 0] &, metricResidualBorn, {3}]
    ];
    lcCourantQ = And @@ Flatten[
        Map[TrueQ[qSimp[#] == 0] &, courantResidualLC, {3}]
    ];
    lcMetricQ = And @@ Flatten[
        Map[TrueQ[qSimp[#] == 0] &, metricResidualLC, {3}]
    ];

    (* ------------------------------------------------------------ *)
    (* Generalised Riemann tensor of the Born-selected LC connection *)
    (* ------------------------------------------------------------ *)

    geometryChecks = <|
        "BornCourant"->bornCourantQ,"BornMetric"->bornMetricQ,
        "BornTorsionPure"->bornTorsionPureQ,"LCCourant"->lcCourantQ,
        "LCMetric"->lcMetricQ,"LCTorsionFree"->lcTorsionFreeQ,
        "BornTorsionSkew12"->ZeroTensorQ[torsionBornLow+Transpose[torsionBornLow,{2,1,3}]],
        "BornTorsionSkew23"->ZeroTensorQ[torsionBornLow+Transpose[torsionBornLow,{1,3,2}]]
    |>;
    bornK=ArrayFlatten[{{ConstantArray[0,{d,d}],psiInv},{psi,ConstantArray[0,{d,d}]}}];
    bornKResidual=Table[qSimp[
        anchorD[A,bornK[[C,B]]] + Total[Flatten[Table[
            bornOmega[[A,E,C]] bornK[[E,B]]-bornK[[C,E]] bornOmega[[A,B,E]],{E,2 d}]]]],
        {A,2 d},{B,2 d},{C,2 d}];
    AssociateTo[geometryChecks,"BornExchangeMapCompatible"->ZeroTensorQ[bornKResidual]];
    If[!And@@Values[geometryChecks],Return[Failure["ConnectionChecksNotProved",geometryChecks]]];
    curvResult=CurvatureFromFrame[name,x,lcOmega,etaDD,geDD,anchor];
    If[FailureQ[curvResult],Return[curvResult]];

    Join[<|
        "Name" -> name,
        "Coordinates" -> x,
        "Dimension" -> d,
        "gDD" -> gd,
        "gUU" -> gu,
        "bDD" -> bd,
        "HTwistDDD" -> hTwistRaw,
        "dbDDD" -> db,
        "HPreferredDDD" -> hPref,
        "HPreferredUp" -> hUp,
        "GammaLC" -> gamma,
        "Psi" -> psi,
        "PsiInverse" -> psiInv,
        "JPlusRaw" -> jPlusRaw,
        "JMinusRaw" -> jMinusRaw,

        "EtaDD" -> etaDD,
        "EtaUU" -> etaUU,
        "GEDD" -> geDD,
        "GEUU" -> geUU,
        "Anchor" -> anchor,

        "BornOmega" -> bornOmega,
        "BornOmegaLow" -> bornOmegaLow,
        "BracketLow" -> bracketLow,
        "BornTorsionLow" -> torsionBornLow,
        "BornTorsionVector" -> torsionBornVec,

        "LCOmega" -> lcOmega,
        "LCOmegaLow" -> lcOmegaLow,
        "LCTorsionLow" -> torsionLCLow,

        "BornTorsionPureQ" -> bornTorsionPureQ,
        "LCTorsionFreeQ" -> lcTorsionFreeQ,
        "BornCourantCompatibleQ" -> bornCourantQ,
        "BornMetricCompatibleQ" -> bornMetricQ,
        "LCCourantCompatibleQ" -> lcCourantQ,
        "LCMetricCompatibleQ" -> lcMetricQ,

        "BornK" -> bornK,
        "ConnectionChecks" -> geometryChecks,
        "BornDivergenceCoefficients" -> Table[qSimp[Total[Flatten[Table[bornOmega[[A,B,A]],{A,2 d}]]]],{B,2 d}],
        "LCDivergenceCoefficients" -> Table[qSimp[Total[Flatten[Table[lcOmega[[A,B,A]],{A,2 d}]]]],{B,2 d}]
    |>,curvResult]
];


(* ================================================================ *)
(* 6. Build ORIGINAL and T-DUAL Born geometries                     *)
(* ================================================================ *)

Print[""];
Print["============================================================"];
Print["Building ORIGINAL Born-selected geometry"];
Print["============================================================"];

OriginalGeometry = BuildBornGeometry[
    "Original",
    coords,
    gOriginalDD,
    bOriginalDD,
    HTwistOriginalDDD,
    PsiOriginal
];
If[FailureQ[OriginalGeometry],Return[OriginalGeometry]];


Print[""];
Print["============================================================"];
Print["Building T-DUAL transported-Born geometry"];
Print["============================================================"];

TDualGeometry = BuildBornGeometry[
    "T-Dual",
    coords,
    gTDualDD,
    bTDualDD,
    HTwistTDualDDD,
    PsiTDual
];
If[FailureQ[TDualGeometry],Return[TDualGeometry]];



(* ================================================================ *)
(* 7. Summary and T-duality scalar checks                           *)
(* ================================================================ *)

ClearAll[PrintGeometrySummary];
PrintGeometrySummary[geom_Association] := Module[{nm = geom["Name"]},
    Print[""];
    Print["---------------- ",nm," ----------------"];
    If[KeyExistsQ[geom,"ConnectionChecks"],Print[geom["ConnectionChecks"]]];
    If[KeyExistsQ[geom,"RawCurvatureChecks"],Print[geom["RawCurvatureChecks"]]];
    Print["Generalised Ricci scalar = "];
    Print[geom["RScalar"]];
    Print["Generalised Kretschmann scalar = "];
    Print[geom["KGen"]];
];

PrintGeometrySummary[OriginalGeometry];
PrintGeometrySummary[TDualGeometry];

RScalarDifference = qSimp[
    TDualGeometry["RScalar"] -
    OriginalGeometry["RScalar"]
];

KGenDifference = qSimp[
    TDualGeometry["KGen"] -
    OriginalGeometry["KGen"]
];

Print[""];
Print["T-duality scalar comparison"];
Print["R_TDual - R_Original = ", RScalarDifference];
Print["K_TDual - K_Original = ", KGenDifference];


(* ================================================================ *)
(* 8. Readable component-display utilities                          *)
(* ================================================================ *)

ClearAll[coordName, basisName];
coordName[geom_Association, i_Integer] :=
    ToString[geom["Coordinates"][[i]], InputForm];

basisName[geom_Association, A_Integer] := Module[
    {d = geom["Dimension"]},

    If[
        KeyExistsQ[geom, "BasisNames"],
        Return[geom["BasisNames"][[A]]]
    ];

    "e[" <>
    coordName[geom, spacetimeIndex[A, d]] <>
    "," <>
    If[sectorSign[A, d] == 1, "+", "-"] <>
    "]"
];


(* ---------------------------------------------------------------- *)
(* A single robust printer for all component tables.                *)
(*                                                                  *)
(* This deliberately does NOT use Flatten on nested Table output.   *)
(* The previous version could leave empty nested lists behind,      *)
(* producing large blank regions in Grid.                           *)
(* ---------------------------------------------------------------- *)

ClearAll[PrintComponentTable];
PrintComponentTable[
    title_String,
    headers_List,
    rows_List
] := Module[
    {cleanRows},

    cleanRows = DeleteDuplicates[
        Select[
            rows,
            ListQ[#] && Length[#] == Length[headers] &
        ]
    ];

    If[
        Length[cleanRows] == 0,
        Print[title, ": all components are zero."],
        (
            Print[
                title,
                ": ",
                Length[cleanRows],
                " non-zero unique components."
            ];

            Print@Grid[
                Prepend[cleanRows, headers],
                Frame -> All,
                Alignment -> Left,
                Spacings -> {1, 0.6}
            ];
        )
    ];
];


(* ---------------------------------------------------------------- *)
(* Connection coefficients                                         *)
(*                                                                  *)
(* A connection has no general symmetry in its three displayed      *)
(* indices, so each ordered non-zero coefficient is genuinely       *)
(* distinct.                                                        *)
(* ---------------------------------------------------------------- *)

ClearAll[ShowNonZeroConnection];
ShowNonZeroConnection[
    geom_Association,
    which_String : "LC"
] := Module[
    {omega, d, harvested, rows},

    d = geom["Dimension"];

    omega = Switch[
        which,
        "Born", geom["BornOmega"],
        "LC", geom["LCOmega"],
        _,
        Print[
            "Unknown connection type. Use \"Born\" or \"LC\"."
        ];
        Return[$Failed]
    ];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[omega[[A, B, C]]]},
                If[
                    nonZeroQ[val],
                    Sow[
                        {
                            basisName[geom, A],
                            basisName[geom, B],
                            basisName[geom, C],
                            val
                        }
                    ]
                ]
            ],
            {A, 2 d}, {B, 2 d}, {C, 2 d}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        geom["Name"] <> " " <> which <> " connection",
        {
            "first argument",
            "second argument",
            "output basis vector",
            "coefficient"
        },
        rows
    ];
];


(* ---------------------------------------------------------------- *)
(* Intrinsic Born torsion                                           *)
(*                                                                  *)
(* T_DB is a three-form.  Therefore only A < B < C are printed.     *)
(* This removes all permutations/sign duplicates.                   *)
(* ---------------------------------------------------------------- *)

ClearAll[ShowNonZeroBornTorsion];
ShowNonZeroBornTorsion[
    geom_Association
] := Module[
    {t, d, harvested, rows},

    t = geom["BornTorsionLow"];
    d = geom["Dimension"];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[t[[A, B, C]]]},
                If[
                    nonZeroQ[val],
                    Sow[
                        {
                            basisName[geom, A],
                            basisName[geom, B],
                            basisName[geom, C],
                            val
                        }
                    ]
                ]
            ],
            {A, 1, 2 d - 2},
            {B, A + 1, 2 d - 1},
            {C, B + 1, 2 d}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        geom["Name"] <> " intrinsic Born torsion",
        {"A", "B", "C", "T_DB(A,B,C)"},
        rows
    ];
];


(* ---------------------------------------------------------------- *)
(* Full Ricci curvature                                             *)
(*                                                                  *)
(* No symmetry is imposed here: only exact duplicate rows are       *)
(* removed.                                                         *)
(* ---------------------------------------------------------------- *)

ClearAll[ShowNonZeroRicci];
ShowNonZeroRicci[
    geom_Association
] := Module[
    {ric, d, harvested, rows},

    ric = geom["RicFull"];
    d = geom["Dimension"];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[ric[[A, B]]]},
                If[
                    nonZeroQ[val],
                    Sow[
                        {
                            basisName[geom, A],
                            basisName[geom, B],
                            val,
                            If[
                                sectorSign[A, d] !=
                                    sectorSign[B, d],
                                "generalised Ricci Rc",
                                "pure Ricci block"
                            ]
                        }
                    ]
                ]
            ],
            {A, 2 d}, {B, 2 d}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        geom["Name"] <>
            " full generalised Ricci curvature",
        {"A", "B", "Ric^D(A,B)", "block"},
        rows
    ];
];


(* ---------------------------------------------------------------- *)
(* Mixed Generalised Ricci tensor                                   *)
(*                                                                  *)
(* Only the non-zero mixed entries already stored in RcGen are      *)
(* printed.  We do not impose an additional symmetry which need not *)
(* hold in general.                                                  *)
(* ---------------------------------------------------------------- *)

ClearAll[ShowNonZeroGeneralisedRicci];
ShowNonZeroGeneralisedRicci[
    geom_Association
] := Module[
    {rc, d, harvested, rows},

    rc = geom["RcGen"];
    d = geom["Dimension"];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[rc[[A, B]]]},
                If[
                    nonZeroQ[val],
                    Sow[
                        {
                            basisName[geom, A],
                            basisName[geom, B],
                            val
                        }
                    ]
                ]
            ],
            {A, 2 d}, {B, 2 d}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        geom["Name"] <>
            " mixed generalised Ricci tensor",
        {"A", "B", "Rc(A,B)"},
        rows
    ];
];


(* ---------------------------------------------------------------- *)
(* Independent Generalised Riemann components                       *)
(*                                                                  *)
(* The Riemann tensor is represented on antisymmetric index pairs.  *)
(* We use i <= j, so antisymmetry in each pair and pair-exchange     *)
(* symmetry are already quotiented out.                             *)
(* ---------------------------------------------------------------- *)

ClearAll[ShowNonZeroIndependentRiemann];
ShowNonZeroIndependentRiemann[
    geom_Association
] := Module[
    {pairsLocal, rmPair, harvested, rows, p, q},

    pairsLocal = geom["Pairs"];
    rmPair = geom["RmPair"];

    harvested = Reap[
        Do[
            p = pairsLocal[[i]];
            q = pairsLocal[[j]];

            With[
                {val = qSimp[rmPair[[i, j]]]},
                If[
                    nonZeroQ[val],
                    Sow[
                        {
                            basisName[geom, p[[1]]] <>
                                " wedge " <>
                                basisName[geom, p[[2]]],

                            basisName[geom, q[[1]]] <>
                                " wedge " <>
                                basisName[geom, q[[2]]],

                            val
                        }
                    ]
                ]
            ],
            {i, Length[pairsLocal]},
            {j, i, Length[pairsLocal]}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        geom["Name"] <>
            " generalised Riemann tensor",
        {
            "first antisymmetric pair",
            "second antisymmetric pair",
            "Rm"
        },
        rows
    ];
];


ClearAll[ShowAllNonZeroComponents];
ShowAllNonZeroComponents[
    geom_Association
] := (
    ShowNonZeroConnection[geom, "Born"];
    ShowNonZeroBornTorsion[geom];
    ShowNonZeroConnection[geom, "LC"];
    ShowNonZeroRicci[geom];
    ShowNonZeroGeneralisedRicci[geom];
    ShowNonZeroIndependentRiemann[geom];
);



(* ================================================================ *)
(* 9. Side-by-side component counts                                 *)
(* ================================================================ *)

ClearAll[CountNonZero];
CountNonZero[array_] := Count[
    Flatten[array],
    x_ /; nonZeroQ[x]
];

ClearAll[PrintComponentCounts];
PrintComponentCounts[] := Module[{},
    Print[""];
    Print["================ COMPONENT COUNTS ================"];

    Print[
        "Original Born connection: ",
        CountNonZero[OriginalGeometry["BornOmega"]]
    ];
    Print[
        "T-dual Born connection: ",
        CountNonZero[TDualGeometry["BornOmega"]]
    ];

    Print[
        "Original Born-selected LC connection: ",
        CountNonZero[OriginalGeometry["LCOmega"]]
    ];
    Print[
        "T-dual Born-selected LC connection: ",
        CountNonZero[TDualGeometry["LCOmega"]]
    ];

    Print[
        "Original full Ricci curvature: ",
        CountNonZero[OriginalGeometry["RicFull"]]
    ];
    Print[
        "T-dual full Ricci curvature: ",
        CountNonZero[TDualGeometry["RicFull"]]
    ];

    Print[
        "Original mixed generalised Ricci tensor: ",
        CountNonZero[OriginalGeometry["RcGen"]]
    ];
    Print[
        "T-dual mixed generalised Ricci tensor: ",
        CountNonZero[TDualGeometry["RcGen"]]
    ];

    Print[
        "Original full Riemann array: ",
        CountNonZero[OriginalGeometry["RmArray"]]
    ];
    Print[
        "T-dual full Riemann array: ",
        CountNonZero[TDualGeometry["RmArray"]]
    ];
];

PrintComponentCounts[];


(* ================================================================ *)
(* 10. Optional tensor-transport helpers                            *)
(* ================================================================ *)

(*
   TFrame maps the ORIGINAL preferred eigenframe to the T-DUAL preferred
   eigenframe, blockwise through TPlus and TMinus.

   If a covariant rank-two tensor A is transported by T-duality, its
   components in the dual preferred frame are

       AHat = TFrame^{-T} . A . TFrame^{-1}.

   This gives a useful Ricci-level check without introducing any other
   generalised connection.
*)

TFrameInverse = qMap[Inverse[TFrame], 2];

ClearAll[TransportCovariant2];
TransportCovariant2[array_] := qMap[
    Transpose[TFrameInverse] . array . TFrameInverse,
    2
];

RicOriginalTransported = TransportCovariant2[
    OriginalGeometry["RicFull"]
];

RcOriginalTransported = TransportCovariant2[
    OriginalGeometry["RcGen"]
];

RicTransportDifference = qMap[
    TDualGeometry["RicFull"] - RicOriginalTransported,
    2
];

RcTransportDifference = qMap[
    TDualGeometry["RcGen"] - RcOriginalTransported,
    2
];

Print[
    "Full Ricci curvature transport check = 0: ",
    ZeroTensorQ[RicTransportDifference]
];

Print[
    "Mixed generalised Ricci transport check = 0: ",
    ZeroTensorQ[RcTransportDifference]
];

(* On-demand transported Riemann component.  We deliberately do not build
   the full rank-four transformed array automatically because that can be
   unnecessarily expensive. *)
ClearAll[TransportedOriginalRmComponent];
TransportedOriginalRmComponent[
    p_Integer, q_Integer, s_Integer, t_Integer
] := TransportedOriginalRmComponent[p, q, s, t] = qSimp[
    Total[Flatten[Table[
        TFrameInverse[[A, p]]
        TFrameInverse[[B, q]]
        TFrameInverse[[C, s]]
        TFrameInverse[[D, t]]
        OriginalGeometry["RmArray"][[A, B, C, D]],
        {A, genDim}, {B, genDim}, {C, genDim}, {D, genDim}
    ]]]
];

ClearAll[CheckTDualRmComponent];
CheckTDualRmComponent[
    p_Integer, q_Integer, s_Integer, t_Integer
] := qSimp[
    TDualGeometry["RmArray"][[p, q, s, t]] -
    TransportedOriginalRmComponent[p, q, s, t]
];


(* ================================================================ *)
(* 11. Specialisation and regularity helpers                        *)
(* ================================================================ *)

ClearAll[SpecialiseF];
SpecialiseF[expr_,fExpr_] := FullSimplify[
    expr /. {Derivative[k_Integer][f][r] :> D[fExpr,{r,k}], f[r]->fExpr},
    Assumptions -> (r>0 && 0<th<Pi && mBH>0 && Element[qCharge,Reals])];

ClearAll[HorizonLimit];
HorizonLimit[expr_,rStar_,direction_:"FromAbove",extraAssumptions_:True] :=
    Block[{$Assumptions=True}, TimeConstrained[
        FullSimplify[Limit[expr,r->rStar,Direction->direction,
            Assumptions->extraAssumptions],Assumptions->extraAssumptions],
        30, Missing["LimitTimedOut"]]];

ClearAll[LaurentData];
LaurentData[expr_, rStar_, order_Integer : 2] :=
    Normal@Series[expr, {r, rStar, order}];

(* ================================================================ *)
(* 12. Transported horizon-regular T-dual frame                     *)
(* ================================================================ *)

(*
   WHY THIS SECTION IS NEEDED
   --------------------------
   The T-dual metric/b-field parametrisation and the preferred tangent
   eigenframe can be singular at f(r)=0.  Individual components of the
   T-dual Generalised Riemann and Ricci tensors can therefore contain 1/f
   even when the underlying transported tensor is regular.

   The natural regular frame is not guessed.  It is obtained by transporting
   the ORIGINAL generalised-metric eigenframe with the factorised T-duality
   map itself.

   Recall that the duality construction already produced

       TGen . JPlusOriginalRaw  = JPlusCanonicalTDualRaw . TPlus,
       TGen . JMinusOriginalRaw = JMinusCanonicalTDualRaw . TMinus.

   Therefore the block matrix

       Creg = diag(TPlus,TMinus)

   gives the coefficients of the TRANSPORTED original frame in the T-dual
   preferred eigenframe.

   The columns of TransportedRegularRawFrame below are the actual transported
   generators; no further rescaling or implicit change of frame is used.

   Tensor components transform covariantly with Creg.

   Connection coefficients DO NOT transform tensorially.  The code below
   includes the derivative/inhomegeneous term.  Consequently, if the Born
   connection really is transported by T-duality, its coefficients in this
   transported frame should agree directly with the ORIGINAL coefficients.
*)

ClearAll[FrameChangeCovariant2];
FrameChangeCovariant2[array_, c_] := qMap[
    Transpose[c] . array . c,
    2
];

ClearAll[FrameChangeCovariant3];
FrameChangeCovariant3[array_, c_] := Module[
    {m, t1, t2, t3},

    m = Length[c];

    t1 = Table[
        qSimp[
            Total[Flatten[Table[
                c[[A, a]] array[[A, B, Cc]],
                {A, m}
            ]]]
        ],
        {a, m}, {B, m}, {Cc, m}
    ];

    t2 = Table[
        qSimp[
            Total[Flatten[Table[
                c[[B, b]] t1[[a, B, Cc]],
                {B, m}
            ]]]
        ],
        {a, m}, {b, m}, {Cc, m}
    ];

    t3 = Table[
        qSimp[
            Total[Flatten[Table[
                c[[Cc, cc]] t2[[a, b, Cc]],
                {Cc, m}
            ]]]
        ],
        {a, m}, {b, m}, {cc, m}
    ];

    t3
];

ClearAll[FrameChangeConnection];
FrameChangeConnection[
    omega_,
    oldAnchor_,
    c_,
    x_List
] := Module[
    {m, d, cInv, newAnchor},

    m = Length[c];
    d = Length[x];

    cInv = qMap[Inverse[c], 2];

    (* If F_a = e_A C^A{}_a, then rho(F_a)=C^A{}_a rho(e_A). *)
    newAnchor = qMap[
        Transpose[c] . oldAnchor,
        2
    ];

    Table[
        qSimp[
            Total[Flatten[Table[
                cInv[[cc, Dd]] (
                    (* Inhomogeneous frame-derivative term. *)
                    Total[Flatten[Table[
                        newAnchor[[a, mu]]
                            D[c[[Dd, b]], x[[mu]]],
                        {mu, d}
                    ]]]
                    +
                    (* Tensorial-looking connection term. *)
                    Total[Flatten[Table[
                        c[[A, a]]
                        c[[B, b]]
                        omega[[A, B, Dd]],
                        {A, m}, {B, m}
                    ]]]
                ),
                {Dd, m}
            ]]]
        ],
        {a, m}, {b, m}, {cc, m}
    ]
];

ClearAll[BuildPairFrameChange];
BuildPairFrameChange[c_, pairsLocal_] := Table[
    Module[
        {
            A = pairsLocal[[i, 1]],
            B = pairsLocal[[i, 2]],
            a = pairsLocal[[j, 1]],
            b = pairsLocal[[j, 2]]
        },
        qSimp[
            c[[A, a]] c[[B, b]]
            -
            c[[A, b]] c[[B, a]]
        ]
    ],
    {i, Length[pairsLocal]},
    {j, Length[pairsLocal]}
];

ClearAll[BuildPairMetricFromInverse];
BuildPairMetricFromInverse[metricUU_, pairsLocal_] := Table[
    Module[
        {
            a = pairsLocal[[i, 1]],
            b = pairsLocal[[i, 2]],
            c = pairsLocal[[j, 1]],
            d = pairsLocal[[j, 2]]
        },
        qSimp[
            metricUU[[a, c]] metricUU[[b, d]]
            -
            metricUU[[a, d]] metricUU[[b, c]]
        ]
    ],
    {i, Length[pairsLocal]},
    {j, Length[pairsLocal]}
];

(* --------------------------------------------------------------- *)
(* 12.1 The transported frame itself                               *)
(* --------------------------------------------------------------- *)

CRegularTDual = qMap[TFrame, 2];
CRegularTDualInverse = qMap[Inverse[CRegularTDual], 2];

(* These columns are the transported original eigenframe written directly
   in the fixed RAW coordinate Courant basis

       {partial_mu, dx^mu}.

   On invariant sections of the local untwisted model, TGen intertwines the
   Courant structures. The following raw bundle frame is the cleanest
   object to inspect when asking whether the FRAME itself extends to f=0.
*)
TransportedRegularPlusRaw = qMap[
    TDualData["JPlusTDualRaw"],
    2
];

TransportedRegularMinusRaw = qMap[
    TDualData["JMinusTDualRaw"],
    2
];

TransportedRegularRawFrame = qMap[
    Join[
        TransportedRegularPlusRaw,
        TransportedRegularMinusRaw,
        2
    ],
    2
];

TransportedRegularRawFrameDet = qSimp[
    Det[TransportedRegularRawFrame]
];

regularBasisNames = Table[
    "T(" <> basisName[OriginalGeometry, A] <> ")",
    {A, genDim}
];

rawCourantBasisNames = Join[
    Table[
        "partial_" <>
        ToString[coords[[mu]], InputForm],
        {mu, n}
    ],
    Table[
        "d" <>
        ToString[coords[[mu]], InputForm],
        {mu, n}
    ]
];

ClearAll[ShowTransportedRegularRawFrame];
ShowTransportedRegularRawFrame[] := Module[
    {rows},

    rows = MapThread[
        Prepend,
        {
            TransportedRegularRawFrame,
            rawCourantBasisNames
        }
    ];

    Print[
        "Transported original eigenframe in the raw T-dual coordinate ",
        "Courant basis:"
    ];

    Print@Grid[
        Prepend[
            rows,
            Prepend[
                regularBasisNames,
                "raw coordinate basis"
            ]
        ],
        Frame -> All,
        Alignment -> Left
    ];

    Print[
        "Determinant of transported raw frame = ",
        TransportedRegularRawFrameDet
    ];
];


(* --------------------------------------------------------------- *)
(* 12.2 Generalised metric, Courant pairing and anchor              *)
(*      in the transported regular frame                            *)
(* --------------------------------------------------------------- *)

EtaRegularTDualDD = FrameChangeCovariant2[
    TDualGeometry["EtaDD"],
    CRegularTDual
];

GERegularTDualDD = FrameChangeCovariant2[
    TDualGeometry["GEDD"],
    CRegularTDual
];

EtaRegularTDualUU = qMap[
    Inverse[EtaRegularTDualDD],
    2
];

GERegularTDualUU = qMap[
    Inverse[GERegularTDualDD],
    2
];

AnchorRegularTDual = qMap[
    Transpose[CRegularTDual] .
        TDualGeometry["Anchor"],
    2
];


(* --------------------------------------------------------------- *)
(* 12.3 Born and Born-selected LC connections                      *)
(* --------------------------------------------------------------- *)

Print[""];
Print[
    "Transforming the T-dual Born connection to the transported ",
    "regular frame..."
];

BornOmegaRegularTDual = FrameChangeConnection[
    TDualGeometry["BornOmega"],
    TDualGeometry["Anchor"],
    CRegularTDual,
    coords
];

Print[
    "Transforming the T-dual Born-selected LC connection to the ",
    "transported regular frame..."
];

LCOmegaRegularTDual = FrameChangeConnection[
    TDualGeometry["LCOmega"],
    TDualGeometry["Anchor"],
    CRegularTDual,
    coords
];

BornOmegaLowRegularTDual = Table[
    qSimp[
        Total[Flatten[Table[
            BornOmegaRegularTDual[[A, B, E]]
                EtaRegularTDualDD[[E, C]],
            {E, genDim}
        ]]]
    ],
    {A, genDim}, {B, genDim}, {C, genDim}
];

LCOmegaLowRegularTDual = Table[
    qSimp[
        Total[Flatten[Table[
            LCOmegaRegularTDual[[A, B, E]]
                EtaRegularTDualDD[[E, C]],
            {E, genDim}
        ]]]
    ],
    {A, genDim}, {B, genDim}, {C, genDim}
];


(* --------------------------------------------------------------- *)
(* 12.4 Intrinsic Born torsion and Ricci tensors                   *)
(* --------------------------------------------------------------- *)

BornTorsionRegularTDual = FrameChangeCovariant3[
    TDualGeometry["BornTorsionLow"],
    CRegularTDual
];

RicRegularTDual = FrameChangeCovariant2[
    TDualGeometry["RicFull"],
    CRegularTDual
];

RcRegularTDual = FrameChangeCovariant2[
    TDualGeometry["RcGen"],
    CRegularTDual
];


(* --------------------------------------------------------------- *)
(* 12.5 Full Generalised Riemann tensor in the regular frame        *)
(* --------------------------------------------------------------- *)

(*
   The full rank-four transformation is performed efficiently on
   antisymmetric pairs.

   If P is the induced frame change on Lambda^2 E, then

       R_regular = P^T . R_preferred . P.

   This is mathematically identical to applying CRegularTDual to all four
   covariant Riemann indices, but is much faster.
*)

PairFrameRegularTDual = BuildPairFrameChange[
    CRegularTDual,
    TDualGeometry["Pairs"]
];

RmPairRegularTDual = qMap[
    Transpose[PairFrameRegularTDual] .
        TDualGeometry["RmPair"] .
        PairFrameRegularTDual,
    2
];

PairMetricRegularTDualUU = BuildPairMetricFromInverse[
    GERegularTDualUU,
    TDualGeometry["Pairs"]
];

RScalarRegularTDual = qSimp[
    1/2 Total[Flatten[Table[
        GERegularTDualUU[[A, B]]
            RicRegularTDual[[A, B]],
        {A, genDim}, {B, genDim}
    ]]]
];

KGenRegularTDual = qSimp[
    4 Tr[
        PairMetricRegularTDualUU .
        RmPairRegularTDual .
        PairMetricRegularTDualUU .
        Transpose[RmPairRegularTDual]
    ]
];


(* --------------------------------------------------------------- *)
(* 12.6 Direct transport/equivariance checks                        *)
(* --------------------------------------------------------------- *)

(*
   The transported regular frame is literally {T(e_A)}.  Therefore the
   T-dual Born connection and all transported curvature tensors should have
   exactly the SAME COMPONENT ARRAYS in this frame as the original objects
   have in {e_A}.

   This is stronger than equality of scalar contractions.
*)

BornConnectionRegularMinusOriginal = qMap[
    BornOmegaRegularTDual -
        OriginalGeometry["BornOmega"],
    3
];

LCConnectionRegularMinusOriginal = qMap[
    LCOmegaRegularTDual -
        OriginalGeometry["LCOmega"],
    3
];

BornTorsionRegularMinusOriginal = qMap[
    BornTorsionRegularTDual -
        OriginalGeometry["BornTorsionLow"],
    3
];

RicRegularMinusOriginal = qMap[
    RicRegularTDual -
        OriginalGeometry["RicFull"],
    2
];

RcRegularMinusOriginal = qMap[
    RcRegularTDual -
        OriginalGeometry["RcGen"],
    2
];

RmPairRegularMinusOriginal = qMap[
    RmPairRegularTDual -
        OriginalGeometry["RmPair"],
    2
];

RScalarRegularMinusOriginal = qSimp[
    RScalarRegularTDual -
        OriginalGeometry["RScalar"]
];

KGenRegularMinusOriginal = qSimp[
    KGenRegularTDual -
        OriginalGeometry["KGen"]
];

Print[""];
Print["============= TRANSPORTED REGULAR-FRAME CHECKS ============="];

Print[
    "Born connection in regular T-dual frame = original Born connection: ",
    ZeroTensorQ[BornConnectionRegularMinusOriginal]
];

Print[
    "Born-selected LC connection in regular T-dual frame = original LC: ",
    ZeroTensorQ[LCConnectionRegularMinusOriginal]
];

Print[
    "Intrinsic Born torsion in regular T-dual frame = original torsion: ",
    ZeroTensorQ[BornTorsionRegularMinusOriginal]
];

Print[
    "Full Ricci curvature in regular T-dual frame = original Ricci: ",
    ZeroTensorQ[RicRegularMinusOriginal]
];

Print[
    "Generalised Ricci tensor in regular T-dual frame = original Rc: ",
    ZeroTensorQ[RcRegularMinusOriginal]
];

Print[
    "Generalised Riemann pair matrix in regular T-dual frame = original Rm: ",
    ZeroTensorQ[RmPairRegularMinusOriginal]
];

Print[
    "Regular-frame Generalised Ricci scalar - original = ",
    RScalarRegularMinusOriginal
];

Print[
    "Regular-frame Generalised Kretschmann - original = ",
    KGenRegularMinusOriginal
];


(* --------------------------------------------------------------- *)
(* 12.7 Association for the existing component-display routines    *)
(* --------------------------------------------------------------- *)

Print["Recomputing curvature directly from the regular-frame connection and anchor..."];
RegularCurvatureDirect=CurvatureFromFrame["T-Dual regular (direct)",coords,
    LCOmegaRegularTDual,EtaRegularTDualDD,GERegularTDualDD,AnchorRegularTDual];
If[FailureQ[RegularCurvatureDirect],Return[RegularCurvatureDirect]];
RegularBracketLow=RawFrameBracketLow[TransportedRegularRawFrame,coords,HTwistTDualDDD];
RegularBornTorsionDirect=TorsionFromLow[BornOmegaLowRegularTDual,RegularBracketLow];
RegularLCTorsionDirect=TorsionFromLow[LCOmegaLowRegularTDual,RegularBracketLow];
TDualRegularGeometry = Join[
    KeyTake[TDualGeometry,{"Coordinates","Dimension"}],
    <|"Name"->"T-Dual transported regular frame","FrameKind"->"Transported original eigenframe",
      "BasisNames"->regularBasisNames,
      "EtaDD"->EtaRegularTDualDD,"EtaUU"->EtaRegularTDualUU,
      "GEDD"->GERegularTDualDD,"GEUU"->GERegularTDualUU,"Anchor"->AnchorRegularTDual,
      "JPlusRaw"->TransportedRegularPlusRaw,"JMinusRaw"->TransportedRegularMinusRaw,
      "RawFrame"->TransportedRegularRawFrame,"RawFrameDeterminant"->TransportedRegularRawFrameDet,
      "BornOmega"->BornOmegaRegularTDual,"BornOmegaLow"->BornOmegaLowRegularTDual,
      "BornTorsionLow"->RegularBornTorsionDirect,
      "BornTorsionVector"->Table[qSimp[Total[Flatten[Table[EtaRegularTDualUU[[C,Dd]] RegularBornTorsionDirect[[A,B,Dd]],{Dd,genDim}]]]],{A,genDim},{B,genDim},{C,genDim}],
      "BracketLow"->RegularBracketLow,
      "LCOmega"->LCOmegaRegularTDual,"LCOmegaLow"->LCOmegaLowRegularTDual,
      "LCTorsionLow"->RegularLCTorsionDirect,
      "BornDivergenceCoefficients"->Table[qSimp[Total[Flatten[Table[BornOmegaRegularTDual[[A,B,A]],{A,genDim}]]]],{B,genDim}],
      "LCDivergenceCoefficients"->Table[qSimp[Total[Flatten[Table[LCOmegaRegularTDual[[A,B,A]],{A,genDim}]]]],{B,genDim}],
      "RmPairViaTensorTransport"->RmPairRegularTDual
    |>,RegularCurvatureDirect
];

(* --------------------------------------------------------------- *)
(* 12.8 Convenience regularity diagnostics                         *)
(* --------------------------------------------------------------- *)

ClearAll[RegularFrameHorizonValue];
RegularFrameHorizonValue[
    expr_,
    fExpr_,
    rStar_,
    direction_ : "FromAbove",
    extraAssumptions_ : True
] := HorizonLimit[
    SpecialiseF[expr, fExpr],
    rStar,
    direction,
    extraAssumptions
];

ClearAll[finiteLimitQ];
finiteLimitQ[expr_] := LimitStatus[expr,True];

ClearAll[ShowRegularRicciHorizonLimits];
ShowRegularRicciHorizonLimits[
    fExpr_,
    rStar_,
    direction_ : "FromAbove",
    extraAssumptions_ : True
] := Module[
    {ric, harvested, rows},

    ric = TDualRegularGeometry["RicFull"];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[ric[[A, B]]]},
                If[
                    nonZeroQ[val],
                    With[
                        {
                            component =
                                SpecialiseF[val, fExpr]
                        },
                        With[
                            {
                                lim =
                                    RegularFrameHorizonValue[
                                        val,
                                        fExpr,
                                        rStar,
                                        direction,
                                        extraAssumptions
                                    ]
                            },
                            Sow[
                                {
                                    basisName[
                                        TDualRegularGeometry,
                                        A
                                    ],
                                    basisName[
                                        TDualRegularGeometry,
                                        B
                                    ],
                                    component,
                                    lim,
                                    finiteLimitQ[lim]
                                }
                            ]
                        ]
                    ]
                ]
            ],
            {A, genDim}, {B, genDim}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        "T-Dual regular-frame Ricci horizon limits",
        {
            "A",
            "B",
            "component",
            "horizon limit",
            "finite?"
        },
        rows
    ];
];

ClearAll[ShowRegularRiemannHorizonLimits];
ShowRegularRiemannHorizonLimits[
    fExpr_,
    rStar_,
    direction_ : "FromAbove",
    extraAssumptions_ : True
] := Module[
    {pairsLocal, rmPair, harvested, rows, p, q},

    pairsLocal = TDualRegularGeometry["Pairs"];
    rmPair = TDualRegularGeometry["RmPair"];

    harvested = Reap[
        Do[
            p = pairsLocal[[i]];
            q = pairsLocal[[j]];

            With[
                {val = qSimp[rmPair[[i, j]]]},
                If[
                    nonZeroQ[val],
                    With[
                        {
                            component =
                                SpecialiseF[val, fExpr]
                        },
                        With[
                            {
                                lim =
                                    RegularFrameHorizonValue[
                                        val,
                                        fExpr,
                                        rStar,
                                        direction,
                                        extraAssumptions
                                    ]
                            },
                            Sow[
                                {
                                    basisName[
                                        TDualRegularGeometry,
                                        p[[1]]
                                    ] <>
                                    " wedge " <>
                                    basisName[
                                        TDualRegularGeometry,
                                        p[[2]]
                                    ],

                                    basisName[
                                        TDualRegularGeometry,
                                        q[[1]]
                                    ] <>
                                    " wedge " <>
                                    basisName[
                                        TDualRegularGeometry,
                                        q[[2]]
                                    ],

                                    component,
                                    lim,
                                    finiteLimitQ[lim]
                                }
                            ]
                        ]
                    ]
                ]
            ],
            {i, Length[pairsLocal]},
            {j, i, Length[pairsLocal]}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        "T-Dual regular-frame Riemann horizon limits",
        {
            "first pair",
            "second pair",
            "component",
            "horizon limit",
            "finite?"
        },
        rows
    ];
];

ClearAll[ShowRegularLCConnectionHorizonLimits];
ShowRegularLCConnectionHorizonLimits[
    fExpr_,
    rStar_,
    direction_ : "FromAbove",
    extraAssumptions_ : True
] := Module[
    {omega, harvested, rows},

    omega = TDualRegularGeometry["LCOmega"];

    harvested = Reap[
        Do[
            With[
                {val = qSimp[omega[[A, B, C]]]},
                If[
                    nonZeroQ[val],
                    With[
                        {
                            component =
                                SpecialiseF[val, fExpr]
                        },
                        With[
                            {
                                lim =
                                    RegularFrameHorizonValue[
                                        val,
                                        fExpr,
                                        rStar,
                                        direction,
                                        extraAssumptions
                                    ]
                            },
                            Sow[
                                {
                                    basisName[
                                        TDualRegularGeometry,
                                        A
                                    ],
                                    basisName[
                                        TDualRegularGeometry,
                                        B
                                    ],
                                    basisName[
                                        TDualRegularGeometry,
                                        C
                                    ],
                                    component,
                                    lim,
                                    finiteLimitQ[lim]
                                }
                            ]
                        ]
                    ]
                ]
            ],
            {A, genDim},
            {B, genDim},
            {C, genDim}
        ]
    ][[2]];

    rows = If[
        harvested === {},
        {},
        First[harvested]
    ];

    PrintComponentTable[
        "T-Dual regular-frame LC connection horizon limits",
        {
            "first argument",
            "second argument",
            "output",
            "coefficient",
            "horizon limit",
            "finite?"
        },
        rows
    ];
];


(* Parameter conventions are copied from the supplied examples. *)
exampleF=If[choice=="Schwarzschild",1-2 mBH/r,-2 mBH/r+qCharge^2/r^2];
exampleHorizon=If[choice=="Schwarzschild",2 mBH,qCharge^2/(2 mBH)];
parameterAssumptions=mBH>0 && Element[qCharge,Reals] && 0<th<Pi &&
    If[choice=="PlanarRN",qCharge!=0,True];
checks=Join[
    <|"Finite-limit checker accepts constants"->(LimitStatus[1,mBH>0]=="Finite"),
      "Finite-limit checker rejects infinity"->(LimitStatus[Infinity,True]=="NonFinite"),
      "Finite-limit checker preserves uncertainty"->(LimitStatus[Missing["Unresolved"],True]=="Inconclusive")|>,
    Association@KeyValueMap[("Original: "<>#1)->#2&,OriginalGeometry["ConnectionChecks"]],
    Association@KeyValueMap[("Dual: "<>#1)->#2&,TDualGeometry["ConnectionChecks"]],
    Association@KeyValueMap[("Original curvature: "<>#1)->#2&,OriginalGeometry["RawCurvatureChecks"]],
    Association@KeyValueMap[("Dual curvature: "<>#1)->#2&,TDualGeometry["RawCurvatureChecks"]],
    Association@KeyValueMap[("Regular curvature: "<>#1)->#2&,TDualRegularGeometry["RawCurvatureChecks"]],
    <|"Born connection transport"->ZeroTensorQ[BornConnectionRegularMinusOriginal],
      "LC connection transport"->ZeroTensorQ[LCConnectionRegularMinusOriginal],
      "Born torsion transport"->ZeroTensorQ[BornTorsionRegularMinusOriginal],
      "Regular torsion recomputation"->ZeroTensorQ[RegularBornTorsionDirect-BornTorsionRegularTDual],
      "Regular LC torsion vanishes"->ZeroTensorQ[RegularLCTorsionDirect],
      "Direct regular Riemann = tensor transform"->ZeroTensorQ[TDualRegularGeometry["RmPair"]-RmPairRegularTDual],
      "Full regular Riemann = original"->ZeroTensorQ[TDualRegularGeometry["RmArray"]-OriginalGeometry["RmArray"]],
      "Direct regular Ricci = tensor transform"->ZeroTensorQ[TDualRegularGeometry["RicFull"]-RicRegularTDual],
      "Direct regular Rc = tensor transform"->ZeroTensorQ[TDualRegularGeometry["RcGen"]-RcRegularTDual],
      "Regular Ricci = original"->ZeroTensorQ[TDualRegularGeometry["RicFull"]-OriginalGeometry["RicFull"]],
      "Regular Rc = original"->ZeroTensorQ[TDualRegularGeometry["RcGen"]-OriginalGeometry["RcGen"]],
      "Scalar transport"->ZeroTensorQ[{RScalarDifference}],
      "Kretschmann transport"->ZeroTensorQ[{KGenDifference}],
      "Direct regular scalar = original"->ZeroTensorQ[{TDualRegularGeometry["RScalar"]-OriginalGeometry["RScalar"]}],
      "Direct regular Kretschmann = original"->ZeroTensorQ[{TDualRegularGeometry["KGen"]-OriginalGeometry["KGen"]}],
      "Regular generalised metric = original"->ZeroTensorQ[GERegularTDualDD-OriginalGeometry["GEDD"]],
      "Regular Courant pairing = original"->ZeroTensorQ[EtaRegularTDualDD-OriginalGeometry["EtaDD"]],
      "Regular raw anchor matches frame"->ZeroTensorQ[AnchorRegularTDual-Transpose[TransportedRegularRawFrame[[1;;n,All]]]],
      "Original torsion removal preserves divergence"->ZeroTensorQ[OriginalGeometry["BornDivergenceCoefficients"]-OriginalGeometry["LCDivergenceCoefficients"]],
      "Dual torsion removal preserves divergence"->ZeroTensorQ[TDualGeometry["BornDivergenceCoefficients"]-TDualGeometry["LCDivergenceCoefficients"]],
      "Divergence transport in invariant frame"->ZeroTensorQ[TDualRegularGeometry["LCDivergenceCoefficients"]-OriginalGeometry["LCDivergenceCoefficients"]]
    |>];
(* Original anchor-induced zero-flux connection has metric divergence.
   The dual trace is also compared directly with the Buscher dilaton shift. *)
metricDivOriginal=Table[qSimp[Total[Flatten[Table[OriginalGeometry["GammaLC"][[mu,mu,nu]],{mu,n}]]]],{nu,n}];
metricDivDual=Table[qSimp[Total[Flatten[Table[TDualGeometry["GammaLC"][[mu,mu,nu]],{mu,n}]]]],{nu,n}];
dPhiDual=Table[qSimp[-D[gOriginalDD[[dualIndex,dualIndex]],coords[[nu]]]/(2 gOriginalDD[[dualIndex,dualIndex]])],{nu,n}];
AssociateTo[checks,
 "Original metric divergence"->ZeroTensorQ[OriginalGeometry["LCDivergenceCoefficients"]-Join[metricDivOriginal,metricDivOriginal]],
 "Dual Buscher divergence"->ZeroTensorQ[TDualGeometry["LCDivergenceCoefficients"]-Join[metricDivDual-2 dPhiDual,metricDivDual-2 dPhiDual]]];
Print["Computing the independent ordinary Riemann contraction..."];
ordinaryK=OrdinaryKretschmann[OriginalGeometry];
AssociateTo[checks,"Generalised Kretschmann = 4 ordinary Kretschmann"->ZeroTensorQ[{OriginalGeometry["KGen"]-4 ordinaryK}]];
result=<|"Example"->choice,"TransverseModel"->transverseModel,
  "fExpression"->exampleF,"Horizon"->exampleHorizon,"ParameterAssumptions"->parameterAssumptions,
  "Original"->OriginalGeometry,"DualPreferred"->TDualGeometry,"DualRegular"->TDualRegularGeometry,
  "SpecialisedOriginal"->Map[SpecialiseF[#,exampleF]&,OriginalGeometry],
  "SpecialisedDualPreferred"->Map[SpecialiseF[#,exampleF]&,TDualGeometry],
  "SpecialisedDualRegular"->Map[SpecialiseF[#,exampleF]&,TDualRegularGeometry],
  "OrdinaryKretschmann"->SpecialiseF[ordinaryK,exampleF],
  "Checks"->checks,"AllChecksProved"->(And@@Values[checks]),
  "Scope"->"Local invariant untwisted EF examples; transported Born structure. No global horizon extension is asserted."|>;
Print[ShowValidation[result]];
result

];
End[];
EndPackage[];
