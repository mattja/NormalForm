(* ::Package:: *)

(* :Title: NormalForm *)

(* :Copyright:
    Copyright 2015-2016 Matthew J. Aburn
    This program is free software: you can redistribute it and/or modify 
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or 
    (at your option) any later version. See <http://www.gnu.org/licenses/>.
*)

(* :Summary:
    Find the normal form transformation of an n-dimensional dynamical system 
    at a local bifurcation point.
    Analyze the effect of noise on oscillations at a Hopf bifurcation.
*)

BeginPackage["NormalForm`"]

NormalFormTransformation::usage = 
"NormalFormTransformation[rhs, {x1,...,xn}, {u1,...,un}, m] transforms the\
 dynamical system with right hand side rhs (expressed in original variables\
 {xi}) to a simpler system (normal form to order m) in the new variables\
 {ui}. Returns a pair {newrhs, trans} where newrhs is the transformed\
 system and trans is a smooth invertible coordinate transformation that\
 maps rhs to newrhs.\
 Options: \
 Verbose: whether to print working at each step (default False).\
 BifurcationParameters: which symbols in `rhs` should be interpreted as\
   small bifurcation parameters (default Global`\[Epsilon])\
 AsymptoticScaling: relative scaling of variables and bifurcation\
   parameters in the asymptotic limit, used when truncating power series.\
   e.g. {u1, u2, u3, Sqrt[Global`\[Epsilon]]} means to assume\
   O(\[Epsilon]) == O(u_i^2) which suits a Hopf or Pitchfork bifurcation.\
 Extended: whether to compute the normal form of the extended system,\
   that is with phase space extended with dimensions for the (rescaled)\
   bifurcation parameters and their equations \\dot{\\alpha} == 0,\
   thus finding a transformation dependent on the bifurcation parameters.\
   If False (the default) then the normal form will be found with\
   respect to the dynamical variables only.";

MultiSeries::usage = 
"MultiSeries[v, {x1,...,xn}, m] generates a multivariate power series\
 expansion for the vector field v (in variables {xi}) about the origin,\
 to order x^m";

MultiSeriesData::usage = 
"MultiSeriesData[{x1,...,xn}, {e1,...,en}, \[Delta], seriesdata] represents a\
 multivariate power series in the variables {xi^ei} about the origin.";

MultiSeriesQ::usage = "MultiSeriesQ[expr] tests whether expr is a MultiSeries";

MultiSeriesFieldQ::usage = 
"MultiSeriesFieldQ[expr] tests whether expr is a vector field with each\
 entry a MultiSeries.";

TransformContravariant::usage = 
"TransformContravariant[U, R] applies the near-identity coordinate\
 transformation U to transform the contravariant vector field R(u).\
 Both U and R should be given in the form of MultiSeries";

TransformNoisyHopf::usage = 
"TransformNoisyHopf[rhs, {x1,...,xn}, {\[Sigma]1,...,\[Sigma]n},\
 {\[Xi]1,...\[Xi]n}, {r,\[Theta]}, {new\[Xi]1,new\[Xi]2}] takes the stochastic\
 dynamical system with right hand side rhs (expressed in variables {xi},\
 small noise parameters {\[Sigma]i} and Langevin noise symbols {\[Xi]i}\
 with Stratonovich interpretation of any multiplicative noise) and\
 transforms it to a simple circular 2 dimensional Hopf normal form system\
 (expressed in new polar variables {r, \[Theta]} and new Langevin noise\
 symbols {new\[Xi]1, new\[Xi]2}).\
 Currently it is assumed that the linear part of the system has already\
 been transformed to Jordan real form, with Hopf in first two variables.\
 Returns:\
   2D polar vector field representing transformed system in a standard form\
 Options:\
 Verbose: whether to print working at each step (default False)\
 BifurcationParameters: which symbols in `rhs` should be interpreted as\
   small bifurcation parameters (default Global`\[Epsilon])\
 AsymptoticScaling: relative scaling of variables and bifurcation\
   parameters in the asymptotic limit, used when truncating power series.\
   e.g. {u1, u2, u3, Sqrt[Global`\[Epsilon]]} means to assume\
   O(\[Epsilon]) == O(u_i^2) which suits a Hopf or Pitchfork bifurcation.\
 Extended: whether to compute the normal form of the extended system,\
   that is with phase space extended with dimensions for the (rescaled)\
   bifurcation parameters and their equations \\dot{\\alpha} == 0,\
   thus finding a transformation dependent on the bifurcation parameters.\
   If False (the default) then the normal form will be found with\
   respect to the dynamical variables only.\
 Average: whether to average around the cycle. (default True)\
 Rescale: whether to linearly rescale the radial variable to make the\
   coefficient of the R^3 term -1. (default True)";

ToPolar::usage = 
"ToPolar[v, {x1,...,xn}, {r, \[Theta]}] transforms a n-dimensional flow,\
 represented as a contravariant vector field v from Cartesian to cylindrical \
 coordinates. (The first two variables {x1, x2} are changed to polar variables \
 {r, \[Theta]}). Each component of v should be an expression in the {xi}";

ToCartesian::usage = 
"ToCartesian[v, {r, \[Theta]}, {x1, x2,...}] transforms a n-dimensional flow, \
 represented as a contravariant vector field v from cylindrical to Cartesian \
 coordinates. (The first two components of v are assumed to be in polar \
 coordinates {r, \[Theta]}).";

BalanceMatrix::usage =
"BalanceMatrix[A] returns the pair {T, B} where T is a similarity\
 transformation and B is the transformed matrix, B = T^-1.A.T, such that\
 B is as close to symmetric as possible. This is used to improve an\
 ill-conditioned matrix A, allowing eigenvalues and eigenvectors to be\
 computed more precisely from matrix B. Ref: Parlett and Reinsch (1969)";

Arrange::usage =
"Arrange[expr, subexpressions] partly factorizes the expression by collecting\
 any subexpressions that are listed. The purpose of this function is to allow\
 printing the expression in a more readable form.";

ArrangePolar::usage =
"ArrangePolar[expr, subexpressions, angleVars] partly factorizes the expression\
 by collecting subexpressions, as well as collecting Sin(n x), Cos(n x) for any\
 symbols x in the list angleVars, for integers n = 1..6. The purpose of this\
 function is to allow printing the expression in a more readable form.";


Begin["`Private`"]

(* Define some utility functions *)

(* Function for verbose print output during computation. This can be turned on
   by giving the option Verbose->True in NormalFormTransformation[] *)
verbose = False;
dPrint[x___] := Block[{$Context="NormalForm`Private`"}, Print @@ {x}] /; verbose
dPrint[x___] := Null /; !verbose

(* Filter for output of numbers.
   If "exactOutput" is set, then the verbose output will print numbers without 
   numerical approximation. Otherwise will output numbers in decimal form with 
   4 significant figures, zeroing values smaller than 10^-10 *)
exactOutput = True;
NN[x_] := Identity[x] /; exactOutput
NN[x_] := N[x, 4] /; !exactOutput

(* Utility functions to collect expressions for easy to read output. *)
Arrange[expr_, syms_List] := Collect[Normal[expr], Flatten[syms, 1]];

(* Also collect trig expressions such as Sin[n\[Theta]], for all n=1..m. *)
ArrangePolar[expr_, syms_List, angleVars_List] :=
    Module[{trigExprs, m=6},
        trigExprs = Flatten@Outer[Map[#1, {#2 #3}]&, 
                                  {Sin, Cos}, Range[m], angleVars];
        Collect[TrigReduce[Normal[expr]],
                Flatten[syms, 1] ~Join~ trigExprs]
    ];

(* ShowIt function to aid debugging, suggested by Leonid Shiffrin *) 
SetAttributes[ShowIt, HoldAll];
ShowIt[expr_] :=
    Module[{y},
        Print[ToString[Unevaluated[expr]], " = ", y = expr]; 
        y
    ]

(* Project vector into a subspace, given an orthogonal basis for the subspace *)
Project[u_, orthogbasis_] := Plus@@Map[x\[Function]Projection[u,x],orthogbasis]

(* Balancing of matrices prior to eigenvalue computation: *)

(* Returns a pair {m,P} where P is a similarity transformation permuting rows 
   and columns of A so that any rows/columns that directly give eigenvalues 
   of A are moved to the top left of the transformed matrix P.A.Transpose[P].
   m is the number of rows/columns thus moved. *)
prPermute[A_?MatrixQ] :=
    Module[{n, m, A0, i, P, B, subP},
        n = Length[A];
        A0 = A - DiagonalMatrix[Diagonal[A]];
        i = FirstPosition[A0, Table[0, {n}], {0}][[1]];
        If[i == 0,
            i = FirstPosition[Transpose[A0], Table[0, {n}], {0}][[1]];
            If[i == 0, Return[{0, IdentityMatrix[n]}]] (* nothing to be done *)
        ] ;
        P = Permute[IdentityMatrix[n], Cycles[{Range[i]}]];
        B = P.A.Transpose[P];
        (* having removed row and column i, recursively treat the 
           remaining submatrix: *)
        {m, subP} = prPermute[B[[2 ;;, 2 ;;]]];
        subP = PadLeft[subP, {n, n}];
        subP[[1, 1]] = 1;
        Return[{m + 1, subP.P}]
    ]

(* BalanceMatrix[A] returns the pair {T, B} where T is a similarity 
   transformation and B is the transformed matrix, B = T^-1.A.T, such that 
   B is as close to symmetric as possible. This is used to improve an 
   ill-conditioned matrix A, allowing eigenvalues and eigenvectors to be 
   computed more precisely from matrix B. Ref: Parlett and Reinsch (1969) *)
BalanceMatrix[A_?MatrixQ, p_:2] :=
    Module[{m, P, Bwhole, B, N, n, D, Dinc, Dcum, Dwhole, R, C, f, fcond, 
            \[Gamma] = 1},
        {m, P} = prPermute[A];
        Bwhole = P.A.Transpose[P];
        B = Bwhole[[m + 1 ;;, m + 1 ;;]];
        N = Length[Bwhole];
        n = Length[B];
        Dcum = Table[1, {n}];
        Dinc = Table[0, {n}];
        While[Dinc != Table[1, {n}],
            R = Map[Norm[#, p] &, B];
            C = Map[Norm[#, p] &, Transpose[B]];
            f = 2^(Ceiling[(Log2[R/C] + 1)/2] - 1);
            fcond[R_, C_, f_] :=
                If[(C f)^p + (R/f)^p < \[Gamma] (C^p + R^p), f, 1];
            Dinc = MapThread[fcond, {R, C, f}];
            Dcum = Dinc Dcum;
            B = DiagonalMatrix[1/Dinc].B.DiagonalMatrix[Dinc];
        ];
        Dwhole = Table[1, {m}]~Join~Dcum;
        Bwhole = DiagonalMatrix[1/Dwhole].Bwhole.DiagonalMatrix[Dwhole];
        {Transpose[P].DiagonalMatrix[Dwhole], Bwhole}
    ] /; NumberQ[p] && NonNegative[p] || p === \[Infinity]

(* Define some useful tools for differential operators and polynomial spaces *)

SymbolListQ[expr_] := VectorQ[expr, Head[#]===Symbol&]

(* check whether x is a polynomial vector field in variables u *)
PolyFieldQ[x_, u_] := VectorQ[x, PolynomialQ[#, u]&]

(* turn a polynomial expression into a function so we 
   can apply the "Derivative" operator to it *)
polyToFunc[poly_, u_] := Function@@{u, poly} /; PolynomialQ[poly, u]

(* turn a polynomial vector field expression into a function *)
polyNToFunc[polyfield_, u_] :=
    Function[{x}, Table[polyToFunc[polyfield[[i]],u]@@x, {i,Length[u]}]]

(* define an algebra of differential operators: *)
ApplyD[A_, expr_, u_] := A expr /; FreeQ[A,Derivative]
ApplyD[L1_+L2_, expr_, u_] := ApplyD[L1,expr,u] + ApplyD[L2,expr,u]
ApplyD[A_ L_, expr_, u_] := A ApplyD[L,expr,u] /; FreeQ[A,Derivative]
ApplyD[A:HoldPattern[D[__]&], expr_, u_] := A[expr]
ApplyD[L1_**L2_, expr_, u_] := Expand[ApplyD[L1, ApplyD[L2, expr, u], u]]
ApplyD[L1_^n_Integer, expr_, u_] := Nest[Expand[ApplyD[L1,#1,u]]&,expr,n] /; n>1
(* tell it how to apply differential operators to polynomial expressions *)
ApplyD[L_, poly_, u_] := L[polyToFunc[poly,u]]@@u /; PolynomialQ[poly,u]
(* and how to apply vector operators to vector field expressions *)
ApplyND[L_, vec_List, u_] := Plus@@Thread[ApplyD[L,vec,u], List, 2]

(* basis table of first order differential operators for n variables *)
diff1[n_Integer?Positive] :=
    Array[i\[Function](Derivative@@ReplacePart[Table[0,{n}], i->1]), n]

(* construct basis for vector space of mth order homogeneous 
   polynomials in n variables *)

products[m_, u_] := Flatten[Outer@@({Times}~Join~Table[u, {i, m}])]

basis1D[m_, u_] := Union/@Gather[products[m, u]] // Flatten

(* construct dual basis (using differential operator tools we defined above) *)
dpoly1D[m_, u_] := Plus@@@Gather[products[m, u]]/Factorial[m]

dualbasis1D[m_, u_] :=
    Module[{n=Length[u]},
        dpoly1D[m, u] /. Thread[u -> diff1[n]] /. Times->NonCommutativeMultiply
    ]

(* Finally, from these ingredients we make the basis and dual basis for the
   space of homogeneous polynomial vector fields of mth order in n variables: *)

rnbasis[n_Integer?Positive] := IdentityMatrix[n]

PolyFieldBasis[m_, u_] := 
    Module[{n=Length[u]},
        Flatten[Transpose[Outer[Times, basis1D[m,u], rnbasis[n]]], 1]
    ]

PolyFieldDualBasis[m_,u_] :=
    Module[{n=Length[u]},
        Flatten[Transpose[Outer[Times, dualbasis1D[m,u], rnbasis[n]]], 1]
    ]

(* define homological operator (Lie bracket of vector fields A.u and q(u) *)
L[A_?MatrixQ, q_, u_?SymbolListQ] := D[q,{u,1}].A.u - A.q /; PolyFieldQ[q, u]


(* Convert n-dimensional flow (which is represented as a vector field)
   from Cartesian to cylindrical coordinates.
   The first two dimensions are put in polar form. *)
ToPolar[field_?VectorQ, u_?SymbolListQ, {r_Symbol, \[Theta]_Symbol}] :=
    Module[{udotpolar, thetadot, rdot, therest},
        udotpolar = field /. {u[[1]]->r Cos[\[Theta]],u[[2]]->r Sin[\[Theta]]};
        thetadot = (Cos[\[Theta]] udotpolar[[2]] - 
                    Sin[\[Theta]]udotpolar[[1]]) / r;
        rdot = (Cos[\[Theta]]udotpolar[[1]] + 
                Sin[\[Theta]]udotpolar[[2]]);
        therest = udotpolar[[3;;]];
        {rdot,thetadot}~Join~therest
    ]

(* Convert n-dimensional flow (which is represented as a vector field)
   from cylindrical to Cartesian coordinates.
   The first two dimensions of the input are assumed to be in polar form. *)
ToCartesian[field_?VectorQ, {r_Symbol, \[Theta]_Symbol}, u_?SymbolListQ] :=
    Module[{u1dot, u2dot, therest},
        u1dot = (Cos[\[Theta]]field[[1]] - r Sin[\[Theta]]field[[2]] /.
            {r->Sqrt[u[[1]]^2+u[[2]]^2], \[Theta]->ArcTan[u[[1]],u[[2]]]}) // 
            TrigExpand;
        u2dot = (Sin[\[Theta]]field[[1]]+r Cos[\[Theta]]field[[2]] /.
            {r->Sqrt[u[[1]]^2+u[[2]]^2], \[Theta]->ArcTan[u[[1]],u[[2]]]}) // 
            TrigExpand;
        therest = (field[[3;;]] /.
            {r->Sqrt[u[[1]]^2+u[[2]]^2], \[Theta]->ArcTan[u[[1]],u[[2]]]}) //
            TrigExpand;
        {u1dot,u2dot}~Join~therest
    ]


(* Mathematica only has proper built-in support for univariate power series. 
  We define some useful tools for multivariate power series about 0.
  This will be much faster than using polynomial expressions, as we will
  automatically truncate terms beyond maxOrder during operations.
  TODO: allow series expansion about points other than the origin *)

(* declare these symbols in NormalForm`Private` module scope for internal use *)
\[Delta]; \[Alpha]; r; \[Theta]; \[Phi];
$Assumptions = $Assumptions && r>=0 && \[Theta]\[Element]Reals &&
               \[Phi]\[Element]Reals && \[Alpha]\[Element]Reals;

(* normalize the list of symbol powers given in a MultiSeries expression *)
processVars[vars_] :=
    Module[{newvars, symbols, exponents},
        If[Head[vars]==List, 
            newvars = Flatten[vars]
        , 
            newvars = {vars}
        ];
        symbols = Replace[newvars, (xxx_Symbol)^(y_) \[Rule] xxx, {1}];
        exponents = Replace[newvars, {Except[(x_)^(y_)] \[Rule] 1,
                                      (x_Symbol)^(y_) \[Rule] y}, {1}];
        {symbols, exponents}
    ];

unionVars[syms1_, exps1_, syms2_, exps2_] :=
    Module[{pairs, newsyms, newexps},
        pairs = Union[Thread[syms1->exps1],Thread[syms2->exps2]];
        newsyms = Keys[pairs];
        newexps = Values[pairs];
        If[Length[newsyms] > CountDistinct[newsyms],
            Print["Don't yet support combining different limits ",
                   syms1^exps1, " and ", syms2^exps2];
            Abort[];
        ];
        {newsyms, newexps}
    ];

unionScaling[U_?MultiSeriesFieldQ, R_?MultiSeriesFieldQ] :=
    Module[{usyms, uexps, rsyms, rexps, newsyms, newexps},
        usyms = U[[1]][[1]];
        uexps = U[[1]][[2]];
        rsyms = R[[1]][[1]];
        rexps = R[[1]][[2]];
        {newsyms, newexps} = unionVars[usyms, uexps, rsyms, rexps];
        newsyms^newexps
    ];

(* Approximate a scalar expression locally to origin with multivariate series *)
MultiSeries[expr_, vars_, maxOrder_Integer?NonNegative] :=
    Module[{symbols, exponents},
        {symbols, exponents} = processVars[vars];
        MultiSeriesData[
            symbols, exponents, \[Delta],
            Series[expr /. Thread[symbols -> symbols*\[Delta]^(1/exponents)],
                   {\[Delta], 0, maxOrder}] + O[\[Delta]]^(maxOrder+1)
        ]
    ];

(* Approximate a field or matrix locally to origin with multivariate series *)
MultiSeries[array_?ArrayQ, vars_, maxOrder_Integer?NonNegative] :=
    Thread[Unevaluated[MultiSeries[array, vars, maxOrder]], List, 1] /;
        ArrayQ[array]

(* test whether an expression is a MultiSeries or MultiSeriesField *)
MultiSeriesQ[expr_] := Head[expr]===MultiSeriesData;
MultiSeriesFieldQ[expr_] := VectorQ[expr, MultiSeriesQ];

(* allow using 'Normal[]' to truncate a series (or series field or matrix) 
   to a polynomial (or polynomial field or matrix): *)
Unprotect[Normal];

Normal[MultiSeriesData[symbols_, exponents_, \[Delta]_, seriesdata_]] := 
    Normal[seriesdata] /. \[Delta]->1 // Expand

Normal[seriesArray_?ArrayQ] :=
    Normal /@ seriesArray /; ArrayQ[seriesArray, _, Head[#]===MultiSeriesData&]

Protect[Normal];

(* define some operators on multivariate series *)
Unprotect[Plus, Times, Power];
(* currently symbols and exponents of the two series must be the same.
   TODO implement operations where these are not the same *)

MultiSeriesData[vars1_, exp1_, \[Delta]1_, seriesdata1_] + 
        MultiSeriesData[vars2_, exp2_, \[Delta]2_, seriesdata2_] :=
    Module[{newvars, newexps},
        {newvars, newexps} = unionVars[vars1, exp1, vars2, exp2];
        MultiSeriesData[newvars, newexps, \[Delta]1, seriesdata1+seriesdata2]
    ] /; \[Delta]1===\[Delta]2

MultiSeriesData[vars1_, exp1_, \[Delta]1_, seriesdata1_] MultiSeriesData[
        vars2_, exp2_, \[Delta]2_, seriesdata2_] :=
    Module[{newvars, newexps},
        {newvars, newexps} = unionVars[vars1, exp1, vars2, exp2];
        MultiSeriesData[newvars, newexps, \[Delta]1, seriesdata1 seriesdata2]
    ] /; \[Delta]1===\[Delta]2

x_?NumberQ MultiSeriesData[vars_, exponents_, \[Delta]_, seriesdata_] :=
    MultiSeriesData[vars, exponents, \[Delta], x seriesdata]

MultiSeriesData[vars_, exponents_, \[Delta]_, seriesdata_]^0 := 1

MultiSeriesData[vars_, exponents_, \[Delta]_, seriesdata_]^n_Integer := 
    MultiSeriesData[vars, exponents, \[Delta], seriesdata^n]

Protect[Plus, Times, Power];

(* calculate matrix power while not keeping terms beyond O[u]^maxOrder *)
TruncatingMatrixPower[x_, n_Integer, u_?SymbolListQ, maxOrder_Integer] := 
    Nest[Dot[x,#]+O[u]^(maxOrder+1)&, IdentityMatrix[Length[x]], n] /; 
        MatrixQ[x] && Length[x]==Length[x[[1]]] && n>=0

(* TODO define arithmetic operations between series and ordinary polynomials *)

(* Define an asymptotic notation for input of multivariate series *)
Unprotect[O];
O[vars_List] :=
    Module[{symbols, exponents},
        {symbols, exponents} = processVars[vars];
        MultiSeriesData[symbols, exponents, \[Delta], O[\[Delta]]]
    ]
Protect[O];

(* tell mathematica how to output series *)
Format[MultiSeriesData[symbols_, exponents_, \[Delta]_, seriesdata_]] :=
    Expand[Normal[seriesdata] /. \[Delta]->1] + 
    Superscript["O[" <> ToString[symbols^exponents, StandardForm] <> "]",
                seriesdata[[5]]/seriesdata[[6]]];

(* Now define how to substitute one series into another series *)

(* First redefining the "/." operator *)
(* TODO extend to correctly cover more cases *)
Unprotect[ReplaceAll];

ReplaceAll[MultiSeriesData[syms1_, exps1_, \[Delta]1_, seriesdata1_],
           x1_->MultiSeriesData[syms2_, exps2_, \[Delta]2_, seriesdata2_]] := 
    Module[{ex},
        ex = exps1[[FirstPosition[syms1, x1, {1}, 1][[1]]]];
        MultiSeriesData[
            Union[Complement[syms1, {x1}], syms2],
            exps1,
            \[Delta]1,
            seriesdata1 /. x1->(Normal[seriesdata2]/(\[Delta]1^(1/ex)))
        ]
    ] /; MemberQ[syms1, x1] && exps1===exps2 && \[Delta]1===\[Delta]2

ReplaceAll[MultiSeriesData[syms1_, exp1_, \[Delta]_, seriesdata1_],
           ruleList:List[(_Symbol->_MultiSeriesData)..]] :=
    Module[{symbolsToReplace, exs, ruleRHSsymbols, newSymbolList, newSubRules},
        symbolsToReplace = ruleList[[All,1]];
        exs = exp1[[Map[Position[syms1, #, {1}, 1][[1]][[1]]&, 
                        symbolsToReplace]]];
        ruleRHSsymbols = Union@@ruleList[[All,2]][[All,1]];
        newSymbolList =
            Union[Complement[syms1,symbolsToReplace],ruleRHSsymbols];
        newSubRules =
            Table[symbolsToReplace[[i]] \[Rule]
                    (Normal[ruleList[[i]][[2]][[4]]]/(\[Delta]^(1/exs[[i]]))),
                  {i, Length[ruleList]}];
        MultiSeriesData[newSymbolList, exp1, \[Delta], 
                        seriesdata1 /. newSubRules]
    ]

ReplaceAll[MultiSeriesData[syms_, exps_, \[Delta]_, seriesdata_], 
           x_->y_?NumberQ] :=
    Module[{pos},
        If[Length[syms]==1,
            Expand[Normal[seriesdata] /. \[Delta]->1] /. x->y
        ,
            pos = FirstPosition[syms, x, {1}, 1][[1]];
            MultiSeriesData[Drop[syms, {pos}], Drop[exps, {pos}], \[Delta],
                            seriesdata /. x->y]
        ]
    ] /; MemberQ[syms, x]

ReplaceAll[ms_MultiSeriesData, ruleList:List[(_Symbol->y_?NumberQ)..]] :=
    Fold[ReplaceAll, ms, ruleList]

ReplaceAll[msf_?MultiSeriesFieldQ, x_Symbol->y_] :=
    Map[ReplaceAll[#, x->y]&, msf]

ReplaceAll[msf_?MultiSeriesFieldQ, ruleList:List[(_Symbol->y_)..]] :=
    Map[ReplaceAll[#, ruleList]&, msf]

Protect[ReplaceAll];

(* TODO validate that all the series in ruleList above use 
   the same \[Delta] as the original series.
   TODO enable substituting an ordinary polynomial expression into a series *)

(* Define transforming a series vector field "field1" by substituting 
   another series vector field "field2" into it *)
SubstituteField[field1_?MultiSeriesFieldQ, field2_?MultiSeriesFieldQ] :=
    Module[{syms=field1[[1]][[1]], exps=field1[[1]][[2]], 
            \[Delta]1=field1[[1]][[3]], m, newfield2},
        (* TODO validate that all series on LHS use the same list of symbols.
           In particular having a common ordering of symbols is crucial. *)
        (* max order of field1: *)
        m = (#[[4]][[5]]/#[[4]][[6]]& /@ field1 // Max) - 1;
        If[Length[syms] < Length[field2],
            dPrint["dimension mismatch in substitution"];
            Abort[]];
        newfield2 = field2;
        (* use identity transformation for any extra variables in field1 *)
        If[Length[syms] > Length[field2],
           newfield2 = field2~Join~Table[MultiSeries[syms[[i]], syms^exps, m],
                                         {i, Length[field2]+1, Length[syms]}]
        ];
        (* Next line uses our definition of /. (ReplaceAll) given above *)
        Table[field1[[i]] /. Thread[syms->newfield2], {i,Length[field1]}]
    ]

(* Compose a list of transformations, each expressed as a series field *)
(* TODO we currently assume here all transformations do not change the list 
   of symbols, i.e. they map u_i to a expression of u_j's, not v_j's *)
ComposeTransformations[transformationsList:List[List[_MultiSeriesData..]..]] :=
    Fold[SubstituteField, transformationsList]

(* Extract terms of homogenous order m from a series or a series field *)

OrderTerms[ms_MultiSeriesData, m_Integer?NonNegative] :=
    Module[{nmin=ms[[4]][[4]]},
        If[m < nmin || m > (Length[ms[[4]][[3]]] + nmin - 1),
            0
        ,
            ms[[4]][[3]][[m-nmin+1]]
        ]
    ]

OrderTerms[msf_?MultiSeriesFieldQ, m_Integer?NonNegative] :=
    Map[OrderTerms[#, m]&, msf]

(* Find sum of terms of homogenous order from min to max inclusive*)

OrderTerms[ms_MultiSeriesData, 
           {min_Integer?NonNegative, max_Integer?NonNegative}] :=
    Sum[OrderTerms[ms, i], {i, min, max}]

OrderTerms[msf_?MultiSeriesFieldQ,
           {min_Integer?NonNegative, max_Integer?NonNegative}] :=
    Sum[OrderTerms[msf, i], {i, min, max}]


(* Explicit direct algorithm (slow and memory hungry) *)
semisimpleAlgorithm1z[rm_?VectorQ, LL_?ArrayQ, m_Integer, u_?VectorQ] :=
    Module[{basis, dualbasis, basissize, image, kernel, orthogToIm, K, Q, sm},
        basis = PolyFieldBasis[m, u];
        dualbasis = PolyFieldDualBasis[m, u];
        basissize = Length[basis];
        (* orthonormal basis for column space of LL: *)
        image = Select[Orthogonalize[Transpose[LL]],!AllTrue[#,PossibleZeroQ]&];
        (* Murdock, Theorem 2.1.3: as LL is semisimple, its kernel and 
           image are linearly independent, so span the whole space *)
        (* As our chosen complement of Image(LL), use Kernel(LL).
           But note in general these two subspaces are not orthogonal 
           (with respect to the inner product induced by our basis) *)
        dPrint["Finding basis for kernel of L..."];
        kernel = NullSpace[LL];
        dPrint["Finding basis for subspace orthogonal to Image(L)..."];
        orthogToIm =
            Orthogonalize[kernel - Map[Project[#,image]&,kernel]] // Simplify;
        (* Want to find that vector that lies in ker LL but has the same
           projection on orthogToIm as our existing terms. That vector 
           will give the terms we retain in the transformed version *)
        dPrint["Solving linear system to find transformed terms."];
        If[Length[kernel]!=0,
            K = Transpose[kernel];
            Q = orthogToIm;
            sm = K.LinearSolve[Q.K,Q.rm];
        , 
            sm = Table[0, {basissize}];
        ];
        sm
    ]


(* Faster direct algorithm from Murdock chapter 2.1 *)
semisimpleAlgorithm1a[rm_?VectorQ, LL_?ArrayQ, m_Integer, u_?VectorQ] :=
    Module[{TTbal, LLbal, uniqueNonzeroEigenvals, ident, n, sI, sLL, tab,
            projector, matrixNorm, sm},
        {TTbal, LLbal} = BalanceMatrix[LL];
        uniqueNonzeroEigenvals = 
            Select[DeleteDuplicates[Chop@Simplify[Eigenvalues[LLbal]], 
                                    Abs[#1-#2]<10^-10&], 
                   !#===0&];
        dPrint["uniqueNonzeroEigenvals == ", uniqueNonzeroEigenvals];
        dPrint["number of uniqueNonzeroEigenvals == ",
               Length[uniqueNonzeroEigenvals]];
        sLL = SparseArray[LL];
        (* Find spectral projection for eigenvalue zero. will map vectors 
           to the kernel of LL: *) 
        (**
        poly[x_] := Product[1-x/uniqueNonzeroEigenvals[[i]],
                            {i,Length[uniqueNonzeroEigenvals]}];
        projector = MatrixFunction[poly, sLL];
        **)
        n = Length[LL];
        sI = SparseArray[Band[{1,1}]->1, {n,n}];
        sS[\[Lambda]] := SparseArray[Band[{1,1}]->1/\[Lambda], {n,n}];
        mnext[M_, \[Lambda]_] := Dot[M, sI - sS[\[Lambda]].sLL];
        projector = Fold[mnext, sI, uniqueNonzeroEigenvals];
        (**
        (* Optional: check that it does... *)
        matrixNorm = Norm[LL.projector];
        dPrint["matrixNorm == ", matrixNorm];
        If[matrixNorm < 10^-6,
            dPrint["Passed check: projector does map to kernel of LL."];
        ,
            dPrint["Failed check: projector does not map to kernel of LL: ",
                   "norm of LL.projector == ", matrixNorm];
            Abort[];
        ];
        **)
        sm = (projector.rm // Simplify)
    ]


(* Explicit direct algorithm for non-semisimple systems *)
innerProductAlgorithm1z[rm_List, LL_, A_, m_, u_] :=
    Module[{basis, dualbasis, basissize, LLstar, image, kernel, orthogToIm, 
            K, Q, sm},
        basis = PolyFieldBasis[m, u];
        dualbasis = PolyFieldDualBasis[m, u];
        basissize = Length[basis];
        (* orthonormal basis for column space of LL: *)
        image = Select[Orthogonalize[Transpose[LL]],!AllTrue[#,PossibleZeroQ]&];
        dPrint["Finding representation of Lstar operator in this basis"];
        LLstar = Outer[ApplyND[#1, #2, u]&, dualbasis,
                       Map[L[ConjugateTranspose[A], #, u]&, basis], 1];
        (* As our chosen complement of Image(LL), use Kernel(LLstar)  *)
        dPrint["Finding basis for kernel of Lstar..."];
        kernel = Orthogonalize[NullSpace[LLstar]];
        dPrint["Checking complementary..."];
        (* verify that kernel(LLstar) is a complement to image(LL): *) 
        If[Length[image]+Length[kernel]==basissize && 
                MatrixRank[Join[image,kernel]]==basissize,
            dPrint["Passed. Image(LL), Kernel(LLstar) are complementary."];
        ,
            dPrint["Failed. Image(LL) and Kernel(LLstar) not complementary."];
            Abort[];
        ];
        dPrint["Finding basis for subspace orthogonal to Image(L)..."];
        orthogToIm =
            Orthogonalize[kernel - Map[Project[#,image]&, kernel]] // Simplify;
        (* Find that vector that lies in ker(LLstar) but has the same 
           projection on orthogToIm as our existing terms. That vector 
           will give the terms we retain in the transformed version *)
        dPrint["Solving linear system to find transformed terms."];
        If[Length[kernel]!=0,
            K = Transpose[kernel];
            Q = orthogToIm;
            sm = K.LinearSolve[Q.K,Q.rm];
        ,
            sm = Table[0,{basissize}];
        ];
        sm
    ]


InvSum[x_, maxOrder_Integer?Positive] :=
    Normal[Series[(1+x)^-1, {x, 0, maxOrder}]]

(* Apply the near-identity coordinate transformation  u -> U(y)
   to transform the contravariant field R(u)  *)
TransformContravariant[U_?MultiSeriesFieldQ, R_?MultiSeriesFieldQ] :=
    Module[{u, asympScaling, maxOrder, n, T, DT, factor1, factor2},
        u = U[[1]][[1]];
        asympScaling = unionScaling[U, R];
        maxOrder = R[[1]][[4]][[5]]/R[[1]][[4]][[6]] - 1;
        n = Length[u];
        T = Normal[U] - u;
        DT = D[T, {u, 1}];
        factor1 = InvSum[x, maxOrder] /.
            {1->IdentityMatrix[n], x->DT, Power->MatrixPower} //
            MultiSeries[#, asympScaling, maxOrder]&;
        factor2 = SubstituteField[R, U];
        factor1.factor2
    ]


(* Defines one iteration, simplifying terms of order m.
   Args:
       {R, U} where 
         R is our system so far,
         U is cumulative transformation used so far (from original system to R).
       m: the order of terms we are currently simplifying
       u: list of phase space variables, e.g. {u1, u2, u3}
       maxOrder: neglect terms of higher than this order

   Assumption: 
       R is already in normal form up to order m-1

   Returns:
       {S, Ucum} where 
       S is the transformed system now in normal form up to order m,
       Ucum is the new cumulative transformation (from original system to S).
*)
simplifyOrder[{R_?MultiSeriesFieldQ, U_?MultiSeriesFieldQ}, 
              m_Integer?(#>=2&),
              u_?SymbolListQ, 
              maxOrder_Integer?(#>=2&)] := 
    Module[{n, basis, basissize, dualbasis, rm, A, sm, tm, LL, image, Sm, 
            factor1, factor2, transformedSys, Tm, DTm, Um, S, Ucum, newrhs},
        dPrint["\n-----------------------------------------------------------"];
        dPrint["Now trying to simplify terms of order ", m, "."];
        n = Length[u];
        basis = PolyFieldBasis[m, u];
        basissize = Length[basis];
        dualbasis = PolyFieldDualBasis[m, u];
        (* Express existing mth order terms as a vector using our basis: *)
        rm = Map[ApplyND[#, OrderTerms[R, m], u]&, dualbasis];
        dPrint["Current order ", m, " terms: ",
              OrderTerms[R, m]//Expand//NN//MatrixForm];
        dPrint["Finding representation of homological operator L ",
               "in this basis."];
        A = D[OrderTerms[R, 1], {u}];
        LL = Outer[ApplyND[#1, #2, u]&, dualbasis, Map[L[A,#,u]&,basis], 1];
        (* dPrint["LL = ", LL//MatrixForm]; *)
        (* dPrint["condition number of LL is ",
               SingularValueList[LL] // Max[#]/Min[#]&]; *)
        If[MatrixRank[LL]==basissize,
            dPrint["Image of L is whole space, so can remove all order ",
                   m, " terms."];
            sm = Table[0, {basissize}];
        ,
            If[AllTrue[Eigenvectors[LL], !AllTrue[#,PossibleZeroQ]&],
                (* then there are no zero vectors in the list *)
                dPrint["Representation of L is semisimple. ",
                       "Choosing semisimple style."];
                sm = semisimpleAlgorithm1z[rm, LL, m, u];
                (* sm = semisimpleAlgorithm1a[rm, LL, m, u]; *)
            ,
                dPrint["Representation of L is NOT semisimple. ",
                       "Choosing inner product style."];
                sm = innerProductAlgorithm1z[rm, LL, A, m, u];
            ];
        ];
        sm = sm // Chop;
        Sm = sm.basis;
        dPrint["Transformed order ", m, " terms: ", Sm//NN//MatrixForm];
        dPrint["Solving linear system to find the required transformation..."];
        tm = LinearSolve[LL, rm-sm] // Simplify // Chop;
        (* nonlinear part of incremental transformation (expression in u): *)
        Tm = tm.basis // Expand; 
        (* incremental transformation (MultiSeries in u): *)
        Um = MultiSeries[u+Tm, u, maxOrder];
        (* composition of all transformations so far (MultiSeries in u): *)
        Ucum = ComposeTransformations[{U, Um}] // Simplify // Chop;
        dPrint["Incremental transformation: ",
               u//MatrixForm," \[Rule] ",Um//Normal//NN//MatrixForm];
        dPrint["Cumulative transformation so far: ",
               u//MatrixForm," \[Rule] ",Ucum//Normal//NN//MatrixForm];
        (* Find new transformed system S  *) 
        dPrint["Computing equations of motion in the new variables..."];
        transformedSys = TransformContravariant[Um, R];
        (* verify that the transformed system has order m terms equal to Sm: *)
        If[Norm[(OrderTerms[transformedSys, m] - Sm) /.Thread[u->0.5]] < 10^-6,
            dPrint["Passed. Transformation gives expected order ", m," terms."];
        ,
            dPrint["Failed. Transformation does not give the expected order ",
                   m, " terms. Got ", OrderTerms[transformedSys,m]//MatrixForm];
            Abort[];
        ];
        (* Orders 0 to m-1 are unchanged from R, and Sm may have more 
           precision than the order m terms of transformedSys so we will take 
           orders 0 to m-1 directly from R and order m directly from Sm   *)
        newrhs = OrderTerms[R, {0, m-1}] + Sm + 
                 OrderTerms[transformedSys, {m+1, maxOrder}];
        S = MultiSeries[newrhs, u, maxOrder] // Simplify // Chop;
        dPrint["New system: ", OverDot/@u//MatrixForm, " = ",
               Normal[S]//NN//MatrixForm, " + ",
               Superscript["O[|u|]", maxOrder+1]];
        {S, Ucum}
    ]


Options[NormalFormTransformation] =
    {Verbose->False,
     BifurcationParameters->{Global`\[Epsilon]},
     AsymptoticScaling->{Sqrt[Global`\[Epsilon]]},
     Extended->False};

(* Compute the normal form of system, to a specified order. 
    Args:
        rhs: the right hand side of the system (a vector expression in the xi)
        vars: list of phase space variables used in the original system, 
              e.g. {x1, x2, x3}
        newvars: list of new variable names to use in the transformed system
              e.g. {u1, u2, u3}
        maxOrder: compute normal form up to terms of this order

    Returns:
        {newrhs, trans} where
          newrhs is the transformed system in the new variables,
          trans is the normal form transformation (expressed as a list of 
            rules x -> f(u) mapping old variables to new).

    Options:
        Verbose: whether to print all working (default False)
        BifurcationParameters: which symbols in `rhs` should be interpreted as
          small bifurcation parameters (default Global`\[Epsilon])
        AsymptoticScaling: relative scaling of variables and bifurcation
          parameters in the asymptotic limit, used when truncating power series.
          e.g. {u1, u2, u3, Sqrt[Global`\[Epsilon]]} means to assume
          O(\[Epsilon]) == O(u_i^2) which suits a Hopf or Pitchfork bifurcation.
        Extended: whether to compute the normal form of the extended system,
          that is with phase space extended with dimensions for the (rescaled)
          bifurcation parameters and their equations \dot{\alpha} == 0,
          thus finding a transformation dependent on the bifurcation parameters.
          If False (the default) then the normal form will be found with
          respect to the dynamical variables only.
*)        
NormalFormTransformation[rhs_?VectorQ, 
                         vars_?SymbolListQ, 
                         newvars_?SymbolListQ,
                         maxOrder_Integer?Positive,
                         OptionsPattern[]] :=
    Module[{n, u, ua, syms, exponents, epsscale, bifParams, asympScaling, RHS,
            RHSseries, startSeries, S, U, identityTrans, newrhs, trans, nbif},
        verbose = OptionValue[Verbose];
        bifParams = OptionValue[BifurcationParameters];
        If[Head[bifParams]=!=List, bifParams={bifParams}];
        nbif = Length[bifParams];
        extended = OptionValue[Extended];
        asympScaling = OptionValue[AsymptoticScaling];
        asympScaling = Select[vars, FreeQ[asympScaling, #]&]~Join~asympScaling;
        dPrint["Deterministic system using asymptotic scaling ", asympScaling]; 
        If[Length[vars] != Length[newvars],
            Print["Number of new variables must match number of old variables"];
            Abort[];
        ];
        n = Length[vars]; (* dimension of phase space *)
        dPrint["Dimension of phase space is ", n];
        (* Internally use the symbols u1..un for the new variables. Use private
           context when generating symbols, to avoid clash with global names *)
        Block[{$Context="NormalForm`Private`"},
            u = Table[Symbol["u"<>ToString[i]], {i,n}];
            \[Alpha] = Table[Symbol["\[Alpha]"<>ToString[i]], {i,nbif}];
        ];
        RHS = rhs /. Thread[vars->u];
        asympScaling = asympScaling /. Thread[vars->u];
        If[extended,
            (* want normal form of system extended with bifurcation params *)
            ua = u~Join~\[Alpha];
            (* rescale bif params to symbols \[Alpha] of same order as u_i *) 
            {syms, exponents} = processVars[asympScaling];
            epsscale = Extract[exponents,
                               Position[syms, _?(MemberQ[bifParams, #]&)]];
            RHS = (RHS /. Thread[bifParams->\[Alpha]^(1/epsscale)]) ~Join~ {0};
            asympScaling = ua;
        ,
            (* else will compute normal form with respect to variables u only *)
            ua = u;
        ];
        dPrint["RHS: ", RHS//MatrixForm];
        (* First approximate the system locally to origin with power series *)
        RHSseries = MultiSeries[RHS, asympScaling, maxOrder] //Simplify//Chop;
        dPrint["Series approximation to the original deterministic system:\n",
              RHSseries // NN // MatrixForm];
        If[extended,
            startSeries = RHSseries;
        ,
            (* else transform based on system exactly at bifurcation point *)
            startSeries = RHSseries /. Thread[bifParams->0];
        ];
        (* Before starting, cumulative transformation is the Identity: u->u *)
        identityTrans = MultiSeries[ua, ua, maxOrder];
        (* Now invoke the main algorithm. Iteratively simplify terms at each
           order, from 2nd order to maxOrder: *)
        {S, U} = Fold[simplifyOrder[#1, #2, ua, maxOrder]&,
                      {startSeries, identityTrans},
                      Range[2, maxOrder]];
        If[extended,
            (* then de-augment *)
            newrhs = (Normal[S][[1;;-2]] /.
                      Thread[\[Alpha]->bifParams^epsscale]) // Simplify // Chop;
            trans = Thread[vars->(Normal[U][[1;;-2]] /.
                                 Thread[\[Alpha]->bifParams^epsscale])];
        ,
            (* else transform original system including bifurcation parameters*)
            newrhs = Normal[TransformContravariant[U, RHSseries]] //
                       Simplify // Chop;
            trans = Thread[vars->Normal[U]];
        ];
        {newrhs, trans} /. Thread[u->newvars]
    ]


(* Take truncated power series field expression in cylindrical coordinates that
   may involve integer powers >= -2 and remove any terms of order greater than 
   O[polarScaling]^maxOrder *)
truncatePolar[field_, polarScaling_List, maxOrder_Integer?NonNegative] :=
    Module[{n=Length[field], maxOrders},
        (* theta equation retains terms of order one less than other eqns: *)
        maxOrders = Table[maxOrder, {n}];
        maxOrders[[2]] = maxOrder - 1; 
        Table[Normal@MultiSeries[field[[i]], polarScaling, maxOrders[[i]]], 
              {i, 1, n}] // Simplify // Chop
    ];


Options[TransformNoisyHopf] =
    {Verbose->False,
     BifurcationParameters->{Global`\[Epsilon]},
     AsymptoticScaling->{Sqrt[Global`\[Epsilon]], Global`\[Sigma]},
     MaxOrder->3,
     Extended->False,
     Average->True,
     Rescale->True};

(* TransformNoisyHopf[rhs, {x1,...,xn}, {\[Sigma]1,...,\[Sigma]n}, 
   {\[Xi]1,...\[Xi]n}, {r,\[Theta]}, {new\[Xi]1,new\[Xi]2}] takes the stochastic
   dynamical system with right hand side rhs (expressed in variables {xi}, 
   small noise parameters {\[Sigma]i} and Langevin noise symbols {\[Xi]i} with 
   Stratonovich interpretation of any multiplicative noise) and transforms it to
   a simple circular 2 dimensional Hopf normal form system (expressed in polar 
   variables {r,\[Theta]} and new Langevin noise symbols {new\[Xi]1,new\[Xi]2}).

   Returns:
       2D polar vector field representing transformed system in a standard form

   Options:
       Verbose: whether to print all working (default False)
       BifurcationParameters: which symbols in `rhs` should be interpreted as
         small bifurcation parameters (default Global`\[Epsilon])
       AsymptoticScaling: relative scaling of variables and bifurcation
         parameters in the asymptotic limit, used when truncating power series.
         e.g. {u1, u2, u3, Sqrt[Global`\[Epsilon]]} means to assume
         O(\[Epsilon]) == O(u_i^2) which suits a Hopf or Pitchfork bifurcation.
       Extended: whether to compute the normal form of the extended system,
         that is with phase space extended with dimensions for the (rescaled)
         bifurcation parameters and their equations \dot{\alpha} == 0,
         thus finding a transformation dependent on the bifurcation parameters.
         If False (the default) then the normal form will be found with
         respect to the dynamical variables only.
       Average: whether to average around the cycle. (default True)
       Rescale: whether to linearly rescale the radial variable to make the
         coefficient of the R^3 term -1. (default True)

   TODO: currently it is assumed that the linear part of the system has already
   been transformed to Jordan real form, with Hopf in first two variables. 
   Should automate that instead. *)
TransformNoisyHopf[rhs_?VectorQ, 
                   vars_?SymbolListQ, 
                   \[Sigma]_?SymbolListQ,
                   \[Xi]_?SymbolListQ,
                   {r_Symbol, \[Theta]_Symbol},
                   {new\[Xi]1_Symbol, new\[Xi]2_Symbol},
                   OptionsPattern[]] :=
    Module[{maxOrder, bifParams, n, asympScaling, deterministicScaling, A,
            nDpolarScaling, polarScaling, polarVars, cylVars, centerEigs,
            deterministicRhs, newrhs, trans, transformedCartesian, fullPolar,
            extended, average, rescale, result},
        verbose = OptionValue[Verbose];
        maxOrder = OptionValue[MaxOrder];
        extended = OptionValue[Extended];
        average = OptionValue[Average];
        rescale = OptionValue[Rescale];
        bifParams = OptionValue[BifurcationParameters];
        If[Head[bifParams]=!=List, bifParams={bifParams}];
        n = Length[vars]; (* dimension of phase space *)
        Block[{$Context="NormalForm`Private`"},
            u = Table[Symbol["u"<>ToString[i]], {i,n}];
        ];

        (* which asymptotic limit to use when truncating series: *)
        (* e.g. {u1, ..., un, eps^(1/2), sigma^(1/2)} *)
        asympScaling = OptionValue[AsymptoticScaling];
        asympScaling = Select[vars, FreeQ[asympScaling, #]&]~Join~asympScaling;
        asympScaling = asympScaling /. Global`\[Epsilon]->bifParams[[1]] /.
                                       Global`\[Sigma]->\[Sigma] //
                                       Flatten;
        dPrint["Stochastic system using asymptotic scaling ", asympScaling]; 
        (* e.g. {u1, ..., un, eps^(1/2} *)
        deterministicScaling = 
            Select[asympScaling, And@@Through[Thread[FreeQ[\[Sigma]]][#]]&];
        asympScaling = asympScaling /. Thread[vars->u];
        (* e.g. {r, u3, ..., un, eps^(1/2), sigma^(1/2)} *)
        nDpolarScaling =
            {r}~Join~Select[asympScaling, (FreeQ[#,u[[1]]]&&FreeQ[#,u[[2]]])&];
        (* e.g. {r, eps^(1/2), sigma^(1/2)} *)
        polarScaling =
            {r}~Join~Select[asympScaling, And@@Through[Thread[FreeQ[u]][#]]&];

        $Assumptions = $Assumptions && And@@Thread[\[Sigma]>=0] &&
                       r>=0 && Element[\[Theta], Reals] &&
                       And@@Map[Element[#, Reals]&, bifParams];
        polarVars = {r, \[Theta]}; 
        cylVars = polarVars~Join~u[[3;;]];

        deterministicRhs = rhs /. Thread[\[Xi]->0];
        {newrhs, trans} = 
            NormalFormTransformation[deterministicRhs, vars, u, maxOrder, 
                                     Verbose->verbose, 
                                     BifurcationParameters->bifParams,
                                     AsymptoticScaling->deterministicScaling,
                                     Extended->extended];

        (* validate that we are at Hopf point in first 2 variables *)
        A = D[newrhs, {u}] /. Thread[u->0] /. Thread[bifParams->0];
        centerEigs = Eigenvalues[A[[1;;2,1;;2]]];
        If[!PossibleZeroQ[Abs[centerEigs[[1]] - Conjugate[centerEigs[[2]]]]] ||
                !PossibleZeroQ[Re[centerEigs[[1]]]],
            Print["First two variables not at Hopf bifurcation: eigenvalues=",
                  centerEigs];
            Abort[]; 
        ];

        (* Utility functions to collect terms for easy-to-read output: *)
        Ar[expr_] := Arrange[expr, Flatten@Outer[Times, \[Xi], \[Sigma]]];
        ArPolar[expr_] := ArrangePolar[
                expr, Flatten@Outer[Times, \[Xi], \[Sigma]], {\[Theta]}];

        dPrint["\n-----------------------------------------------------------"];
        dPrint["Transformed deterministic system, polar coordinates: ",
               OverDot /@ cylVars // MatrixForm, " = ", 
               ToPolar[newrhs, u, polarVars] //
                   NN // TrigReduce // Simplify // Chop // MatrixForm, " + ", 
               Superscript["O[r]", maxOrder+1]];

        (* Having found the desired change of coordinates, 
           Now look at the full system in these new coordinates, including 
           noise and bifurcation parameter *)
        transformedCartesian = 
            TransformContravariant[
                MultiSeries[trans[[All,2]], u, maxOrder],
                MultiSeries[rhs /. Thread[vars->u], asympScaling, maxOrder]
            ];
        (* TODO make versions of ToPolar and ToCartesian for series fields *)
        dPrint["transformed stochastic system, cartesian:"];
        dPrint[OverDot/@u//MatrixForm, " = ",
               transformedCartesian // Ar // NN // Chop // MatrixForm];
        transformedPolar = truncatePolar[ToPolar[Normal[transformedCartesian],
                u, polarVars], nDpolarScaling, maxOrder];
        dPrint["transformed stochastic system, polar:"];
        dPrint[OverDot/@cylVars//MatrixForm, " = ",
               transformedPolar // ArPolar // NN // Simplify // MatrixForm];

        (* Use an averaging approximation of the stable variables to reduce
           the system to two dimensions: *)
        reducedCartesian =
            NoisyHopfReduce2D[Normal[transformedCartesian], u, \[Sigma], \[Xi],
                              BifurcationParameters->bifParams, 
                              MaxOrder->maxOrder];
        dPrint["Two-dimensional reduced system, cartesian: "];
        dPrint[OverDot/@u[[1;;2]]//MatrixForm, " = ",
               reducedCartesian // Ar // NN // Chop // MatrixForm];
        result = truncatePolar[ToPolar[Normal[reducedCartesian],
                                       u, polarVars], nDpolarScaling, maxOrder];
        dPrint["Two-dimensional reduced system, polar: "];
        dPrint[OverDot /@ polarVars // MatrixForm, " = ",
               result // ArPolar // NN // MatrixForm];

        If[average,
            (* Use averaging of the diffusion process around the cycle to find
               a circularly symmetrical approximation for longer time scales. *)
            result = AverageCycle[result, polarVars, polarScaling,
                                  \[Sigma], \[Xi], {new\[Xi]1, new\[Xi]2},
                                  BifurcationParameters->bifParams,
                                  MaxOrder->maxOrder];
        ];

        If[rescale,
            (* rescale r so that the coefficient of the r^3 term is -1 *)
            result = HopfRescaleRadius[result, polarVars];
        ];
        result
    ]


Options[NoisyHopfReduce2D] =
    {BifurcationParameters->{Global`\[Epsilon]},
     MaxOrder->3};

(* Reduce the noisy Hopf normal form system to two dimensions (the center 
   eigenspace variables u1 and u2) by approximating the stochastic coupling 
   terms by zero. It is assumed the deterministic terms are already in Hopf
   normal form and are expressed in cartesian coordinates {u1, u2, ...} with
   the center eigenspace variables listed first.
        
   TODO: currently assumes the frequency of stable eigenvalues are not close
   to resonance with the Hopf frequency. Should test for this explicitly.

   TODO: also assumes there is a non-zero additive noise term to which these
   zero-mean coupling terms are added. Should test this explicitly. *)
NoisyHopfReduce2D[rhs_?VectorQ,
                  u_?SymbolListQ,
                  \[Sigma]_?SymbolListQ,
                  \[Xi]_?SymbolListQ,
                  OptionsPattern[]] :=
    Module[{maxOrder, bifParams, cVars, sVars, A, reStableEigenvals,
            weakestStableEigenvalue, result},
        maxOrder = OptionValue[MaxOrder];
        bifParams = OptionValue[BifurcationParameters];
        cVars = u[[1;;2]]; (* center eigenspace variables *)
        sVars = u[[3;;]]; (* stable eigenspace variables  *)
        (* Check criterion for approximate decoupling of the center variables *)
        A = D[rhs[[3;;]], {sVars}] /. Thread[u->0]
                                   /. Thread[bifParams->0] /. Thread[\[Xi]->0];
        reStableEigenvals = Diagonal[A];
        weakestStableEigenvalue = Max[reStableEigenvals];
        If[Re[weakestStableEigenvalue]>0,
            Print["Error. There is an unstable eigenspace. ",
                  "System will not de-couple."];
            Abort[];
        ];
        dPrint["Assuming now that \[Sigma]^2/",
               NN[(-Re[weakestStableEigenvalue])^(3/2)], " << 1  so that the ",
               "oscillations approximately decouple from the stable variables"];
        (* Approximate by omitting the stochastic coupling terms *)
        result = rhs[[1;;2]] /. Thread[sVars->0] // Simplify
    ]


Options[AverageCycle] =
    {BifurcationParameters->{Global`\[Epsilon]},
     MaxOrder->3};

(* Remove angle dependency from a 2D polar SDE system by averaging around the
   cycle. The averaging is done on the corresponding Fokker Planck equation
   in a rotating frame, then a choice of SDE is made that corresponds to the
   averaged system. 
   
   Note: this approximation will destroy fidelity of system behaviour for
   phenomena that occur on a shorter time scale than one period. The benefit
   is a much simpler system. *)
AverageCycle[rhs_?VectorQ,
             polarVars_?SymbolListQ,
             polarScaling_List,
             \[Sigma]_?SymbolListQ,
             \[Xi]_?SymbolListQ,
             {new\[Xi]1_Symbol, new\[Xi]2_Symbol},
             OptionsPattern[]] :=
    Module[{maxOrder, bifParams, r, \[Omega], rotPolar, stratonovichRHS,
            stratonovichDrift, noiseTerms, G, stratToIto, itoDrift, itoRHS,
            diffusion, diffusionAv, driftAv, Gav, noiseTermsAv, itoToStrat,
            stratonovichDriftAv, avPolarSys},
        maxOrder = OptionValue[MaxOrder];
        bifParams = OptionValue[BifurcationParameters];
        {r, \[Theta]} = polarVars;

        (* Utility functions to collect terms for easy-to-read output: *)
        Ar[expr_] := Arrange[expr, Flatten@Outer[Times, \[Xi], \[Sigma]]];
        ArPolar[expr_] := ArrangePolar[
                expr, Flatten@Outer[Times, \[Xi], \[Sigma]], {\[Phi]}];

        (* Find angular frequency at Hopf point: *)
        \[Omega] = (rhs[[2]] /. Thread[Join[u,\[Xi],bifParams]->0] // Expand)
                /. r->0 // FullSimplify // Chop;
        dPrint["\[Omega] == ", \[Omega]];
        (* Change to rotating frame:*)
        rotPolar = rhs - {0,\[Omega]} /. \[Theta]->\[Phi]+\[Omega] t;
        dPrint["Reduced system in rotating frame  \[Phi]=\[Theta]-\[Omega]t:"];
        dPrint[OverDot/@{r,\[Phi]} // MatrixForm, " = ",
               rotPolar // ArPolar // NN // MatrixForm];

        (* Convert the above result to Ito form *)
        stratonovichRHS = rotPolar;
        (* the deterministic part: *)
        stratonovichDrift = stratonovichRHS /. Thread[\[Xi]->0];
        (* the noise part, including additive and multiplicative noise terms: *)
        noiseTerms = stratonovichRHS - stratonovichDrift;
        G = D[noiseTerms, {\[Xi]}];
        dPrint["Matrix of noise coefficients G = ", G//NN//MatrixForm]; 
        (* transpose is needed to apply trace to correct ranks of the tensor: *)
        stratToIto = 1/2 Tr[Transpose[D[G,{{r,\[Phi]}}].G,{3,1,2}],Plus,2];
        dPrint["Stratonovich to Ito drift adjustment: ",
               stratToIto//NN//MatrixForm];
        itoDrift = stratonovichDrift + stratToIto // Expand // Chop;
        itoRHS = itoDrift + noiseTerms;
        dPrint["Ito form:"];
        dPrint[OverDot /@ {r, \[Phi]} // MatrixForm, " = ",
               itoRHS // ArPolar // NN // MatrixForm];
        dPrint["----------Removing insignificant terms... ----------"];
        dPrint[OverDot /@ {r, \[Phi]} // MatrixForm, " = ",
               truncatePolar[itoRHS, polarScaling, maxOrder]//NN//MatrixForm];

        (* Convert to Fokker-Planck equation and average around the cycle *)
        diffusion = G.Transpose[G] // Expand // Chop;
        dPrint["Fokker-Planck diffusion matrix D = ",diffusion//MatrixForm];
        dPrint["Fokker-Planck drift: ",itoDrift//MatrixForm];
        dPrint["----------Removing insignificant terms... ----------"];
        diffusionT = Normal@MultiSeries[diffusion, polarScaling, maxOrder*2];
        driftT = truncatePolar[itoDrift, polarScaling, maxOrder];
        dPrint["Fokker-Planck diffusion matrix D = ",diffusionT//MatrixForm];
        dPrint["Fokker-Planck drift: ",driftT//MatrixForm];
        (* TODO might need a bit of validation before integrating *)
        diffusionAv = Integrate[diffusion,{\[Phi],-\[Pi],\[Pi]}]/(2\[Pi]) //
            FullSimplify // Chop;
        driftAv = Integrate[itoDrift,{\[Phi],-\[Pi],\[Pi]}]/(2\[Pi]) //
            FullSimplify // Chop;
        dPrint["----------After averaging: --------------"];
        dPrint["Averaged Fokker-Planck diffusion: ",diffusionAv//MatrixForm];
        dPrint["Averaged Fokker-Planck drift: ",driftAv//MatrixForm];
        dPrint["----------Removing insignificant terms... ----------"];
        diffusionAv = 
            Normal@MultiSeries[diffusionAv, polarScaling, maxOrder*2];
        dPrint["Averaged Fokker-Planck diffusion: ",diffusionAv//MatrixForm];
        driftAv = truncatePolar[driftAv, polarScaling, maxOrder];
        dPrint["Averaged Fokker-Planck drift: ",driftAv//MatrixForm];
         
        (* Choose an Ito SDE corresponding to the Fokker-Planck equation: *)
        Gav = Transpose[CholeskyDecomposition[diffusionAv]] // Simplify // Chop;
        noiseTermsAv = Gav.{new\[Xi]1, new\[Xi]2};
        dPrint["Averaged normal form system in Ito form: "];
        dPrint[OverDot/@{r,\[Phi]}//MatrixForm, " = ",
               driftAv+noiseTermsAv // Ar // NN // MatrixForm];

        dPrint["----------Removing insignificant terms... ----------"];
        noiseTermsAv = truncatePolar[noiseTermsAv, polarScaling, maxOrder];
        dPrint["Averaged normal form system in Ito form: "];
        dPrint[OverDot/@{r,\[Phi]}// MatrixForm, " = ",
               driftAv+noiseTermsAv // Ar // NN // MatrixForm];

        (* Convert SDE to Stratonovich form *)
        (* transpose is needed to apply trace to correct ranks of the tensor: *)
        itoToStrat = 
            -1/2 Tr[Transpose[D[Gav,{{r,\[Phi]}}].Gav,{3,1,2}],Plus,2] // 
                Simplify // Chop;
        dPrint["Ito to Stratonovich drift adjustment: ",
               itoToStrat // NN // MatrixForm];
        itoToStrat = truncatePolar[itoToStrat, polarScaling, maxOrder];
        dPrint["truncated Ito to Stratonovich drift adjustment: ",
               itoToStrat // NN // MatrixForm];
        stratonovichDriftAv = driftAv + itoToStrat;
        avPolarSys = stratonovichDriftAv + noiseTermsAv;
        dPrint["Averaged normal form system in Stratonovich form: "];
        dPrint[OverDot/@{r,\[Phi]}//MatrixForm, " = ",
               avPolarSys // Ar // NN // MatrixForm];
        dPrint["----------Removing insignificant terms... ----------"];
        avPolarSys = truncatePolar[avPolarSys, polarScaling, maxOrder];
        dPrint["Averaged normal form system in Stratonovich form: "];
        dPrint[OverDot/@{r,\[Phi]}//MatrixForm, " = ",
               avPolarSys // Ar // NN // MatrixForm];

        (* change back to non-rotating frame: *)
        result = avPolarSys + {0, \[Omega]}
    ]


(* Rescale r so that the coefficient of the r^3 term is -1
   Assumes the system is expressed in polar or cylindrical coordinates *)
HopfRescaleRadius[rhs_?VectorQ, cylVars_?SymbolListQ] :=
    Module[{r, coeff, result},
        r = cylVars[[1]];
        coeff = SeriesCoefficient[rhs[[1]], {r, 0, 3}];
        If[PossibleZeroQ[coeff],
            Print["Handling degenerate Hopf bifurcation not yet implemented."];
            Abort[];
        ];
        (* TODO: for degenerate Hopf, rescale so coefficient of r^5 is -1 *)
        result = rhs /. r -> r/Sqrt[-coeff];
        result[[1]] = result[[1]] * Sqrt[-coeff];
        result // Simplify
    ]


End[]

EndPackage[]
