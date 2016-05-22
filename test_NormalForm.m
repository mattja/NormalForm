<<NormalForm`;


(* define some variables and noise parameters: *)
x = {x1,x2,x3};
u = {u1,u2,u3};
\[Sigma] = {\[Sigma]1,\[Sigma]2,\[Sigma]3};
\[Xi] = {\[Xi]1,\[Xi]2,\[Xi]3};

(* define our original system right-hand-side: *)
rhs={\[Epsilon] x1-\[Omega] x2-x3 x1+\[Sigma]1 \[Xi]1,
     \[Omega] x1+\[Epsilon] x2-x3 x2,
     -\[Lambda] x3+4 (x1^2+x2^2)+\[Sigma]3 \[Xi]3};

$Assumptions = \[Omega]\[Element]Reals && \[Lambda]>0;
\[Omega]=9;
\[Lambda]=2;

(* First test transforming a noise-free system *)
noiseFree = rhs /. Thread[\[Sigma]->0];
(* Compute normal form to third order: *)
{newrhs, trans} = NormalFormTransformation[noiseFree, x, u, 3];

Print["transformed system: ", OverDot/@u//MatrixForm," = ",newrhs//MatrixForm];
Print["transformation: ", trans // MatrixForm];


(* Now analyze the full system including noise. Output is a Stratonovich SDE: *)
result = TransformNoisyHopf[rhs, x, \[Sigma], \[Xi],
                            {r, \[Theta]}, {\[Xi]A, \[Xi]B}];

Print[OverDot/@{r,\[Theta]}//MatrixForm, " = ", result//MatrixForm];
