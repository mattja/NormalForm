<<NormalForm`;

x={x1,x2,x3};
 
ms1=MultiSeries[Cos[x1] Cos[\[Epsilon]] , {x,Sqrt[\[Epsilon]]}, 7]


ms1//Normal


ms1/.\[Epsilon]->0
ms1/.\[Epsilon]->0/.x1->0/.x2->0
ms1/.\[Epsilon]->0/.Thread[x->0]


ms1
NormalForm`Private`OrderTerms[ms1,{0,4}]
 

msf1=MultiSeries[{Cos[x1] Cos[\[Epsilon]],x3^2,x1+x2}, {x,Sqrt[\[Epsilon]]}, 5];
msf1//MatrixForm


msf1//Normal//MatrixForm
msf1/.\[Epsilon]->0//MatrixForm
msf1/.Thread[x->0]//MatrixForm


msf1//MatrixForm
NormalForm`Private`OrderTerms[msf1, {0,3}]//MatrixForm
 

ms1
ms2 = MultiSeries[x1^2-4x2^3+7, {x,Sqrt[\[Epsilon]]},7]


ms1+ms2
ms1 ms2
 

msf1//MatrixForm
msf2=MultiSeries[{2x1,x2,x2^2},{x,Sqrt[\[Epsilon]]},7];
msf2//MatrixForm
NormalForm`Private`SubstituteField[msf1,msf2]//MatrixForm
