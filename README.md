(*Module-1*)

(*----------Basic code to Solve Quantum annealing problems in Mathematica----------*)

(*Quantum Annealing in Mathematica:Test the code for n>=3 Qubits, n>=2 Qubits needs some modifications in logic: The code will be used to study various types of feature selections and quantum annealing problems...........*)

\[Sigma]x = {{0, 1}}, {{1, 0}};

\[Sigma]y = {{0, -I}}, {{I, 0}};

\[Sigma]z = {{1, 0}}, {{0, -1}};

id2 = IdentityMatrix[2];

(*code for Initial Hamiltoanin which is prepared in ground state*)

Clear[iNitalHamiltonain];
iNitalHamiltonain[n_] := Module[{l1, l2, d, i2, mYexpression1},
  If[n == 0 || n == 1, Print["Provide (n\[GreaterEqual]2) Qubits"],
   l1 = Table[1, {i, 1, n}];
   l2 = Table[\[Sigma], {i, 1, n}];
   d = DiagonalMatrix[l1*l2];
   i2 = d /. {0 -> id};
   mYexpression1 = Total[Apply[KroneckerProduct, i2, {1}]];(*Testing: 
   one can check here the correct tensor product expression*)
   (-1)*mYexpression1 /. {\[Sigma] -> \[Sigma]x, id -> id2}]
  ]

iNitalHamiltonain[3] // MatrixForm;

(*code for the Magnetific field term \!\(
\*SubscriptBox[\(\[Sum]\), \(i\)]\(
\*SubscriptBox[\(h\), \(i\)]
\*SubscriptBox[
SuperscriptBox[\(\[Sigma]\), \(z\)], \(i\)]\)\)*)

Clear[mAgnetifcField];
mAgnetifcField[n_] := Module[{i1, l1, l2, mYexpression1, mYexpression2},
  SeedRandom[12345];
  If[n == 0 || n == 1, Print["Provide (n\[GreaterEqual]2) Qubits"],
   l1 = Table[RandomReal[{0, 1}], {i, 1, n}];
   (*List l1 is the list of Magnetic filed strengths*)
   l2 = Table[\[Sigma], {i, 1, n}];
   i1 = DiagonalMatrix[l2] /. {0 -> id};
   mYexpression1 = Total[l1*Apply[KroneckerProduct, i1, {1}]];
   (*Testing: one can check here the correct tensor product expression*)
   mYexpression2 = mYexpression1 /. {\[Sigma] -> \[Sigma]z, id -> id2}]
  ]

mAgnetifcField[3] // MatrixForm;

(*code for the Interaction term \!\(
\*SubscriptBox[
SuperscriptBox[\(\[Sum]\), \(n\)], \(\(\(<\)\(i\)\), j, i \[NotEqual] j > \)]\(
\*SubscriptBox[\(J\), \(ij\)]
\*SubscriptBox[
SuperscriptBox[\(\[Sigma]\), \(z\)], \(i\)]\ \[TensorProduct]
\*SubscriptBox[
SuperscriptBox[\(\[Sigma]\), \(z\)], \(j\)]\)\)*)

Clear[interaction];
interaction[n_] := Module[{i, j, l1, l2, l3, eXpression1, eXpression2},
  SeedRandom[4567];
  l1 = Table[id, {i, 1, n}];
  If[n == 0 || n == 1 , Print["Provide (n\[GreaterEqual]2) Qubits"],
   If[n == 2, 
    Subscript[J, 
     12]*(KroneckerProduct[\[Sigma]z, id2] + KroneckerProduct[id2, \[Sigma]z]),
    (*l2=Table[If[i<j,Subscript[J,i,j],0],{i,1,n},{j,1,n}];*)
    (*Basically l2 is an upper-traingular matrix with Cofficients Subscript[J,
     ij]*)
    (*Lets create l2 as a random upper-
    traingular matrix along with Diagonalwith SeedRandom[12345] fixed: 
    In that sence our Mutual Information matrix is fixed*)
    l2 = Table[
      If[i <= j, RandomVariate[NormalDistribution[0, 1]], 0], {i, 1, n}, {j, 
       1, n}];
    l3 = Table[
      If[i < j, 
       Apply[KroneckerProduct, ReplacePart[l1, {{i}}, {{j}} -> \[Sigma]]], 
       0], {i, 1, n}, {j, 1, n}];
    eXpression1 = Total[Total[l2*l3]];
    (*Testing: one can check here the correct tensor product expression*)
    eXpression2 = eXpression1 /. {\[Sigma] -> \[Sigma]z, id -> id2 }]]
  ]

interaction[3] // MatrixForm;

Clear[constraint];
constraint[n_] := Module[{l1, l2, l3, l4, l5, l6, l7},
  l1 = DiagonalMatrix[Table[\[Sigma], {i, 1, n}]] /. {0 -> id};
  l2 = Apply[KroneckerProduct, l1, {1}];
  l3 = Total[l2];(*\!\(
\*SubscriptBox[
SuperscriptBox[\(\[Sum]\), \(n\)], \(i\)]
\*SubscriptBox[\(\[Sigma]\), \(i\)]\)*)
  l4 = Table[id, {i, 1, n}];
  l5 = k*Apply[KroneckerProduct, 
     l4];(*k(I\[TensorProduct]I\[TensorProduct]I......)*)
  l6 = \[Alpha]*(l3 - l5)^2;
  l7 = l6 /. {\[Sigma] -> \[Sigma]z, id -> id2}
  ]

constraint[3] // MatrixForm;

(*code for Problem Hamiltonain Subscript[H^1, p]=Subscript[H, p]+Subscript[\
\[Alpha]H, constraint]*)

Clear[pRoblemHamiltonain];
pRoblemHamiltonain[n_] := Module[{},
  mAgnetifcField[n] + interaction[n] + constraint[n](*With Constraint*)
  ]

pRoblemHamiltonain[3] // MatrixForm;
(*Final Hamiltinan H(s) is the function of paremeters s whcih is used to \
switch the initial Hamiltonian to problem hamiltonian*)


(*Quantum annealing Hamiltonain H(s)*)

Clear[qAnnealing];
qAnnealing[n_] := Module[{},
  (1 - s)*iNitalHamiltonain[n] + s*pRoblemHamiltonain[n]
  ]

qAnnealing[3] // MatrixForm;

(*Now find out the differnce between two initial eiganvalues of the Quantum \
Annealing Hamiltonain between 0\[LessEqual]s\[LessEqual]1*)

Clear[eIgenqAnnealing];
eIgenqAnnealing[n_, k1_] := Module[{e1, e2, dim, diff},
  e1 = Table[
    Table[Min[Eigenvalues[qAnnealing[n] /. {k -> k1}]], {s, 0, 1, 
      0.01}], {\[Alpha], 0, 1, 0.2}];
  dim = Dimensions[e1][[1]];
  e2 = Table[
    Table[RankedMin[Eigenvalues[qAnnealing[n] /. {k -> k1}], 2], {s, 0, 1, 
      0.01}], {\[Alpha], 0, 1, 0.2}];
  Table[ListLinePlot[{e1[[i]], e2[[i]]}, DataRange -> {0, 1}, Frame -> True, 
    FrameLabel -> {"s", 
      "(\!\(\*SubscriptBox[\(E\), \(1\)]\),\!\(\*SubscriptBox[\(E\), \
\(0\)]\))"}, PlotLabel -> Row@{"\[Alpha]=", N[((i - 1)/((dim) - 1))]}], {i, 1,
     dim}](*Here we Plot the Ist and 2nd eigenvalues spectrum*)
  ]

(*Results:......for different values of K features and with 0\[LessEqual]\
\[Alpha]\[LessEqual]1*)

eIgenqAnnealing[3, 5]

{\!\(\*

{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1lAtQVGUUx+/dXUARfJIUpWCAgikyOCqi8VdkVDQTUEuUKBUwywQMSjFm
RHmoOCOiKDmtj1htzBFkErJUvg/lIQYCAq6CArvAsjwWTIXKu4/cvd+embs7
e/fec87//zvnm7Y1LjxGwnFczJvL/G0N3vI5iu5zKlihSk5DfN5kRckdO/b7
JCI8fuu4lW5Hb0xN7w+7cRb5IymrnUKs/yuw68rijFIHO3rmR3NcxgW307hd
b8ueL4Bfz+Dcibm2dCDMfKcIwcIqHxJhy96/jtrSiNGVU22pvTI//kP73yGP
zrrm1WXD8v2B874rpvVetqHekcfuDHvdAvEJ3Dk63oblL0VxcFtJ1nwbqnjz
tjKfQl6jXbzbIGP1ynBxvm84KZfRuZa4i2WH90cnHZWx+uWQ/1NiOLVeRu8O
e72pUIG1vo6V06fIWD9VcEmWZ8zQSGm4peA9VM/xnPVToZT1V42N18fcz9wr
papkc8L7qMq4Ltctk7J+azCp2H1iw1gpTTC3F1+L3KP6u/NaJKz/B3hg25nu
eklCRRr12OcUEpOzW8L01KP3QFFhDiR0paVAA0bu7P7Xw1HC9DVgeYjM9eNW
npqrRR57CD5VtWHSrzzT24iYEMeZyXt5apGjbIT/8SVZh0J4pr8JfFFtU7AL
T13N7UxtRmjcEmVxP8f8aMZ/j4lafZuj22PN8QhBQ34PG7M55o8SJ73dnY9H
c9SczalAiUcb1la6L+SYX49xs6ZZkTuWoyPmdMOPQcvOfK2uNhHRvyewb3bN
PrjHRAItCVswPu223QxvExH9bEHCgbohTauRWNpLb8X5gbGr2rONRPT3KZRl
XeN8VhpJbY05nsLncixe8kYi+v0MXmvWtYeWGshblgbbsOC7xpGvUgxE9L8N
NzMXJG2HgVjSRbajKiXBZpfMQEQeHTh3ovLeLzV6IurpwOyrK14F5umJyKcD
M18GrInarieivg4on79K9Q3QE5GXCjXLmqPax+uZXhUmjwpW5/QLROSngjTq
2Pz4+wLTr0ZrzBh53lWBiDzVqHiyaJNXrsD8UOP7cYtK5qUKROTbiSDPt6fX
JwjMn05Iz5UVmWIFIuLpROLovy+WfyEwv7ogvFd+xvlzgYj8u8B5L63r3SIw
/7pQUBbb4LlDIOI8dMOt0mPbtUSB+dmNuqDezB1pAhHnoxutE5oUs/ME5q8G
uftNZx8VCkScFw1MW9WaxX8JzG8N0Hgh2r9PIOL89CB05qLACw565n8PxruY
EkP99EScpx5onyZsdtusZzy0aH/3/Z6BDCsPLY5G7Ok8UWzlocWm85PGDGus
PLQ4qFvo+MLFwHj04rXzN2m7wgyMRy+q8zLVa48YGI9eXJJIP0mpMDAefQgr
XH2lXWIk2RYeffjBduOsLUFGxqMP4Uke+3UHjYxHP16fXpiVUGWd136s8zGF
tTqYGI9+7FGEPXtnvYnxGMCn/KlpU+QmxmMAHyk+0zb3mBiPAfhfXN7m4sKx
/dShYv23O52Xcuw80WF6/jaHli+t+6qD54TEpMjjHBV5DKIp5eXPeX9a93cQ
fsO+dYe6OHb+DCLAzT9uxjie7fMQ0lPd58QF8FTkMYTOKt3gplie7fcQXjRc
Jd05PDuvnuNI1Acqe8rT/wFSxzMx
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtMU3cUxhtl6mAYHwWn4oPJeBhcxpjOKOtnFAJuVQEXJxHkMcH5nDjd
eCybm0qA4HyLgqkzFMUwga1qIUP/DgSRxwSFFRCQR2lLKdDbe1sQcZmtpye5
ue9zvu/7Hfe4b8LjJ4lEovjXh/VsrzfXRpYqLgruSTkGpggaf/yO/f4cKhRV
szvnGlnpwuODYaUyHN/ReN7Ly/5eDrfby6qLPzaynEvWuoH4i79qDq+1f18E
zwyDf0qokRnCrE/+QNzL5v7y7fb/b2FptGyZZJ+ROaryDnzqqETszGPL/0u1
9ytDVvC1Ra8yjcwn8mSF2bscIYV6h9U59v73EDBaxitvGJn89d+qvPvITL4w
/G2Zfd7f0PSl6fc/MjJ/W1UiQxraUNBmn/8AoTeS97jrjazS7P16QhV+3Dsr
rXXcruchrvgtr2lw4li4bWANLsolD0QLONL3CJ5CsfHABxzrSbE2rMWk1Qr5
gjUc6a1HV/gn4inhHEu0yjvQAKdj2Z/77OBI/z+Yn1qxNO07jr2h0Qh9/5Po
+Rkc+WnEFwucYrW5HAuxDWhCZYxHgq6II39N+D2lJ8atgmPWaZEnn8D1s2bh
SAtHfp+itSFo/+wBjtnsqJ7ioqdfaeMER/6bEdm9Luz2DBNbZJWzsAWH57fh
roeJ8mhBSJzeQb3SxHYmWOtfrK8fO+GxwUT5qHBw68T7ybEmZu0mLlIhUfJ9
/LPDJsqrFdKCbvfATBOzWNuZW7F6lzKyUGai/NpwVqH4aJrCxCS2hu24lX6/
LvShifJsx16lXJryzMRs8o4/w6ozp+uSRkyUbwd+DvBdj8k8a6i3Vgf0G11K
mCtPeXfi8eD+Tr0Pz1xsArvgcH3H7ZIAnvLvwvo7K4JFm3hmaxf5HHtmFTvU
xPDEoxs6petMSyJPfrqxJPdcYdIvPPHphqzGu1Ryhid/3XCunrxr5VWeePXA
b7Z86tZinvz24OqFtrYz5Tzx60G71815vTU8+e9FaEHW+Lpmnnj2QrQ2/3x5
F0959CLqnLpui44nvn0IOpKn9uB4yqcP+7KkA2te8MS7Dx0f3r1WLxIoLzX0
L6Ok3FSB+KtxKmqioGu6QPmpsWFz+/RaF4H2oR+RklyjwU2gPPtxa/fBjUc9
BNqPfizOvuxc6ytQvhoUPq3e/Gq5QPuiQc0hn1db1giUtwY7pxTCWSrQ/mgR
2OufHhYhUP5a7OK8czZ9LdA+adF5QSL4JgvEQ4f0PGXytCyBeOhwlqvaNv6b
QDx0mJPuGSAuFYiHDn0dK47ENgnEYwCb7slqtQaBeAwgYNGq09cdzcRjANPS
FlddXmomHnrsa4o/yqRmdsrGQw8li6h2SDQTDz3ES36YEZ1tJh6DyHlucGLM
TDwGMeXFzXlivZl4DOKO208tUlcL8TAguzMi6ctAC/EwYOHNSwHvHbIQDwME
5ZWrufkW4jEEZYQ0/89WC/EYQoR5YmiL8yjxGIJWpdEmrRslHsMIn1uiejt1
lHgMozlujvqlYpR4DOMv323vrh0eJR4jaBzTxDT6jBGPEZR8JTuamzBGPEYw
dfuw+JJ8jHgY8daEbDfrG2P/A3nuaGU=
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}}, {[{True, True}},
FrameLabel->{{}}
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}}, {{Automatic, Automatic}},
GridLines->{None, None},
GridLinesStyle->Directive[
GrayLevel[0.5, 0.4]],
ImagePadding->All,
Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
(Part[{{Identity, Identity}}, {{Identity, Identity}}, 1, 2][#]& )[
Part[#, 1]], 
(Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Part[{{Identity, Identity}}, {{Identity, Identity}}, 1, 2][#]& )[
Part[#, 1]], 
(Part[{{Identity, Identity}}, {{Identity, Identity}}, 2, 2][#]& )[
Part[#, 2]]}& )}},
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 0}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.02]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhUNkKyJLWFWIQB00tlqtC7UzeH5tq61VBqR2RkxdkMUFY5xq
GcLYVgUErMCA2IJWDSiCYgAHLVTtY9EqECxUMAiCELZACEmlIIvD1OT/e2fe
e/O2e8853/UOPbQlnM/j8cLfHsbr/2VmOltXxDorNnTJ4rBiMsWj+JYluz+L
vTlJiZL9FhVlwnhtUNlF+Gy3U3ktMGfvr2BqIjH1uoZfkZ1lrALcFxpKdAoz
9r0Cye+PfVAXw6sYCjI+KcHGc0dL1qdPc/T/UjyLwoHC1W84G1Wu1N/mNyTb
FxWdGpjgaL9yyGrcl2TJX3MicWrV6MJ7kG/2jcvZPcrR/n9gla1dQoRohLvy
9m9VbgWWvq6XFowbODqvEprIS06hF3XcclNVYyK2eEfStgE2/wHmJxjKQ3x6
uerRhW8nPER+YtPmA+OdTM8jVMnr7ubdec5tMQ18jI1rXwaLMxuZvhosSuOL
It2ruS6ZsWEtpr3y/LvyspheJXLid43e5pfjsFGetB4jTTUrr81Ugup/ggtr
RUtKbzeD0mjArvwzjyYD20H9NGC4sjdyQZQan5sGNGKdUNIi8u0D9dcI2UFp
bNOzARiniVP/xsiNhhVHv9eB+n2KKq8Oq8IUA0x2VE8R+1l+Tve2V6D+m1CQ
EHPZwvtfzDPKETYj3loStLJnFDSPZuz4+ppXZu5rREYY6xnEquiz63dOgOaj
wrytgQH73aZg7OasUCEk/LyNl/INaF4tSJbcuJxwbBpjxnajLVgTHWXW4cIj
NL/nUI5/e0mi5ZE1poatyC8NxUnOjNA8WzFZ9jROlM4nJnnxbbDdWpcj2TOD
0HxfwAEl4p3LzUm90lgvkLb7BytrvgWhebcjwObd00f+siAuJoEdOCyYFVZ8
3pLQ/Duguj105c8IK2JqJ36J1YvTxiuXWhPKoxPKm725hRPWhPrpRJTILTGt
8h1C+XRCanF3fnSiDaH+OrFC75S9L2Amoby6EPBm2w2pky3z2wVdgFXlWZUt
ofy64HLaR9uePYv5V2PJd8e7g7+xI5SnGsne0/bTQnuWhxolY7/Edb+0J5Rv
N+z0FScsAh1YPt3wI8fj7hc5EIqnG4oc+QsPe0eWVw9OOL+S+UocCeXfA6m6
eFG70pHl14N1BmHGJ+8JCN2HXjQ8NPiJkwQsz14UHTsXI+oXELofvViUGpx5
81Mnlm8fgl55Xh2UOxG6L33wXfnxJfW0E8u7D+25p7b+vN2Z0P3ph1dM1BGb
MmeWfz/8hzcf93d2IXSf+mF5XdayTOrCeGjg4Bh5XVPnwnhosKzRcTx0gSvj
oUGbLiRFfsKV8dBgRqYV72q7K+MxgIzf3S4e/MiN8RjAOeuxg1PpbozHAJ5k
9HgE6NwYj0EcymzWh21wJ2kmHoMoLM9WrJG7Mx6DeMBNJD6fdGc8tFBk3h35
8KvZjIcWSt/FsRtuzmY8tAizrT4lsJzDeAzhp/MFDik75jAeQ+jf5ON3684c
xmMICy3NHx6xm8t46FC7KyzpcfhcxkOHe8pV6RfuzWU8dOBHu25qE3gwHsPo
EWzJ37fXg/EYxoGs1gzf+x6MxzC4f1q/0Dp4Mh56zPcb23RmjyfjoceUz5e1
6lJPxkOPH2X1YTXmQsbDAP3JX6sWBwvJf6Go/3g=
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtMk3cUxWsHKBXl1UIpFGuQiZsmDLOonXD+ZioqQXkJOphkBMFZGeJU
BKbDDNh8VJhTo0QmBmRzSgERUedWmBplPAYTBBSklEdpeakIbCIy2u/uJl++
973nnN+dHxUftJ3P4/G2Tx/G8//FXT9XpwhVvp3JaVBvDpeVvOin+1NQOZav
TOzXqW+6pvcH3vwRGenahLjBLnqfjzMV4yvChjvU2eeMdRlZEn9rpbCFvlch
zm7EO3dDvXog0PikBP4J4qptbpX0/3Ukb9IVKOVKtaA5b7e3oBzLv7G9qCi9
Da7fLaS8/Pbt2egaLIrI/GPU4w5mv70WmVDQBK7/71jnO7ba27Id+dN/N+dV
ILZ7ybzgJVpw8yrxwS5Lr6DaHiw11V1YWm9uOButBzf/HgKqLf5aMDKAu6Me
0xPuY7yar0he+Rycngdwzv9shdnECwSZBj6EXuRp71k2QvqqUDSVXfl25yg6
k40N/4QwpPjnBJdx0luD/Zn1bUer/kGCUd7uWkRWR8v997wm/XWQlWgO3XF4
A45GPUKm5kS135gkP/XIEeSoyoKmsM40oAF+uQcXv7uPxzh/DZgRo65q9JnB
jNMiMv9G8NWMeIUFn3F+H+FBaeWx4Wo+M9lpfgTfXfnesZnvMM5/I9yjH4e2
BpixeUY5rk3Yucaw1d/GnHF5NMFvh+v+6lpzFhtjrMeItP46NeKIBePyaUb4
m3MGy49nMmM3oaoZSf9alLe/nsm4vFrgc0+8S1s8i40Z2422wK1NaO623ZJx
+bVi49NQ7S+OAuZjavgEYYH3b6Y/FDAuzye4JphIVSfOZiZ56U/h8TI1Jtrd
inH5tuH0fEV5aoMVq60xVhskC3MbZV/NYVze7TjWauYW6D6XiUwCn8FCdntS
WjuXcfk/Q0Sc8vIPX1ozU7uIDgwktSyItLJhHA8Nbll0GUK22ZAfDQrycmaV
FNkwjo8GSt6bVUk8W/KnQXzDAnlpgC3jeHXCa5KXEJ1rS347UZu7ddXJYVvi
1wksT+5Y5WNH/rXQNPltOXDcjnhq8V3Ysyz5EzvKQ4tLhpC1yoX2xLcLE23v
z923157y6YJ1Tl3oWIU98e7Cta59O0RzhJRXN3gf+cY83iIk/t14tbd4k3e+
kPLrhkY2MWPjsJD2oQeJRSK5lVxEefZg/Lfvz6SkiWg/enB2quXchToR5duL
xYfOX00RO9C+9ML9aJrMPsqB8u7F+Zr3ln5+xYH2RwenIoXt4VcOlL8OvoOy
p596O9I+6RD0MuInfoYj8eiDymztkT11jsSjD8H3b50udhATjz6k9Q7q1dvE
xKMPQgSX5hWIiYceXh2iFZFDYuKhR42iRDf6oRPx0OPYhmWRioNOxMMA8yTH
9Lv3nFiWiYcB68sOZ/OtJMTDAFVG/ZYlwRLi0Y+MT/Si9dkS4tGPRp1kJLxT
Qjz6sXp+lVOchzPxGMD1G1JBerwz8RiApLCNp7rhTDwGECWYlfxq0pl4DKKC
798Uu8aFeAziRMAXUmelC/EYhFdYot69yYV4DCHmwMVll6RS4jGEDCu/8V9j
pMRjCOHH5eUni6XEYxitZZ6FmyekxGMYbheuBHv6uhKP6fcbN1xYe8qVeDzH
4nWL9hRqXdl/62/oTg==
"]]}}, {}, {}, {{}}, {{}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}}, {{True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}}, {{Automatic, Automatic}},
GridLines->{None, None},
GridLinesStyle->Directive[
GrayLevel[0.5, 0.4]],
ImagePadding->All,
Method->{"CoordinatesToolOptions" -> {"DisplayFunction" -> ({
(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
Part[#, 1]], 
(Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
Part[#, 2]]}& ), "CopiedValueFunction" -> ({
(Part[{{Identity, Identity}, {Identity, Identity}}, 1, 2][#]& )[
Part[#, 1]], 
(Part[{{Identity, Identity}, {Identity, Identity}}, 2, 2][#]& )[
Part[#, 2]]}& )}},
Quit
