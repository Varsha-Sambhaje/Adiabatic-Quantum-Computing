(*Module-1*)

(*----------Basic code to Solve Quantum annealing problems in Mathematica----------*)

(*Quantum Annealing in Mathematica:Test the code for n>=3 Qubits, n>=2 Qubits needs some modifications in logic: The code will be used to study various types of feature selections and quantum annealing problems...........*)

\[Sigma]x = {{0, 1}}, {{1, 0}};

\[Sigma]y = {{0, -I}, {I, 0}};

\[Sigma]z = {{1, 0}, {0, -1}};

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
       Apply[KroneckerProduct, ReplacePart[l1, {{i}, {j}} -> \[Sigma]]], 
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
GraphicsBox[{{}, {{}, {}, 
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
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
GraphicsBox[{{}, {{}, {}, 
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
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.2`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 7.721364768762253}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kntMU3cUx6kM5ekLqBYE+qBQMuJrOgVx3x9oxI1ohAkObMam4yECgmaK
+PpjiooYILiJhhmYuGGIRJSHzoXiEAUjBqwKBcQWLM+WMhG0OOps79lJ7r25
r3O+3+/niHbsiYibYWVlFffxMF//L57lbFt/yKUiVJN5HFbjG3K6ymzo/hym
d5bYbEiwrr/leWIk/NYlpKe+F/MCefS+FPkLk0wZMCkuXjDXVSzI0YV8uXhK
wX1fgVUHspVdaycVunDzk0oERRScrN39WsH9XwWlESnOu/UK+/bLaWvta1Ef
15SxdXk/9buN1NUJE+kRPQo/ee7fE7K/UFdh4xhY2Er96xB5uGG0J/OqovTj
3+2X62E9rmza3d0Abt5dRKU4Ofzb+ByfWaoB9+3vRh6q0ICbfw8ZUv2Wjq4B
NEzIPk5oRLDrrvVSvR6cngdIz7jLRK//QYRlYBPyztvm79W8AaevGe+cAqs3
Nb6FJtPc8CGErdZxiZemwOl9hKaYorrElGmkm+WltUCVlfbjPoEV4/Q/Ruec
82G1B3iMo9GKzrABSX3bDMb5acXA76aQZt9P2EbLgDYcuB643y7ThnH+2pDm
1pB9s2kmM0+T5z5BwpF7NZPOtozzq8S61DWTpu12zGKnXYllD4vkE8X2jPP/
FCtXFcit+hyYl1mO5zO8Oz4+c73YiXF5PMOGU6+De7+dzRLizfUcMx7nGqcK
5zAun3ZcabSVRiXOZeZuLhXtyB4vj7Cznse4vDrw1bLug6VF89ikud1EB1ab
jn6zYuV8xuWnwk9jDtaNj+azLywNOxF/5q3h+53OjMuzEwFqLY9vdGYWeSe6
sFrwJG/krAvj8u3GvpLwvTqRK2t5ZK5uOA0588TVrozL+wXsFH9Kfg7lM1eL
wB48V4m/3qziU/49OCpfWhyZtIBZ2slfoka75PaNqQXEQw1Hg2o0OXsh+VEj
Q2cdcEYgID5qVP5mu8alTED+1FCcNoTO/tyNeGmQejM6/liDG/nVQF0QsDVl
izvx08BTGJP9otud/Pfi5dQtPEtcRDx7EdRWkxz9ZhHl0Ys1ecLHqcc8iG8f
PO44DwrsPSmfPvDVCtvYc57Euw8uuplhIZ5elNcrVI1NFz74w4v4v0Ki1nfb
wYVCyu8Vlh0JPtsnF9I+aNE4ni/yLhFSnlr4xw7fCNIKaT+02DEtalnsJ6J8
+9EcYzP8PllE+9KPrs6AvVeuiyjvftwR8ev834hofwbwXccefsEqMeU/gOHg
2GpVppj2aQAJrsb7vDox8RjEU/sPuxx5EuIxiLJTJ69PrZMQj0EMdhQVt2RJ
iMcgZm96uiKrWUI8hpBU+UIudfQmHkPITXUwlW/2Jh5DKK5SjrvlexOPYcTP
Cpu7T+nN8iw8hpF0rHBRtauUeAxDPu8HlWablHiMoEKGD8YLUuIxgksvZ8HU
JSUeIyhwDw0xePgQDx1i2k6Xt8T6EA8dtJXV/IslPsRDh18ORwsi+3yIhx4h
y40hH7x9iYceLbFRHhfifYmHHqKaKn9JmS/xGMWk0H3s1yFf4jGKnTn5GrtP
ZcRjFO1Z2/cnJMuIhwHXhpZ41V6TEQ8DsipTNhpHZcTDgBv+m2MWL/UjHmPg
l7rcjkr3Y/8Bw3HxTw==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQ1HUUxxdLRCGxDTk8gAXk2F2WRcsyk+9PUTwGQhzTYgiFAk1EIKxR
nLAQRFE0RxtESlJQJJQFyVEJLw4BEeSQo+VwQdhlYQ81VvBIYndfb+Y///u9
7/f7ebzwmHURkzgcTsTEoT//X8brJ7f2WBWu7E1IRkrB/HtVLQq6P4Hm7aet
1xx7dOuafcpw0LXTOHctKNlvzwN6n4v+VaW7ohZH3jqVqa98qHx4Q3krGmD8
vhAJR2IUi717oArSPynGsmNOc6b+KIfx/z8RYXJYG1+jxrT2nNgl067i8oWl
KzbtfgZjv+uoWy95rRrXwSPkaLnOvQxLn0nOh+16Qf1vYl+J1yfs8WvkTvzd
nnMbl3ZOlaQ6cphx3h2s3qA+Lqw1YQsMVYGNMT2S1qi3mHF+JcTqYlmd2WRW
oXOfmFCF3CIVN/h3U2bUU409eyd71s83Y+sMA2swKTGwKv72VGbUV4vW5Mgz
v602Z70J+ob3oLb8Z0pWgwUz6r2PU5zNI38FTmdxenmx9Yjb9tInsN6SGfU3
wDs3e6AxZwYz0miEldn2gcRN7zKjn0bYaLIsW+24bJVhQBO4Y6nakWYu+WtC
43fBnDdp7zH9tJCjzfh3Yenxd5ZZkd8WHLBXV/uMWTGDnfYWjC7UpP1xaSb5
fwjB/gd5W8KtmYNejn0rvG16OtOsbSiPVkgrfJ+63rNhWyL11Qb7vYu6/X+w
pXzakT30ocxcbMf03awK23HjlfJ8XJ8d5dWBvEtPt6WdmMWe69vpOsD1f8AN
85tN+f0NUWns2vHR2czH0FCKA3Ntb0bnz6E8pZBrvnlREjyXGeSldGKay6Gy
DnN7yrcLqq1rdL1l9qz+vr66EJYUIn0Y7UB5dyM9LXQ82NKRzTQI7MGz9HPF
lVscKf8eXDcL9fW47cgM7UIeIWKrycdptjziIYPnSI5iMJZHfmRI/r4m1q+W
R3xk0JVnNuXwnMifDI2tGX0mCU7EqxdfpVaFhDU7kd9e7Hwe3l3OdyZ+vfhS
DaXrPmfy34eh3CJteqcz8eyDIEoUMLbAhfLog4XsYkDkYRfi+xizigpEHf0u
lM/E/Zg579Ml84j3Y5ypk6yp+WUe5dWPkivKaj/NPOLfj6vxdy7W+rlSfv04
cpZrGpTtSvswAGQclHaNulKeAxgVJAmi17rRfgzAoa9OZ5LvRvnK8W3ohvlZ
Ju60L3Js1Na0fRTsTnnLoUya0ia97E77o0CBKM/xJ3MPyl+BpaL15wRfe9A+
KXDYTBAkLfMgHoMQBpvOSJ/JJx6DCMosv++7g088BhHP2R335i6feAxic8Fl
RZmDgHgoYdtywnbvLgHxUKLtpal6eZOAeChRU90SMJ0vJB5DqHSrcu5MErKf
DTyGMPh52caCTiHxGALW3W1OfN+TeAwje1yY+Fm6J/EYxrIrKb5eck/iMYwd
C25YWEBEPFQ4dJBXOZwhIh4qDMzZtr7hiYh4qDA+4+aFktVexEMN05XXL2ad
9SIeakx2ifXf/8qLeKhR7byodBFfTDw0sBHKDxRvEBMPDYaT3j7I3ycmHhqs
kJwszpWIiYcWRSvlLxy6xMRDC7/jMV/8auZNPLTI141U233gTTyeoNZz55KT
Yd7sP+3M9og=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.4`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 17.381229248870827`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhYksRUQligQxEtlCAAMJxtYR6fkptlC3KkhxKopYl46MEnUE
xXYKFmzdAPdlOo4acBQGWqyU0KmEzZFWo8gWRFkTdhDFRtAq1uTd3pn33rzt
3nPOd902JkRsnmBhYbH5/WG6/l8889m2bL9jQVhHchrSQutXT86xovtTiHPt
CtMm88rUrukDq9QXEV51NXb+qzca7n02Os/Gq1KWjmounDfVdRxefPGTE/rn
Gu77AgxYB8Vd2dWnGVxlelKIokJEMac2+v8m2ofbfgt3qNbY6VTKYLtiBB9E
QlR8Cbh+Jcha1XQnzasJPjGZFUbJnyifpDwrvGkA178UraPi6iPOQ8h+/7dO
VYbS/ugTs1Qj4OaVo3zH2lBcfol55qrEUYWT0v7Qa3DzqzAmXq07EzOOSqPk
/YTbCJzd9lPzFB7j9NxBcWxUnlXPBBZhHliNLwxbN6QVWzFO31/wiR7Ls0i1
YR3JpoZ/I/orq2NTF9syTu89THSe57KAZ8d2muQptXhiL1jyj3oS4/TfR5G9
TdKN+MmMo1GD3E2KvMaZUxnnpwYlmcfbWuIcWLh5wENcHpzr06zgM87fQwRd
WV9rbTONmabFZNai07nK94eGaYzzW4eMq8b0narpzGxHV4cSx/CDQwmOjPNf
j4/qVsbyg2YwkUmOawMudIfwjVZOjMujAStbLC8WaJ3Y1i2maoRr1IMFEacF
lI8OK3b3K0fWOjNTN8cCHRzDBGsuuc2kvJogtJbmxHfPZC9N7YxNkK//rG5j
rgvl9wguLp3BJ7fPYh+bGzajhB++cKJMSHk2I0m237/muZCZ5aU/RtG/Yzkj
N2ZTvk/QPTFlU/JuV6a9Z6onaDnvkKScJ6K8W3DVMiq0ZUTEZpgFtuLuLu0S
w5o5lH8rPq0taxSq5zBzu5g2vDigyZ4lcCMe7Ui7lNek3+NGftrxKuFh7N56
N+LTju9tgo2PA93JXztus3wv/nF34tWBrzOlkS7D7uS3A5M9Y9WvlnkQvw7w
/PfFF+Z6kP9OjNt/EBdq60k8O6H4dduOws2elEcn2osly95WeBJfPbZvsDwn
meNF+ehxWtPbqPjOi3jrwf9GnejR7EV5GfDO3fDA+KGY+Bswnz8szj8ppvwM
mKC4n7/kmZj2oQuG6dL0B8u8Kc8uqEbP3Vp03Zv2owsWhlNrTlpLKN9uvN10
7vP6OAntSzeOjbxOtCiVUN7d+PGM9g9nFx/anx5URWTrZyf6UP49yLuyssih
1of2qQfeqzP6nkt9iUcvcnUqWcUhX+LRiy8njcoPdPkSj15kKc8eDQzxIx69
eK0q1dX97Ec8+vCyJUy1ZcyPePTBVibaOxA5l3j0YcjjetPGX+YSj34sivQI
0dpJWZaZRz8SU8e9/bZIiUc/hL8Hjn9bLiUeAwi+eyygUuhPPAagf/TGfTzJ
n3gMYE9XzUL/On/iMYhW0bv6SP8A4jGI44J1SxMOBRCPQYiFj16kGAKIxxDO
ZERnVYlkxGMIge7LHRYvlxGPIVRXpIo0yTLi8RRB+U6d86/JiMdTCHYUJF9r
kBGPp1ArUuwElnLiMYx9vG2FqTI58RjGlEpeRu86OfEYxoo9ObeWHpETj2ew
jzgck6eWs/8AsfDlZA==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgdUk3cUxQFHLUdAEZXDEowgAQkJYSWO+6dqcKAIndRYBQ2IoqKtiOPU
USnSiqCUWu0QBayiUgdVXAFcCEIdVMKGCIGQQEQsII5ak7y+c77z7ffuvb/n
ErUuXGZmYmIie3foz/+X8bqnaKtNfrByy26UPf95WIe8le5/QOT108r0ifeL
Cp2StGGFv4GJ7RYdcyyC8X0OUsKi975Ob8LhQ/o6CYdVG7bzYjUwfp+Pf6M8
y1/k9KIrTP/kHCwjgiecHBig/wtwMFLkN8h7A3NFdvx080v4M7ay7Uc3U2bs
dxkxQ7OnifKGMK407Uaf+zV4RhzPXMUZzoz95Qi1OQVOxgiW8+5vRXYx5Ph8
3MKX5sw4rwTVcb37yiIsmNBQN/Et9+j4leetmHH+LXQlVjVsMhvNbva5v5tw
G6o7plkjudbMqKcUQdvjsoIXjGHhhoF3kXhq4EzCOhvSV4bmzLie5rSxTLlF
37AckslJ/SX540hvBQa9XV5/WDGerdfLi6/EEnnq4jtqW9L/F5ITYuxnDrVj
RhoP4Dx+5/c6J3vy8wA7vK6oOgId2BzDgIfYseCMZm64I/l7iOLLrV95rnZi
+mnStEcoGLfz7ZFvJpDfKtzLSrJqWODMDHYUVaj0tuovrnUm/39DcGHh2d0y
FzZBL8fpMXw2e3805ZkL5fEY6hSJ8Nq2iSwmWl/VSG3Z/iJgBIfyUeDql7El
xzI4TN/NJl+BC/POC986TaK8ajDDLjotNG8S69e366sB58EX8kw/V8qvFu2x
n0RWFbuyGYaGdTg4ra/u/RA3yrMOObWeEWKFGzPIS6pHCtehURY1mfJtwDLv
e0Gp3ZNZZYW+GpBsXxF0LtGd8m6ERanFtkdDuGysQWATVqB707N9XMq/CX8c
yDthaefBDO2kzWiM+yzeI9eDeLRgP3/xxNl8T/LTgtu6PcOWXvUkPi243/ji
40TJFPLXgoKj+yT7H04hXkpMvzXfLE/qRX6V6LWPP3Ojw4v4KcHfFR9dv4FH
/p+g1Lpc8s8bHvF8grOBl5ZbpHhTHk8g+VoWM8acT3xbMcpr6vysUD7l0wq2
cvpyr0w+8W7F7+vdCq7U8SmvNixdWxA0x1lA/Nsgba20rJYJKL825HOOOq44
JaB9UCEopmhdb4+A8lShwdV3+E5/H9oPFeqnHlJabfOhfNuRPc/d7EiJD+1L
OxIulsXw3hNS3u2oyXSxlocIaX860FJ72GThASHl34E33fE+TQoh7VMHqocU
Hl/r6Es81MgRrZC+jfIlHmpoSkLD00/4Eg815nuzZBedL/FQY66pdvC80I94
dMJ07+CJWZv9iEcnpItGZlTL/YhHJ8xrrl9cOdSfeGhQ65pr/WquP0s38NBg
491Pc1PT/ImHBuyX1WucH/sTDy3qJYmxF+wCiIcWNRttf5IsCyAeWpxWawZq
cwOIRxd0IRnJa7QBxKMLts+vzDYVBBKPLrhHjhZkJgQSj26Ul4bM4V4LJB7d
uLnE47vrpiLi0Q2z42v7woJFxEOHWcdUe9r3ioiHDnb9UTO3PhIRDx0OvSrh
jLIVE4+nKKxQuecuEROPp1j8weUwcbaYeDzFGC3v1/tqMfHowcGXzFzGm8r+
Awjb5os=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.6`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 27.18122924887083}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg1Q03UYxxGykB3WcKDIeBmiwBDcGCwt5ftTzF48vBOBvMA0O/IlUqDQ
jmFKCATXKWmgxmkcIF5nzLNMSMlN0QQdiCJs40UZbIONwRBDEfOMbU/P3f//
v//b83y/388j2Lo7PtXZyckpdfqwXf+vGfazq1LGk7+ryz6I+Bbd8rlnXOj+
R6TXKNeOr3mpqPfLH15ffwqeRTzvjPanCsf7aigNQ9HXe8cUP52w1S9QLu7y
Ep3W0/dymDkJDasL2hSW9bYn5xGVPN83JaIBjv8v4GOzJdbteC/c1FXpK9zq
cKH85OPZlSY4+v2JvNYjU5WccYSmHL42EdIAjja2J4b/FI7+V3Dv7T+kcP0X
1dN/q6uUyCh4Z8p7pxNzzLuKvAme6uVVZyaxVyP2BOdrst1nMsf86zDt2ZTz
LO411jgRMj3hBlzERXH1B2cxh56bKPsnwtj8G4fF2wc2gW04scpX684c+pqh
O3BpRt7k60yXbWt4CwfS23UaVy5z6FXh6HiCqdrDg2XY5KW3oHVNLt/iPYc5
9Lfi8rJOucWPxxw02rC/Zv+pXoEnc/hpQ6BSV2MN9GLv2QfcRX/rvDfWBc4l
f3fxu+ejLf4B85htWsrhe8gN1D/N53uT33YUKvWHKubOZ3Y76na48YJqi7k+
5P8+6kq+nrXJjc/8bXL8OpCVWP0w3NmX8ugAd9/mt7jPfNm2z2zVCVe5UuJt
9aN81Dg7tTk8Se/PbN14cjWOWGsPfJ4WQHlpcIx7qizBGMCe2NpNaPBm1XWZ
ZbOA8tPC3LOtaalWwGLsDbuw5U7i8Zj4QMqzCyUV54853w5kdnn53ZhVt3J2
UewCyrcHtWvSTR2XF7AWla16ENwScsMiCaK8e1EQVve47WwQ87QLfICXM5sO
5i9YSPk/wBzV6bXu5QuZvV3KQzzyv6Xf6bGIePTBKtJcqihaRH76kGHmcs47
BROfPqQVfWg6uTeY/PVhhVARu2M0mHjpMC5sM/BSQ8ivDneUabLK7hDip0Of
OvPonPhQ8t+Pbim/ckdTKPHsR7s5yKUmRkh59OO74CvPmy8Iie8Avpi9oaFT
GEb5DCAxLGhrS0UY8R7AQtNfk3KvxZSXHieuxclyvl9M/PX4oFDbIXUOp/z0
+EqebHywN5z2wYCh5Jv7skbCKU8Dlmc+2fJ8awTthwEaVd77GZoIytcI4bdH
jeq4JbQvRmT1FvcuaVxCeRsxeDPFsJwvov0ZBOfSVEB3kojyH8Ru99LBzB9E
tE+DSCjfJXtFJSIeQyjPSPApeVVMPIaQmdjp7LVSTDyGkMSR7SyTiYnH0HTe
Y4Xci2LiYcKZ4WWHCsfExMMEk6mzdlIYSTxMGDlj4aWmRhIPM/jF262qnyNZ
iZ2HGbntgm2irkjiYUZaw2TVYZ6EeAxDsP3xbfM6CfEYRuvARpdVRRLiMYy/
k26llTZKiIcFiQUfhRpeSIiHBbmtVRvFS6OIhwWf1PbMy86MIh4j6HpuLlD8
GkU8RuAzknVlxmAU8RhB/vATLRNEE49RBGxvHs9JjiYeo9CeWya4WBpNPEax
QnMkx3InmnhYIbkf5hPgJiUeVnyq+NJj/Wop8bDi3K4Xe/Z/IyUeY+DEF6ec
rZey/wCqVPbg
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglM03cUx7HTbJ6IIIcgh04t0I2j5SrH91fEAXEYQGEeOEGHKENFxBEQ
iE7UxDlFRSZO3RxHhk5QUIsMh0zRqeA4JgXUUigthRYKBjS4gLPt20v++d/v
fb/fz3PYvCsynmNkZBT//tCd/y/D9XDtPrOy4O6MHMSqt15NiO+i+zwcK827
9bVwK6psD6kjqi5iJNZH/K5WCsP7IqzfXhAj5Q/hXIGuSmFjlbSlYtVr+r4M
EwvHXbavmYAmQvfkOtJ2GLEsSw4z/H8DIyFi997UaWyGpDDZf4YYPwVOXm6o
+4gZ+t3Gk9PXloVPmcUcY078OcatQc/BuSEhHsbM0P8PGOccFPzmYsKK3v8t
KbyL7enhwRl35zHDvDoEH7zwtCLMjPH1dQ9Xzre9U0jm0/z7SLCofJW7yYLd
G+O+n1APL2mYylZhSXoeIqgqrW1g2wIWqR/4Fzw3OjqK1Nak7xEcPFvOxSUt
ZN0ZuoaP4Xt7zdQDGlvS24DJxw8PKzzs2W6dvORGXJyW2e42Zk/6n6Kud/dn
ghsOzECjCf2lcVPGUxaRnybMmTVqetRtMQvRD2hGCe/oA612MflrRnt/6Ihr
2cdMNy3mRAteNXpjddIS8tuKI/uGvtrgtJTp7Uha0ScwCwtTLSX//yA46qQf
r2QZs9PJsX0GZeyo9u0WLuXxDBXRcZW/OziyhK26akOq0OSb5C5HykeCoJl+
nQsuODFdN7MyCd5+ZxNYs96Z8mrHTNcfO6Mseey1rt1YO7RqkyLVMx7l1wHv
PTnpe09/wgL0DTtxLT/ffCL8U8qzE3Jups3+OS5ML+/Qc5ivm2f+7okL5fsC
KYVF/neCXFljg65e4HzB5Lr0QlfK+yVk8bzZAo4bm68XKEVDcWL2YKwb5S9F
5dXRv4tr3Zi+XUwXPEXVdhtt3YmHDMUB40dMs9zJjwzj8lq7R8/diY8Mld9u
+DdLyCd/MiQ23OK5F/CJVzc2RV++p3jDJ7/dSNoveno2WkD8unHq+7G1K28K
yH8PwpNlOydMPYhnD2qOqazKUzwojx7szVZujm32IL5yCOrrI0xcPSkfOSrM
Lw3WHfck3nI8fN3plzLoSXn1onEhf8Wiz72Ify+i7Dtmtlz2ovx6wa3OPHNg
ujftgwK5xjFK123elKcCUudwjuyBN+2HAqXjEs3xJT6UrxLVKodS/xwf2hcl
yneYCTU9PpS3Eq8ORP1yTiSk/elDaFBEb8jPQsq/Dy/bwHkzKaR96sPRU06T
RRt9iYcKZsXDHZE1vsRDBT5v5Q9G1n7EQ4Xyl01eZel+xEOFztDxmg3tfsSj
H9b1dtzpXv7Eox/LOYczbp3xJx79OJQQfXPLqD/xGEDmkpKOuasDWK6exwB+
tUhT3bkeQDwGsDb1dk/iXBAPNbIjwh5Z7ALxUMPkg6kX7jeCeKhx6Up4ubc5
Ix4aWK3f+eGJFYx4aMCd7faFIpURDw2ER/ov+hYx4jEIS7FcerKVEY9BOFsm
WKg4IuIxCOMZ4pAAdxHxGMJye9M9eXEi4jGE6kRx/kCuiHgMYUT8opLdFREP
LfY3Vz3O14qIhxY2ldnPNbaBxEOLmDtfKgNXBRKPYVyXH1afzQpk/wGb0+5H

"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.8`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 36.981229248870825`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgs0lHkYxik51KZY2tSSUDNqjGHcBu3zP9lOupyOpNpaq42hs0vbsO1F
97ZxTmdKOKuS7SZRVE6caLKYUhqELarBVm4jlzGErtu0rZl59z3n+77z3d73
eZ7fOzdqW1jMBBMTk5jxQ3/9v0wNZ4ubO20Ll3bukKKoNcSltnAC3WfgmGnm
u+2K9wq5Y7Jmtfw0RMtqc5VvRxXG9+exPqXILbSqR5F1Ql/5yF+5fExRVkff
F0IVE7jh4NsGDK7WPylCIqvsKtGqYfz/Gqp2bdfFdbzAZFWOZNHk6+Bslen+
sngDY78bUK9pidz/jw5uEalVr7jlCFfOjv8YacqM/SsRmhIc6OBuxs6P/63K
uQlp5bbNPw+YM+O8WzicnDvw3QlLJjTUbUw+zh0bCfyEGeffgf0D+5OlzVbs
9ivu+IRq1B9r3W9VNJ0Z9SiRULExy+2lNQszDKzBoaKLzy4LP2VGfbUwydNd
nySxZZ079A3rUJS2tCU9344Z9dbDWV5Zn90xgyXo5UkaME9dab/Tbibpb4R0
t/2bDSH2zEjjPg6VWikSk2aRn/toWLHHcezibBZiGPAA4c7D1hNVn5O/B1gW
fFtePdGR6adFpDbhFN/BfKPHHPLbDPcZ1QVanhMz2FE1g7PK4+XVW07k/yH6
b3RdMF03l83Ry3F8hPXx0dH9/XMpj0eoK/9VIdntzLbE6usxRoo1kjPTXSgf
FSzajop35bgwfTfbQhU2Ne1RTfV1pbxaMJWTF7ehxpW91rd71QIhTycSb5xH
+bUi/+7rTq/BeewLQ8M21Jg7HK/dPZ/ybEN82SMebxqHGeQl/43mQblT5FkO
5fsEQR5ZOVGeXNZQr68nyOEnZQdVcSnvpzjpuLdmKMyN2RkEPoN3qVvqzm43
yv8ZpKcVop4fFzBDu4h2WDjNb3c3W0g8OuAfz9u6PmMh+elA+sdvLke78ohP
B8SuccHh13jkrwNHUrzLFnzpTrw6kRiUvKuv2Z38duLqFumTI9F84tcJTmOv
7ZwxPvnvQm6ddnnWbx7Eswtf3bWpMbEUUB5duCedNFoQKiC+3dhcXTFrZaaA
8umGddrEYnW7gHh3Y0ZkutNPHE/KSw2Z65maDz94En81eE3fmu8t9aT81KiS
JQneffCkfejB2TvizPglXpRnD2xWPpa1Hvai/ejBHRubUPbQi/J9Dn5mEC97
tpD25TkE2ec36aKElPdzRHWZc8IKhLQ/vTDj1/95bkRI+fdi7MLi8CF/b9qn
Xgi2N7n47PMmHn2Ysoz/9S9Kb+LRB8nifdwSKx/i0YcjK36/PrTWh3j0Yd33
7zxdT/kQj340eXko1qp9iEc/zM6USw8s9CUe/ZBI1+RdSfQlHgOwa9MseXjD
l6UZeAwgVvf0wBsTP+IxgEb793s+C/EjHhpYaHPCvVP9iIcG5/qeOq167Ec8
NAgejf1X7OBPPAbhe09pniT2Jx6DyENMmOySP/EYBMcn7W3mqD/x0CLtSsa0
XJGIeGgx0qS8XLhPRDy0kDgeVZcoRcRjCLcKZHVlVgHEYwjzq8QJFWsDiMcQ
VpfM1FScDCAew7C0bA8p7w4gHsMonrToD/mCQOIxDO61qaPFCYHE4wWmhMki
LskD2X9PkO+T
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
       AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQ1HUUxxEVZBMvVhkxkMNVqI1TDmXZ74+FBdZSLkE0lEkyjBgEbUYB
C0skpVAqHFwX0wZSpEAwFDATRREBkcPWxcVwl4Xl2OVQU2fUydzl9Wb+87/f
+36/n+ewdUfUNlMTE5Ntbw7D+f+aup5syOJWhqozc2Ba/XimdrOC7gsRc7ri
eVxsB+rsDugi635CcVGvR1bxKKbel2LtUetE3y3PcFxqqLPQDPqaSVa8pu8r
sW3wbtkZ5XSmjzQ8qcYmrxLrMItZbOr/GpyS/enT5jSbcRQlaQGcWoR+H7JV
uWgem+pXj6ih5oZn3AXMJf5I41Pny2hP7ZGnW3LZVP8rOMmta4mYtoiVvvlb
UXIVCvyxMvWRNc27hkwfpTS7bzHzMtZ1FHX842TTsoTm3wD3q8PutdW27PpT
5zcTmlAg33Wu7thS0tOM2EIm9WmwZ1HGgbfw+3XzpkPVDqSvBcmp3PHUEkem
zjQ0bEWR1SubF4VOpPc2vp3hGumfu4ylG+SltUM7YfoZdvNI/x1YX955bd72
5WyKRid6s1ura+NWkJ9OvIz51cxf4szCjAO68OJc4ppTq1zIXxd+3PGdud7l
HWaYFn+kG+XofOlg8y75vYsnxXqIOXxmtKO4C8+KdUlxL/jk/y8IOTzp5tH3
2FKDHDs52s6vL41VulIecvAyUuaKW91Y0ieGuocwt/oFCze7Uz4KxBy22Tdd
7s4M3biVCuQ7RoeNf+BBefWg2b67W37Dgz0ztHvaA2FO455LAk/K7z6asiO+
Ka7xZEJjQyVCOQcD9/K9KE8lvo5QqTeWejGjvAO9sL6xX+b99krK9wFiwj+S
zS1cydpvG+oBtDXzLIbf8qa8/0aLyw/mDfu92UKjwD4IyiYuFL70pvz7cC2g
ye3TXT7M2C7+IaTSmnyBzod4qFAT2zw4J9GX/Kgg6W18X6X0JT4qjDildVVF
+ZE/FdqS7ufua/UjXmpwXLMzw0WryK8a5bzcettLq4ifGk33+8J1HqvJfz82
8LeH1Z1dTTz7ERpvWZHj4E959CPcX5UbIfUnvhp8aefet2S+gPLR4MK/Jy5p
DwqItwbT+K5O1SYBlNcAPIY32GXtCSD+A5jtbnU+eDKA8huAbHJAY7ldSPsw
iOgHI3X3Hgopz0HEWKq8T24A7ccgLGRuCUkdoHy18HbMKC5yYrQvWpzmfZGY
lcAoby0WrE0xS5Ax2p8hiJc9yhUpGOU/hI9bjmt5VoG0T0PIc3m43CI8kHgM
o3+mrUSfF0g8hsGaT6zruBlIPIYh4xX5nTcVEY9hCMv4s44KRcRjBL2m6Vd2
Z4qIxwiCopPiN10UEY8R3BHN0Agei4jHKGrmuEUvdQ1iBUYeo9iZ0VgxLTmI
eIyi9Uz+E80vQcRDh2Tdet5NdRDx0OG3vd1BZbbBxEOHFF5teN7GYOKhh1tf
uSTlaDDx0KNty1rPdV3BxEOP/hIfjrulmHiMgZ/8vGu+REw8xhB4anHekxwx
8RiDXae9l/yqmHiMo2pfXvvFV2LiMY6Mgjsbj/mFEI9x7D+U0JPxeQjxmEDR
RLnkw6oQ4jGBWwcXVwr0IcRjAqt/nmtm5xxKPCZhvWYs8nViKPsPqwrh6g==

"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "1.`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 46.78122924887082}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\)}

GraphicsGrid[{{\!\(\*
GraphicsBox[{{}, {{}, {}, 
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
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtQjXkYxk9329hNN5trmFkyTa1DJa08bItipFqGxmG1Lg2rkokotnU5
i6YmG6t1qzilNQjjllz+XzptagtJV6c6p3O/O1QSw3bO/7/vzPd9893e93me
3zv555S4TfY8Hm/T8GG9/l92tvMILtOrfLEs4xBy6jJOvYxxYfcnUGx8Yzkm
cuYqJgr1sRWFmBu6YPeMQSf2vgRjfWvyyqOduNOnrHUJQ8qTpc5ljuz7cvjI
9iwMsnPkDLHWJzdQGLE1cNZaB/b/Lbi1Xtr8udKec20TbQ93vYv6+JDlf461
Z/3uIV93+o4l046bLsh73O/3APFunnr3Hh7r/wh5PkP/vo3gcSXDf7eJOOhW
udz2X/mZ0HlViE/aedz0zScyy1bVcDt4aNyMoY+EzhcjIKjopaXpA6nu9xue
UIO1Wed4s68OEaqnFg2X+bcdc96TONvAJzhROeNtQvIgofrqMIoXc3/VindE
lmFtWA/ufKOvJXyAUL0NSPGMifvOv5+kWuVtb0SJeszVheP7CNX/FE2ivXcm
u78llMZz9CnmT5B/8YZQP88xa8nhK0XOFhJpG9CEHPPu4BSX18xfE2RTPna0
+piIdZog7wVWBi5N/8nfwPw2o3JN34cbi3TEZqetGd8HpCe6/qJh/l/ioP2j
ORf/UhFfq5yJLTBs4A5XP1OwPFowc8tFfdEoOUncbK1W1LzP/zF9jYzl04Y7
hXVbnBJ6iLWbV3kbxs2dH2/0lbC82tFd4cl37GsnA9Z2/e24Vq5JOtXVwvLr
QPZUzT295AWZZ2vYiTEfI5ctUzxjeQ7fb4gtqwqrJzZ5wlfoWGcJDX4oZvlK
sIA7+Xt/6yPS2GAtCQoubOpKyr7F8u5CHnezXjitlHjbBHZjkcNoly5+Msu/
G53CiSHi+8WwtRP0QBER03tUcAOUhxQrqr+eUjiyEtSPFNzsbMeI0CpQPlKc
/btAnLiwBtSfFEnNO7d6fFsHykuGffkjQ84PNID6lSG4ZnzW/ejnoPxkeJho
/0emxwtQ/72wGI47LKhtBuXZi2L3q173klpA2/WirF3wKse5DZSvHBe8w1av
z28HzUeObdLbVU3unaB45FjEX+20Zzg3mpcCqRv37VIbJaD8FUgwRonKlnaD
5qfApM53keFFPaD7oMS62uLr+WlS0DyVEE66OVnvJQPdDyUe3LXX9VyXgear
woF5yp7SyF7QfVEhKuFJd7KkFzRvFRr5It2WbXLQ/VHj18xrK6oH5Sx/NY4E
5klL9itA90mNf9KDtMEuSsZDg9y+kcIzR5WMhwYHxB6XBkeoGA8NvN/sX5Is
VDEeGvg8dc4a91nFeGixKmTpmQm71IyHFsGLW5cf06kZDy3Gh84V5wo0jIcO
WZ2i64ENGhyz8dChsfbKzuw5WsZDh6LImpl3S7SMhx7Xztm9qvtSx3jocc7x
SGRrmo7x0MP1cmm0uUPHeBhw9qvz9VPD9YyHAaP91t3MLdQzHgaE/KCZFvRJ
z3gYUUD69k4XGBgPI9xcUqP2VhgYDyO2RwfsCPM0Mh4mFO+4FZS6zch4mJCe
tZ4fIDYyHia4RVemZYw1MR5m+FXwf9uYYmI8zIjRhsUNPDYxHmb4pEmSwrzN
jMdrzDdUyMI2m/Ef76YwDQ==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQk1cUhVmUJaDASKcKMyDV1l2Ljop19DztaERcJsSKrYzVDIuGFlxa
t1LFVqq4lkWsLdIKqQNaGbQoEZEfxAURlRhDEAJkAxISErRNSluWkrzXO/PP
v997zvlumCglOt7NxcUlfuRwnP8vet3HfRVYwtccPArXjEeHPky2sPscTK6y
v67uMHPSkHSTQJoP2bl0t283mNh7CcYcOsDLbjByP15wVDFmhy4z9vMN7PsS
pD1foah42MWZBY4n1yG23pvetqqT/V8GlXZ4oeiZjuMpC3cu4ZVDm/q25+oY
Let3Gwu/PvZulk7NTYs9e882tRKy9CK5uKud9a/CyTiSEuOr4iQjfysLq9HC
j0mXR7xi82ogLsguqElq4uY5qxbBFu3ggETO5t+HsfKT4U36Rq7WNnVkwgMo
IiPPLRM2MD2PEPPx+PxZxodctHNgHXpFH3CGjBqm7zH4n+Vu+nWHlNMcdDSs
h0eYy6R2yVWmtwEZJpXGHneC2+WQt/MpUiuGvf0u54Hqf4a71pdhLZnXQWk0
YvNB/9Bkyx1QP40QisbkSifWYpVzgAyqnivCiyvrQP3JUOC6v/Qd8VM4psWe
fYGbyQm+UStloH7lMC32ExUr5XDaGTl3JQTt7djRBOr/JS5ezVwgGWxGqENO
iAI/aDriS7NbQfNQoLSrTlg5sx2JCY5qwrMpOcVVcWrQfJTov5eYKripgaNb
YIkSa9bVB4m9dKB5NWNr/S/hG7foYXe0szUj27b3/djyTtD8XsH/WrT19rhu
LHU2bMFQ4I3Uwj0G0Dxb0LokZdJypRFOeemtqBxyv/X7UhNovirMuOQzjXfF
jKcNjlLhRcz8b6ImWEDzbkP4+fzg3aeseMspsB3n1k//6daRPpZ/O57njTX6
4zWc7WI7MNN+qn3+0GvGQ43yFEW8H/eG+VFDaAkJPHrkD8ZHjfMFVxT7VvzJ
/Knh+50ut4FnY7w0EJbK9CKZjfnV4LhfWlrQBTvjp0Hnyy3aF9v+Yv61uBQx
eDVpRj/jqUVRbY7nA1s/y0OLgbX3bVz134yvDik7D59Yeuoflo8OeTmSoMmb
/mW8dXg4x7D5o8kDLC89tjQ275f2DTD+eqj55UHz7w6y/PRYu65GXJYxxPah
E6d7Rb5hG4dZnp2Y21B3fF6AC6H70Ym8/Dm2RftcCM23C7Oi108NbnMhdF+6
YOCHjn++3JXQvLvgsT9x5boiV0L3pxuRexbdyBzjRmj+3eCFT6vL2u1G6D51
I6t6gSBK6UYoDwOubSRbby52J5SHASVjpUue/OxOKA8DCo9HCI65jyKUhwGX
937upU8YRSgPI7a35Fa0Ph5FKA8juPiyY9tnjiaUhxFnjCFzj54ZTSiPHpzf
PdA3pW80+d7JoweqDacXrRF4EMqjB4qAJ6mG6x6E8jAhYmBh05sAT0J5mGA/
EFz06S5PQnmYULetZWJYoyehPMz4oswkeW+2F6E8zKjYdnlC7EkvQnmYsS+5
2a+424tQHr2o/y1bbl3uTSiPXtTuuL3H96I349GLlLhIk9zmzXhYoJonvjNu
LY/xsCD8VdSXWQU8xsOCpDRBW4Cdx3hY4V8tWC3k+zAeVpgDjiSF5PowHlYk
Jh/OjNT6MB59qLJLV9fM8CX/AdTMH4A=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.2`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 3.521364768762253}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtMU2cYhksLtZarRqdLNlScigQnTplzRt/fKHEqGGFephYiCsVtMjrn
BS9LXIY4FEVF3ZwiTIjRGZ2ABopupwgMJ1Kpclux2httaQtUoRQRcLT/vy85
5+Tcvu993+ebsjUtLpnP4/GSRw739f/y8pxFigPjbi7X7c9A6JxWUYpOyO7P
oNThlPLqfBTlwYdtseWXUPAgLjb1njd7X4TggozsmDsCxa/n3XUNpSXPti+5
w2ff38S2U0dibXe9FPZY95NimPzIK8cDHvv/NjIbbw27yt5y4pZC2SJxGU5P
5Sfs7hjiaD859rUNmrhpg9xMSc59Z+g9iDZKxct3DHC0/1/QVK/P/fFeP1c0
8ndLoQJ2pTTKNd7F0XmVWHD1miQ83cnN9VQV6uveGbqt7+Ho/GokXf/+Qti6
V1yVM3RkQg0SP5idMuuxg+mphSYneWdtQRcX5xn4ALLmkIk9t2xM3z8I21n/
BkoLp9vvbvgQpWdNFe/2tzO9j3AgOj2jd5aB+9YtT1aPCbbChOg0LdOvRP+G
GN6KvW0cpdGA7b9tTPef1sz8NECkXDvkY27gPvMMUOF0hTz+bnot86dC1Zuy
U6GnKjj3NEnOE8zLlUZd155nfp9CUVd+rLL8d3jstDyFMDn/wwt7FaD+G5EQ
cs5/35SHmOSWE9wEw8HHhzKaR/p68mgCNybV+ufxZqRI3dWMJbeFsoGP20Dz
acHwhMAtNU9ewN1t3M0WbM6b/1YVoQfNqxX52QqLeKcRfe52zlYsyyzZtfKG
CTS/f/FeXviYOr0Fiz0N1fh89pri6rE20DzVQO58645FnfDIO9yGzrSQANfW
btB8n+HJXK+pDysdqH/krmdomvJ3o/zxS9C8NZgzpy8sse0VxnsEPkdL9MXK
JmMPaP7PMcOxKTrc3gtPO8kLnNzGL0h/6QTloYVz9+W79b19zI8WftnBqkV9
LlA+WtgqVsSonf3MnxayopwXpT2vQXnpUKya9YO6e4D51UGbwDsYb30Dyk+H
pVqDbJVhkPnXo2z3jfgS9RAoTz3mCuWX8xuGWR561Bz2E02ueQvK14C+AaFu
koRHaD4GfP2J8VOfJh6heAw4cyD+UlSMF6F5GXHUPm6tqdqLUP5GJF3rTbQv
5BOanxHDeSuzJSV8QvehHfnhuesXhgoIzbMdtVqL4HiegND9aEe22u/ndWO9
Cc3XhPK0ftUvmd6E7osJQfuN8zYNeBOatwlVt24MXkz1IXR/zNjTvmdyktaH
0PzNeD8/cXVJnJDQfTJDIA8zZ1ULCeVhgeD8V/32yFGE8rBAOX16pObKKEJ5
WDBvV5R98wQRoTwsyBlKGvzyiIhQHh3QqFZVj3KJCOXRAePEiYGR0tGE8ujA
2WXmAUfjaEJ5WBEt3ti6YKmYnPTwsEIVUcXzLxYTysOK5zMiC6TBvoTysOF+
WHrakmO+jIcN/pLhL865fBkPG+TKmXzJNj/Gw47ryIrIVPoxHnbkvNbwfBb4
Mx52/FRRbKy67M94dCJ3Yankqm8A49EJU9IfqSe/C2A8OpFlUuvXqAMYjy7c
X9ycV4pAxqMLERuOZu0pDGQ8urB6bHDIFX4Q49GNrrWa0bKPghiPbhwqSonr
3hLEeHRjs/mbWq8TQYyHA+usgpln5EHkPyV3H7U=
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg9Q03UYxudg/BkwhhJaiURghyEh4amQ8HxLUDvDxPLsjCTl1FPxkDRR
OAGnk8gMD/+FJ0ogcJqiMyIQZOo0iRwxEDZhwMbGBm5j6sRZpMX2/fbe/e73
/32f5/m8weszVm7gcjicDROH4/x/0etH0hz/miXa7APYuGz+tp/vmtn9MXxZ
N34t+fSItH6G2JRcfwYlr2VxPHMN7P05BMsKro+n66SnShx1HpKZostks4Z9
X4PI9/Mkd4Z6pOZkxxMJzvSm7pBN6WL/1yJfItP9Vt8m5Ssrtsfxf8W0uK29
T3JlrF8DjqgzEGC7Ip2VUnRrLKwJUcs7TKU3ikH7N6P//Lsz1lgacG7ib2XF
DfB2FmqD0u6CzruJS3nXl0f92Y5oZ8nw0ZKxZdL7XaDzb0MRY+u+c6sHsrGw
iQl3cPhyjmt3yQConrsISxd5fhAxiJXOgS1wOTdLflWkB9X3O9ThA/YCuQHa
bEfDVjwoyF/vJRgB1XsPmbxtInGCCZkOedvlWHQw0L4w08L0t+GweNG8Xces
oDTacTEq5UKM8hHz0w5RqOZpw8BjLHUOUGDvlpYDPbonzJ8CtSP3z9frbHBM
SynqgFFw6ZvtA0+Z306U16XOCVSOwWlH2Yks/Se1qtZnzP99aGKfL226ZkeQ
Q86MLtzKrFynqXrO8ujCscjVbSlFf2HTRkd1o6N7imrpzr9ZPkoczfb/vHHV
OBzd/GuUMDTb+1qj/2F5qbC7XlR5WPACzxztxlQI0lXfftXwguX3AN6bmsLE
114i3tmwB9+1Wn4cOvQvy7MHsfsDdoZFcIhTnrgXjZUZy0eLOYTmq0Zqojr0
YzuHyO85So3xuT/UrFozidC8++ASfKDYu2kSecUpsB828cW03OlcQvPvR5ki
qfSXHC5xtksZwMbOb880PeASykODXfG8qtJ5LoT60YDbtU++vtiFUD4ahJBG
XoDFhVB/Grz0zJ98c7Eroby0iNZHFWw760qoXy30p+tkoXZXQvlp4c41PLYm
8Qj1P4j673NcVBU8QnkOYvqbRVeGnvMIbTcI4+zju4KT3Ajlq0PVVf3Jo2Vu
LB8dot668XqCzY1QPDpYc/vaYhPdWV56ZPRtTsg74U4ofz3Sq8TKaUZ3lp8e
MwPnXxLM9yB0H4bg9l5n3aaDHizPIXDl/bVvd3kQuh9DqG5Pj10b4snyNSCh
kSvxyPQkdF8m7s2+X8xp9mR5G/AHvNJ6+XxC98cIhWqF3Hc1n+VvRNJcXXFn
OZ/QfTKigRN/dvYon/EYRkl3Vu8bMV6MxzD0fdyYuv1ejMcw4oO+bn4i92I8
hrG4QXtCPdWb8RjBlois0Ox13ozHCJrDV0ztuODNeIygrOlKocnmzXg8xOLA
FSFtC33IESePh1A/ViWKxD6Mx0MkNOQn+rX5MB4mpCoKJeIAAeNhQpjiusS0
VsB4mPBOSVZecrWA8TAj6vZZbotVwHiYEb6v5XTaAl/Gw4yTR8sj4vb5Mh4W
RO85vmZHqy/jYYH/Tx0+AXwh42HB7vY9ccsihYzHKMSlC0KefSpkPEZRsGhl
UUyOkPEYRad+sjGuXMh4WPHVyazwyFYh42HFZ2WzTy2wCRkPK+yHtpr3Bvox
Ho9QXFOkmvKhH/kPa40Gwg==
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.4`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 4.321364768762253}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\)}, {\!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1knlQU2cUxQlLCAQMIbgVxEERrXUbqQuZlPNBWxD+cKuM2oJTrUGnOggD
2BEqtgw4zqCgRahYKlawFKHUrRVBeRFwo0RBsGGNBAgEwmaEyFKgJt/XO/Pe
m7fde875XY99R7bLLS0sLOTvDtP1/+KZzwJFvEtxoCYuCRd97gXtieGz+/M4
9XjfTMqktaLEPVm/reQSopZnzHY7b8Xe56EyWsb7VWapuJhlqgKsip7aWvfG
gn1fjJTD8hV3I2e4/m2mJzcwb8u1nX7OUxz9/zaW3vXoW1Q9wdmrciM/sr+D
OT862MZljnG0313UP8+9lR5r5N4PTasYXXYPJWG7T/fKRzjavxx5WdkBC8IN
XN67v1W5CiwoDN89EjPM5j1AkWPYaP6iAc7bXJWY43FdEB/Ty+ZXIfLg/pOZ
bVqucnTZuwkPkaeW7fgitIPpeQx/dWpCrrqN224e+ARd5aLDwctVTN9T+Gco
5vll1XKaOFPDagjqQ1XXtlcxvTVItV6yw7eqmIsyyYtU4sFfsr+fq66A6n+G
Hwrb/G5qFaA0apEabznurakB9VOL0oTv4nJCGrDJPKAO3VrPkFWLm0H91eGp
F+k/0vwKpmmhaS8Q/Cb/a2l5B6jfesjsvspwTdbCbEdVD48y522b/XWg/hsg
8lVskBn7sNAkx/0lUst4CQG/DIDm8RJTpW5hnR7DOBBuqn/QcmKsxvbIa9B8
VKh+L7NyfYkBpm4uxSo8Ol6Qkz71BjSvRkhEqauAURhN7UYbMTn9uX7ncSNo
fk0Y113eOHHnLXzNDZsR/iQogwyNgebZjONeNsqPF0/ALC+5BfbBP3u67pgE
zbcVfzacTGz9/l8oa0zViuF4Q3p20RRo3m0olF0ejGiYxmyzQDVcDxStDh+f
Yfmr8UFylbD6jAUxtwt9BcO93+qSFvAI5dGODztrz+kLeYT6aYfi2OllfB9L
Qvm0o98rZHVrlSWh/tohlg1bx26xIpSXBvNLshpUjVaE+tWgqWdTtnCvNaH8
NCiLcF87V2dNqP8OlEqFQbwIG0J5dmDRLS6m1mBDaLsOrJMaVqQc5RPKtxOj
ncLaDeN8QvPpROBM0NXWOFtC8XRCWb5PfWLSltC8urBHV3phabyAUP5dUE76
qJvGBYTm1wXP20+vZn9jR+g+aBG49Ex97IgdoXlqUTa7ZeOhSHtC90OL+TL+
+US9PaH5dmNNwVTaA7mQ0H3pRkvgs8qVr4SE5t2NtgHHmzU7HQjdnx54yzc8
yX3uwPLvwR9202ElAY6E7lMP5t6/dFt035Hx0MFvibexeO0sxkOHnL3CdT/l
z2I8dHhdJD7U5CpiPHQ4d0N6LCJNxHj0wjcsMX01z4nx6IXmdyvxLqkT49GL
L0Oqp9dGOzEefbi6a31WRaETOWvm0Ye6634iuy4nxqMPhoTDiQI3MeOhR8KW
m/zKz8SMhx67vE5dkaaIGQ89to50bI6qEDMe/ai3Dnh7cFzMePRjrKvx2OI1
zoxHP4zSpPwr4c6MxwCU0y6f6LKdGY8BZCTK+cYXzozHAJLvpFxQCiSMxyCS
3B0uRPpKGI9BxKqcMtujJYzHIKI+zZu7sEDCeAxh8OFgxkq1hPEYgu/00QRH
iQvjMYT84P3eikAXxmMYjzQVwsBvXch/tukJVA==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtUk3UYxidaCog4NhiMiZZ60vRoZueUeHn+pAkmmOIF7MBSiIuFstIQ
hbwFGil57Xi85CkhBLSBRxRRY0NBYOIEcRsq4jbkto2BLhBBNLfv7T3nO9/9
fZ/n+b3vRCWGxjjxeLyYN4f9/H9x192KFKE80LA1DeGnCnNcJprp/giW5KZt
Wx/Xqrjsl25edvkUtLr5ythKA73Pxi6J5UVbxgPF8WP2yoO6f/GsfEUdfS9H
ZWn8wkDpdYVlmf3JebQ6bVwC/S/0fxE+9++7t2nyP3DRZcnmuhSjM0mrurBe
Da5fCVZELG/f1q3F5Ij913smXcO6VYNrZVuawPUvRd+/je9aO4zIfvO3LkuJ
gcLtF3t3tIKbV4bk0B1KvZMJMx11AxPC7zwrlHWCm18OmX/jrYKYbtzomfRm
QgV+l0dWmJ2fkZ5KTCiZ5ZHylw2hjoFVyCgW2io+6SF91RgoctNMK++FYau9
oQqap2F8QVAf6a1Bf0jptVvlL/CdXZ7sNh4tydso9x8g/WpM3FMtGJX/EhyN
Wphrxq6eLXhFfmpR/2dyeHrSawQ5BtRB/zTXJPmGxzh/dZh7aeqVUK8hzD4t
Yv9deK/htzkrhzDObz2UxjMX4uKcmMOOrh4FMWNVB0YOZZz/e/DNm/NeXsFQ
NtYux0+DwI/0vHtLhzEuDw3ios35M7uHsbhYe2kRKzuR2pH5FuPy0SE6pavY
8/23mb2bUK6D7fkHm3TlbzMurwaUdn62L1g6nPXa2/U0wC9EdPBg73DG5Xcf
VUPdlGWZI9g8R8MH2Pw4YbB3vDPj8nyAmckL7y4qcWYOeekPMW3lrvV1IS6M
y7cRfI/TZ7IMLux2jb0a4RKtGKHZ5Mq4vB9hdXVkbfLwkczTIbAJh8LPvfrj
2EjG5d8ES6xeKp3ixhztIh4j0lt6SHnNjXE89NiZkSjRhowiP3qESldtONc0
inF89ChZmmn7ItGd/Onh0Szyq3ntzjheBqRUR5V4Bo0mvwbMs+apZh8YTfwM
mHF0en5Ew2jyb0TYV3u/TB3HJ55GxAhWbj8Vz6c8jMh9XhR1s5BPfJtRZVP1
9/bxKZ9mqC6YrB8GeBDvZsgXnD2ZmuFBeT1BgGu1RFvnQfyf4PjH9wMDxALK
7wnKXdVXFFEC2ocWtG1Z/DD0rIDybEHAzVsPB2wC2o8WiHlHFFfnCCnfVvwQ
nCw9nC6kfWlFwgRZ/E9qIeXdisN3++P3iTxpf9oQl71CJV/jSfm3Ye+GHxtM
eZ60T21InMovZzZP4tGOgrO27MI5XsSjHS2ra1P9d3sRj3bM35A4q+mOF/Fo
x2CV4dIJHxHx6EDSYFrZ99Ei4tGBeKFBFP23iHh0YH7O2oSE5yLiYYL+pTn2
1wBvdsDBw4SKnuU5VXu9iYcJrV+fLxNrvYmHGeNLji5KG+dDPMzIOT2jeti3
PsTDDFnSp5UnLvoQDwtMQ9wrg3hi4mFB6ubpi1wXi4mHBSfdg2uMv4mJRycW
WCpuq/Vi4tEJzdUQlXqKL/HohPvP6+YZk3yJhxW5ruacEdd9iYcVmdqwnQFu
EuJhRVHSguB94RLi0YXIqqw97VkS4tEF+TbXiyu7JMSjCxqjNazefwzx6Mb4
4tqeNbvHsP8AS1cHGQ==
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.6`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 7.381229248870826}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg1Qk3Ucx8fkxUHoBs+gpEioTsB8OzHtDL9/s4RyXiKaEtyhFioHCvh2
KiBcpoinDpQMQckTEkVbSgSI6NClTHCBgkwQYRtsIC8DBYZAQG7/f7+753nu
efv9vt/v5+exKXp1OJ/H44W/OczX/8vKcp5cHsfJ/LX7f8KGCmV8vsCW3afj
Zuz0uj99rctL3A91BZZk47BmuO7nBD57n4t5t7bvHNTwyjPPmOsylL5zg44U
j8vp9zKc+uFkY0b9qLw70PzkOgwPVfmbpg3L6f+FKBop3lh9wCS3V+fE+NkX
Q5lyyy2YNyCn/W7gK9VtlSnrpdw7VHp30KsMY/djBVy2kfW/DZNVTP652Z3y
3Dd/q3PK0fBrwvyzej2bdwdD7nZ5OXKtfL6lFNjidSpz285nbP7fELy6cH7b
hVq5YtDrzYR7EMrOqtdOVTI9Ffi09aCPJL5AvtoyUAl/0+WQa+pcUH0PICnJ
9VgZpIB2v7lhJWxqIh2eS2tA9T5E2EXblR+vUiPWLC9GhYl5v8n03c9B9f+D
SZVHDUJfHSiNGuyNC+S1bdSD+qnBO3u8N1QmdSDAMuARWnSXjknTukD9PcLa
01EyY7oR5mmh0scovR4RUHmvD9RvLVJjd1wVLn4Fix11La4YEsbT8/tB/ddh
9pqPFGe5QbxvluP+BG/zbZ0i9plA83iC5QHHHRY2DGHLZnPVo8lb40d8h0Hz
UeMwt2G09OgIzN04mRrLHEy5jU2joHk9RYiS11ntMwaTud3gU/BOJE4od42z
/BqwdeJp9cvSCSyxNGzEDN8FTqUHeITm2YiIa71ti92tiEXeoWcY3CteLyqz
IjTfJtSrow5K1vOJ6qG5mhBwcOaV8Vd8QvN+jl3hQSsWHZtExBaBzUheMOEh
+NCa0PybcfzbqkVJpdbE0i60BZdSE3uufGNDKA8NLp52HTnTakOoHw0kRZ/H
r9tjSygfDaRR2//ts7Mj1J8GEbcid+3IsCOUlxZ3bcrWts+YTKhfLcIjpdLv
iiYTyk8L0ev49PplAkL96/CuYiQsvEZAKE8d+hdlxQlD7VkeOqiSna21BntC
+bbCL2KKrCXGgeXTCg+dqpIbcSAUTyu638uWHPnxLZZXG7i5blVL7R0J5d8G
jf+K5IA0R5ZfG+4WnojOcZ1C6D7ocXxs8wfB56awPPWYu8BlKNpzKqH7ocde
metQz8WpLF8DCtJv7I4SCwndFwNSKnbPclsnZHkb0O9jU3EoQ0jo/rRjbJX1
A1mDkOXfjvrsqnMZ00SE7lM7rhJF+pchIsajA7KcvNbiLBHj0YF9fHney2ci
xqMDgpO89gE3J8ajA6GNs9IUIU6MxwuU329OCstyYjxewM2UFF/Z6MR4vEDV
468X2k5zZjw6sdXkuZoLdiapFh6d8AqJDzX94sx4dOL7hhlPfq93Zjy6sCXK
beZnYo7x6ELzX3NuXwjiGI8ulM+va9GlcYxHNwozSmvHqjnGoxthBftqBhzF
jEc3hMu9eA9WiBmPHqzx70uKSxEzHj2YpcwMd6wQMx49CJ7gihKtXRgPI+7c
NEbXLnVhPIz4IvPxYUGiC+NhxIjkdZtnmQvj0YvkpWmR04ddGI9ecL1z+vmf
uDIevTgv8pFU7XBlPPowzJmW7P/DlfwH3hwIHg==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhQOI4IIosoUlJBAVVBaHymK15zfVurRBQLQgWNBWU6UgaFVE
C2Jh7JRRxK1D7TCC2ipYVFzqGqBqXQbQRJPYAEogQBLWaUVxHKQm7/bOvHn7
veec74rWbopdZ83j8da9P8zn/4u7HqjZ6Vy1SJedD/uhrYveFBjp/jAi1wQL
P2HtNVcFBd0xV0vxsswq4uhqLb0/ia5H59LTJA01P5eY6wwUFQ+V4fMr6Psq
dJTFr7K+L0dPjPnJBfibcDozVgnu/0uImBvnrSxtxljNiYx5Y/9AxOg8XlZU
O7h+12Al9E2RBBgQkFT056D/TWhfFqry23rA9ZdjeZNEIpMP4OT7vzUnalEd
Ll8xb/gfmleHuqx9zzVBgwi11G10OGRPnRX3mubfgbV087TA9De4Pej/fsJd
rB6efvmX796SnnvobUy8wr4fRqxl4H3ETD3r1pA7Qvoe4DCrKxu8xmO6bHPD
h7jkJLQRnLdinN56nOPvTz9bbs0yzfIyGhByfkqyqdiGcfobscvtYqQgZxTj
aDzGD5qHMYUyW8b5eYw9IbmqzKjRbLFlgAJa27qzE0LtGOdPgSXL1SuPudgz
87SkIiUqBasapK/sGef3CQKyZJnxqjHMYkfzBA6FM21N1WMZ5/8prB2l1YFF
45iPWY5AhXu1YV6fbhzPuDxUuNB0d9f6BQ5Mtt5cauyvinhw3HsC4/LRQL05
qcXt1QRm7uZcpYFh5ZuKzgZHxuX1DP2GlssHYyeyV+Z2g8+wM9HR9PLWRMbl
9zdet0w84RcwiX1kaahF7tr0vGmHJ1GeWlSqBIH2I5OYRV5BEyo3RH5Tv8GJ
8m1GxZz4iOynTqyh3lzNmF17QOmKyZR3C06FDql+PTOZuVgEPkfQ7kuLg5yd
Kf/ncBtdXnk+x5lZ2iW9QOzgMftZRmfi0YoZ17NGLix3IT+tGLe0UTFb7kJ8
WiFJOZQm93clf62wCr64TXrIlXjp8FZ/epl+2JX86qBSh4vyv3YjfjoEIk4X
+MSN/LdhyY5lCW3z3IlnGzLTE2aUn3anPNrQWi3TpE3mE992GIuKaxbm8Cmf
dije9c2YbuQT73aMktlFecZ5UF566BubPubXeBB/PRZ6/+Qhnu5J+emxwK7v
8YdHPGkfOmCTEJ6QwvOiPDuwLXXBseJUL9qPDvAf8bY9UntRvp0o3/mili/x
pn3pRNBMT2nG796UdyfmvG7sVLoLaH+6sO+37hXz8wWUfxfyXjgk3+gX0D51
4UZpxB2W6EM8DNB9sX2h4i8f4mFARvTl6DC+kHgYUHKtdemYZULiYYAtX7FJ
ly8kHkY4MetS+XUh8TDiW5+Ri8cHhMTDiPrCPQf2ThURDxOChg86bUkSsQMW
HibEd5XwvzwoIh4mrAnbXvD5fRHx6Mb81OKw6Hci4tGN3XfbHKQf+BKPbpRq
bXqlG32JRw/CBkKvxB73JR49iPwq+bNEtS/x6EHyOvEh2Xg/4tGLtDB1+naJ
H/Hoxang3Ds/ZvkRj16UVOp2lFX5EY8+HLHL2npD70c8+uCg3VL6zENMPPoQ
knOzeShaTDz6scex191rr5h49EPXEzNXcktMPPqh1HSEpP4rJh4DGKqL0h0N
mML+Az5E5cA=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.8`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 10.581229248870827`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\), \!\(\*
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhYMRCQgIITEYFIsVC0pdIl1cz0+lBbUdCmKHKjqtFh0Hq9ji
4E7L5rhQESimiE4VnApWJIIKVgigDgoS2ZRV1hASCCQFwzKilOT/e2fee/O2
e88533XZsT8gZAqHwwmZPIzX/8vMdOYVHxVk+3QciYFozCVj1RVzdp+Mo6nH
csYOcovznWP7/PMvY7Wi9SPfCDP2PgOqQdfv/1wzIU/9w1iZ4KyfGvqB4o2c
fp+NRfGWN5Qpo3Ktv/GJDBbZjS2uCQY5/T8PieGb757OHZRb1aeHrbG6h4Er
q1WJIzrWrwDzratHDU965e7B50oNbg/g5ZOlb6vtZv2L4GMhs53N75BnTP5d
n16MzyWeJ5+PN7B5JSjNPSs+61YlX26qh2iyHc16daKIzX+EDLW0JOhMDB4a
3CYnPEZe+uYIQ0AJqJ4yZH04yCsNr0KAaeATlB+PypdZNYDqe4oL7/1SvGF/
GzqOGBuWo+DwA73dpi5Qvc+wQ+/S90mjCgeM8sIqkWa554tzG3tB9SuQuanm
7Z4b/aA0qmA/dLP03TE9qJ8qtLgp/5HMHoSvaUA1LBu2B9TmDoH6q0aj/4Qi
yNsA47TgczVw7Q/cJ1EMg/qtxXwxWV/mPwqTnfpa8J0uOkU8HwP1XweprZ8u
1ecN5hrlOL/ArZ36nNj74yyPF7BKGJLFuL3D7l3Geol9dXec2s9PsHzq8avF
cNaPRznE2E2QXY9Hu71kTYvMCM2rARve95UGNpmRYWM7QwO2bKkZGo6bQmh+
jZg+7bOU5mVcstbUsAlbdTNyxM1cQvNsQs/SYm151FRikhfbjHHd7/Fm7uaE
5tuCqbe9omoqzUnlM2O1gCx4EPn1gWmE5v0KbncihPEOFkRoEtgKnrdNiTTP
gtD8W+GpjeZEB/KIqV1wGyyu++3cMsQjlEc7pANzUz3OWzI/7agrTI6dttiK
UD7tKI8cWPf6qRXz146k8Io1ZiHTCeXVgcJKruTTienMbwfS/PJIutSaUH4d
kCyIM2xcZsP8d+LyGS//VU9tCOXZiUz1vZrD39myPDrhu8I+csaILaF8u/D6
75MabvwMlk8Xklq/mrPE3o5QPF3Qzg25mLjVjuWlRFDmrqu3rtkRyl8JriR/
xSmdHctPibLoSMxbYU/oPnRjb9vC0ONR9izPbvTerPS7VGFP6H5040ZK4PYY
AZ/lq8LypqQez218QvdFhW1pEmfZNT7LW4UX0oOOb/r5hO5PD/79ra9a8LED
y78Hmri3HpwTDmyfehCm/eF10WMHxkONkMBTLQE2AsZDjdVXy2IKAwWMhxoN
LhmHOGkCxkONHGvN0jldAsZDg6yXtwXChULGQ4PxmNAjmjAh46EBd+Xd6xfv
CRmPXmy1/rLe452QJJh49GKVe8WSy94zGY9eOEsOK/SnZzIefRiRjZW6Vs9k
PPpgfmjMbq1IxHj0wenC85yV20SMhxa8AkHC7HQR46HFt3uXJSrVIsZDi1Cu
MiFpsSPj0Q9eUlHwwnBHxqMf3cmXCjMLHBmPyfd/iX7mTzgyHgPY3KzyDPGe
xXgMoG5iY+7VU7MYjwF8s27vhQrFLMZDh+KEzOh2BzHjocP9EYNDZ5CY8dCh
Nmxea9UlMeOhh4ev+083O8XkP0zMAFE=
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
         AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg1Q03UYx0FKFAUVGW/jbYBsQ8mTkV508v1NU/A1JQkTUOFEEk8BPU0h
UrskonkgUamnGIpk8tYSDihUcBiI4LFAZxOQgYwxxgQDNfElt//Tc/e///vz
fL/fzyOISwqPn2RhYRH/5jCd/y/ueuRamkNZqCb1K9xRthdv0g/QfR4UO7++
H+3Xc63a4+jQ+up8XJCcDtgRr6T3hZBJJZq46xevnTppql8Q8/lyiWXYDXDf
lyE99nFBlVoFw3rTEzmk4b+FGFf1gvu/Au/InEPle3WwUZ1PXmxThU5vn9aV
iUZw/WoQLh0xJLs+hjg6+/q4qBYOP3c8nWwYo/5XkaAdSzQUP0Xhm79V5+vw
p+239f1Rz2lePeo2ZQoPv3wBibkUqI255MvPfU3zG8A7nq44YGvJFOOiNxNu
gB9r3JGzehLj9DQiXx2UrsqwYuHmgU1oXKsQN9a+xTh9NyFxrHRqM77NNKmm
hs1YM9rlG+VuzTi9LSjBaN3NsCksxSQvuRXJuft5+1OmMk7/bZyId3W9/KMN
42i0YSL5WHHLH9MY56cNv7fky627p7Mw8wAlEtWln5x7Zcs4f0o0n10X3+E+
g5mmRWf/hWDZeMOA+0zG+W1HYgkzbimZycx2VO0o1cd5bw+exTj/Hfii/tOD
05pmMU+THI870AR0DW6OsKc87mBX0Mwju3rtWcJ2U91FiEv4iSVJsykfFeY/
S5usm5jNTN0cylRIbNq5ZnumA+V1D+nSCEWDA489MbUbv4fYbWNVUwt4lN/f
OLuw5fvgAEcWYm6oRsya8tdRNY6Upxo1RWp9yjInZpZ39D6KhCOVh5ROlG8n
Cn71WJwZ48xaW0zViYi9iqFjg86UdxfmObfcyt3nwnhmgd3419HtuxOWrpR/
N7aMLyj/6ZgrM7eLfoC6VUs/LHHhE48e7Llisav2Ap/89CCmKKBDucCN+PTA
OShQPHzFjfz1YJ1Vs9FupTvx0mDUfzhr0V138quBrMd2VUKcB/HTwD6rPiHf
6EH+e3HYc9C6K9WTePbiI/vKguNWXpRHL8ryFbnVK7yIbx/qy/dE9mV7UT59
iNxXODTrrhfx7kNTUFboB24Cyush/qkO25oaJyD+D5Fz8Kio4qKA8nsIFD7I
GTUKaB/6cVkRmRn4rjfl2Q953MUX+9O8aT/6kVJ1XH213pvy1aKG72pnM8WH
9kWLreqGLzeu9aG8tThUWOdzKc+H9mcAnh12na/UPpT/ALIPFv7wscCX9mkA
vBkZIfIEX+Khgzan+6pdmS/x0EF2xt4qacyXeOjgFND2RBk8h3jo4LCJn7Ho
yBziMYjWoo1nzjbOIR6DcOFHiG3s/IjHIDZmrXb6bIMf8dDj9nLjBu0pP5Zj
5qGHPn1rW6TGj3joIUlZcuCWUEg8hiCrLn1fultIPIawoGKhTU2FkHgMYfo6
Q3vghJB4GNAXFPlNmVREPAyIOnyAPy9TRDwM8HJpSiu+LSIew5DJA0/P44mJ
xzC0mU5J5VFi4jGM954tG5CcExMPI5apnj+t0YmJhxGV17vypPP9iYcR21a/
kjfv8ycejxDorlwaUetPPB5BYBcAzaS5xOMRLELCTu5eMZd4jGCzqDL0ZfZc
9h+btO2Q
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "1.`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 13.781229248870826`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}]\)}}]

\!\(\*
GraphicsBox[{{}, {{InsetBox[
GraphicsBox[{{}, {{}, {}, 
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
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 0}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.02]}},
Ticks->{Automatic, Automatic}], {193.5, -118.66252583997979`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}], InsetBox[
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtQjXkYxk9329hNN5trmFkyTa1DJa08bItipFqGxmG1Lg2rkokotnU5
i6YmG6t1qzilNQjjllz+XzptagtJV6c6p3O/O1QSw3bO/7/vzPd9893e93me
3zv555S4TfY8Hm/T8GG9/l92tvMILtOrfLEs4xBy6jJOvYxxYfcnUGx8Yzkm
cuYqJgr1sRWFmBu6YPeMQSf2vgRjfWvyyqOduNOnrHUJQ8qTpc5ljuz7cvjI
9iwMsnPkDLHWJzdQGLE1cNZaB/b/Lbi1Xtr8udKec20TbQ93vYv6+JDlf461
Z/3uIV93+o4l046bLsh73O/3APFunnr3Hh7r/wh5PkP/vo3gcSXDf7eJOOhW
udz2X/mZ0HlViE/aedz0zScyy1bVcDt4aNyMoY+EzhcjIKjopaXpA6nu9xue
UIO1Wed4s68OEaqnFg2X+bcdc96TONvAJzhROeNtQvIgofrqMIoXc3/VindE
lmFtWA/ufKOvJXyAUL0NSPGMifvOv5+kWuVtb0SJeszVheP7CNX/FE2ivXcm
u78llMZz9CnmT5B/8YZQP88xa8nhK0XOFhJpG9CEHPPu4BSX18xfE2RTPna0
+piIdZog7wVWBi5N/8nfwPw2o3JN34cbi3TEZqetGd8HpCe6/qJh/l/ioP2j
ORf/UhFfq5yJLTBs4A5XP1OwPFowc8tFfdEoOUncbK1W1LzP/zF9jYzl04Y7
hXVbnBJ6iLWbV3kbxs2dH2/0lbC82tFd4cl37GsnA9Z2/e24Vq5JOtXVwvLr
QPZUzT295AWZZ2vYiTEfI5ctUzxjeQ7fb4gtqwqrJzZ5wlfoWGcJDX4oZvlK
sIA7+Xt/6yPS2GAtCQoubOpKyr7F8u5CHnezXjitlHjbBHZjkcNoly5+Msu/
G53CiSHi+8WwtRP0QBER03tUcAOUhxQrqr+eUjiyEtSPFNzsbMeI0CpQPlKc
/btAnLiwBtSfFEnNO7d6fFsHykuGffkjQ84PNID6lSG4ZnzW/ejnoPxkeJho
/0emxwtQ/72wGI47LKhtBuXZi2L3q173klpA2/WirF3wKse5DZSvHBe8w1av
z28HzUeObdLbVU3unaB45FjEX+20Zzg3mpcCqRv37VIbJaD8FUgwRonKlnaD
5qfApM53keFFPaD7oMS62uLr+WlS0DyVEE66OVnvJQPdDyUe3LXX9VyXgear
woF5yp7SyF7QfVEhKuFJd7KkFzRvFRr5It2WbXLQ/VHj18xrK6oH5Sx/NY4E
5klL9itA90mNf9KDtMEuSsZDg9y+kcIzR5WMhwYHxB6XBkeoGA8NvN/sX5Is
VDEeGvg8dc4a91nFeGixKmTpmQm71IyHFsGLW5cf06kZDy3Gh84V5wo0jIcO
WZ2i64ENGhyz8dChsfbKzuw5WsZDh6LImpl3S7SMhx7Xztm9qvtSx3jocc7x
SGRrmo7x0MP1cmm0uUPHeBhw9qvz9VPD9YyHAaP91t3MLdQzHgaE/KCZFvRJ
z3gYUUD69k4XGBgPI9xcUqP2VhgYDyO2RwfsCPM0Mh4mFO+4FZS6zch4mJCe
tZ4fIDYyHia4RVemZYw1MR5m+FXwf9uYYmI8zIjRhsUNPDYxHmb4pEmSwrzN
jMdrzDdUyMI2m/Ef76YwDQ==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQk1cUhVmUJaDASKcKMyDV1l2Ljop19DztaERcJsSKrYzVDIuGFlxa
t1LFVqq4lkWsLdIKqQNaGbQoEZEfxAURlRhDEAJkAxISErRNSluWkrzXO/PP
v997zvlumCglOt7NxcUlfuRwnP8vet3HfRVYwtccPArXjEeHPky2sPscTK6y
v67uMHPSkHSTQJoP2bl0t283mNh7CcYcOsDLbjByP15wVDFmhy4z9vMN7PsS
pD1foah42MWZBY4n1yG23pvetqqT/V8GlXZ4oeiZjuMpC3cu4ZVDm/q25+oY
Let3Gwu/PvZulk7NTYs9e882tRKy9CK5uKud9a/CyTiSEuOr4iQjfysLq9HC
j0mXR7xi82ogLsguqElq4uY5qxbBFu3ggETO5t+HsfKT4U36Rq7WNnVkwgMo
IiPPLRM2MD2PEPPx+PxZxodctHNgHXpFH3CGjBqm7zH4n+Vu+nWHlNMcdDSs
h0eYy6R2yVWmtwEZJpXGHneC2+WQt/MpUiuGvf0u54Hqf4a71pdhLZnXQWk0
YvNB/9Bkyx1QP40QisbkSifWYpVzgAyqnivCiyvrQP3JUOC6v/Qd8VM4psWe
fYGbyQm+UStloH7lMC32ExUr5XDaGTl3JQTt7djRBOr/JS5ezVwgGWxGqENO
iAI/aDriS7NbQfNQoLSrTlg5sx2JCY5qwrMpOcVVcWrQfJTov5eYKripgaNb
YIkSa9bVB4m9dKB5NWNr/S/hG7foYXe0szUj27b3/djyTtD8XsH/WrT19rhu
LHU2bMFQ4I3Uwj0G0Dxb0LokZdJypRFOeemtqBxyv/X7UhNovirMuOQzjXfF
jKcNjlLhRcz8b6ImWEDzbkP4+fzg3aeseMspsB3n1k//6daRPpZ/O57njTX6
4zWc7WI7MNN+qn3+0GvGQ43yFEW8H/eG+VFDaAkJPHrkD8ZHjfMFVxT7VvzJ
/Knh+50ut4FnY7w0EJbK9CKZjfnV4LhfWlrQBTvjp0Hnyy3aF9v+Yv61uBQx
eDVpRj/jqUVRbY7nA1s/y0OLgbX3bVz134yvDik7D59Yeuoflo8OeTmSoMmb
/mW8dXg4x7D5o8kDLC89tjQ275f2DTD+eqj55UHz7w6y/PRYu65GXJYxxPah
E6d7Rb5hG4dZnp2Y21B3fF6AC6H70Ym8/Dm2RftcCM23C7Oi108NbnMhdF+6
YOCHjn++3JXQvLvgsT9x5boiV0L3pxuRexbdyBzjRmj+3eCFT6vL2u1G6D51
I6t6gSBK6UYoDwOubSRbby52J5SHASVjpUue/OxOKA8DCo9HCI65jyKUhwGX
937upU8YRSgPI7a35Fa0Ph5FKA8juPiyY9tnjiaUhxFnjCFzj54ZTSiPHpzf
PdA3pW80+d7JoweqDacXrRF4EMqjB4qAJ6mG6x6E8jAhYmBh05sAT0J5mGA/
EFz06S5PQnmYULetZWJYoyehPMz4oswkeW+2F6E8zKjYdnlC7EkvQnmYsS+5
2a+424tQHr2o/y1bbl3uTSiPXtTuuL3H96I349GLlLhIk9zmzXhYoJonvjNu
LY/xsCD8VdSXWQU8xsOCpDRBW4Cdx3hY4V8tWC3k+zAeVpgDjiSF5PowHlYk
Jh/OjNT6MB59qLJLV9fM8CX/AdTMH4A=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.2`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 3.521364768762253}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}], {580.5, -118.66252583997979`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}], InsetBox[
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtMU2cYhksLtZarRqdLNlScigQnTplzRt/fKHEqGGFephYiCsVtMjrn
BS9LXIY4FEVF3ZwiTIjRGZ2ABopupwgMJ1Kpclux2httaQtUoRQRcLT/vy85
5+Tcvu993+ebsjUtLpnP4/GSRw739f/y8pxFigPjbi7X7c9A6JxWUYpOyO7P
oNThlPLqfBTlwYdtseWXUPAgLjb1njd7X4TggozsmDsCxa/n3XUNpSXPti+5
w2ff38S2U0dibXe9FPZY95NimPzIK8cDHvv/NjIbbw27yt5y4pZC2SJxGU5P
5Sfs7hjiaD859rUNmrhpg9xMSc59Z+g9iDZKxct3DHC0/1/QVK/P/fFeP1c0
8ndLoQJ2pTTKNd7F0XmVWHD1miQ83cnN9VQV6uveGbqt7+Ho/GokXf/+Qti6
V1yVM3RkQg0SP5idMuuxg+mphSYneWdtQRcX5xn4ALLmkIk9t2xM3z8I21n/
BkoLp9vvbvgQpWdNFe/2tzO9j3AgOj2jd5aB+9YtT1aPCbbChOg0LdOvRP+G
GN6KvW0cpdGA7b9tTPef1sz8NECkXDvkY27gPvMMUOF0hTz+bnot86dC1Zuy
U6GnKjj3NEnOE8zLlUZd155nfp9CUVd+rLL8d3jstDyFMDn/wwt7FaD+G5EQ
cs5/35SHmOSWE9wEw8HHhzKaR/p68mgCNybV+ufxZqRI3dWMJbeFsoGP20Dz
acHwhMAtNU9ewN1t3M0WbM6b/1YVoQfNqxX52QqLeKcRfe52zlYsyyzZtfKG
CTS/f/FeXviYOr0Fiz0N1fh89pri6rE20DzVQO58645FnfDIO9yGzrSQANfW
btB8n+HJXK+pDysdqH/krmdomvJ3o/zxS9C8NZgzpy8sse0VxnsEPkdL9MXK
JmMPaP7PMcOxKTrc3gtPO8kLnNzGL0h/6QTloYVz9+W79b19zI8WftnBqkV9
LlA+WtgqVsSonf3MnxayopwXpT2vQXnpUKya9YO6e4D51UGbwDsYb30Dyk+H
pVqDbJVhkPnXo2z3jfgS9RAoTz3mCuWX8xuGWR561Bz2E02ueQvK14C+AaFu
koRHaD4GfP2J8VOfJh6heAw4cyD+UlSMF6F5GXHUPm6tqdqLUP5GJF3rTbQv
5BOanxHDeSuzJSV8QvehHfnhuesXhgoIzbMdtVqL4HiegND9aEe22u/ndWO9
Cc3XhPK0ftUvmd6E7osJQfuN8zYNeBOatwlVt24MXkz1IXR/zNjTvmdyktaH
0PzNeD8/cXVJnJDQfTJDIA8zZ1ULCeVhgeD8V/32yFGE8rBAOX16pObKKEJ5
WDBvV5R98wQRoTwsyBlKGvzyiIhQHh3QqFZVj3KJCOXRAePEiYGR0tGE8ujA
2WXmAUfjaEJ5WBEt3ti6YKmYnPTwsEIVUcXzLxYTysOK5zMiC6TBvoTysOF+
WHrakmO+jIcN/pLhL865fBkPG+TKmXzJNj/Gw47ryIrIVPoxHnbkvNbwfBb4
Mx52/FRRbKy67M94dCJ3Yankqm8A49EJU9IfqSe/C2A8OpFlUuvXqAMYjy7c
X9ycV4pAxqMLERuOZu0pDGQ8urB6bHDIFX4Q49GNrrWa0bKPghiPbhwqSonr
3hLEeHRjs/mbWq8TQYyHA+usgpln5EHkPyV3H7U=
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg9Q03UYxudg/BkwhhJaiURghyEh4amQ8HxLUDvDxPLsjCTl1FPxkDRR
OAGnk8gMD/+FJ0ogcJqiMyIQZOo0iRwxEDZhwMbGBm5j6sRZpMX2/fbe/e73
/32f5/m8weszVm7gcjicDROH4/x/0etH0hz/miXa7APYuGz+tp/vmtn9MXxZ
N34t+fSItH6G2JRcfwYlr2VxPHMN7P05BMsKro+n66SnShx1HpKZostks4Z9
X4PI9/Mkd4Z6pOZkxxMJzvSm7pBN6WL/1yJfItP9Vt8m5Ssrtsfxf8W0uK29
T3JlrF8DjqgzEGC7Ip2VUnRrLKwJUcs7TKU3ikH7N6P//Lsz1lgacG7ib2XF
DfB2FmqD0u6CzruJS3nXl0f92Y5oZ8nw0ZKxZdL7XaDzb0MRY+u+c6sHsrGw
iQl3cPhyjmt3yQConrsISxd5fhAxiJXOgS1wOTdLflWkB9X3O9ThA/YCuQHa
bEfDVjwoyF/vJRgB1XsPmbxtInGCCZkOedvlWHQw0L4w08L0t+GweNG8Xces
oDTacTEq5UKM8hHz0w5RqOZpw8BjLHUOUGDvlpYDPbonzJ8CtSP3z9frbHBM
SynqgFFw6ZvtA0+Z306U16XOCVSOwWlH2Yks/Se1qtZnzP99aGKfL226ZkeQ
Q86MLtzKrFynqXrO8ujCscjVbSlFf2HTRkd1o6N7imrpzr9ZPkoczfb/vHHV
OBzd/GuUMDTb+1qj/2F5qbC7XlR5WPACzxztxlQI0lXfftXwguX3AN6bmsLE
114i3tmwB9+1Wn4cOvQvy7MHsfsDdoZFcIhTnrgXjZUZy0eLOYTmq0Zqojr0
YzuHyO85So3xuT/UrFozidC8++ASfKDYu2kSecUpsB828cW03OlcQvPvR5ki
qfSXHC5xtksZwMbOb880PeASykODXfG8qtJ5LoT60YDbtU++vtiFUD4ahJBG
XoDFhVB/Grz0zJ98c7Eroby0iNZHFWw760qoXy30p+tkoXZXQvlp4c41PLYm
8Qj1P4j673NcVBU8QnkOYvqbRVeGnvMIbTcI4+zju4KT3Ajlq0PVVf3Jo2Vu
LB8dot668XqCzY1QPDpYc/vaYhPdWV56ZPRtTsg74U4ofz3Sq8TKaUZ3lp8e
MwPnXxLM9yB0H4bg9l5n3aaDHizPIXDl/bVvd3kQuh9DqG5Pj10b4snyNSCh
kSvxyPQkdF8m7s2+X8xp9mR5G/AHvNJ6+XxC98cIhWqF3Hc1n+VvRNJcXXFn
OZ/QfTKigRN/dvYon/EYRkl3Vu8bMV6MxzD0fdyYuv1ejMcw4oO+bn4i92I8
hrG4QXtCPdWb8RjBlois0Ox13ozHCJrDV0ztuODNeIygrOlKocnmzXg8xOLA
FSFtC33IESePh1A/ViWKxD6Mx0MkNOQn+rX5MB4mpCoKJeIAAeNhQpjiusS0
VsB4mPBOSVZecrWA8TAj6vZZbotVwHiYEb6v5XTaAl/Gw4yTR8sj4vb5Mh4W
RO85vmZHqy/jYYH/Tx0+AXwh42HB7vY9ccsihYzHKMSlC0KefSpkPEZRsGhl
UUyOkPEYRad+sjGuXMh4WPHVyazwyFYh42HFZ2WzTy2wCRkPK+yHtpr3Bvox
Ho9QXFOkmvKhH/kPa40Gwg==
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.4`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 4.321364768762253}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}], {967.5, -118.66252583997979`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}]}, {InsetBox[
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1knlQU2cUxQlLCAQMIbgVxEERrXUbqQuZlPNBWxD+cKuM2oJTrUGnOggD
2BEqtgw4zqCgRahYKlawFKHUrRVBeRFwo0RBsGGNBAgEwmaEyFKgJt/XO/Pe
m7fde875XY99R7bLLS0sLOTvDtP1/+KZzwJFvEtxoCYuCRd97gXtieGz+/M4
9XjfTMqktaLEPVm/reQSopZnzHY7b8Xe56EyWsb7VWapuJhlqgKsip7aWvfG
gn1fjJTD8hV3I2e4/m2mJzcwb8u1nX7OUxz9/zaW3vXoW1Q9wdmrciM/sr+D
OT862MZljnG0313UP8+9lR5r5N4PTasYXXYPJWG7T/fKRzjavxx5WdkBC8IN
XN67v1W5CiwoDN89EjPM5j1AkWPYaP6iAc7bXJWY43FdEB/Ty+ZXIfLg/pOZ
bVqucnTZuwkPkaeW7fgitIPpeQx/dWpCrrqN224e+ARd5aLDwctVTN9T+Gco
5vll1XKaOFPDagjqQ1XXtlcxvTVItV6yw7eqmIsyyYtU4sFfsr+fq66A6n+G
Hwrb/G5qFaA0apEabznurakB9VOL0oTv4nJCGrDJPKAO3VrPkFWLm0H91eGp
F+k/0vwKpmmhaS8Q/Cb/a2l5B6jfesjsvspwTdbCbEdVD48y522b/XWg/hsg
8lVskBn7sNAkx/0lUst4CQG/DIDm8RJTpW5hnR7DOBBuqn/QcmKsxvbIa9B8
VKh+L7NyfYkBpm4uxSo8Ol6Qkz71BjSvRkhEqauAURhN7UYbMTn9uX7ncSNo
fk0Y113eOHHnLXzNDZsR/iQogwyNgebZjONeNsqPF0/ALC+5BfbBP3u67pgE
zbcVfzacTGz9/l8oa0zViuF4Q3p20RRo3m0olF0ejGiYxmyzQDVcDxStDh+f
Yfmr8UFylbD6jAUxtwt9BcO93+qSFvAI5dGODztrz+kLeYT6aYfi2OllfB9L
Qvm0o98rZHVrlSWh/tohlg1bx26xIpSXBvNLshpUjVaE+tWgqWdTtnCvNaH8
NCiLcF87V2dNqP8OlEqFQbwIG0J5dmDRLS6m1mBDaLsOrJMaVqQc5RPKtxOj
ncLaDeN8QvPpROBM0NXWOFtC8XRCWb5PfWLSltC8urBHV3phabyAUP5dUE76
qJvGBYTm1wXP20+vZn9jR+g+aBG49Ex97IgdoXlqUTa7ZeOhSHtC90OL+TL+
+US9PaH5dmNNwVTaA7mQ0H3pRkvgs8qVr4SE5t2NtgHHmzU7HQjdnx54yzc8
yX3uwPLvwR9202ElAY6E7lMP5t6/dFt035Hx0MFvibexeO0sxkOHnL3CdT/l
z2I8dHhdJD7U5CpiPHQ4d0N6LCJNxHj0wjcsMX01z4nx6IXmdyvxLqkT49GL
L0Oqp9dGOzEefbi6a31WRaETOWvm0Ye6634iuy4nxqMPhoTDiQI3MeOhR8KW
m/zKz8SMhx67vE5dkaaIGQ89to50bI6qEDMe/ai3Dnh7cFzMePRjrKvx2OI1
zoxHP4zSpPwr4c6MxwCU0y6f6LKdGY8BZCTK+cYXzozHAJLvpFxQCiSMxyCS
3B0uRPpKGI9BxKqcMtujJYzHIKI+zZu7sEDCeAxh8OFgxkq1hPEYgu/00QRH
iQvjMYT84P3eikAXxmMYjzQVwsBvXch/tukJVA==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kgtUk3UYxidaCog4NhiMiZZ60vRoZueUeHn+pAkmmOIF7MBSiIuFstIQ
hbwFGil57Xi85CkhBLSBRxRRY0NBYOIEcRsq4jbkto2BLhBBNLfv7T3nO9/9
fZ/n+b3vRCWGxjjxeLyYN4f9/H9x192KFKE80LA1DeGnCnNcJprp/giW5KZt
Wx/Xqrjsl25edvkUtLr5ythKA73Pxi6J5UVbxgPF8WP2yoO6f/GsfEUdfS9H
ZWn8wkDpdYVlmf3JebQ6bVwC/S/0fxE+9++7t2nyP3DRZcnmuhSjM0mrurBe
Da5fCVZELG/f1q3F5Ij913smXcO6VYNrZVuawPUvRd+/je9aO4zIfvO3LkuJ
gcLtF3t3tIKbV4bk0B1KvZMJMx11AxPC7zwrlHWCm18OmX/jrYKYbtzomfRm
QgV+l0dWmJ2fkZ5KTCiZ5ZHylw2hjoFVyCgW2io+6SF91RgoctNMK++FYau9
oQqap2F8QVAf6a1Bf0jptVvlL/CdXZ7sNh4tydso9x8g/WpM3FMtGJX/EhyN
Wphrxq6eLXhFfmpR/2dyeHrSawQ5BtRB/zTXJPmGxzh/dZh7aeqVUK8hzD4t
Yv9deK/htzkrhzDObz2UxjMX4uKcmMOOrh4FMWNVB0YOZZz/e/DNm/NeXsFQ
NtYux0+DwI/0vHtLhzEuDw3ios35M7uHsbhYe2kRKzuR2pH5FuPy0SE6pavY
8/23mb2bUK6D7fkHm3TlbzMurwaUdn62L1g6nPXa2/U0wC9EdPBg73DG5Xcf
VUPdlGWZI9g8R8MH2Pw4YbB3vDPj8nyAmckL7y4qcWYOeekPMW3lrvV1IS6M
y7cRfI/TZ7IMLux2jb0a4RKtGKHZ5Mq4vB9hdXVkbfLwkczTIbAJh8LPvfrj
2EjG5d8ES6xeKp3ixhztIh4j0lt6SHnNjXE89NiZkSjRhowiP3qESldtONc0
inF89ChZmmn7ItGd/Onh0Szyq3ntzjheBqRUR5V4Bo0mvwbMs+apZh8YTfwM
mHF0en5Ew2jyb0TYV3u/TB3HJ55GxAhWbj8Vz6c8jMh9XhR1s5BPfJtRZVP1
9/bxKZ9mqC6YrB8GeBDvZsgXnD2ZmuFBeT1BgGu1RFvnQfyf4PjH9wMDxALK
7wnKXdVXFFEC2ocWtG1Z/DD0rIDybEHAzVsPB2wC2o8WiHlHFFfnCCnfVvwQ
nCw9nC6kfWlFwgRZ/E9qIeXdisN3++P3iTxpf9oQl71CJV/jSfm3Ye+GHxtM
eZ60T21InMovZzZP4tGOgrO27MI5XsSjHS2ra1P9d3sRj3bM35A4q+mOF/Fo
x2CV4dIJHxHx6EDSYFrZ99Ei4tGBeKFBFP23iHh0YH7O2oSE5yLiYYL+pTn2
1wBvdsDBw4SKnuU5VXu9iYcJrV+fLxNrvYmHGeNLji5KG+dDPMzIOT2jeti3
PsTDDFnSp5UnLvoQDwtMQ9wrg3hi4mFB6ubpi1wXi4mHBSfdg2uMv4mJRycW
WCpuq/Vi4tEJzdUQlXqKL/HohPvP6+YZk3yJhxW5ruacEdd9iYcVmdqwnQFu
EuJhRVHSguB94RLi0YXIqqw97VkS4tEF+TbXiyu7JMSjCxqjNazefwzx6Mb4
4tqeNbvHsP8AS1cHGQ==
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.6`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 7.381229248870826}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}], {193.5, -355.98757751993935`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}], InsetBox[
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg1Qk3Ucx8fkxUHoBs+gpEioTsB8OzHtDL9/s4RyXiKaEtyhFioHCvh2
KiBcpoinDpQMQckTEkVbSgSI6NClTHCBgkwQYRtsIC8DBYZAQG7/f7+753nu
efv9vt/v5+exKXp1OJ/H44W/OczX/8vKcp5cHsfJ/LX7f8KGCmV8vsCW3afj
Zuz0uj99rctL3A91BZZk47BmuO7nBD57n4t5t7bvHNTwyjPPmOsylL5zg44U
j8vp9zKc+uFkY0b9qLw70PzkOgwPVfmbpg3L6f+FKBop3lh9wCS3V+fE+NkX
Q5lyyy2YNyCn/W7gK9VtlSnrpdw7VHp30KsMY/djBVy2kfW/DZNVTP652Z3y
3Dd/q3PK0fBrwvyzej2bdwdD7nZ5OXKtfL6lFNjidSpz285nbP7fELy6cH7b
hVq5YtDrzYR7EMrOqtdOVTI9Ffi09aCPJL5AvtoyUAl/0+WQa+pcUH0PICnJ
9VgZpIB2v7lhJWxqIh2eS2tA9T5E2EXblR+vUiPWLC9GhYl5v8n03c9B9f+D
SZVHDUJfHSiNGuyNC+S1bdSD+qnBO3u8N1QmdSDAMuARWnSXjknTukD9PcLa
01EyY7oR5mmh0scovR4RUHmvD9RvLVJjd1wVLn4Fix11La4YEsbT8/tB/ddh
9pqPFGe5QbxvluP+BG/zbZ0i9plA83iC5QHHHRY2DGHLZnPVo8lb40d8h0Hz
UeMwt2G09OgIzN04mRrLHEy5jU2joHk9RYiS11ntMwaTud3gU/BOJE4od42z
/BqwdeJp9cvSCSyxNGzEDN8FTqUHeITm2YiIa71ti92tiEXeoWcY3CteLyqz
IjTfJtSrow5K1vOJ6qG5mhBwcOaV8Vd8QvN+jl3hQSsWHZtExBaBzUheMOEh
+NCa0PybcfzbqkVJpdbE0i60BZdSE3uufGNDKA8NLp52HTnTakOoHw0kRZ/H
r9tjSygfDaRR2//ts7Mj1J8GEbcid+3IsCOUlxZ3bcrWts+YTKhfLcIjpdLv
iiYTyk8L0ev49PplAkL96/CuYiQsvEZAKE8d+hdlxQlD7VkeOqiSna21BntC
+bbCL2KKrCXGgeXTCg+dqpIbcSAUTyu638uWHPnxLZZXG7i5blVL7R0J5d8G
jf+K5IA0R5ZfG+4WnojOcZ1C6D7ocXxs8wfB56awPPWYu8BlKNpzKqH7ocde
metQz8WpLF8DCtJv7I4SCwndFwNSKnbPclsnZHkb0O9jU3EoQ0jo/rRjbJX1
A1mDkOXfjvrsqnMZ00SE7lM7rhJF+pchIsajA7KcvNbiLBHj0YF9fHney2ci
xqMDgpO89gE3J8ajA6GNs9IUIU6MxwuU329OCstyYjxewM2UFF/Z6MR4vEDV
468X2k5zZjw6sdXkuZoLdiapFh6d8AqJDzX94sx4dOL7hhlPfq93Zjy6sCXK
beZnYo7x6ELzX3NuXwjiGI8ulM+va9GlcYxHNwozSmvHqjnGoxthBftqBhzF
jEc3hMu9eA9WiBmPHqzx70uKSxEzHj2YpcwMd6wQMx49CJ7gihKtXRgPI+7c
NEbXLnVhPIz4IvPxYUGiC+NhxIjkdZtnmQvj0YvkpWmR04ddGI9ecL1z+vmf
uDIevTgv8pFU7XBlPPowzJmW7P/DlfwH3hwIHg==
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhQOI4IIosoUlJBAVVBaHymK15zfVurRBQLQgWNBWU6UgaFVE
C2Jh7JRRxK1D7TCC2ipYVFzqGqBqXQbQRJPYAEogQBLWaUVxHKQm7/bOvHn7
veec74rWbopdZ83j8da9P8zn/4u7HqjZ6Vy1SJedD/uhrYveFBjp/jAi1wQL
P2HtNVcFBd0xV0vxsswq4uhqLb0/ia5H59LTJA01P5eY6wwUFQ+V4fMr6Psq
dJTFr7K+L0dPjPnJBfibcDozVgnu/0uImBvnrSxtxljNiYx5Y/9AxOg8XlZU
O7h+12Al9E2RBBgQkFT056D/TWhfFqry23rA9ZdjeZNEIpMP4OT7vzUnalEd
Ll8xb/gfmleHuqx9zzVBgwi11G10OGRPnRX3mubfgbV087TA9De4Pej/fsJd
rB6efvmX796SnnvobUy8wr4fRqxl4H3ETD3r1pA7Qvoe4DCrKxu8xmO6bHPD
h7jkJLQRnLdinN56nOPvTz9bbs0yzfIyGhByfkqyqdiGcfobscvtYqQgZxTj
aDzGD5qHMYUyW8b5eYw9IbmqzKjRbLFlgAJa27qzE0LtGOdPgSXL1SuPudgz
87SkIiUqBasapK/sGef3CQKyZJnxqjHMYkfzBA6FM21N1WMZ5/8prB2l1YFF
45iPWY5AhXu1YV6fbhzPuDxUuNB0d9f6BQ5Mtt5cauyvinhw3HsC4/LRQL05
qcXt1QRm7uZcpYFh5ZuKzgZHxuX1DP2GlssHYyeyV+Z2g8+wM9HR9PLWRMbl
9zdet0w84RcwiX1kaahF7tr0vGmHJ1GeWlSqBIH2I5OYRV5BEyo3RH5Tv8GJ
8m1GxZz4iOynTqyh3lzNmF17QOmKyZR3C06FDql+PTOZuVgEPkfQ7kuLg5yd
Kf/ncBtdXnk+x5lZ2iW9QOzgMftZRmfi0YoZ17NGLix3IT+tGLe0UTFb7kJ8
WiFJOZQm93clf62wCr64TXrIlXjp8FZ/epl+2JX86qBSh4vyv3YjfjoEIk4X
+MSN/LdhyY5lCW3z3IlnGzLTE2aUn3anPNrQWi3TpE3mE992GIuKaxbm8Cmf
dije9c2YbuQT73aMktlFecZ5UF566BubPubXeBB/PRZ6/+Qhnu5J+emxwK7v
8YdHPGkfOmCTEJ6QwvOiPDuwLXXBseJUL9qPDvAf8bY9UntRvp0o3/mili/x
pn3pRNBMT2nG796UdyfmvG7sVLoLaH+6sO+37hXz8wWUfxfyXjgk3+gX0D51
4UZpxB2W6EM8DNB9sX2h4i8f4mFARvTl6DC+kHgYUHKtdemYZULiYYAtX7FJ
ly8kHkY4MetS+XUh8TDiW5+Ri8cHhMTDiPrCPQf2ThURDxOChg86bUkSsQMW
HibEd5XwvzwoIh4mrAnbXvD5fRHx6Mb81OKw6Hci4tGN3XfbHKQf+BKPbpRq
bXqlG32JRw/CBkKvxB73JR49iPwq+bNEtS/x6EHyOvEh2Xg/4tGLtDB1+naJ
H/Hoxang3Ds/ZvkRj16UVOp2lFX5EY8+HLHL2npD70c8+uCg3VL6zENMPPoQ
knOzeShaTDz6scex191rr5h49EPXEzNXcktMPPqh1HSEpP4rJh4DGKqL0h0N
mML+Az5E5cA=
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "0.8`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 10.581229248870827`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}], {580.5, -355.98757751993935`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}], InsetBox[
GraphicsBox[{{}, {{}, {}, 
{RGBColor[0.368417, 0.506779, 0.709798], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kglQU2cUhYMRCQgIITEYFIsVC0pdIl1cz0+lBbUdCmKHKjqtFh0Hq9ji
4E7L5rhQESimiE4VnApWJIIKVgigDgoS2ZRV1hASCCQFwzKilOT/e2fee/O2
e88533XZsT8gZAqHwwmZPIzX/8vMdOYVHxVk+3QciYFozCVj1RVzdp+Mo6nH
csYOcovznWP7/PMvY7Wi9SPfCDP2PgOqQdfv/1wzIU/9w1iZ4KyfGvqB4o2c
fp+NRfGWN5Qpo3Ktv/GJDBbZjS2uCQY5/T8PieGb757OHZRb1aeHrbG6h4Er
q1WJIzrWrwDzratHDU965e7B50oNbg/g5ZOlb6vtZv2L4GMhs53N75BnTP5d
n16MzyWeJ5+PN7B5JSjNPSs+61YlX26qh2iyHc16daKIzX+EDLW0JOhMDB4a
3CYnPEZe+uYIQ0AJqJ4yZH04yCsNr0KAaeATlB+PypdZNYDqe4oL7/1SvGF/
GzqOGBuWo+DwA73dpi5Qvc+wQ+/S90mjCgeM8sIqkWa554tzG3tB9SuQuanm
7Z4b/aA0qmA/dLP03TE9qJ8qtLgp/5HMHoSvaUA1LBu2B9TmDoH6q0aj/4Qi
yNsA47TgczVw7Q/cJ1EMg/qtxXwxWV/mPwqTnfpa8J0uOkU8HwP1XweprZ8u
1ecN5hrlOL/ArZ36nNj74yyPF7BKGJLFuL3D7l3Geol9dXec2s9PsHzq8avF
cNaPRznE2E2QXY9Hu71kTYvMCM2rARve95UGNpmRYWM7QwO2bKkZGo6bQmh+
jZg+7bOU5mVcstbUsAlbdTNyxM1cQvNsQs/SYm151FRikhfbjHHd7/Fm7uaE
5tuCqbe9omoqzUnlM2O1gCx4EPn1gWmE5v0KbncihPEOFkRoEtgKnrdNiTTP
gtD8W+GpjeZEB/KIqV1wGyyu++3cMsQjlEc7pANzUz3OWzI/7agrTI6dttiK
UD7tKI8cWPf6qRXz146k8Io1ZiHTCeXVgcJKruTTienMbwfS/PJIutSaUH4d
kCyIM2xcZsP8d+LyGS//VU9tCOXZiUz1vZrD39myPDrhu8I+csaILaF8u/D6
75MabvwMlk8Xklq/mrPE3o5QPF3Qzg25mLjVjuWlRFDmrqu3rtkRyl8JriR/
xSmdHctPibLoSMxbYU/oPnRjb9vC0ONR9izPbvTerPS7VGFP6H5040ZK4PYY
AZ/lq8LypqQez218QvdFhW1pEmfZNT7LW4UX0oOOb/r5hO5PD/79ra9a8LED
y78Hmri3HpwTDmyfehCm/eF10WMHxkONkMBTLQE2AsZDjdVXy2IKAwWMhxoN
LhmHOGkCxkONHGvN0jldAsZDg6yXtwXChULGQ4PxmNAjmjAh46EBd+Xd6xfv
CRmPXmy1/rLe452QJJh49GKVe8WSy94zGY9eOEsOK/SnZzIefRiRjZW6Vs9k
PPpgfmjMbq1IxHj0wenC85yV20SMhxa8AkHC7HQR46HFt3uXJSrVIsZDi1Cu
MiFpsSPj0Q9eUlHwwnBHxqMf3cmXCjMLHBmPyfd/iX7mTzgyHgPY3KzyDPGe
xXgMoG5iY+7VU7MYjwF8s27vhQrFLMZDh+KEzOh2BzHjocP9EYNDZ5CY8dCh
Nmxea9UlMeOhh4ev+083O8XkP0zMAFE=
"]]}, 
{RGBColor[0.880722, 0.611041, 0.142051], PointSize[0.01388888888888889], 
           AbsoluteThickness[1.6], LineBox[CompressedData["
1:eJw1kg1Q03UYx0FKFAUVGW/jbYBsQ8mTkV508v1NU/A1JQkTUOFEEk8BPU0h
UrskonkgUamnGIpk8tYSDihUcBiI4LFAZxOQgYwxxgQDNfElt//Tc/e///vz
fL/fzyOISwqPn2RhYRH/5jCd/y/ueuRamkNZqCb1K9xRthdv0g/QfR4UO7++
H+3Xc63a4+jQ+up8XJCcDtgRr6T3hZBJJZq46xevnTppql8Q8/lyiWXYDXDf
lyE99nFBlVoFw3rTEzmk4b+FGFf1gvu/Au/InEPle3WwUZ1PXmxThU5vn9aV
iUZw/WoQLh0xJLs+hjg6+/q4qBYOP3c8nWwYo/5XkaAdSzQUP0Xhm79V5+vw
p+239f1Rz2lePeo2ZQoPv3wBibkUqI255MvPfU3zG8A7nq44YGvJFOOiNxNu
gB9r3JGzehLj9DQiXx2UrsqwYuHmgU1oXKsQN9a+xTh9NyFxrHRqM77NNKmm
hs1YM9rlG+VuzTi9LSjBaN3NsCksxSQvuRXJuft5+1OmMk7/bZyId3W9/KMN
42i0YSL5WHHLH9MY56cNv7fky627p7Mw8wAlEtWln5x7Zcs4f0o0n10X3+E+
g5mmRWf/hWDZeMOA+0zG+W1HYgkzbimZycx2VO0o1cd5bw+exTj/Hfii/tOD
05pmMU+THI870AR0DW6OsKc87mBX0Mwju3rtWcJ2U91FiEv4iSVJsykfFeY/
S5usm5jNTN0cylRIbNq5ZnumA+V1D+nSCEWDA489MbUbv4fYbWNVUwt4lN/f
OLuw5fvgAEcWYm6oRsya8tdRNY6Upxo1RWp9yjInZpZ39D6KhCOVh5ROlG8n
Cn71WJwZ48xaW0zViYi9iqFjg86UdxfmObfcyt3nwnhmgd3419HtuxOWrpR/
N7aMLyj/6ZgrM7eLfoC6VUs/LHHhE48e7Llisav2Ap/89CCmKKBDucCN+PTA
OShQPHzFjfz1YJ1Vs9FupTvx0mDUfzhr0V138quBrMd2VUKcB/HTwD6rPiHf
6EH+e3HYc9C6K9WTePbiI/vKguNWXpRHL8ryFbnVK7yIbx/qy/dE9mV7UT59
iNxXODTrrhfx7kNTUFboB24Cyush/qkO25oaJyD+D5Fz8Kio4qKA8nsIFD7I
GTUKaB/6cVkRmRn4rjfl2Q953MUX+9O8aT/6kVJ1XH213pvy1aKG72pnM8WH
9kWLreqGLzeu9aG8tThUWOdzKc+H9mcAnh12na/UPpT/ALIPFv7wscCX9mkA
vBkZIfIEX+Khgzan+6pdmS/x0EF2xt4qacyXeOjgFND2RBk8h3jo4LCJn7Ho
yBziMYjWoo1nzjbOIR6DcOFHiG3s/IjHIDZmrXb6bIMf8dDj9nLjBu0pP5Zj
5qGHPn1rW6TGj3joIUlZcuCWUEg8hiCrLn1fultIPIawoGKhTU2FkHgMYfo6
Q3vghJB4GNAXFPlNmVREPAyIOnyAPy9TRDwM8HJpSiu+LSIew5DJA0/P44mJ
xzC0mU5J5VFi4jGM954tG5CcExMPI5apnj+t0YmJhxGV17vypPP9iYcR21a/
kjfv8ycejxDorlwaUetPPB5BYBcAzaS5xOMRLELCTu5eMZd4jGCzqDL0ZfZc
9h+btO2Q
"]]}}, {}, {}, {{}, {}}},
AspectRatio->0.6180339887498948,
Axes->{True, True},
AxesLabel->{None, None},
AxesOrigin->{0, 0},
DisplayFunction->Identity,
Frame->{{True, True}, {True, True}},
FrameLabel->{{
FormBox["\"(\\!\\(\\*SubscriptBox[\\(E\\), \\(1\\)]\\),\\!\\(\\*SubscriptBox[\
\\(E\\), \\(0\\)]\\))\"", TraditionalForm], None}, {
FormBox["\"s\"", TraditionalForm], None}},
FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
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
PlotLabel->FormBox[
TemplateBox[{"\"\[Alpha]=\"", "1.`"}, "RowDefault"], TraditionalForm],
PlotRange->{{0, 1.}, {-3.0000000000000004`, 13.781229248870826`}},
PlotRangeClipping->True,
PlotRangePadding->{{
Scaled[0.02], 
Scaled[0.02]}, {
Scaled[0.05], 
Scaled[0.05]}},
Ticks->{Automatic, Automatic}], {967.5, -355.98757751993935`}, 
       ImageScaled[{0.5, 0.5}], {360., 222.4922359499621}]}}, {}},
ContentSelectable->True,
ImageMargins->0.,
ImageSize->{605.9999999999982, Automatic},
PlotRangePadding->{6, 5}]\)
Quit
