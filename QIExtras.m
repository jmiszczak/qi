(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package header*)


(* Mathematica Package *)

BeginPackage["QIExtras`", { "QI`"}]
(* Exported symbols added here with SymbolName::usage *)  
Unprotect@@Names["QIExtras`*"]
Clear@@Names["QIExtras`*" ]


(* ::Section:: *)
(*Public definitions*)


(* ::Subsection::Closed:: *)
(*Common matrices*)


KroneckerDeltaMatrix::usage = "KroneckerDeltaMatrix[i,j,d] returns d\[Cross]d matrix with 1 at position (i,j) and zeros elsewhere.";

PauliMatrices::usage = "Predefined list of Pauli matrices {\!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\)}. Use Map[MatrixForm[#]&,PauliMatrices] to get this list in more readible form.";

GellMannMatrices::usage = "List of Gell-Mann matrices. Use Map[MatrixForm[#]&,GellMannMatrices] to get this list in more readable form.";

UpperTriangularOnes::usage = "UpperTriangularOnes[k,dim] returns not normal matirx of dimension dim with 1<k<dim-1 bands of ones over the diagonal.";

UpperBandOnes::usage = "UpperBandOnes[k,dim] returns not normal matirx of dimension dim with bands at position k of ones over the diagonal.";


Swap::usage="Swap[d] returns permutation operator \!\(\*UnderoverscriptBox[RowBox[{\" \", \"\[Sum]\"}], RowBox[{\"i\", \"=\", \"0\"}], RowBox[{\"n\", \"-\", \"1\"}]]\)\!\(\*UnderoverscriptBox[RowBox[{\" \", \"\[Sum]\"}], RowBox[{\"j\", \"=\", \"0\"}], RowBox[{\"n\", \"-\", \"1\"}]]\)|i\[RightAngleBracket]\[LeftAngleBracket]j|\[CircleTimes]|j\[RightAngleBracket]\[LeftAngleBracket]i\[VerticalSeparator] acting on d-dimensional, d=\!\(\*SuperscriptBox[\"n\", \"2\"]\) space and exchanging two \[Sqrt]d-dimensional subsystems.";

QFT::usage = "QFT[n,method] - quantum Fourier transform of dimension n. This function accepts second optional argument, which specifies method used in calculation. Parameter method can be equal to 'Symbolic', which is default, or 'Numerical'. The second option makes this function much faster.";

GeneralizedPauliX::usage = "Generalized Pauli matrix X. See also: \!\(\[Sigma]\_x\)";

GeneralizedPauliZ::usage = "Generalized Pauli matrix Z. See also: \!\(\[Sigma]\_z\)";


Ket::usage = "Ket[i,d] returns |i\[RightAngleBracket] in d-dimensional Hilbert space. See also: StateVector for a different parametrization.";

Ketbra::usage = "This function can be used in two ways. Ketbra[i,j,d] returns \[VerticalSeparator]i\[RightAngleBracket]\[LeftAngleBracket]j\[VerticalSeparator] acting on d-dimensional space. See also: Proj. Ketbra[v1,v2] returns the appropriate operator for vectors v1 and v2..";

KetFromDigits::usage = "<f>KetFromDigits</f>[<v>str</v>,<v>bs</v>] - ket vector labeled by a number given as a string <v>str</v> in the base <v>bs</v>.\n
<f>KetFromDigits</f>[<v>ls</v>,<v>bs</v>] - ket vector labeled by a number with digits, in the base <v>bs</v>, provided in the list <v>ls</v>.";

MaxMix::usage = "MaxMix[n] - the maximally mixed state in a n-dimensional space of density matrices.";

MaxEnt::usage = "MaxEnt[n] - maximally entangled state in n dimensional vector space. Note that n must be perfect square.";

WernerState::usage = "WernerState[d,p] - generalized Werner state d\[Cross]d-dimensional system with the mixing parameter p\[Element][0,1]. This state is defined as p Proj[MaxEnt[d]] + (1-p) MaxMix[d],. See also: MaxEnt, MaxMix..";

IsotropicState::usage = "IsotropicState[d,p] - isotropic state of dimensions d\[Cross]d with a parameter p\[Element][0,1]. This state is defined as p Proj[MaxEnt[d]] + (1-p)/(\!\(\*SuperscriptBox[\"d\", \"2\"]\)-1)(I-Proj[MaxEnt[d]]]). This family of states is invariant for the operation of the form U\[CircleTimes]\!\(\*SuperscriptBox[\"U\", \"\[Star]\"]\).";



ReshufflePermutation::usage = "ReshufflePermutation[dim1,dim2] - permutation matrix equivalent to the reshuffling operation on dim1\[Cross]dim2-dimensional system. See also: Reshuffle.";

ReshufflePermutationPrim::usage = "ReshufflePermutation2[dim1,dim2] - permutation matrix equivalent to the alternative reshuffling operation on dim1\[Cross]dim2-dimensional system. See also: Reshuffle.";


ReshufflePrim::usage = "\
ReshufflePrim[\[Rho], drows, dcols] for a matrix of dimensions (drows[[1]]\[Times]drows[[2]])\[Times](dcols[[1]]\[Times]dcols[[2]]) returns reshuffled (prim) matrix with dimensions \
(dcols[[2]]\[Times]drows[[2]])\[Times](dcols[[1]]\[Times]drows[[1]]), \
drows and dcols parameters can be ommited for a square matrix.";


IdentityChannel::usage = "IdentityChannel[n,\[Rho]] - apply the identity operation to a n-dimensional density matrix \[Rho].";


TransposeChannel::usage = "TransposeChannel[n,\[Rho]] - apply the transposition operation to a n-dimensional density matrix \[Rho]. Note that this operations is not completely positive.";


DepolarizingChannel::usage = "DepolarizingChannel[n,p,\[Rho]] - apply the completely depolarizing channel with parameter p acting to a n-dimensional input state \[Rho]. See also: QubitDepolarizingKraus, HolevoWernerChannel.";


HolevoWernerChannel::usage = "HolevoWernerChannel[n,p,\[Rho]] - apply the Holeve-Werner channel, also known as transpose-depolarizing channel, with parameter p acting to a n-dimensional input state \[Rho]. See also: DepolarizingChannel.";

GeneralizedPauliKraus::usage = "GeneralizedPauliKraus[d,P] - list of Kraus operators for d-dimensional generalized Pauli channel with the d-dimesnional matrix of parameters P. See: M. Hayashi, Quantum Information An Introduction, Springer 2006, Example 5.8, p. 126.";

ExtendKraus::usage = "ExtendKraus[ch,n] - produces n-fold tensor products of Kraus operators from the list ch.";

Concurrence4::usage = "Concurrence4[\[Rho]] returns quantum concurrence of a density matrix \[Rho] representing a state of two-qubit system. This function uses Chop to provide numerical results.";

EntanglementOfFormation4::usage = "EntanglementOfFormation4[\[Rho]] returns entanglement of formation of density matrix  \[Rho] representing a state of two-qubit system."

Negativity::usage = "Negativity[\[Rho],{m,n}] returns the absolute value of the sum of negative eigenvalues of the density matrix \[Rho]\[Element]\!\(\*SubscriptBox[\"\[DoubleStruckCapitalM]\", 
RowBox[{\"m\", \"\[Cross]\", \"n\"}]]\) after their partial transposition with respect to the first subsystem. See: G. Vidal, R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002) DOI[10.1103/PhysRevA.65.032314].";



(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


QubitBitflipChannel::usage = "BitflipChannel[p,\[Rho]] applies bif-flip channel to the input state \[Rho]. See also: QubitBitflipKraus.";


QubitPhaseflipChannel::usage  = "QubitPhaseflipChannel[p,\[Rho]] applies phase-flip channel to the input state \[Rho]. See also: QubitPhaseflipKraus.";


QubitBitphaseflipChannel::usage  = "QubitBitphaseflipChannel[p,\[Rho]] applies bit-phase-flip channel to the input state \[Rho]. See also: QubitPhaseflipKraus.";


QubitDepolarizingKraus::usage = "Kraus operators of the depolarizing channel for one qubit. Note that it gives maximally mixed state for p=1.";


QubitDecayKraus::usage = "Kraus operators of the decay channel, also know as amplitude damping, for one qubit.";


QubitPhaseKraus::usage = "Kraus operators for one qubit phase damping channel.";


QubitBitflipKraus::usage = "Kraus operators for one qubit bit-flip channel.";


QubitPhaseflipKraus::usage = "Kraus operators for one qubit phase-flip channel.";


QubitBitphaseflipKraus::usage = "Kraus operators for one qubit bit-phase-flip channel.";


QubitDynamicalMatrix::usage = "QubitDynamicalMatrix[\!\(\*SubscriptBox[\"\[Kappa]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Kappa]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Kappa]\", \"z\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"z\"]\)] returns parametrization of one-qubit dynamical matrix. See: I. Bengtsson, K. Zyczkowski, Geometry of Quantum States, Chapter 10, Eg.(10.81).";


QubitDaviesSuperoperator::usage = "QubitDaviesSuperoperator[a,c,p] returns a superoperator matrix for one-qubit Davies channel with parameters a and c and the stationary state (p,1-p).";


(* ::Subsection::Closed:: *)
(*One-qutrit channels*)


QutritSpontaneousEmissionKraus::usage="QutritSpontaneousEmissionKraus[A1,A2,t] Kraus operators for qutrit spontaneous emission channel with parameters A1, A2, t >= 0. See: A. Checinska, K. Wodkiewicz, Noisy Qutrit Channels, arXiv:quant-ph/0610127v2.";


(* ::Subsection::Closed:: *)
(*Entropy*)


Log0::usage = "Log0[x] is equal to Log[2,x] for x>0 and 1 for x=0.";


QuantumEntropy::usage = "QuantumEntropy[m] - von Neuman entropy for the matrix m.";


QuantumChannelEntropy::usage = "QuantumChannelEntropy[ch] - von Neuman entropy of the quantum channel calculated as a von Neuman entropy for the image of this channel in Jamiolkowski isomorphism. See also: Jamiolkowski, Superoperator.";



(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


VandermondeMatrix::usage = "VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}] - Vandermonde matrix for variables (\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)).";


ProdSum::usage = "ProdSum[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] gives \!\(\*SubsuperscriptBox[\"\[Product]\", 
RowBox[{\"i\", \"<\", \"j\"}], \"n\"]\)\!\(\*SubscriptBox[\"x\", \"i\"]\)+\!\(\*SubscriptBox[\"x\", \"j\"]\).";


ProdDiff2::usage = "ProdDiff2[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] is equivalent to Det[VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}]\!\(\*SuperscriptBox[\"]\", \"2\"]\) and gives a discriminant of the polynomial with roots {\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}.";


ProbBuresNorm::usage = "ProbBNorm[n] - Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Bures distance.";


ProbBures::usage = "ProbBures[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},\[Delta]] - Joint probability distribution of eigenvalues \[Lambda] = {\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}of a matrix according to Bures distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: ''Indicator''";


ProbHSNorm::usage = "Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Hilbert-Schmidt distribution.";


ProbHS::usage = "ProbHS[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},\[Delta]] Probability distribution of eigenvalues \[Lambda] = {\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)} of a matrix according to Hilbert-Schmidt distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: ''Indicator''";


RandomProductKet::usage = "RandomProductKet[{dim1,dim2,...,dimN}] - random pure state (ket vector) of the tensor product form with dimensions of subspaces specified dim1, dim2,...,dimN.";


RandomProductNumericalRange::usage = "RandomLocalNumericalRange[M,{dim1,dim2,...,dimN},n] returns n points from the product numerical range of the matrix M with respect to division specified as {dim1,dim2,...,dimN}. Note that dim1\[Cross]dim2\[Cross]...\[Cross]dimN must be equal to the dimension of matrix M.";


RandomMaximallyEntangledNumericalRange::usage = "RandomMaximallyEntangledNumericalRange[M,n] returns n points from the maximally entangled numerical range of the matrix M with respect to division Sqrt[dim[M]]\[Cross]Sqrt[dim[M]].";


NumericalRangeBound::usage = "NumericalRangeBound[A,dx] - bound of numerical range of matrix A calculated with given step dx. Default value of dx is 0.01. Ref: Carl C. Cowen, Elad Harel, An Effective Algorithm for Computing the Numerical Range. Technical report, Dep. of Math. Purdue University, 1995.";


IntegrateSU2::usage = "IntegrateSU2[f,U] - gives the integral \[Integral]f dU, where dU is Haar measure on the group SU(2) - special unitary matrices of size 2. 
\nIntegrateSU2[f,U,V,...] - gives the multiple integral  \[Integral]f dU dV dW ... on the group SU(2). 
\nExample: Integration squares of absolute values of elements of ranom unitary matrix: IntegrateSU2[Abs[U\!\(\*SuperscriptBox[\(]\), \(2\)]\),U]"


Unitary3::usage = "<f>Unitary3</f>[<v>\[Alpha]</v>,<v>\[Beta]</v>,<v>\[Gamma]</v>,<v>\[Tau]</v>,<v>a</v>,<v>b</v>,<v>c</v>,<v>ph</v>] returns the Euler parametrization of <v>U</v>(3).";

Unitary4Canonical::usage = "<f>Unitary4Canonical</f>[<v>a1</v>,<v>a2</v>,<v>a3</v>] returns the parametrization of non-local unitary matrices for two qubits. \
See: B. Kraus, J.I. Cirac, Phys. Rev. A 63, 062309 (2001), quant-ph/0011050v1.";

Mub::usage ="<f>Mub</f>[<v>p,m</v>] for prime number <v>p</v> and positive integer <v>m</v> returns <v>p^m+1</v> mutually unbaised bases of <v>p^m</v> dimensional Hilbert space." 


(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"] (* Begin Private Context *) 

QIDocRep = {"<v>" -> "\!\(\*StyleBox[\"" , "</v>" -> "\", \"TI\"]\)", "<f>"->"\!\(\*StyleBox[\"", "</f>" -> "\", \"Input\"]\)", "<s>" -> "", "</s>" -> ""} 
(MessageName[Evaluate[ToExpression[#]], "usage"] = StringReplace[MessageName[Evaluate[ToExpression[#]], "usage"],QIDocRep])& /@ Names["QIExtras`*"];


(* ::Subsection:: *)
(*Internal functions (without usage strings)*)


qiExtrasAuthors = "Jaroslaw Miszczak <miszczak[at]iitis[dot]pl>, Piotr Gawron <gawron[at]iitis[dot]pl>, Zbigniew Puchala <z.puchala[at]iitis[dot]pl>";

qiExtrasLicense = "GPLv3 <http://www.gnu.org/licenses/gpl.html>";

qiExtrasHistory = {
	{"0.0.1", "13/10/2011", "Zbyszek,Jarek,Piotr", "Some functions moved from QI to QIExtras."},
	{"0.0.6", "17/12/2011", "Jarek", "Fixed Negativity"},
	{"0.0.7", "23/12/2011", "Zbyszek", "Mubs"},
    {"0.0.8", "23/01/2012", "Jarek", "Versioning added and KetFromDigits improved."},
    {"0.0.9", "02/02/2012", "Gawron", "Concurrence4 fixed, thanks Maciej Demianowicz."},
	{"0.0.10", "04/02/2012", "Jarek", "BaseVectors and BaseMatrices moved from QIExtras to QI"}
};  

qiExtrasVersion = Last[qiExtrasHistory][[1]];

qiExtrasLastModification = Last[qiExtrasHistory][[2]];

qiExtrasAbout = "QIExtrtas is a package of functions for Mathematica computer algebra system, .";



(* ::Subsection::Closed:: *)
(*Common matrices*)


KroneckerDeltaMatrix[m_,n_,dim_]:=Block[{mtx},
    mtx=Table[0,{dim},{dim}];
    mtx[[m,n]]=1;
    mtx
];


PauliMatrices = GeneralizedPauliMatrices[2];


GellMannMatrices = GeneralizedPauliMatrices[3];


UpperTriangularOnes[rn_,dim_]:=Table[If[i<j&&j-i<rn+1,1,0],{i,1,dim},{j,1,dim}];


UpperBandOnes[bandNo_,dim_]:=Table[If[i<j&&j-i==bandNo,1,0],{i,1,dim},{j,1,dim}];


Swap[dim_]:=Plus@@Flatten[Table[KroneckerProduct[Ketbra[i,j,Sqrt[dim]],Ketbra[j,i,Sqrt[dim]]],{i,0,Sqrt[dim]-1},{j,0,Sqrt[dim]-1}],1];

QFT[n_,method_:"Symbolic"]:=Block[{\[Omega]},
    If [method=="Numerical",\[Omega]=N[Exp[2 \[Pi] I/n]],\[Omega]=Exp[2 \[Pi] I/n]];
    Table[\[Omega]^(i*k) ,{i,1,n},{k,1,n}]
];

GeneralizedPauliX[d_]:=Sum[Ketbra[Mod[j-1,d],j,d],{j,0,d-1}];

GeneralizedPauliZ[d_]:=DiagonalMatrix[Table[Exp[2\[Pi]*I*j/d],{j,0,d-1}]];


Ket[label_,dim_]:=Block[{vec},
    If[label<dim,
        vec =Table[0,{dim}];
        vec[[label+1]]=1;
        vec,
        (* else *)
        Message[Ket::argerr,label,dim];
    ]
];
Ket::argerr = "Requested index `1` not smaller than dimension of vector: `2`.";


Ketbra[i_,j_,dim_]:=KroneckerProduct[Ket[i,dim],Ket[j,dim]];
Ketbra[v1_?VectorQ,v2_?VectorQ]:={v1}\[ConjugateTranspose]\[CircleTimes]{v2};


KetFromDigits[l_List,b_:2]:=Ket[FromDigits[l,b],b^Length[l]];


KetFromDigits[str_String,b_:2]:=Block[{l=IntegerDigits[FromDigits[str,b],b,Length[Characters[str]]]},KetFromDigits[l,b]];


MaxMix[n_Integer]:=(1/n)*IdentityMatrix[n];


MaxEnt[dim_]:=Block[{subDim=Sqrt[dim]},
  If[IntegerQ[subDim],1/Sqrt[subDim] Plus@@Table[Ket[i,subDim]\[CircleTimes]Ket[i,subDim],{i,0,subDim-1}]]
];


WernerState[dim_,p_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],
    p Proj[MaxEnt[dim]] + (1-p)*MaxMix[dim],
    Message[WernerState::argerr,dim];
]
];
WernerState::argerr = "The first `1` argument is not a perfect square.";


IsotropicState[dim_,p_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],
    p Proj[MaxEnt[dim]] + (1-p)/(subDim^2-1)(IdentityMatrix[dim]-Proj[MaxEnt[dim]]),
    Message[WernerState::argerr,dim];
]
];
IsotropicState::argerr = "The first `1` argument is not a perfect square.";


ReshufflePermutation[dim1_,dim2_]:=Block[{initPos},
    initPos=Flatten[Reshuffle[Partition[Range[dim1*dim1*dim2*dim2],dim1*dim2],{dim1,dim2},{dim1,dim2}]];
    Table[UnitVector[dim1*dim1*dim2*dim2,Position[initPos,i][[1,1]]],{i,1,dim1*dim1*dim2*dim2}]
];


ReshufflePermutationPrim[dim1_,dim2_]:=Block[{initPos},
    initPos=Flatten[ReshufflePrim[Partition[Range[dim1*dim1*dim2*dim2],dim1*dim2],{dim1,dim2},{dim1,dim2}]];
    Table[UnitVector[dim1*dim1*dim2*dim2,Position[initPos,i][[1,1]]],{i,1,dim1*dim1*dim2*dim2}]
];


ReshufflePrim[\[Rho]_]:=Block[{dim},
    dim = Sqrt[Length[\[Rho]]];
    If[And [SquareMatrixQ[\[Rho]] , IntegerQ[dim]] ,
        ReshufflePrim[\[Rho], {dim,dim},{dim,dim}]
    (*else*),
        Message[ReshufflePrim::argerr]
    ]
]
ReshufflePrim::argerr = "ReshufflePrim works only for square matrices of dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\)\[Times]\!\(\*SuperscriptBox[\"d\", \"2\"]\), where d is an Integer, for other dimensions use Reshuffle with extra parameters";

ReshufflePrim[A_,{n_,m_}]:=Flatten[
    Table[Flatten[Part[A,1+i1;;n[[2]]+i1,1+i2;;m[[2]]+i2]\[Transpose]],{i2,0,m[[1]]*m[[2]]-1,m[[2]]},{i1,0,n[[1]]*n[[2]]-1,n[[2]]}]
,1]\[Transpose];



IdentityChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]];


TransposeChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]\[Transpose]];


DepolarizingChannel=Function[{dim,p,\[Rho]},(1-p) IdentityMatrix[dim].\[Rho]+(p) Tr[\[Rho]]MaxMix[dim]];

HolevoWernerChannel=Function[{dim,p,\[Rho]},( p \[Rho]\[Transpose]+ (1-p)Tr[\[Rho]]MaxMix[dim])];


GeneralizedPauliKraus[d_,p_]:= Flatten[Table[Sqrt[ p[[i+1]][[j+1]]] (MatrixPower[GeneralizedPauliX[d],i].MatrixPower[GeneralizedPauliZ[d],j])\[ConjugateTranspose],{i,0,d-1},{j,0,d-1}],1];

ExtendKraus[operators_,n_] := Block[{tpl},tpl=Tuples[operators,n];Table[KroneckerProduct@@tpl[[i]],{i,1,Length[tpl]}]];

Concurrence4[m_]:=Block[{sqrtM=MatrixSqrt[m],evl,R,rhotilde},
	rhotilde=(sy\[CircleTimes]sy).Conjugate[m].(sy\[CircleTimes]sy);
	R=MatrixSqrt[sqrtM.rhotilde.sqrtM];
	evl=Eigenvalues[R];
    Max[ 0 , Re[evl[[1]]-evl[[2]]-evl[[3]]-evl[[4]]] ]
];

EntanglementOfFormation4[rho_] := Block[{h},
  h[x_] := -x*Log0[x] - (1 - x)*Log0[1 - x];
  h[1/2*(1 + Sqrt[1 - Concurrence4[rho]^2])]
];

Negativity[\[Rho]_, {m_, n_}] := Plus@@Abs[Select[Eigenvalues[PartialTranspose[\[Rho], {m, n}, {1}]], # < 0 &]];


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


QubitBitflipChannel=Function[{p,\[Rho]},(1-p) id.\[Rho]+p sx.\[Rho].sx\[ConjugateTranspose]];


QubitPhaseflipChannel=Function[{p,\[Rho]},(1-p) id.\[Rho]+p sz.\[Rho].sz\[ConjugateTranspose]];


QubitBitphaseflipChannel=Function[{dim,p,\[Rho]},(1-p) id.\[Rho]+p sy.\[Rho].sy\[ConjugateTranspose]];


QubitDepolarizingKraus[p_]:={\[Sqrt](1-3p/4)*id,\[Sqrt](p/4)*sx,\[Sqrt](p/4)*sy,\[Sqrt](p/4)*sz};


QubitDecayKraus[p_]:={ {{1,0},{0,\[Sqrt](1-p)}} , {{0,\[Sqrt]p},{0,0}} };


QubitPhaseKraus[p_]:={ {{1,0},{0,\[Sqrt](1-p)}} , {{0,0},{0,\[Sqrt]p}} };


QubitBitflipKraus[p_]:={ \[Sqrt](1-p)*id,\[Sqrt]p*sx};


QubitPhaseflipKraus[p_]:={\[Sqrt](1-p)*id,\[Sqrt]p*sz};


QubitBitphaseflipKraus[p_]:={\[Sqrt](1-p)*id,\[Sqrt]p*sy};


QubitDynamicalMatrix[kx_,ky_,kz_,nx_,ny_,nz_]:= 1/2{
    {1 + nz + kz, 0, kx + I ky, nx + ny},
    {0, 1 - nz + kz, nx - ny, kx + I ky},
    {kx - I ky, nx - ny, 1 - nz - kz, 0},
    {nx + ny, kx - I ky, 0, 1 + nz - kz}
}


QubitDaviesSuperoperator[a_,c_,p_]:={{1 - a,0,0,(a*p)/(1-p) },{0,c,0,0},{0,0,c,0},{a,0,0,1-(a*p)/(1-p)}};


(* ::Subsection::Closed:: *)
(*One-qutrit channels*)


QutritSpontaneousEmissionKraus[A1_,A2_,t_]:={{{1,0,0},{0,Exp[-(A1*t/2)],0},{0,0,Exp[-(A2*t/2)]}},{{0,Sqrt[1-Exp[-(A1*t)]],0},{0,0,0},{0,0,0}},{{0,0,Sqrt[1-Exp[-(A2*t)]]},{0,0,0},{0,0,0}}};


(* ::Subsection::Closed:: *)
(*Entropy*)


Log0[x_]:=If[x==0,0,Log[2,x]];
SetAttributes[Log0,Listable];


\[Eta][x_]:= -x Log[2,x];


\[Eta]2[x_]:= -x Log[2,x] - (1-x)Log[2,1-x];


QuantumEntropy[m_]:=Block[{eigvals,qe},eigvals=Chop[Eigenvalues[m]];
    qe=Sum[eigvals[[i]] Log0[eigvals[[i]]],{i,1,Length[eigvals]}];
    - Chop[qe]
];


QuantumChannelEntropy[ch_List]:=QuantumEntropy[Jamiolkowski[ch]];
QuantumChannelEntropy[fun_Function,dim_Integer]:=QuantumEntropy[Jamiolkowski[fun,dim]];


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


\[Delta][a_,type_:"Dirac"] := Switch[type,
    "Dirac", DiracDelta[a],
    "Indicator", DiscreteDelta[a]
]; 


VandermondeMatrix[l_]:=Table[Table[l[[j]]^i,{i,0,Length[l]-1}],{j,1,Length[l]}];


ProdSum[l_]:=Times@@Flatten[Table[Table[l[[i]] + l[[j]], {i, 1, j - 1}], {j, 2, Length[l]}]];


ProdDiff2[l_]:=Block[{x},Discriminant[Times@@Table[(x-l[[i]]),{i,1,Length[l]}],x]];


ProbBuresNorm[N_]:=2^(N^2-N) Gamma[N^2/2]/(\[Pi]^(N/2) Product[Gamma[j+1],{j,1,N}]);


ProbBures[l_,delta_:"Dirac"]:=FullSimplify[ProbBuresNorm[Length[l]] \[Delta][\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)-1,delta]/\[Sqrt](\!\(
\*SubsuperscriptBox[\(\[Product]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)) Det[VandermondeMatrix[l]]^2/ProdSum[l]];


ProbHSNorm[N_]:=Gamma[N^2]/Product[Gamma[N-j] Gamma[N-j+1],{j,0,N-1}];


ProbHS[l_,delta_:"Dirac"]:=ProbHSNorm[Length[l]]\[Delta][1-(Plus@@l),delta] Det[VandermondeMatrix[l]]^2;

RandomProductKet[n_?ListQ]:=Flatten[Fold[KroneckerProduct[#1,RandomKet[#2]]&,{1},n]];

RandomProductNumericalRange[A_,sys_,noPoints_:1]:=Block[{prod},
    Table[prod=RandomProductKet[sys];Tr[Proj[prod].A],{noPoints}]
];


RandomMaximallyEntangledNumericalRange[A_,noPoints_]:=Block[{ent,dim},
    dim=Dimensions[A][[1]];
    Table[ent=RandomEntangledUnitVector[dim];ent\[Conjugate].A.ent,{noPoints}]
];

NumericalRangeBound[A_?MatrixQ,step_:0.01]:=Block[
    {w,Ath,Hth,m,s,Kth,pKp,ee,rr,mm,sm,mM,sM,e,r},
    w={};
    Table[
    Ath=Exp[I*(-\[Theta])]*A;
    Hth=MatrixRe[Ath];
    {e,r}=Eigensystem[Hth];
    e=Re[e];
    m=Max[e];
    s=Position[e,m];
    If[
        Length[s]==1,(*then*)
        AppendTo[w,ArrayFlatten[Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose]]],
        (*else*)
        Kth=I*(Hth-Ath); 
        pKp=Extract[r,s]\[Conjugate].Kth.Extract[r,s]\[Transpose]; 
        {ee,rr}=Eigensystem[pKp]; 
        ee=Re[ee]; mm=Min[ee]; 
        sm=Position[ee,mm];
        AppendTo[w,ArrayFlatten[Extract[rr,sm]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sm]\[Transpose]]];
        mM=Max[ee];
        sM=Position[ee,mM];AppendTo[w,ArrayFlatten[Extract[rr,sM]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sM]\[Transpose]]]
    (*end if*)
    ]
    ,{\[Theta],0,2\[Pi],step}];
Flatten[w,2]
];

IntegrateSU2[f_,list__]:=Block[{variables,matrices,matricesCT,matricesC,replacement,range,func},
variables=Map[SymbolicVector[#,3]&,{list}];
matrices=Map[Unitary2Euler[0,ArcSin[Sqrt[#[[1]]]],#[[2]],#[[3]]]&, variables];
matricesCT=Map[Unitary2Euler[0,-ArcSin[Sqrt[#[[1]]]],-#[[2]],#[[3]]]&, variables];
matricesC=Map[Unitary2Euler[0,ArcSin[Sqrt[#[[1]]]],-#[[2]],-#[[3]]]&, variables];
replacement = Join[ Map[#[[1]]-> #[[2]]&,{{list},matrices}\[Transpose]], Map[#[[1]]\[ConjugateTranspose]-> #[[2]]&,{{list},matricesCT}\[Transpose]],Map[Conjugate[#[[1]]]-> #[[2]]&,{{list},matricesC}\[Transpose]]];
range=Flatten[Map[{{#[[1]],0,1},{#[[2]],0,2 * \[Pi]},{#[[3]],0,2 *\[Pi]}}&, variables],1];
func = Expand[f/.replacement];
1/(2*\[Pi])^(2 *Length[{list}])* Apply[Integrate[func,##]&,range]
];

Unitary3[al_,be_,ga_,th_,a_,b_,c_,ph_]:=MatrixExp[I *\[Lambda]3*al].MatrixExp[I*\[Lambda]2*be].
    MatrixExp[I*\[Lambda]3*ga].MatrixExp[I*\[Lambda]5*th].MatrixExp[I*\[Lambda]3*a].
    MatrixExp[I*\[Lambda]2*b].MatrixExp[I*\[Lambda]3*c].MatrixExp[I*\[Lambda]8*ph];

Unitary4Canonical[a1_,a2_,a3_]:=MatrixExp[I*a1*KroneckerProduct[sx,sx]+a2*I*KroneckerProduct[sy,sy]+a3*I*KroneckerProduct[sz,sz]];

Needs["FiniteFields`"]

MubElement[0, t_, p_, m_: 1] := UnitVector[p^m, t + 1];
fieldSqrt[\[Omega]_, k_, q_, p_, m_] := 
 Block[{qF, prod = 1, pow2n, qmod2n, km1},
  If[p == 2,
   qF = FromElementCode[GF[p, m], q];
   Product[
    If[qF =!= 0 && qF[[1]][[n + 1]] != 0,
     pow2n = FromElementCode[GF[p, m], 2^n];
     qmod2n = FromElementCode[GF[p, m], Mod[q, 2^n]];
     km1 = FromElementCode[GF[p, m], k - 1];
     I^(ToElementCode[
         ReduceElement[km1*pow2n*pow2n]])*\[Omega]^(ToElementCode[
         ReduceElement[km1*pow2n*qmod2n]])
     ,(*else*)
     1
     ]
    , {n, 0, m - 1}]
   ,(*else*)
   \[Omega]^(ToElementCode[
      ReduceElement[
       FromElementCode[GF[p, m], k - 1] FromElementCode[GF[p, m], 
         q] FromElementCode[GF[p, m], q]/
         FromElementCode[GF[p, m], 2]]])
   ]
  ]
MubElement[k_, j_, p_, m_: 1] := Block[{},
   \[Omega][kk_, pp_] := Exp[2  \[Pi] I/pp]^(kk);
   1/Sqrt[p^m] Plus @@
     Table[
      \[Omega][
        ToElementCode[
         ReduceElement[-FromElementCode[GF[p, m], q] FromElementCode[
            GF[p, m], j]]], p]*fieldSqrt[\[Omega][1, p], k, q, p, m]*
       MubElement[0, q, p, m],
      {q, 0, p^m - 1}
      ]
   ];
Mub[p_, m_: 1] := Block[{},
   Table[
    Table[
     MubElement[b, t, p, m], {t, 0, p^m - 1}
     ]
    , {b, 0, p^m}
    ]
   ];




(* ::Section:: *)
(*Package footer*)


Print["Package QIExtras ", QIExtras`Private`qiExtrasVersion, " (last modification: ", QIExtras`Private`qiExtrasLastModification, ")."];



End[] (* End Private Context *)

Protect@@Names["QIExtras`*"]

EndPackage[]
