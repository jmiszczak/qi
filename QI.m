(* ::Package:: *)

(* ::Section:: *)
(*Package header*)


(* File: QI.m *)
(* Description: Mathematica package for the analysis of quantum states and operations *)
(* Authors: Jaroslaw Miszczak <miszczak@iitis.pl>, Piotr Gawron <gawron@iitis.pl>, Zbigniew Puchala <z.puchala@iitis.pl> *)
(* License: GPLv3 *)


BeginPackage["QI`"];
Unprotect@@Names["QI`*"]
Clear@@Names["QI`*" ]


$PrePrint = If[SquareMatrixQ[#], MatrixForm[#], #]&;


(* ::Section:: *)
(*Help messages*)


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)

SquareMatrixQ::usage = "SquareMatrixQ[A] returns True only if A is a square matrix, and gives False otherwise.";

SymbolicMatrix::usage = "SymbolicMatrix[a,m,n] returns m\[Cross]n-matrix with elements a[i,j], i=1,...,m, j=1,...,n. If the third argument is ommited this function returns square m\[Cross]m matrix. This functions can save you some keystrokes and, thanks to TeXForm function, its results can be easily incorporated in LaTeX documents.";

SymbolicVector::usage = "SymbolicVector[a,n] is equivalent to Matrix[a,n,1] and it returns a vector with m elements a[i],i=1,...,n.";

SymbolicHermitianMatrix::usage = "SymbolicHermitianMatrix[sym,n] produces a n\[Cross]n Hermitian matrix. See also: SymbolicMatrix, SymbolicVector.";

SymbolicBistochasticMatrix::usage = "SymbolicBistochasticMatrix[sym, dim] produces symbolic bistochastic matrix size dim. See also: SymbolicMatrix, SymbolicVector."; 

ComplexToPoint::usage = "ComplexToPoint[z] returns a real and an imaginary parts of a complex number z as a pair of real numbers.";

MatrixSqrt::usage= "MatrixSqrt[A] returns square root for the matrix A.";

MatrixAbs::usage= "MatrixAbs[A] returns absolute value for matrix A defined as MatrixSqrt[A.A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)]. See also: MatrixSqrt.";

MatrixRe::usage = "Hermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A+A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)).";

MatrixIm::usage = "Antyhermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A-A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)).";

Proj::usage = "Proj[{\!\(\*SubscriptBox[\"v\", \"1\"]\),\!\(\*SubscriptBox[\"v\", \"2\"]\),...,\!\(\*SubscriptBox[\"v\", \"n\"]\)}] returns projectors for the vectors in the input list.";

(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)

Fidelity::usage = "Fidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns the quantum fidelity between states \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) calculated using a simplified formula as (\[Sum]\!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\)\!\(\*SuperscriptBox[\")\", \"2\"]\), where \!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\) are the eigenvalues of \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\).";

Superfidelity::usage = "Superfidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] calculates superfidelity between \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) defined as Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] + Sqrt[1-Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)]]Sqrt[1-Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)]]. See: J.A. Miszczak et al., Quantum Information & Computation, Vol.9 No.1&2 (2009)."; 

Subfidelity::usage = "Subfidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns subfidelity between states \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) See: J.A. Miszczak et al., Quantum Information & Computation, Vol.9 No.1&2 (2009).";

TraceNorm::usage = "TraceNorm[A] = \[Sum]\!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\), where \!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\) are the singular values of A. See also: TraceDistance.";

TraceDistance::usage = "TraceDistance[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns the trace distance between matrices \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\), which is defined as \!\(\*FractionBox[\"1\", \"2\"]\)tr|\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)-\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)|.";

GateFidelity::usage = "GateFidelity[U,V] is equivalent to 1/d tr|UV\[ConjugateTranspose]|.";

(* ::Subsection::Closed:: *)
(*Commonly used matrices*)

sx::usage = "Pauli matrix sx";
sy::usage = "Pauli matrix sy";
sz::usage = "Pauli matrix sz";
id::usage = "Identity matrix for one qubit. See also: IdentityMatrix.";
wh::usage = "Hadamard gate for one qubit. See also: QFT.";
cnot::usage = "Controlled not matrix for two qubits.";

(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


(*VectorSchmidtDecomposition::usage = "VectorSchmidtDecomposition[vec,dim] - \ 
Schmidt decomposition of the vector vec in dimd[[1]]\[Cross]dim[[2]]-dimensional \
Hilbert space.";
*)

(*
OperatorSchmidtDecomposition::usage = "OperatorSchmidtDecomposition[mtx,drows,dcols] - \ 
Schmidt decomposition of mtx in the Hilbert-Schmidt space of matrices of dimension \
(drows[[1]]\[Times]drows[[2]])\[Times](dcols[[1]]\[Times]dcols[[2]]) \
into list of 3-tuples: Schmidt number, matrix (drows[[1]]\[Times]dcols[[1]]), \   
matrix (drows[[2]]\[Times]dcols[[2]]).
W = OperatorSchmidtDecomposition[mtx, {r1, r2}, {c1, c2}];
mtx == Sum[W[[i]][[1]]*KroneckerProduct[W[[i]][[2]], W[[i]][[3]]], {i,Length[W]}]";
*)

SchmidtDecomposition::usage = "SchmidtDecomposition[e,dim] - accepts a vector \
or a matrix as a first argument and returns apropriate Schmidt decomposition. \
"; (*TODO: Change SchmidtDecomposition::usage*)


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec::usage = "Vec[A] - vectorization of the matrix A column by column. See also: Res.";

Unvec::usage = "Unvec[v,c] - de-vectorization of the vector into the matrix\
with c columns. If the second parameter is omitted then it is assumed that v\ 
can be mapped into square matrix. See also: Unres, Vec.";

Res::usage = "Res[A] is equivalent to Vec[Transpose[A]]. Reshaping maps\
matrix A  into a vector row by row. Note, that this is different then the\
reshape operation in Matlab or GNU Octave.";

Unres::usage = "Unres[v,c] - de-reshaping of the vector into a matrix with c columns. If the second parameter is omitted then it is assumed that v can be mapped into a square matrix. See also: Unvec, Res.";

Reshuffle::usage = "\
Reshuffle[\[Rho], drows, dcols] for a matrix of dimensions (drows[[1]]\[Times]drows[[2]])\[Times](dcols[[1]]\[Times]dcols[[2]]) returns reshuffled matrix with dimensions \
(drows[[1]]\[Times]dcols[[1]])\[Times](drows[[2]]\[Times]dcols[[2]]), \
drows and dcols parameters can be ommited for a square matrix.";

(* ::Subsection::Closed:: *)
(*Parametrizations*)

Unitary2::usage = "Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]] returns a parametrization of U(2).";

Unitary2Euler::usage = "Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]] returns the Euler parametrization of U(2).";

SpecialUnitary2::usage ="SpecialUnitary2[\[Beta],\[Gamma],\[Delta]] returns a parametrization of SU(2). This is equivalent to Unitary2[0,\[Beta],\[Gamma],\[Delta]].";

Unitary3::usage = "Unitary3[\[Alpha],\[Beta],\[Gamma],\[Tau],a,b,c,ph] returns the Euler parametrization of U(3).";

Unitary4Canonical::usage = "Parametrization of non-local unitary matrices for two qubits. See: B. Kraus, J.I. Cirac, Phys. Rev. A 63, 062309 (2001), quant-ph/0011050v1.";

StateVector::usage = "StateVector[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\),\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"n\", \"+\", \"1\"}]]\),...,\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"2\", \" \", \"n\"}]]\)}] returns pure n+1-dimensional pure state (ket vector) constructed form probability distribution parametrize by numbers {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)} and phases {\!\(\*SubscriptBox[\"\[Phi]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \"n\"]\)}. See also: SymbolicVector.";

(* ::Subsection::Closed:: *)
(*One-qubit states*)

QubitKet::usage = "QubitKet[\[Alpha],\[Beta]] parametrization of the pure state (as a state vector) for one qubit as (Cos[\[Alpha]] Exp[i\[Beta]], Sin[\[Alpha]]). This is equivalent to StateVector[{\[Alpha],\[Beta]}]. See also: QubitPureState, StateVector.";

QubitPureState::usage = "QubitPureState[\[Alpha],\[Beta]] - a parametrization of the pure state as a density matrix for one qubit. This is just a alias for Proj[QubitKet[\[Alpha],\[Beta]]]. See also: QubitKet.";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


ApplyKraus::usage = "ApplyKraus[ck,\[Rho]] - apply channel ck, given as a list of Kraus operators, to the input state \[Rho]. See also: ApplyUnitary, ApplyChannel.";

ChannelToMatrix::usage = "ChannelToMatrix[E,d] returns matrix representation of a channel E acting on d-dimensional state space. First argument should be a pure function E such that E[\[Rho]] transforms input state according to the channel definition.";

ApplyChannel::usage = "ApplayChannel[f,\[Rho]] - apply channel f, given as a pure function, to the input state \[Rho]. See also: ApplyUnitary, ApplyKraus."

Superoperator::usage = "Superoperator[kl] returns matrix representation of quantum channel given as a list of Kraus operators. Superoperator[fun,dim] is just am alternative name for ChannelToMatrix[fun,dim] and returns matrix representation of quantum channel, given as a pure function, acting on dim-dimensional space. So Superoperator[DepolarizingChannel[2,p,#]&,2] and Superoperator[QubitDepolarizingKraus[p]] returns the same matrix. See also: ChannelToMatrix.";

DynamicalMatrix::usage = "Dynamical matrix of quantum channel given as a list of Kraus operators (DynamicalMatrix[ch]) or as a function fun action on dim-dimensional space (DynamicalMatrix[fun,dim]). See also: Superoperator, ChannelToMatrix.";

Jamiolkowski::usage = "Jamiolkowski[K] gives the image of the Jamiolkowski isomorphism for the channel given as the list of Karus operators K. Jamiolkowski[fun,dim] gives the image of the Jamiolkowski isomorphism for the channel given as a function fun action on dim-dimensional space. See also: Superoperator, ChannelToMatrix, DynamicalMatrix.";

TPChannelQ::usage = "Performs some checks on Kraus operators. Use this if you want to check if they represent quantum channel.";

SuperoperatorToKraus::usage = "Finds Kraus operators for a given super operator";

ProductSuperoperator::usage = "ProductSuperoperator[\[CapitalPsi],\[CapitalPhi]] computes a product superoperator of superoperatos \[CapitalPsi] and \[CapitalPhi].";


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)

PartialTranspose::usage = "PartialTranspose[\[Rho],dim,sys] - Returns the partial transpose, according to systems sys, of density matrix \[Rho] composed of subsystems of dimensions dims.";

PartialTrace::usage = "PartialTrace[\[Rho],dim,sys] - Returns the partial trace, according to systems sys, of density matrix \[Rho] composed of subsystems of dimensions dim.";

(* ::Subsection::Closed:: *)
(*Entanglement*)


(* ::Subsection::Closed:: *)
(*Random states and operations*)


RandomSimplex::usage = "RandomSimplex[d] generates a point on a d-dimensional simplex according to the uniform distibution.";


RandomKet::usage = "RandomKet[d] - ...  random ket vector in d-dimensional space. d may be a list of integers in this case ket will be in product form d[[1]]\[CircleTimes]...\[CircleTimes]d[[k]]. See: T. Radtke, S. Fritzsche, Comp. Phys. Comm., Vol. 179, No. 9, p. 647-664.";

RandomDynamicalMatrix::usage = "RandomDynamicalMatrix[d,k] returns dynamical matrix of operation acting on d-dimensional states with k eigenvalues equal to 0. Thanks to Wojtek Bruzda. see Random Quantum Operations DOI[10.1016/j.physleta.2008.11.043]";

GinibreMatrix::usage = "GinibreMatrix[m,n] returns complex matrix of dimension m\[Cross]n with normal distribution of real and imaginary parts.";

RandomSpecialUnitary::usage = "Random special unitary matrix. See RandomUnitary";

RandomUnitary::usage = "Random unitary matrix using QR decomposition. F. Mezzadri, See: NOTICES of the AMS, Vol. 54 (2007), 592-604";

RandomOrthogonal::usage = "Random orthogonal matrix using QR decomposition. F. Mezzadri, See: NOTICES of the AMS, Vol. 54 (2007), 592-604"

RandomState::usage = "RandomState[d,dist] - random density matrix of dimension d. Argument dist can be ''HS'' (default value) or ''Bures'' or an integer K. ''HS'' gives uniform distribution with respect to the Hilbert-Schmidt measure. ''Bures'' gives random state distributed according to Bures measure. If dist is given as an integer K, the state is generated with respect to induced measure with an ancilla system od dimension K.";


(* ::Subsection::Closed:: *)
(*Bloch Representation*)


StateToBloch::usage = "StateToBloch[A] - for a square matrix A returns a vector of coefficients obtained from expansion on normed generalized Pauli matrices. See also: GeneralizedPauliMatrices.";

BlochToState::usage = "BlochToState[v] - returns a matrix of appropriate dimension from Bloch vector, i.e. coefficients treated as coefficients from expansion on normalized generalized Pauli matrices. See also: GeneralizedPauliMatrices.";

(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"];


qiAuthors = "Jaroslaw Miszczak <miszczak[at]iitis[dot]pl>, Piotr Gawron <gawron[at]iitis[dot]pl>, Zbigniew Puchala <z.puchala[at]iitis[dot]pl>";

qiLicense = "GPLv3 <http://www.gnu.org/licenses/gpl.html>";

qiHistory = {
	{"0.1.0", "05/06/2009", "Jarek", "Initial version"}, 
	{"0.1.1", "09/06/2009", "Zbyszek", "Fixed \[Eta] and \[Eta]2 functions, fixed problem with protected symbols."}, 
	{"0.1.2", "18/06/2009", "Gawron", "Added quantum channel parametrization for one qubit."}, 
	{"0.1.3", "03/07/2009", "Zbyszek", "Added alternative reshuffling."},
	{"0.1.4", "", "Jarek", "Changed default print output."}, 
	{"0.2.0", "09/07/2009", "Jarek and Zbyszek", "Documentation generator added."},
	{"0.2.1", "", "Jarek and Garwon", "Changed QubitGeneralState function."},
	{"0.2.2", "04/08/2009", "Zbyszek", "Added reshuffling permutation and product of superoperators."},
	{"0.2.3", "13/08/2009", "Jarek", "Minor update in documentation."},
	{"0.2.4", "14/09/2009", "Jarek and Zbyszek", "Fixed \[Eta] function."},
	{"0.2.5", "27/08/2009", "", "Fixed Werner state and added IsotropicState."},
	{"0.2.6", "06/10/2009", "Jarek", "Spelling improvements."},
	{"0.2.7", "02/11/2009", "Jarek", "Improved '\[CircleTimes]' usage, added ApplyUnitary."},
	{"0.2.8", "03/22.1009", "", "Fixed small problem with MaxEnt."},
	{"0.2.9", "", "Jarek", "Some code cleanups, fixed SchmidtDecomposition."},
	{"0.3.0", "06/11/2009", "", "Added OperatorSchmidtDecomposition."},
	{"0.3.1", "12/11/2009", "", "Added Concurrence4, fixed Cnot problem, some code cleanups."},
	{"0.3.2", "17/11/2009", "Jarek", "Documentation improvements, SquareMatrixQ predicate."},
	{"0.3.3", "18/11/2009", "Jarek", "All qi* functions and constants moved to the QI`Private context, added list of QI`* names."},
	{"0.3.4", "19/11/2009", "", "Added Negativity, fixed SchmidtDecomposition."},
	{"0.3.5", "21/11/2009", "Jarek", "SchmidtDecomposition now accepts vectos as well as matrices."},
	{"0.3.6", "24/11/2009", "Jarek", "Minor update in qiNames."},
	{"0.3.7", "04/12/2009", "", "OperatorSchmidtDecomposition fixed."},
	{"0.3.8", "04/01/2010", "", "Added local vars in RandomState."},
	{"0.3.9", "06/01/2010", "", "Added error message in Ket."},
	{"0.3.10", "19/01/2010", "Jarek", "Improved Davies map."},
	{"0.3.11", "26/01/2010", "Jarek", "Improved simplex generation algorithm, added some function for random vectors."},
	{"0.3.12", "04/03/2010", "Jarek", "Changed parameter in Swap gate."},
	{"0.3.13", "26/03/2010", "", "Fixed bug with state parametrization."},
	{"0.3.14", "25/05/2010", "Jarek", "Fixed numerical bug in Concurrence4 - Chop function added."},
	{"0.3.15", "07/06/2010", "", "RandomMaximallyEntangledNumericalRange added."},
	{"0.3.16", "20/06/2010", "", "Alternative version of Ketbra function."},
	{"0.3.17", "11/07/2010", "Jarek", "Name changed for Davies map."},
	{"0.3.18", "11/08/2010", "Jarek", "Fiexd bug in GeneralizedPauliKraus function reported by Fatih Ozaydin and one syntax error."},
	{"0.3.19", "13/09/2010", "Gawron", "Fixed bug in QubitDecayKraus and QubitDepolarizingKraus, QubitBitflipKraus, QubitPhaseflipKraus, QubitBitphaseflipKraus."},
	{"0.3.20", "04/10/2010", "Jarek", "Fixed inconsistency in QubitDepolarizingKraus and DepolarizingChannel."},
	{"0.3.21", "08/10/2010", "Jarek", "Added ReshufflePermutation2 and fixed Reshuffle2, qiHistory now stores commiter name."},
	{"0.3.22", "09/11/2010", "Gawron", "RandomUnitary -> RandomUniatryEuler, new RandomUnitary based on QR decomposition."},
	{"0.3.23", "09/11/2010", "Jarek", "Fixed package loading. New function RandomUnitaryQR and modified RandomUnitary. Improved usage messages."},
	{"0.3.24", "17/11/2010", "Zbyszek", "New function SuperoperatorToKraus."},
	{"0.3.25", "18/11/2010", "Zbyszek", "Renaming of Reshuffling functions."},
	{"0.3.26", "22/12/2010", "Jarek", "Two new functions: UpperBandOnes and UpperTriangularOnes."},
	{"0.3.27", "14/02/2011", "Jarek", "Removed PartialTrace from the RandomState function and added induced measures."},
	{"0.3.28", "16/03/2011", "Jarek", "Minor update in documentation and GateFidelity function added."},
	{"0.3.29", "31/03/2011", "Gawron", "List of names to protect is automatically generated now. RandomUnitaryEuler removed."},
	{"0.3.30", "01/04/2011", "Gawron", "New PartialTrace, all other PartialTrace* functions removed as obsolete."},
	{"0.3.31", "04/04/2011", "Gawron", "New PartialTranspose, all other PartialTranspose* functions removed as obsolete. Reshuffle and ReshufflePrim cleaned up."},
	{"0.3.32", "05/04/2011", "Gawron", "RandomSimplex changed."},
	{"0.3.33", "08/04/2011", "Gawron, Zbyszek", "*SchmidtDecomposition changed."},
	{"0.3.34", "29/04/2011", "Gawron, Zbyszek", "ProdSum fixed."},
	{"0.3.35", "11/05/2011", "Gawron, Zbyszek", "SchmidtDecomposition fixed, RandomKet enhenced."},
    {"0.3.36", "18/05/2011", "Zbyszek", "Added functions: Unitary2Euler, IntegrateSU2, RandomOrthogonal."},
    {"0.3.37", "07/07/2011", "Gawron, Jarek", "Added function: SymbolicBistochasticMatrtix."},
	{"0.3.38", "05/08/2011", "Zbyszek, Jarek", "Added HyperlinkToString and DOIToString functions."},
	{"0.4.0",  "13/10/2011", "Zbyszek, Jarek, Gawron", "Big changes"}
};

qiVersion = Last[qiHistory][[1]];

qiLastModification = Last[qiHistory][[2]];

qiAbout = "QI is a package of functions for Mathematica computer algebra system, which implements \
number of functions used in the analysis of quantum states and quantum operations. In contrast to \
many available packages for symbolic and numerical simulation of quantum computation presented \
package is focused on geometrical aspects of quantum information theory.";

(* ::Subsection:: *)
(*Miscellaneous functions*)

HyperlinkToString::usage = "HyperlinkToString[text,link] creates a link labeled with text to the given URL and returns it as a Mathematica string.";

HyperlinkToString[text_,link_]:="\!\(\*ButtonBox[StyleBox[\""<>text<>"\", \"SR\"],Active->True,BaseStyle->\"Link\",ButtonData->\""<>link<>"\"]\)";

DOIToString::usage = "DOIToString[text,doi] creates a link labeled with text to the given DOI and returns it as a Mathematica string.";

DOIToString[text_,doi_]:="\!\(\*ButtonBox[StyleBox[\""<>text<>"\", \"SR\"],Active->True,BaseStyle->\"Link\",ButtonData->\"http://dx.doi.org/"<>doi<>"\"]\)";

(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


(*x_?MatrixQ \[CircleTimes] y_?MatrixQ := KroneckerProduct[x,y];*)
(*x_?VectorQ \[CircleTimes] y_?VectorQ := Flatten[KroneckerProduct[x,y]];*)
(*x_ \[CircleTimes] y_ \[CircleTimes] z_ := (x \[CircleTimes] y) \[CircleTimes] z;*)

CircleTimes[x_?MatrixQ,y_?MatrixQ] := KroneckerProduct[x,y];
CircleTimes[x_?VectorQ,y_?VectorQ] := Flatten[KroneckerProduct[x,y]];
CircleTimes[a_,b__] := CircleTimes[a, CircleTimes[b]]

SquareMatrixQ[A_]:= Block[{dims=Dimensions[A]},
	(Length[dims]==2 )&&(dims[[1]]==dims[[2]])
];

SymbolicMatrix[sym_,d1_,d2_:0] := Which[
	d2==0, Table[Subscript[sym, i,j], {i,1,d1},{j,1,d1}],
	d2==1, Table[Subscript[sym, i], {i,1,d1}],
	True, Table[Subscript[sym, i,j], {i,1,d1},{j,1,d2}] 
];

SymbolicVector[sym_,d1_]:= SymbolicMatrix[sym,d1,1];

SymbolicHermitianMatrix[sym_,d_]:=Block[{mtx},
	mtx=Table[0,{d},{d}];
	Table[mtx[[i,j]]=Subscript[sym, i,j],{i,1,d},{j,1,i}];
	mtx=mtx+mtx\[ConjugateTranspose];
	Table[mtx[[i,i]]=Subscript[sym, i,i],{i,1,d}];
	mtx
];

SymbolicBistochasticMatrix[sym_, dim_] := Block[{f},
  f = Map[Append[#, 1 - Apply[Plus, #]] &, #] &;
  f[f[SymbolicMatrix[sym, dim - 1]]\[Transpose]]\[Transpose]
];

ComplexToPoint[z_]:={Re[z],Im[z]};
SetAttributes[ComplexToPoint,Listable];

MatrixRe[A_?SquareMatrixQ]:=(A+A\[ConjugateTranspose])/2;

MatrixIm[A_?SquareMatrixQ]:=(A-A\[ConjugateTranspose])/2;

Proj[v_]:=Table[v[[i]]Conjugate[v[[j]]],{i,1,Length[v]},{j,1,Length[v]}];

(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)

MatrixSqrt[m_?SquareMatrixQ]:=MatrixPower[m,1/2];

MatrixAbs[a_?SquareMatrixQ]:=MatrixSqrt[a.(a\[ConjugateTranspose])];

Fidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=(Plus@@(Sqrt[Eigenvalues[a.b]]))^2;

Superfidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=Tr[a.b]+Sqrt[(1-Tr[a.a])]*Sqrt[(1-Tr[b.b])];

Subfidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=Block[{prod = a.b}, Tr[prod] +Sqrt[2]Sqrt[(Tr[prod]*Tr[prod]-Tr[prod.prod])]];

TraceNorm[a_?SquareMatrixQ]:=Plus@@SingularValueList[a];

TraceDistance[a_?SquareMatrixQ,b_?SquareMatrixQ]:=1/2*TraceNorm[a-b];

GateFidelity[mU_?SquareMatrixQ,mV_?SquareMatrixQ]:=Block[{dimU=Dimensions[mU][[1]]},
	If[dimU==Dimensions[mV][[1]],
		1/dimU Abs[Tr[mU.ConjugateTranspose[mV]]],
		Message[GateFidelity::argerr]
	]
];
GateFidelity::argerr = "Both matrices have to be of the same dimension.";


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)

sx = {{0,1},{1,0}};
sy = {{0,-I},{I,0}};
sz = {{1,0},{0,-1}};
id = {{1,0},{0,1}};
wh = {{1,1},{1,-1}};
cnot = {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


VectorSchmidtDecomposition[vec_?VectorQ, d_?ListQ] := 
  Block[{mtx, u, w, v, vals, snum = Min[d[[1]], d[[2]]]},
   	mtx = Partition[vec, d[[2]]];
   	{u, w, v} = SingularValueDecomposition[mtx];
   	vals = Select[Diagonal[w], # != 0 &];
   	snum = Length[vals];
     {vals, u\[Transpose][[1 ;; snum]], 
     v\[Transpose][[1 ;; snum]]}\[Transpose]
   ];

VectorSchmidtDecomposition[vec_?VectorQ] :=
Block[{sqrtDim},
  sqrtDim = Sqrt[Dimensions[vec][[1]]];
  If[IntegerQ[sqrtDim],
   VectorSchmidtDecomposition[vec, {sqrtDim, sqrtDim}]
   ]
  ]
  
OperatorSchmidtDecomposition[op_, {n_?VectorQ, m_?VectorQ}] := 
  Block[{mtx, u, w, v, snum = Min[n[[1]]*m[[1]], n[[2]]*m[[2]]]},
   	mtx = Reshuffle[op, n, m];
   	{u, w, v} = SingularValueDecomposition[mtx];	
   Table[ {w[[i, i]], Partition[u\[Transpose][[i]], m[[1]]], 
     Partition[(v\[Transpose])[[i]]\[Conjugate], m[[2]]]}, {i, 1, 
     snum}]	
   ];
 
OperatorSchmidtDecomposition[op_?SquareMatrixQ, n_?VectorQ] := 
 OperatorSchmidtDecomposition[op, {n, n}]

OperatorSchmidtDecomposition[op_?SquareMatrixQ] := Block[{sqrtDim},
  sqrtDim = Sqrt[Dimensions[op][[1]]];
  If[IntegerQ[sqrtDim],
   OperatorSchmidtDecomposition[
    op, {{sqrtDim, sqrtDim} , {sqrtDim, sqrtDim}}]
   ]
  ]

SchmidtDecomposition[e_,dim___]:=Which[
	MatrixQ[e], OperatorSchmidtDecomposition[e,dim],
	VectorQ[e], VectorSchmidtDecomposition[e,dim],
	True, Message[SchmidtDecomposition::argerr]
];
SchmidtDecomposition::argerr = "First argument should be a vector or a matrix.";


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec[m_]:=Flatten[Transpose[m]]; 

Unvec[v_List,cols_:0]:=Which[
 (cols== 0)&&IntegerQ[\[Sqrt]Length[v]],Transpose[Partition[v,\[Sqrt]Length[v]]],
 Mod[Length[v],cols]==0,Transpose[Partition[v,cols]]
];

Res[m_List]:=Flatten[m]; 

Unres[v_List,cols_:0]:=Which[
 (cols== 0)&&IntegerQ[\[Sqrt]Length[v]],Partition[v,\[Sqrt]Length[v]],
 Mod[Length[v],cols]==0,Partition[v,cols]
];

Reshuffle[\[Rho]_]:=Block[{dim},
	dim = Sqrt[Length[\[Rho]]];
	If[And [SquareMatrixQ[\[Rho]] , IntegerQ[dim]] ,
		Reshuffle[\[Rho], {{dim,dim},{dim,dim}}]
	(*else*),
		Message[Reshuffle::argerr]
	]
]
Reshuffle::argerr = "Reshuffle works only for square matrices of dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\)\[Times]\!\(\*SuperscriptBox[\"d\", \"2\"]\), where d is an Integer, for other dimensions use ReshuffleGeneral";

Reshuffle[A_,{n_,m_}]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n[[2]]+i1,1+i2;;m[[2]]+i2]],{i1,0,n[[1]] n[[2]]-1,n[[2]]},{i2,0,m[[1]]*m[[2]]-1,m[[2]]}]
,1];


(* ::Subsection::Closed:: *)
(*Parametrizations*)

Unitary2[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_]:=
Exp[I \[Alpha]] DiagonalMatrix[{Exp[-I \[Beta]/2],Exp[I \[Beta]/2]}].{{Cos[\[Gamma]/2],-Sin[\[Gamma]/2]},{Sin[\[Gamma]/2],Cos[\[Gamma]/2]}}.DiagonalMatrix[{Exp[-I \[Delta]/2],Exp[I \[Delta]/2]}];

Unitary2Euler[\[Alpha]_,\[Theta]_,\[Phi]_,\[Psi]_]:={{E^(I*\[Alpha] + I*\[Phi])*Cos[\[Theta]], E^(I*\[Alpha] + I*\[Psi])*Sin[\[Theta]]}, {-(E^(I*\[Alpha] - I*\[Psi])*Sin[\[Theta]]), E^(I*\[Alpha] - I*\[Phi])*Cos[\[Theta]]}}

SpecialUnitary2[\[Beta]_,\[Gamma]_,\[Delta]_]:=Unitary2[0,\[Beta],\[Gamma],\[Delta]];

Unitary3[al_,be_,ga_,th_,a_,b_,c_,ph_]:=MatrixExp[I *\[Lambda]3*al].MatrixExp[I*\[Lambda]2*be].
	MatrixExp[I*\[Lambda]3*ga].MatrixExp[I*\[Lambda]5*th].MatrixExp[I*\[Lambda]3*a].
	MatrixExp[I*\[Lambda]2*b].MatrixExp[I*\[Lambda]3*c].MatrixExp[I*\[Lambda]8*ph];

Unitary4Canonical[a1_,a2_,a3_]:=MatrixExp[I*a1*KroneckerProduct[\[Sigma]x,\[Sigma]x]+a2*I*KroneckerProduct[\[Sigma]y,\[Sigma]y]+a3*I*KroneckerProduct[\[Sigma]z,\[Sigma]z]];

StateVector[l_]:=Block[{pr,ph,N,ProbablityVector},
	ProbablityVector[li_]:=Block[{ll,Ni},
		Ni=Length[li]+2;
        ll=Prepend[li,\[Pi]/2];
        Table[Sin[ll[[i-1]]]^2*Product[Cos[ll[[j-1]]]^2,{j,i+1,Ni}],{i,2,Ni}]
    ];
	N=Length[l]/2;
	pr=ProbablityVector[l[[1;;N]]];
	ph=Prepend[Exp[I*l[[N+1;;2*N]]],1];
	FullSimplify[Sqrt[pr]*ph, Assumptions -> Table[0<l[[i]]<\[Pi]/2,{i,1,N}]]
];


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet[\[Alpha]_,\[Beta]_]:={Cos[\[Alpha]], Exp[I*\[Beta]]*Sin[\[Alpha]]};

QubitPureState[\[Alpha]_,\[Beta]_]:=Proj[QubitKet[\[Alpha],\[Beta]]];


(* ::Subsection::Closed:: *)
(*Quantum channels*)


ApplyKraus[ch_,\[Rho]_]:=Sum[ch[[k]].\[Rho].(ch[[k]]\[ConjugateTranspose]),{k,1,Length[ch]}];

ApplyChannel[f_,\[Rho]_] := Map[f,\[Rho],{0}];

ChannelToMatrix[fun_,dim_] := Map[Res[fun[#]] &, BaseMatrices[dim]];       

Superoperator[ch_List] := Sum[ch[[i]]\[CircleTimes]ch[[i]]\[Conjugate],{i,1,Length[ch]}];

Superoperator[fun_,dim_] := ChannelToMatrix[fun,dim];

DynamicalMatrix[ch_List] := Reshuffle[Superoperator[ch]];

DynamicalMatrix[fun_Function,dim_Integer] := Reshuffle[Superoperator[fun,dim]];

Jamiolkowski[ch_List] := 1/Length[ch[[1]]]*DynamicalMatrix[ch];

Jamiolkowski[fun_Function,dim_Integer] := 1/dim*DynamicalMatrix[fun,dim];

TPChannelQ[operators_] := Sum[operators[[i]]\[ConjugateTranspose].operators[[i]],{i,Length[operators]}] == IdentityMatrix[Length[operators[[1]]]];

SuperoperatorToKraus[m_]:=Block[{val,vec}, {val,vec} = Eigensystem[Reshuffle[m]]; Sqrt[val] (Unres[#]&/@vec)];

ProductSuperoperator[m1_,m2_]:=Block[{dim1=Length[m1],dim2=Length[m2],perm},
    perm=ReshufflePermutation[Sqrt[dim1],Sqrt[dim2]];
    perm.(m1\[CircleTimes]m2).perm\[Transpose]
];

(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


ListReshape[list_, shape_] := 
  FlattenAt[Fold[Partition[#1, #2] &, Flatten[list], Reverse[shape]], 
   1];

PartialTrace[\[Rho]_,dim_?ListQ,sys_?ListQ] := Block[
	{offset, keep, dispose, keepdim, disposedim, perm1, perm2, perm, tensor},
	offset=Length[dim];
	keep=Complement[Range[offset], sys];
	dispose=Union[sys];
	perm1=Join[dispose,keep];
	perm2=perm1+offset;
	perm=Ordering[Join[perm1,perm2]];
	tensor=ListReshape[\[Rho], Join[dim,dim]];
	keepdim=Apply[Times, Join[dim, dim][[keep]]];
	disposedim=Apply[Times, Join[dim, dim][[dispose]]];
	tensor=Transpose[tensor,perm];
	tensor=ListReshape[tensor,{disposedim,keepdim,disposedim,keepdim}];
	Sum[tensor[[i,All,i,All]],{i,1,disposedim}]
];

PartialTranspose[\[Rho]_,dim_?ListQ,sys_?ListQ]:=Block[{offset,tensor,perm,idx1,idx2,s,targetsys},
	offset=Length[dim];
	tensor=ListReshape[\[Rho], Join[dim,dim]];
	targetsys=Union[sys];
	perm=Range[offset*2];
	For[s=1, s<=Length[targetsys], s+=1, 
		idx1 = Position[perm, targetsys[[s]]][[1, 1]];
		idx2 = Position[perm, targetsys[[s]] + offset][[1, 1]];
		{perm[[idx1]],perm[[idx2]]}={perm[[idx2]],perm[[idx1]]};
	];
	tensor=Transpose[tensor,InversePermutation[perm]];
	ListReshape[tensor,Dimensions[\[Rho]]]
]


(* ::Subsection::Closed:: *)
(*Entanglement*)


(* ::Subsection::Closed:: *)
(*Random states and operations*)


RandomSimplex[d_]:=Block[{r,r1,r2},
	r=Sort[Table[RandomReal[{0,1}],{i,1,d-1}]];
	r1=Append[r,1];r2=Prepend[r,0];r1-r2
];


RandomKet[n_?IntegerQ]:=Block[{p,ph},
	p=Sqrt[RandomSimplex[n]];
	ph=Exp[I*RandomReal[{0,2\[Pi]},n-1]];
	ph=Prepend[ph,1];
	p*ph
];

RandomKet[{d1_?IntegerQ, d2_?IntegerQ}, a_: "Sep"] := 
  Block[{l, v, U, V, w},
   Switch[a,
    "MaxEnt", 
            l = ConstantArray[1/Sqrt[Min[d1, d2]], Min[d1, d2]];,
    "Sep", 
            l = {1};,
    _List, 
            If[Total[a] == 1 && Length[a] <= Min[d1, d2], l = a, Message[RandomKet::argerr];],
    _, 
            Message[RandomKet::argerr];
    ];
   v = Sum[Sqrt[l[[i]]]*(UnitVector[d1, i]\[CircleTimes]UnitVector[d2, i]), {i, 1,Length[l]}];
   U = RandomUnitary[d1];
   V = RandomUnitary[d2];
   (U\[CircleTimes]V).v
   ];
RandomKet::argerr = "Error"; (*TODO: Write usage for RandomKet*)

RandomDynamicalMatrix[n_,m_:0]:=Block[{X,Y,sY},	
	X=GinibreMatrix[n^2,n^2-m];
	Y=PartialTrace[X.X\[ConjugateTranspose],{n,n},{1}];
	sY=MatrixPower[Y,-1/2];
	KroneckerProduct[IdentityMatrix[n],sY].X.X\[ConjugateTranspose].KroneckerProduct[IdentityMatrix[n],sY]
];

GinibreMatrix[m_,n_]:=RandomReal[NormalDistribution[0,1],{m,n}] + I RandomReal[NormalDistribution[0,1],{m,n}];


RandomSpecialUnitary[dim_]:=Module[{U},
	U=RandomUnitary[dim];
	U/Det[U]
];

RandomUnitary[dim_]:=Module[{q,r,d,ph},
	{q,r}=QRDecomposition[GinibreMatrix[dim,dim]];
	d=Diagonal[r];
	ph=d/Abs[d];
	Transpose[Transpose[q]*ph]
];

RandomOrthogonal[dim_]:=Module[{q,r,d,ph},
    {q,r}=QRDecomposition[RandomReal[NormalDistribution[0,1],{dim,dim}]];
	d=Diagonal[r];
	ph=d/Abs[d];
	Transpose[Transpose[q]*ph]
];

RandomState[d_,dist_:"HS"]:=Block[{A,U},
	Switch[dist,
		"HS",
			A=GinibreMatrix[d,d];
			A=(A.ConjugateTranspose[A]);
			A=Chop[A/Tr[A]],
		"Bures",
			A=GinibreMatrix[d,d];
			U=RandomUnitary[d];
			A=(IdentityMatrix[d]+U).A.A\[ConjugateTranspose].(IdentityMatrix[d]+U)\[ConjugateTranspose];
			Chop[A]\[ConjugateTranspose]/Tr[A],		
		_, 
			If[IntegerQ[dist] && dist >=d,
				A=GinibreMatrix[d,dist];
				A=(A.ConjugateTranspose[A]);
				A=Chop[A/Tr[A]],
				Message[RandomState::argerr,dist]
			]
	]

];
RandomState::argerr = "The second argument should be \"HS\" or \"Bures\" or an integer K>2, mesure \"`1`\" not implemented yet.";

(* ::Subsection::Closed:: *)
(*Bloch Representation*)


StateToBloch[A_]:=Block[{dim},
  dim=Length[A]; 1/Sqrt[2](Tr[A\[ConjugateTranspose].#]&/@GeneralizedPauliMatrices[dim])
];

BlochToState[vec_]:=Block[{dim},
	If[IntegerQ[Sqrt[Length[vec]+1]],
		dim= Sqrt[Length[vec]+1];
		1/dim IdentityMatrix[dim] + vec.GeneralizedPauliMatrices[dim]/Sqrt[2],
		Message[StateFromBlochVector::argerr, vec];
		Beep[];
	]
];
BlochToState::argerr= "Given vector (`1`) is not a Bloch vector of any dimension.";


(* ::Section::Closed:: *)
(*Package footer*)

Print["Package QI version ", QI`Private`qiVersion, " (last modification: ", QI`Private`qiLastModification, ")."];

End[];

Protect@@Names["QI`*"]

EndPackage[];
   