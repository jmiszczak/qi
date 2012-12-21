(* ::Package:: *)

(* ::Section::Closed:: *)
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
(*Commonly used matrices*)


sx::usage = "<v>sx</v> - Pauli matrix <v>sx</v>.";
sy::usage = "<v>sy</v> - Pauli matrix <v>sy</v>.";
sz::usage = "<v>sz</v> - Pauli matrix <v>sz</v>.";
id::usage = "<v>id</v> - Identity matrix for one qubit. See also: IdentityMatrix.";
wh::usage = "<v>wh</v> - Hadamard gate for one qubit.";
cnot::usage = "<v>cnot</v> - Controlled not matrix for two qubits.";
Id::usage = "<f>Id</f>[<v>n</v>] returns an identity matrix of dimension <v>n</v>. This is equivalent to IdentityMatrix[<v>n</v>].";


BaseVectors::usage = "<f>BaseVectors</f>[<v>n</v>] returns a list with the canonical basis in </v>n</v>-dimensional Hilbert space. See also: BaseMatrices.";


BaseMatrices::usage = "<f>BaseMatrices</f>[<v>n</v>] returns a list with the canonical basis in <v>n</v>\[Cross]</v>n</v>-dimensional Hilbert-Schmidt space of matrices. See also: BaseVectors.";


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


SquareMatrixQ::usage = "<f>SquareMatrixQ</f>[<v>A</v>] returns True only if <v>A</v> is a square matrix, and gives False otherwise.";

SymbolicMatrix::usage = "<f>SymbolicMatrix</f>[<v>a,m,n</v>] returns <s>m\[Cross]n</s>-matrix with elements <v>a[i,j], i=1,...,m, j=1,...,n</v>.\
If the third argument is ommited this function returns square <s>m\[Cross]m</s> matrix. This functions can save you some keystrokes and, thanks \
to TeXForm function, its results can be easily incorporated in LaTeX documents.";

SymbolicVector::usage = "<f>SymbolicVector</f>[<v>a,n</v>] returns a vector with <v>n</v> elements <v>a[i],i=1,...,n</v>.";

SymbolicHermitianMatrix::usage = "<f>SymbolicHermitianMatrix</f>[<v>sym,n</v>] produces a <s>n\[Cross]n</s> Hermitian matrix. \
See also: <f>SymbolicMatrix, SymbolicVector</f>.";

SymbolicBistochasticMatrix::usage = "<f>SymbolicBistochasticMatrix</f>[<v>sym, dim</v>] produces symbolic bistochastic matrix size <v>dim</v>. \
See also: <f>SymbolicMatrix, SymbolicVector</f>."; 

ComplexToPoint::usage = "<f>ComplexToPoint</f>[<v>z</v>] returns a real and an imaginary parts of a complex number <v>z</v> as a pair of real numbers.";

MatrixSqrt::usage= "<f>MatrixSqrt</f>[<v>A</v>] returns square root for the matrix <v>A</v>.";

MatrixAbs::usage= "<f>MatrixAbs</f>[<v>A</v>] returns absolute value for matrix <v>A</v> defined as <f>MatrixSqrt</f>[<v>A.A</v><s>\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)</s>]. \
See also: <f>MatrixSqrt</f>.";

MatrixRe::usage = "<f>MatrixRe</f>[<v>A</v>] returns a hermitian part of the matrix <v>A</v> i.e. <v>1/2(A+A</v><s>\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)</s><v>)</v>.";

MatrixIm::usage = "<f>MatrixIm</f>[<v>A</v>] returns an antyhermitian part of the matrix <v>A</v> i.e. <v>1/2(A-A</v><s>\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)</s><v>)</v>.";

Proj::usage = "<f>Proj</f>[<v>v</v>] returns projector of the vector <v>v</v>.";



(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec::usage = "<f>Vec</f>[<v>A</v>] - vectorization of the matrix <v>A</v> column by column. See also: <f>Res</f>.";

Unvec::usage = "<f>Unvec</f>[<v>v,c</v>] - de-vectorization of the vector into the matrix \
with <v>c</v> columns. If the second parameter is omitted then it is assumed that <v>v</v> \
can be mapped into square matrix. See also: <f>Unres, Vec</f>.";

Res::usage = "<f>Res</f>[<v>A</v>] is equivalent to <f>Vec</f>[Transpose[<v>A</v>]]. Reshaping maps\
matrix <v>A</v> into a vector row by row. Note, that this is different then the\
reshape operation in <v>Matlab</v> or <v>GNU Octave</v>.";

Unres::usage = "<f>Unres</f>[<v>v</v>,<v>c</v>] - de-reshaping of the vector into a matrix with <v>c</v> columns. \
If the second parameter is omitted, then it is assumed that <v>v</v> can be mapped into a square matrix. See also: <f>Unvec</f>, <f>Res</f>.";

Reshuffle::usage = "\
<f>Reshuffle</f>[<s>\[Rho]</s>, {<v>drows</v>, <v>dcols</v>}] for a matrix of dimensions <s>(drows[[1]]\[Times]drows[[2]])\[Times](dcols[[1]]\[Times]dcols[[2]])</s> \
returns a reshuffled matrix with dimensions <s>(drows[[1]]\[Times]dcols[[1]])\[Times](drows[[2]]\[Times]dcols[[2]])</s>.
Parameters {<v>drows</v>,<v>dcols</v>} can be ommited for a square matrix of dimension <s>n^2\[Times]n^2</s>.";


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


SchmidtDecomposition::usage = "<f>SchmidtDecomposition</f>[<v>x</v>,<v>dim</v>] - accepts a vector \
or a matrix as a first argument and returns apropriate Schmidt decomposition. The second argument \
is optional and specifies the dimensions of subsystems.
If <v>x</v> is a vector
\t <f>SchmidtDecomposition</f>[<v>x</v>] assumes that <v>x</v> is (<v>n^2</v>)-dimensional,
\t <f>SchmidtDecomposition</f>[<v>x</v>,{<v>n</v>,<v>m</v>}] assumes that the vector is a <v>n \[Times]m</v>-dimensional.
If <v>x</v> is a matrix, this function can be used in three different ways.
\t <f>SchmidtDecomposition</f>[<v>x</v>] assumes that <v>x</v> is (<v>n^2\[Times]n^2</v>)-dimensional, 
\t <f>SchmidtDecomposition</f>[<v>x</v>,{<v>n</v>,<v>m</v>}] assume that the matrix is a <v>n m\[Times]n m</v> square matrix and
\t <f>SchmidtDecomposition</f>[<v>x</v>,{{<v>r1</v>,<v>r2</v>},{<v>c1</v>,<v>c2</v>}}]
For example, for a matrix <v>mtx</v> of dimension <v>r1 r2\[Times] c1 c2</v> one can obtain a Schmidt \
decomposition on <v>r1 c1\[CircleTimes] r2 c2</v> system as 
\t <v>sd</v> = <f>SchmidtDecomposition</f>[<v>mtx</v>, {{<v>r1</v>, <v>r2</v>}, {<v>c1</v>, <v>c2</v>}}];
and reconstruct the original matrix as
\t <v>mtx</v> == Sum[<v>sd</v>[[<v>i</v>,1]]*KroneckerProduct[<v>sd</v>[[<v>i</v>,2]], <v>sd</v>[[<v>i</v>,3]]], {<v>i</v>, Length[<v>sd</v>]}];"; 


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTranspose::usage = "<f>PartialTranspose</f>[<s>\[Rho]</s>,<v>dim</v>,<v>sys</v>] - \
Returns the partial transpose, according to systems <v>sys</v>, of density matrix <s>\[Rho]</s> composed of subsystems of dimensions <v>dim</v>.";

PartialTrace::usage = "<f>PartialTrace</f>[<s>\[Rho]</s>,<v>sys</v>] - \
Returns the partial trace of an operator <s>\[Rho]</s> acting on a bipartite (<s>d\[Times]d</s>)-dimensional system, \
assuming that matrix <s>\[Rho]</s> is (<s>\!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\)</s>) dimensional. \
In this case the system specification can be given as a integer 1 or 2, by the list of integers consisting only 1 or 2 or by the empty list. 
<f>PartialTrace</f>[<s>\[Rho]</s>,<v>dims</v>,<v>sys</v>] - \
Returns the partial trace of an operator <s>\[Rho]</s> acting on a composite system with subsystem dimensions given in the list \
<v>dims</v>. List <v>sys</v> specifies systems to be discarded.";


(* ::Subsection::Closed:: *)
(*Distance measures*)


Fidelity::usage = "<f>Fidelity</f>[<s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)</s>,<s>\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s>] returns the quantum fidelity \
between states <s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)</s> and <s>\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s> calculated using a simplified formula as \
<s>(\[Sum]\!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\)\!\(\*SuperscriptBox[\")\", \"2\"]\)</s>, where <s>\!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\)</s> are the \
eigenvalues of <s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s>.";

Superfidelity::usage = "<f>Superfidelity</f><s>[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)]</s> calculates superfidelity \
between <s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)</s> and <s>\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s> defined as \
<s>tr[\!\(\*SubscriptBox[\(\[Rho]\), \(1\)]\)\!\(\*SubscriptBox[\(\[Rho]\), \(2\)]\)]+\!\(\*SqrtBox[\(\((1 - tr[\*SubsuperscriptBox[\(\[Rho]\), \(1\), \(2\)]])\) \((1 - tr[\*SubsuperscriptBox[\(\[Rho]\), \(2\), \(2\)]])\)\)]\)</s>
See: J.A. Miszczak et al., Quantum Information \& Computation, Vol.9 No.1\&2 (2009)."; 

Subfidelity::usage = "<f>Subfidelity</f><s>[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)]</s> returns subfidelity between \
states <s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)</s> and <s>\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s> \
See: J.A. Miszczak et al., Quantum Information \& Computation, Vol.9 No.1\&2 (2009).";

TraceNorm::usage = "<f>TraceNorm</f>[<v>A</v>] = <s>\[Sum]\!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\)</s>, where <s>\!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\)</s> are the singular values of <v>A</v>. \
See also: <f>TraceDistance</f>.";

TraceDistance::usage = "<f>TraceDistance</f>[<s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s>] \
returns the trace distance between matrices <s>\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)</s> and <s>\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)</s>, \
which is defined as <s>\!\(\*FractionBox[\"1\", \"2\"]\)tr|\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)-\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)|</s>.";

GateFidelity::usage = "<f>GateFidelity</f>[<v>U,V</v>] is equivalent to <v>1/d</v><s>tr|UV\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)</s><v>|</v>.";



(* ::Subsection::Closed:: *)
(*Bloch Representation*)


Lambda1::usage = "<f>Lambda1</f>[<v>i,j,n</v>] generalized Pauli matrix. For example <f>Lambda1</f>[<v>1,2,2</v>] is equal to Pauli <s>\[Sigma]x</s>. See also: <f>GeneralizedPauliMatrices</f>.";

Lambda2::usage = "<f>Lambda2</f>[<v>i,j,n</v>] generalized Pauli matrix. For example <f>Lambda2</f>[<v>1,2,2</v>] is equal to Pauli <s>\[Sigma]y</s>. See also: <f>GeneralizedPauliMatrices</f>.";

Lambda3::usage ="<f>Lambda3</f>[<v>i,n</v>] generalized Pauli matrix. For example <f>Lambda3</f>[<v>2,2</v>] is equal to Pauli <s>\[Sigma]z</s>. See also: <f>GeneralizedPauliMatrices</f>.";

GeneralizedPauliMatrices::usage = "<f>GeneralizedPauliMatrices</f>[<v>n</v>] returns list of generalized Pauli matrices for \
<v>SU(n)</v>. For <v>n=2</v> these are just Pauli matrices and for <v>n=3</v> - Gell-Mann matrices. \
Note that identity matrix is not included in the list.";  (*See also: PauliMatrices, GellMannMatrices, \[Lambda], Lambda1, Lambda2, Lambda3.";*)

StateToBloch::usage = "<f>StateToBloch</f>[<v>A</v>] - for a square matrix <v>A</v> returns a vector of coefficients obtained from expansion on \
normed generalized Pauli matrices. See also: <f>GeneralizedPauliMatrices</f>.";

BlochToState::usage = "<f>BlochToState</f>[<v>v</v>] - \
returns a matrix of appropriate dimension from Bloch vector, i.e. coefficients treated as coefficients from \
expansion on normalized generalized Pauli matrices. See also: <f>GeneralizedPauliMatrices</f>.";


(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2::usage = "<f>Unitary2</f>[<s>\[Alpha]</s>,<s>\[Beta]</s>,<s>\[Gamma]</s>,<s>\[Delta]</s>] returns a parametrization of <v>U</v>(2).";

Unitary2Euler::usage = "<f>Unitary2</f>[<s>\[Alpha]</s>,<s>\[Beta]</s>,<s>\[Gamma]</s>,<s>\[Delta]</s>] returns the Euler parametrization of <v>U</v>(2).";

SpecialUnitary2::usage ="<f>SpecialUnitary2</f>[<s>\[Beta],\[Gamma],\[Delta]</s>] returns a parametrization of <v>SU</v>(2). \
This is equivalent to <f>Unitary2</f>[<s>0,\[Beta],\[Gamma],\[Delta]</s>].";

SpecialUnitary::usage = "<f>SpecialUnitary</f>[<v>d</v>,<v>params</v>] returns the special unitary matrix of size <v>d</v> with Euler parameters given in <v>params</v>. \
<v>params</v> must be a list of length <v>d^2 - 1</v>."

StateVector::usage = "<f>StateVector</f>[{<s>\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\),\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"n\", \"+\", \"1\"}]]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \
RowBox[{\"2\", \" \", \"n\"}]]\)</s>}] returns pure <v>n</v>+1-dimensional pure state (ket vector) constructed form probability distribution parametrize \
by numbers {<s>\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)</s>} and phases \
{<s>\!\(\*SubscriptBox[\"\[Phi]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \"n\"]\)</s>}. \
See also: <f>SymbolicVector</f>.";


QubitKet::usage = "<f>QubitKet</f>[<s>\[Alpha]</s>,<s>\[Beta]</s>] parametrization of the pure state (as a state vector) for one qubit as \
(<s>Cos[\[Alpha]] Exp[i\[Beta]], Sin[\[Alpha]]</s>). This is equivalent to <f>StateVector</f>[{<s>\[Alpha],\[Beta]</s>}]. See also: <f>QubitPureState</f>, <f>StateVector</f>.";

QubitPureState::usage = "<f>QubitPureState</f>[<s>\[Alpha],\[Beta]</s>] - \
a parametrization of the pure state as a density matrix for one qubit. \
This is just a alias for <f>Proj</f>[<f>QubitKet</f>[<s>\[Alpha],\[Beta]</s>]]. See also: <f>QubitKet</f>.";



(* ::Subsection::Closed:: *)
(*Random states*)


RandomSimplex::usage = "<f>RandomSimplex</f>[<v>d</v>] generates a point on a <v>d</v>-dimensional simplex according to the uniform distibution.";

RandomKet::usage = "<f>RandomKet</f>[<v>d</v>] - for integer <v>d</v> returns a random ket vector in <v>d</v>-dimensional space. \
See: T. Radtke, S. Fritzsche, Comp. Phys. Comm., Vol. 179, No. 9, p. 647-664. 
<f>RandomKet</f>[{<v>d1</v>,<v>d2</v>}, <v>type</v>] - for integers <v>d1</v>, <v>d2</v> returns a random ket \
vector in <v>d1 d2</v>-dimensional space with distribution specified by <v>type</v>. 
The parameter <v>type</v> can be:
\t ''Sep'' - in this case function returns separable random vectors,
\t ''MaxEnt'' - in this case function returns maximally entangled random vectors,
\t <v>l</v> list of positive numbers summing up to 1, the length of the list must be less or equal to the Min[<v>d1</v>,<v>d2</v>] \
- in this case function returns random vectors with fixed Schmidt numbers given by <v>l</v>.";

(*<v>d</v> may be a list of integers in \
this case ket will be in product form d[[1]]\[CircleTimes]...\[CircleTimes]d[[k]]. See: T. Radtke, S. Fritzsche, Comp. Phys. Comm., Vol. 179, No. 9, p. 647-664.";
*)

GinibreMatrix::usage = "<f>GinibreMatrix</f>[<v>m</v>,<v>n</v>] returns complex matrix of dimension <s>m\[Cross]n</s> with the standard normal distribution of real and imaginary parts.";

RandomSpecialUnitary::usage = "<f>RandomSpecialUnitary</f>[<v>d</v>] returns a random special unitary matrix of size <v>d</v>. See <f>RandomUnitary</f>.";

RandomUnitary::usage = "<f>RandomUnitary</f>[<v>d</v>] returns a random unitary matrix of size <v>d</v> using <v>QR</v> decomposition. \
See: F. Mezzadri,  NOTICES of the AMS, Vol. 54 (2007), 592-604.";

RandomOrthogonal::usage = "<f>RandomOrthogonal</f>[<v>d</v>] returns a random orthogonal matrix of size <v>d</v> using <v>QR</v> decomposition. \
See: F. Mezzadri, NOTICES of the AMS, Vol. 54 (2007), 592-604."

RandomState::usage = "<f>RandomState</f>[<v>d</v>,<v>dist</v>] - random density matrix of dimension <v>d</v>. \
Argument <v>dist</v> can be ''<v>HS</v>'' (default value), ''<v>Bures</v>'' or an integer <v>K</v>. 
\t ''<v>HS</v>'' - gives uniform distribution with respect to the Hilbert-Schmidt measure. 
\t ''<v>Bures</v>'' - gives a random state distributed according to Bures measure. 
\t Integer <v>K</v> - gives a random state generated with respect to induced measure with an ancilla system od dimension K.";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


RandomDynamicalMatrix::usage = "<f>RandomDynamicalMatrix</f>[<v>d</v>,<v>k</v>] returns dynamical matrix of operation acting \
on <v>d</v>-dimensional states with <v>k</v> eigenvalues equal to 0. Thanks to Wojtek Bruzda. see Random Quantum Operations DOI[10.1016/j.physleta.2008.11.043].";

ApplyKraus::usage = "<f>ApplyKraus</f>[<v>ck</v>,<s>\[Rho]</s>] - apply channel <v>ck</v>, given as a list of Kraus operators, \
to the input state <s>\[Rho]</s>. See also: <f>ApplyUnitary</f>, <f>ApplyChannel</f>.";

ChannelToMatrix::usage = "<f>ChannelToMatrix</f>[<v>E</v>,<v>d</v>] returns matrix representation of a \
channel <v>E</v> acting on <v>d</v>-dimensional state space. First argument should be a pure function <v>E</v> \
such that <v>E</v>[<s>\[Rho]</s>] transforms input state according to the channel definition.";

ApplyChannel::usage = "<f>ApplayChannel</f>[<v>f</v>,<s>\[Rho]</s>] - apply channel <v>f</v>, given as a \
pure function, to the input state <s>\[Rho]</s>. See also: <f>ApplyUnitary</f>, <f>ApplyKraus</f>."

Superoperator::usage = "<f>Superoperator</f>[<v>kl</v>] returns matrix representation of quantum channel \
given as a list of Kraus operators. <f>Superoperator</f>[<v>fun</v>,<v>dim</v>] is just am alternative name \
for <f>ChannelToMatrix</f>[<v>fun</v>,<v>dim</v>] and returns matrix representation of quantum channel, \
given as a pure function, acting on <v>dim</v>-dimensional space. \
See also: <f>ChannelToMatrix</f>.";

DynamicalMatrix::usage = "<f>DynamicalMatrix</f> returns a dynamical matrix of quantum channel
<f>DynamicalMatrix</f>[<v>ch</v>] -  operates on a quantum channel given as a list of Kraus operators.
<f>DynamicalMatrix</f>[<v>fun</v>,<v>dim</v>]  - operates on a a function <v>fun</v> acting on <v>dim</v>-dimensional space. 
See also: <f>Superoperator</f>, <f>ChannelToMatrix</f>.";

Jamiolkowski::usage = "<f>Jamiolkowski</f>[<v>K</v>] gives the image of the Jamiolkowski isomorphism for the channel given as the list of \
Karus operators <v>K</v>. 
<f>Jamiolkowski</f>[<v>fun</v>,<v>dim</v>] gives the image of the Jamiolkowski isomorphism for the channel given \
as a function fun action on <v>dim</v>-dimensional space. See also: <f>Superoperator</f>, <f>ChannelToMatrix</f>, <f>DynamicalMatrix</f>.";

TPChannelQ::usage = "<f>TPChannelQ</f>[<v>ck</v>] performs some checks on Kraus operators <v>ck</v>. Use this if you want to check if they represent quantum channel.";

SuperoperatorToKraus::usage = "<f>SuperoperatorToKraus</f>[<v>m</v>] returns Kraus operators for a given super operator <v>m</v>.";

ProductSuperoperator::usage = "<f>ProductSuperoperator</f>[<s>\[CapitalPsi]</s>,<s>\[CapitalPhi]</s>] \
computes a product superoperator of superoperatos <s>\[CapitalPsi]</s> and <s>\[CapitalPhi]</s>.";


(* ::Section:: *)
(*Private definitions*)


(* ::Subsection:: *)
(*Internal functions (without usage strings)*)


Begin["`Private`"];

QIDocRep = {"<v>" -> "\!\(\*StyleBox[\"" , "</v>" -> "\", \"TI\"]\)", "<f>"->"\!\(\*StyleBox[\"", "</f>" -> "\", \"Input\"]\)", "<s>" -> "", "</s>" -> ""}
(MessageName[Evaluate[ToExpression[#]], "usage"] = StringReplace[MessageName[Evaluate[ToExpression[#]], "usage"],QIDocRep])& /@ Names["QI`*"];

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
	{"0.4.0",  "13/10/2011", "Zbyszek, Jarek, Gawron", "Big changes , some functions moved to QIExtras Package."},
	{"0.4.1",  "14/10/2011", "Zbyszek, Jarek", "Documentation improved."},
	{"0.4.2",  "18/10/2011", "Zbyszek, Jarek", "Partial trace improved."},
	{"0.4.3",  "19/10/2011", "Zbyszek, Jarek", "ProductSuperoperator impoved."},
	{"0.4.31", "25/10/2011", "Zbyszek", "Small changes."},
	{"0.4.32", "16/12/2011", "Zbyszek, Jarek", "Documentation improved."},
	{"0.4.33", "17/12/2011", "Jarek", "Negativity fixed - thanks to Fatih \[CapitalODoubleDot]zayd\[DotlessI]n"},
	{"0.4.34", "11/01/2012", "Zbyszek", "VectorSchmidtDecomposition fixed"},
	{"0.4.35", "26/01/2012", "Gawron", "Reintroduced special unitary prametrization"},
	{"0.4.36", "04/02/2012", "Jarek", "BaseVectors and BaseMatrices moved from QIExtras to QI."},
	{"0.4.37", "21/12/2012", "Jarek", "Fixed SuperoperatorToKraus function - thanks to Vinayak Jagadish."}
};  

qiVersion = Last[qiHistory][[1]];

qiLastModification = Last[qiHistory][[2]];

qiAbout = "QI is a package of functions for Mathematica computer algebra system, which implements \
number of functions used in the analysis of quantum states and quantum operations. In contrast to \
many available packages for symbolic and numerical simulation of quantum computation presented \
package is focused on geometrical aspects of quantum information theory.";



HyperlinkToString::usage = "HyperlinkToString[text,link] creates a link labeled with text to the given URL and returns it as a Mathematica string.";

HyperlinkToString[text_,link_]:="\!\(\*ButtonBox[StyleBox[\""<>text<>"\", \"SR\"],Active->True,BaseStyle->\"Link\",ButtonData->\""<>link<>"\"]\)";

DOIToString::usage = "DOIToString[text,doi] creates a link labeled with text to the given DOI and returns it as a Mathematica string.";

DOIToString[text_,doi_]:="\!\(\*ButtonBox[StyleBox[\""<>text<>"\", \"SR\"],Active->True,BaseStyle->\"Link\",ButtonData->\"http://dx.doi.org/"<>doi<>"\"]\)";



(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


sx = {{0,1},{1,0}};
sy = {{0,-I},{I,0}};
sz = {{1,0},{0,-1}};
id = {{1,0},{0,1}};
wh = {{1,1},{1,-1}};
cnot = {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
Id[n_]:=IdentityMatrix[n];


BaseVectors[n_Integer]:=Table[UnitVector[n,k],{k,1,n}];


BaseMatrices[n_Integer]:=Table[Unres[UnitVector[n^2,k]],{k,1,n^2}];


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


CircleTimes[x_?MatrixQ,y_?MatrixQ] := KroneckerProduct[x,y];
CircleTimes[x_?VectorQ,y_?VectorQ] := Flatten[KroneckerProduct[x,y]];
CircleTimes[a_,b__] := CircleTimes[a, CircleTimes[b]]

SquareMatrixQ[A_]:= Block[{dims=Dimensions[A]},
	(Length[dims]==2 )&&(dims[[1]]==dims[[2]])
];

SymbolicMatrix[sym_,d1_?IntegerQ]:=SymbolicMatrix[sym,d1,d1];
SymbolicMatrix[sym_,d1_?IntegerQ,d2_?IntegerQ] := Table[Subscript[sym, i,j], {i,1,d1}, {j,1,d2}];

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
Reshuffle::argerr = "Reshuffle works only for square matrices of dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\)\[Times]\!\(\*SuperscriptBox[\"d\", \"2\"]\), \
where d is an Integer, for other dimensions use ReshuffleGeneral";

Reshuffle[A_,{n_,m_}]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n[[2]]+i1,1+i2;;m[[2]]+i2]],{i1,0,n[[1]] n[[2]]-1,n[[2]]},{i2,0,m[[1]]*m[[2]]-1,m[[2]]}]
,1];


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


VectorSchmidtDecomposition[vec_?VectorQ, d_?ListQ] := 
  Block[{mtx, u, w, v, vals, snum = Min[d[[1]], d[[2]]]},
   	mtx = Partition[vec, d[[2]]];
   	{u, w, v} = SingularValueDecomposition[mtx];
   	vals = Select[Diagonal[w], # != 0 &];
   	snum = Length[vals];
     {vals, u\[Transpose][[1 ;; snum]], 
     v\[ConjugateTranspose][[1 ;; snum]]}\[Transpose]
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
   	mtx = Reshuffle[op, {n, m}];
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
(*Partial trace and transposition*)


PartialTrace[A_, {}]:=A;
PartialTrace[A_, dim_, {}]:=A;

PartialTrace[A_, s_] := Block[ {d = Length[A], sys},
    If[ IntegerQ[s],
      sys = {s}, (*else*)
      sys = s
    ];  
    If[  IntegerQ[Sqrt[d]],
    If[ Length[sys]==1 ,
      Switch[sys[[1]],
       1, PartialTrace1[A],
       2, PartialTrace2[A],
       _,Print["B1"]; Message[PartialTrace::syserr, sys]
       ], (*else*)
      If[ sys=={1,2} || sys=={2,1},
        Tr[A], (*else*)
        Print["B2"]; Message[PartialTrace::syserr, sys]
      ]
    ]
    ,
    Message[PartialTrace::dimerr, Dimensions[A]]    
    ]
  ];
  


PartialTrace[A_,dim_?VectorQ,s_]:=Block[{sys},
	If[IntegerQ[s], sys={s}, (*else*) sys=s];
    If[Length[dim] == 2 && Length[sys]==1,
	   Switch[sys[[1]],
        1, PartialTrace1[A,dim],
        2, PartialTrace2[A,dim],
        _, Message[PartialTrace::sysspecerr, Length[dim], Select[sys, Or[#>Length[dim],  # < 1]&]]
        ]
        ,(*else*)
        If[Fold[#1 && (1<= #2 <= Length[dim])&,True,sys], 
        	If[Union[sys] == Range[Length[dim]], 
        		Tr[A], (*else*)
                PartialTraceGeneral[A,dim,sys]
        	],
        (*else*)
            Message[PartialTrace::sysspecerr, Length[dim], Select[sys, Or[#>Length[dim],  # < 1]&]]
        ]
    ]
];

PartialTrace::syserr = "The second argument is expected to be 1, 2 or {1,2} (`1` was given)";
PartialTrace::dimerr = "This function expects square matrix of size d^2\[Times]d^2 (the matrix of size `1` was given)";
PartialTrace::sysspecerr = "The system specification is invalid. In the case of `1`-partite systems, it is inpossible to trace out with respect to sub-systems `2`.";

PartialTrace1[X_] := Block[{d = Sqrt[Length[X]]}, Total[X[[1 + d # ;; d + d #, 1 + d # ;; d + d # ]] & /@ Range[0, d - 1]]];
PartialTrace2[X_] := Block[{d = Sqrt[Length[X]]}, Total[X[[1 + # ;; d^2 ;; d, 1 + # ;; d^2 ;; d]] & /@ Range[0, d - 1]]];
PartialTrace1[X_, {d1_, d2_}] := Total[X[[1 + d2 # ;; d2 + d2 #, 1 + d2 # ;; d2 + d2 # ]] & /@ Range[0, d1 - 1]];
PartialTrace2[X_, {d1_, d2_}] := Total[X[[1 + # ;; d1*d2 ;; d2, 1 + # ;; d1*d2 ;; d2]] & /@ Range[0, d2 - 1]];

ListReshape[list_, shape_] := 
  FlattenAt[Fold[Partition[#1, #2] &, Flatten[list], Reverse[shape]], 
   1];

PartialTraceGeneral[A_,dim_?VectorQ,sys_?VectorQ] := Block[
	{offset, keep, dispose, keepdim, disposedim, perm1, perm2, perm, tensor},
	offset=Length[dim];
	keep=Complement[Range[offset], sys];
	dispose=Union[sys];
	perm1=Join[dispose,keep];
	perm2=perm1+offset;
	perm=Ordering[Join[perm1,perm2]];
	tensor=ListReshape[A, Join[dim,dim]];
	keepdim=Apply[Times, Join[dim, dim][[keep]]];
	disposedim=Apply[Times, Join[dim, dim][[dispose]]];
	tensor=Transpose[tensor,perm];
	tensor=ListReshape[tensor,{disposedim,keepdim,disposedim,keepdim}];
	Sum[tensor[[i,All,i,All]],{i,1,disposedim}]
];

PartialTranspose[\[Rho]_,dim_?VectorQ,sys_?VectorQ]:=Block[{offset,tensor,perm,idx1,idx2,s,targetsys},
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
];


(* ::Subsection::Closed:: *)
(*Distance measures*)


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
(*Bloch Representation*)


Lambda1[i_ ,j_,n_]:=Table[KroneckerDelta[j,\[Mu]]KroneckerDelta[i,\[Nu]] + KroneckerDelta[j,\[Nu]]KroneckerDelta[i,\[Mu]] ,{\[Mu],1,n},{\[Nu],1,n}];

Lambda2[i_ ,j_,n_]:=Table[-I(KroneckerDelta[i,\[Mu]]KroneckerDelta[j,\[Nu]] - KroneckerDelta[i,\[Nu]]KroneckerDelta[j,\[Mu]]) ,{\[Mu],1,n},{\[Nu],1,n}];

Lambda3[i_,n_]:=Sqrt[2/(i^2-i)]DiagonalMatrix[Join[Append[Table[1,{i-1}],-(i-1)],Table[0,{n-i}]]];

GeneralizedPauliMatrices[n_]:=Block[{l1,l2,l3,i,j},
    l1=Flatten[Table[Lambda1[i,j,n],{i,1,n},{j,i+1,n}],1];
    l2=Flatten[Table[Lambda2[i,j,n],{i,1,n},{j,i+1,n}],1];
    l3=Table[Lambda3[i,n],{i,2,n}];
    Join[l1,l2,l3]
];


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



(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_]:=
Exp[I \[Alpha]] DiagonalMatrix[{Exp[-I \[Beta]/2],Exp[I \[Beta]/2]}].{{Cos[\[Gamma]/2],-Sin[\[Gamma]/2]},{Sin[\[Gamma]/2],Cos[\[Gamma]/2]}}.DiagonalMatrix[{Exp[-I \[Delta]/2],Exp[I \[Delta]/2]}];

Unitary2Euler[\[Alpha]_,\[Theta]_,\[Phi]_,\[Psi]_]:={{E^(I*\[Alpha] + I*\[Phi])*Cos[\[Theta]], E^(I*\[Alpha] + I*\[Psi])*Sin[\[Theta]]}, {-(E^(I*\[Alpha] - I*\[Psi])*Sin[\[Theta]]), E^(I*\[Alpha] - I*\[Phi])*Cos[\[Theta]]}}

SpecialUnitary2[\[Beta]_,\[Gamma]_,\[Delta]_]:=Unitary2[0,\[Beta],\[Gamma],\[Delta]];

SpecialUnitary[d_, params_] :=  
  Block[{psi, chi, r, s, phi, i, j, u, e, phi0, psi0, chi0, lparams,k = 1}, 
   lparams = Map[FractionalPart[Abs[#]] &, params];
   Do[psi[r, s] = 2*Pi*lparams[[k++]];, {r, 1, d - 1}, {s, r + 1, d}];
   Do[chi[r, s] = 0;, {r, 2, d - 1}, {s, r + 1, d}];
   Do[chi[1, s] = 2*Pi*lparams[[k++]];, {s, 2, d}];
   Do[phi[r, s] = ArcSin[(lparams[[k++]])^(1/(2 r))];, {r, 1, d - 1}, {s, r + 1, d}];
   e=SparseArray[{}, {d,d,d,d}];
   Do[
   	e[[r, s]] = IdentityMatrix[d];
    e[[r, s, r, r]] = Cos[phi0]*Exp[I*psi0];
    e[[r, s, s, s]] = Cos[phi0]*Exp[-I*psi0];
    e[[r, s, r, s]] = Sin[phi0]*Exp[I*chi0];
    e[[r, s, s, r]] = -Sin[phi0]*Exp[-I*chi0];, 
    	{r, 1, d - 1}, 
    		{s, r + 1, d}
    		];
   u = IdentityMatrix[d];
   Do[u = (Normal[e[[r, r + 1]]] /. {phi0 -> phi[d - r, s + 1], 
          psi0 -> psi[d - r, s + 1], 
          chi0 -> chi[d - r, s + 1]}).u;, {s, d - 1, 1, -1}, {r, 
     d - 1, d - s, -1}];
   u];


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


QubitKet[\[Alpha]_,\[Beta]_]:={Cos[\[Alpha]], Exp[I*\[Beta]]*Sin[\[Alpha]]};

QubitPureState[\[Alpha]_,\[Beta]_]:=Proj[QubitKet[\[Alpha],\[Beta]]];


(* ::Subsection::Closed:: *)
(*Random states*)


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
RandomKet::argerr = "Error - the last parameter should be a list of Schmidt numbers or one of the predefined values: ''Sep'' or ''MaxEnt''."; 



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

SuperoperatorToKraus[m_] :=  Block[{val, vec}, {val, vec} = Eigensystem[Reshuffle[m]];   Sqrt[val] (Unvec[#/Norm[#]] & /@ vec)];


ProductSuperoperator[M1_,M2_]:=Block[{prod = M1\[CircleTimes]M2, q1, d1 = Sqrt[Length[M1]], d2 = Sqrt[Length[M2]]},
    q1 = Table[Res[Reshuffle[Unres[prod[[i]], d2 d2], {{d1, d1}, {d2, d2}}]], {i,1, d1*d1*d2*d2}];
    Table[Res[ Reshuffle[ Unres[(q1\[Transpose])[[i]], d2 d2], {{d1, d1}, {d2, d2}}]], {i,1,d1*d1*d2*d2}]\[Transpose]
];

RandomDynamicalMatrix[n_,m_:0]:=Block[{X,Y,sY},	
	X=GinibreMatrix[n^2,n^2-m];
	Y=PartialTrace[X.X\[ConjugateTranspose],{n,n},{1}];
	sY=MatrixPower[Y,-1/2];
	KroneckerProduct[IdentityMatrix[n],sY].X.X\[ConjugateTranspose].KroneckerProduct[IdentityMatrix[n],sY]
];


(* ::Section:: *)
(*Package footer*)


Print["Package QI version ", QI`Private`qiVersion, " (last modification: ", QI`Private`qiLastModification, ")."];

End[];

Protect@@Names["QI`*"]

EndPackage[];
   
