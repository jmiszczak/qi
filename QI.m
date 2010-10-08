(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package header*)


(* File: QI.m *)
(* Description: Mathematica package for the analysis of quantum states and operations *)
(* Authors: Jaroslaw Miszczak <miszczak@iitis.pl>, Piotr Gawron <gawron@iitis.pl>, Zbigniew Puchala <z.puchala@iitis.pl> *)
(* License: GPLv3 *)


BeginPackage["QI`"];


Begin["`Private`"]; (* Start `Private` context *)


qiAuthors = "Jaroslaw Miszczak <miszczak[at]iitis[dot]pl>, Piotr Gawron <gawron[at]iitis[dot]pl>, Zbigniew Puchala <z.puchala[at]iitis[dot]pl>";


qiLicense = "GPLv3 <http://www.gnu.org/licenses/gpl.html>";


qiVersion = "0.3.20";


qiLastModification = "October 4, 2010";


qiHistory = {
	{"0.1.0", "05/06/2009", "Initial version"}, 
	{"0.1.1", "09/06/2009", "Fixed \[Eta] and \[Eta]2 functions, fixed problem with protected symbols."}, 
	{"0.1.2", "18/06/2009", "Added quantum channel parametrization for one qubit."}, 
	{"0.1.3", "03/07/2009", "Added alternative reshuffling."},
	{"0.1.4", "", "Changed default print output."}, 
	{"0.2.0", "09/07/2009", "Documentation generator added."},
	{"0.2.1", "", "Changed QubitGeneralState function."},
	{"0.2.2", "04/08/2009", "Added reshuffling permutation and product of superoperators."},
	{"0.2.3", "13/08/2009", "Minor update in documentation."},
	{"0.2.4", "14/09/2009", "Fixed \[Eta] function."},
	{"0.2.5", "27/08/2009", "Fixed Werner state and added IsotropicState."},
	{"0.2.6", "06/10/2009", "Spelling improvements."},
	{"0.2.7", "02/11/2009", "Improved '\[CircleTimes]' usage, added applyUnitary."},
	{"0.2.8", "03/22.1009", "Fixed small problem with MaxEnt."},
	{"0.2.9", "", "Some code cleanups, fixed SchmidtDecomposition."},
	{"0.3.0", "06/11/2009", "Added OperatorSchmidtDecomposition."},
	{"0.3.1", "12/11/2009", "Added Concurrence4, fixed Cnot problem, some code cleanups."},
	{"0.3.2", "17/11/2009", "Documentation improvements, SquareMatrixQ predicate."},
	{"0.3.3", "18/11/2009", "All qi* functions and constants moved to the QI`Private context, added list of QI`* names."},
	{"0.3.4", "19/11/2009", "Added Negativity, fixed SchmidtDecomposition."},
	{"0.3.5", "21/11/2009", "SchmidtDecomposition now accepts vectos as well as matrices."},
	{"0.3.6", "24/11/2009", "Minor update in qiNames."},
	{"0.3.7", "04/12/2009", "OperatorSchmidtDecomposition fixed."},
	{"0.3.8", "04/01/2010", "Added local vars in RandomState."},
	{"0.3.9", "06/01/2010", "Added error message in Ket."},
	{"0.3.10", "19/01/2010", "Improved Davies map."},
	{"0.3.11", "26/01/2010", "Improved simplex generation algorithm, added some function for random vectors."},
	{"0.3.12", "04/03/2010", "Changed parameter in Swap gate."},
	{"0.3.13", "26/03/2010", "Fixed bug with state parametrization."},
	{"0.3.14", "25/05/2010", "Fixed numerical bug in Concurrence4 - Chop function added."},
	{"0.3.15", "07/06/2010", "RandomMaximallyEntangledNumericalRange added."},
	{"0.3.16", "20/06/2010", "Alternative version of Ketbra function."},
	{"0.3.17", "11/07/2010", "Name changed for Davies map."},
	{"0.3.18", "11/08/2010", "Fiexd bug in GeneralizedPauliKraus function reported by Fatih Ozaydin and one syntax error."},
	{"0.3.19", "13/09/2010", "Fixed bug in QubitDecayKraus and QubitDepolarizingKraus, QubitBitflipKraus, QubitPhaseflipKraus, QubitBitphaseflipKraus."},
	{"0.3.20", "04/10/2010", "Fixed inconsistency in QubitDepolarizingKraus and DepolarizingChannel."}
};


qiAbout = "QI is a package of functions for Mathematica computer algebra system, which implements 
number of functions used in the analysis of quantum states and quantum operations. In contrast to 
many available packages for symbolic and numerical simulation of quantum computation presented 
package is focused on geometrical aspects of quantum information theory.";


qiGenDoc::usage = "Generate documentation for the QI package. This function accepts up to two arguments - file name and output directory. By default the first one is set to 'qi_usage.tex' and the second one to '~/zksi-repo/qi/doc'.";


qiAbout::usage = "Display a short information about the QI package.";


qiVersion::usage = "Display the version of the QI package.";


qiHistory::usage = "Display the history of modifications for the QI package.";


qiConstInfo = ""(*" This is predefined constant."*);


qiNames = {"ApplyChannel","ApplyKraus","ApplyUnitary","BaseMatrices","BaseVectors","BlochVector","ChannelToMatrix","cnot","Commutator","ComplexToPoint","Concurrence4","DepolarizingChannel","DynamicalMatrix","ExpectationValue","ExtendKraus","Fidelity","GellMannMatrices","GeneralizedPauliKraus","GeneralizedPauliMatrices","GeneralizedPauliX","GeneralizedPauliZ","GinibreMatrix","HolevoWernerChannel","id","IdentityChannel","IsotropicState","Jamiolkowski","Ket","Ketbra","KetFromDigits","KroneckerDeltaMatrix","KroneckerSum","Lambda1","Lambda2","Lambda3","Log0","MatrixAbs","MatrixElement","MatrixIm","MatrixRe","MatrixSqrt","MaxEnt","MaxMix","Negativity","NumericalRangeBound","OperatorSchmidtDecomposition","PartialTraceA","PartialTraceB","PartialTraceGeneral","PartialTransposeA","PartialTransposeB","PartialTransposeGeneral","PauliMatrices","ProbablityVector","ProbBures","ProbBuresNorm","ProbHS","ProbHSNorm","ProdDiff2","ProdSum","ProductSuperoperator","Proj","QFT","QuantumChannelEntropy","QuantumEntropy","QubitBitflipChannel","QubitBitflipKraus","QubitBitphaseflipChannel","QubitBitphaseflipKraus","QubitBlochState","QubitDaviesSuperoperator","QubitDecayKraus","QubitDepolarizingKraus","QubitDynamicalMatrix","QubitGeneralState","QubitKet","QubitPhaseflipChannel","QubitPhaseflipKraus","QubitPhaseKraus","QubitPureState","QutritSpontaneousEmissionKraus","RandomComplexUnitVector","RandomDynamicalMatrix","RandomEntangledUnitVector","RandomKet","RandomNormalMatrix","RandomProductKet","RandomProductNumericalRange","RandomRealUnitVector","RandomSimplex","RandomSpecialUnitary","RandomState","RandomUnitary","RandomUnitVector","RandomUnitVectorSchmidt","Res","Reshuffle","Reshuffle2","ReshuffleGeneral","ReshuffleGeneral2","\[AliasDelimiter]Permutation","SchmidtDecomposition","SpecialUnitary2","SquareMatrixQ","StateFromBlochVector","StateVector","Subfidelity","Superfidelity","Superoperator","Swap","sx","sy","SymbolicHermitianMatrix","SymbolicMatrix","SymbolicVector","sz","TPChannelQ","TraceDistance","TraceNorm","TransposeChannel","Unitary2","Unitary3","Unitary4Canonical","Unres","Unvec","VandermondeMatrix","Vec","VectorSchmidtDecomposition","WernerState","wh","\[Delta]","\[Eta]","\[Eta]2","\[Lambda]","\[Lambda]1","\[Lambda]2","\[Lambda]3","\[Lambda]4","\[Lambda]5","\[Lambda]6","\[Lambda]7","\[Lambda]8","\[Sigma]x","\[Sigma]y","\[Sigma]z"};


End[]; (* End of `Private` context *)


Print["Package QI version ", QI`Private`qiVersion, " (last modification: ", QI`Private`qiLastModification, ")."];


Unprotect@@QI`Private`qiNames;
Clear@@QI`Private`qiNames;


$PrePrint = If[SquareMatrixQ[#], MatrixForm[#], #]&;


(* ::Section::Closed:: *)
(*Help messages*)


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


KroneckerSum::usage = "KroneckerSum[A,B] returns the Kronecker sum of matrices A and B defined as A\[CircleTimes]1+1\[CircleTimes]B. Alternative syntax A\[CirclePlus]B for KroneckerSum[A,B] is also provided. See also: KroneckerProduct.";


SquareMatrixQ::usage = "SquareMatrixQ[A] returns True only if A is a square matrix, and gives False otherwise.";


SymbolicMatrix::usage = "SymbolicMatrix[a,m,n] returns m\[Cross]n-matrix with elements a[i,j], i=1,...,m, j=1,...,n. If the third argument is ommited this function returns square m\[Cross]m matrix. This functions can save you some keystrokes and, thanks to TeXForm function, its results can be easily incorporated in LaTeX documents.";


SymbolicVector::usage = "SymbolicVector[a,n] is equivalent to Matrix[a,n,1] and it returns a vector with m elements a[i],i=1,...,n.";


SymbolicHermitianMatrix::usage = "SymbolicHermitianMatrix[sym,n] produces a n\[Cross]n Hermitian matrix. See also: SymbolicMatrix, SymbolicVector.";


ComplexToPoint::usage = "ComplexToPoint[z] returns a real and an imaginary parts of a complex number z as a pair of real numbers.";


MatrixSqrt::usage= "MatrixSqrt[A] returns square root for the matrix A.";


MatrixAbs::usage= "MatrixAbs[A] returns absolute value for matrix A defined as MatrixSqrt[A.A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)]. See also: MatrixSqrt.";


MatrixRe::usage = "Hermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A+A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)).";


MatrixIm::usage = "Antyhermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A-A\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)).";


ExpectationValue::usage = "ExpectationValue[s,A] - accepts a state vector or a density matrix as a first argument and calculates the expectation value for the measurement of A in the state s.";


Commutator::usage = "Commutator[A,B] returns the commutator of matrices A and B i.e. Commutator[A,B] = A.B - B.A.";


(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


Fidelity::usage = "Fidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns the quantum fidelity between states \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) calculated using a simplified formula as (\[Sum]\!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\)\!\(\*SuperscriptBox[\")\", \"2\"]\), where \!\(\*SubscriptBox[\"\[Lambda]\", \"i\"]\) are the eigenvalues of \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\).";


Superfidelity::usage = "Superfidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] calculates superfidelity between \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) defined as Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] + Sqrt[1-Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)]]Sqrt[1-Tr[\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\).\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)]]. See: J.A. Miszczak et al., Quantum Information & Computation, Vol.9 No.1&2 (2009)."; 


Subfidelity::usage = "Subfidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns subfidelity between states \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\) See: J.A. Miszczak et al., Quantum Information & Computation, Vol.9 No.1&2 (2009).";


TraceNorm::usage = "TraceNorm[A] = \[Sum]\!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\), where \!\(\*SubscriptBox[\"\[Sigma]\", \"i\"]\) are the singular values of A. See also: TraceDistance.";


TraceDistance::usage = "TraceDistance[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns the trace distance between matrices \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\), which is defined as \!\(\*FractionBox[\"1\", \"2\"]\)tr|\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)-\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)|.";


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


sx::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\)." <> QI`Private`qiConstInfo;
sy::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\)." <> QI`Private`qiConstInfo;
sz::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\)." <> QI`Private`qiConstInfo;
\[Sigma]x::usage = sx::usage <> " This is an alternative notation for sx.";
\[Sigma]y::usage = sy::usage <> " This is an alternative notation for sy.";
\[Sigma]z::usage = sz::usage <> " This is an alternative notation for sz.";
id::usage = "Identity matrix for one qubit. See also: IdentityMatrix." <> QI`Private`qiConstInfo;
wh::usage = "Hadamard gate for one qubit. See also: QFT." <> QI`Private`qiConstInfo;


\[Lambda]1::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]2::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]3::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]4::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"4\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]5::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"5\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]6::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"6\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]7::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"7\"]\)." <> QI`Private`qiConstInfo;
\[Lambda]8::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"8\"]\)." <> QI`Private`qiConstInfo;


Proj::usage = "Proj[{\!\(\*SubscriptBox[\"v\", \"1\"]\),\!\(\*SubscriptBox[\"v\", \"2\"]\),...,\!\(\*SubscriptBox[\"v\", \"n\"]\)}] returns projectors for the vectors in the input list.";


BaseVectors::usage = "BaseVectors[n] returns a list with the canonical basis in n-dimensional Hilbert space \!\(\*SuperscriptBox[\"\[DoubleStruckCapitalC]\", \"n\"]\). See also: BaseMatrices.";


BaseMatrices::usage = "BaseMatrices[n] returns a list with the canonical basis in n\[Cross]n-dimensional Hilbert-Schmidt space \!\(\*SubscriptBox[\"\[DoubleStruckCapitalM]\", \"n\"]\). See also: BaseVectors.";


KroneckerDeltaMatrix::usage = "KroneckerDeltaMatrix[i,j,d] returns d\[Cross]d matrix with 1 at position (i,j) and zeros elsewhere.";


Lambda1::usage = "Lambda1[i,j,n] generalized Pauli matrix. For example Lambda1[1,2,2] is equal to Pauli \[Sigma]x. See also: GeneralizedPauliMatrices.";


Lambda2::usage = "Lambda2[i,j,n] generalized Pauli matrix. For example Lambda2[1,2,2] is equal to \[Sigma]y. See also: GeneralizedPauliMatrices.";


Lambda3::usage ="Lambda3[i,n] generalized Pauli matrix. For example Lambda3[2,2] is equal to \[Sigma]z. See also: GeneralizedPauliMatrices.";


GeneralizedPauliMatrices::usage = "GeneralizedPauliMatrices[n] returns list of generalized Pauli matrices for SU(n). For n=2 these are just Pauli matrices and for n=3 - Gell-Mann matrices. Note that identity matrix is not included in the list. See also: PauliMatrices, GellMannMatrices, \[Lambda], Lambda1, Lambda2, Lambda3.";


\[Lambda]::usage = "\[Lambda][i,n] is defined as GeneralizedPauliMatrices[n][[i]].";


PauliMatrices::usage = "Predefined list of Pauli matrices {\!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\)}. Use Map[MatrixForm[#]&,PauliMatrices] to get this list in more readible form.";


GellMannMatrices::usage = "List of Gell-Mann matrices. Use Map[MatrixForm[#]&,GellMannMatrices] to get this list in more readable form.";


(* ::Subsection::Closed:: *)
(*Quantum gates*)


Swap::usage="Swap[d] returns permutation operator \!\(\*UnderoverscriptBox[RowBox[{\" \", \"\[Sum]\"}], RowBox[{\"i\", \"=\", \"0\"}], RowBox[{\"n\", \"-\", \"1\"}]]\)\!\(\*UnderoverscriptBox[RowBox[{\" \", \"\[Sum]\"}], RowBox[{\"j\", \"=\", \"0\"}], RowBox[{\"n\", \"-\", \"1\"}]]\)|i\[RightAngleBracket]\[LeftAngleBracket]j|\[CircleTimes]|j\[RightAngleBracket]\[LeftAngleBracket]i\[VerticalSeparator] acting on d-dimensional, d=\!\(\*SuperscriptBox[\"n\", \"2\"]\) space and exchanging two \[Sqrt]d-dimensional subsystems.";


QFT::usage = "QFT[n,method] - quantum Fourier transform of dimension n. This function accepts second optional argument, which specifies method used in calculation. Parameter method can be equal to 'Symbolic', which is default, or 'Numerical'. The second option makes this function much faster.";


cnot::usage = "Controlled not matrix for two qubits."<> QI`Private`qiConstInfo ;


GeneralizedPauliX::usage = "Generalized Pauli matrix X. See also: \!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\)";


GeneralizedPauliZ::usage = "Generalized Pauli matrix Z. See also: \!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\)";


(* ::Subsection::Closed:: *)
(*Special states*)


Ket::usage = "Ket[i,d] returns |i\[RightAngleBracket] in d-dimensional Hilbert space. See also: StateVector for a different parametrization.";


Ketbra::usage = "This function can be used in two ways. Ketbra[i,j,d] returns \[VerticalSeparator]i\[RightAngleBracket]\[LeftAngleBracket]j\[VerticalSeparator] acting on d-dimensional space. See also: Proj. Ketbra[v1,v2] returns the apropritae operator for vectors v1 and v2..";


KetFromDigits::usage = "KetFromDigits[list,base] - ket vector labeled by a list of digits represented in given base.";


MaxMix::usage = "MaxMix[n] - the maximally mixed state in a n-dimensional space of density matrices.";


MaxEnt::usage = "MaxEnt[n] - maximally entangled state in n dimensional vector space. Note that n must be perfect square.";


WernerState::usage = "WernerState[d,p] - generalized Werner state d\[Cross]d-dimensional system with the mixing parameter p\[Element][0,1]. This state is defined as p Proj[MaxEnt[d]] + (1-p) MaxMix[d],. See also: MaxEnt, MaxMix..";


IsotropicState::usage = "IsotropicState[d,p] - isotropic state of dimensions d\[Cross]d with a parameter p\[Element][0,1]. This state is defined as p Proj[MaxEnt[d]] + (1-p)/(\!\(\*SuperscriptBox[\"d\", \"2\"]\)-1)(I-Proj[MaxEnt[d]]]). This family of states is invariant for the operation of the form U\[CircleTimes]\!\(\*SuperscriptBox[\"U\", \"\[Star]\"]\).";


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


VectorSchmidtDecomposition::usage = "VectorSchmidtDecomposition[vec,d1,d2] - Schmidt decomposition of the vector vec in d1\[Cross]d2-dimensional Hilbert space.";


OperatorSchmidtDecomposition::usage = "OperatorSchmidtDecomposition[mtx,d1,d2] - Schmidt decomposition of mtx in the Hilbert-Schmidt space of matrices of dimension d1\[Cross]d2.";


SchmidtDecomposition::usage = "SchmidtDecomposition[e,d1,d2] - accepts a vector or a matrix as a first argument and returns apropriate Schmidt decomposition. See also: VectorSchmidtDecomposition, OperatorSchmidtDecomposition.";


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec::usage = "Vec[A] - vectorization of the matrix A column by column. See also: Res.";


Unvec::usage = "Unvec[v,c] - de-vectorization of the vector into the matrix with c columns. If the second parameter is omitted then it is assumed that v can be mapped into square matrix. See also: Unres, Vec.";


Res::usage = "Res[A] is equivalent to Vec[Transpose[A]]. Reshaping maps matrix A into a vector row by row.";


Unres::usage = "Unres[v,c] - de-reshaping of the vector into a matrix with c columns. If the second parameter is omitted then it is assumed that v can be mapped into a square matrix. See also: Unvec, Res.";


Reshuffle::usage = "Reshuffle[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices. If  the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See also: ReshuffleGeneral, Reshuffle2.";


Reshuffle2::usage = "Alternative definition of the reshuffling operation. Reshuffle2[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices which are transposed versions of standard base matrices. If the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See: See also: ReshuffleGeneral, Reshuffle, BaseMatrices";


ReshuffleGeneral::usage = "ReshuffleGeneral[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix.";


ReshuffleGeneral2::usage = "ReshuffleGeneral2[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix - given by alternative definition of the reshuffling operation.";


MatrixElement::usage = "MatrixElement[n,\[Nu],m,\[Mu],dim,A] - returns the matrix element of a matrix A indexed by two double indices n, \[Nu] and m, \[Mu] of the composite sytem of dimensions dim=dim1*dim2.";


ReshufflePermutation::usage = "ReshufflePermutation[dim1,dim2] - permutation matrix equivalent to the reshuffling operation on dim1\[Cross]dim2-dimensional system. See also: Reshuffling.";


ProductSuperoperator::usage = "ProductSuperoperator[\[CapitalPsi],\[CapitalPhi]] computes a product superoperator of superoperatos \[CapitalPsi] and \[CapitalPhi].";


(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2::usage = "Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]] returns the Euler parametrization of U(2).";


SpecialUnitary2::usage ="SpecialUnitary2[\[Beta],\[Gamma],\[Delta]] returns the Euler parametrization of SU(2). This is equivalent to Unitary2[0,\[Beta],\[Gamma],\[Delta]].";


Unitary3::usage = "Unitary3[\[Alpha],\[Beta],\[Gamma],\[Tau],a,b,c,ph] returns the Euler parametrization of U(3).";


Unitary4Canonical::usage = "Parametrization of non-local unitary matrices for two qubits. See: B. Kraus, J.I. Cirac, Phys. Rev. A 63, 062309 (2001), quant-ph/0011050v1.";


ProbablityVector::usage = "ProbablityVector[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}] returns probability vectors of dimension n+1 parametrize with {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}. See also: StateVector.";


StateVector::usage = "StateVector[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\),\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"n\", \"+\", \"1\"}]]\),...,\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"2\", \" \", \"n\"}]]\)}] returns pure n+1-dimensional pure state (ket vector) constructed form probability distribution parametrize by numbers {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)} and phases {\!\(\*SubscriptBox[\"\[Phi]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \"n\"]\)}. See also: ProbablityVector, SymbolicVector.";


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet::usage = "QubitKet[\[Alpha],\[Beta]] parametrization of the pure state (as a state vector) for one qubit as (Cos[\[Alpha]] Exp[i\[Beta]], Sin[\[Alpha]]). This is equivalent to StateVector[{\[Alpha],\[Beta]}]. See also: QubitPureState, StateVector.";


QubitPureState::usage = "QubitPureState[\[Alpha],\[Beta]] - a parametrization of the pure state as a density matrix for one qubit. This is just a alias for Proj[QubitKet[\[Alpha],\[Beta]]]. See also: QubitKet.";


QubitBlochState::usage = "QubitBlochState[\[Rho]] - a parametrization of the one-qubit mixed state on the Bloch sphere.";


QubitGeneralState::usage = "QubitGeneralState[\[Alpha],\[Beta],\[Gamma],\[Delta],\[Lambda]] - Parametrization of the one-qubit mixed state using rotations and eigenvalues. Returns one-qubits density matrix with eigenvalues \[Lambda] and 1-\[Lambda] rotated as U.diag(\[Lambda],1-\[Lambda]).\!\(\*SuperscriptBox[\"U\", \"\[Dagger]\"]\) with U defined by parameters \[Alpha],\[Beta],\[Gamma] and \[Delta].";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


IdentityChannel::usage = "IdentityChannel[n,\[Rho]] - apply the identity operation to a n-dimensional density matrix \[Rho].";


TransposeChannel::usage = "TransposeChannel[n,\[Rho]] - apply the transposition operation to a n-dimensional density matrix \[Rho]. Note that this operations is not completely positive.";


DepolarizingChannel::usage = "DepolarizingChannel[n,p,\[Rho]] - apply the completely depolarizing channel with parameter p acting to a n-dimensional input state \[Rho]. See also: QubitDepolarizingKraus, HolevoWernerChannel.";


HolevoWernerChannel::usage = "HolevoWernerChannel[n,p,\[Rho]] - apply the Holeve-Werner channel, also known as transpose-depolarizing channel, with parameter p acting to a n-dimensional input state \[Rho]. See also: DepolarizingChannel.";


ChannelToMatrix::usage = "ChannelToMatrix[E,d] returns matrix representation of a channel E acting on d-dimensional state space. First argument should be a pure function E such that E[\[Rho]] transforms input state according to the channel definition.";


GeneralizedPauliKraus::usage = "GeneralizedPauliKraus[d,P] - list of Kraus operators for d-dimensional generalized Pauli channel with the d-dimesnional matrix of parameters P. See: M. Hayashi, Quantum Information An Introduction, Springer 2006, Example 5.8, p. 126.";


ApplyKraus::usage = "ApplyKraus[ck,\[Rho]] - apply channel ck, given as a list of Kraus operators, to the input state \[Rho]. See also: ApplyUnitary, ApplyChannel.";


ApplyUnitary::usage = "ApplyUnitary[U,\[Rho]] - apply unitary a unitary matrix U to the input state \[Rho]. See also: ApplyKraus, ApplyChannel.";


ApplyChannel::usage = "ApplayChannel[f,\[Rho]] - apply channel f, given as a pure function, to the input state \[Rho]. See also: ApplyUnitary, ApplyKraus."


Superoperator::usage = "Superoperator[kl] returns matrix representation of quantum channel given as a list of Kraus operators. Superoperator[fun,dim] is just am alternative name for ChannelToMatrix[fun,dim] and returns matrix representation of quantum channel, given as a pure function, acting on dim-dimensional space. So Superoperator[DepolarizingChannel[2,p,#]&,2] and Superoperator[QubitDepolarizingKraus[p]] returns the same matrix. See also: ChannelToMatrix.";


DynamicalMatrix::usage = "Dynamical matrix of quantum channel given as a list of Kraus operators (DynamicalMatrix[ch]) or as a function fun action on dim-dimensional space (DynamicalMatrix[fun,dim]). See also: Superoperator, ChannelToMatrix.";


Jamiolkowski::usage = "Jamiolkowski[K] gives the image of the Jamiolkowski isomorphism for the channel given as the list of Karus operators K. Jamiolkowski[fun,dim] gives the image of the Jamiolkowski isomorphism for the channel given as a function fun action on dim-dimensional space. See also: Superoperator, ChannelToMatrix, DynamicalMatrix.";


TPChannelQ::usage = "Performs some checks on Kraus operators. Use this if you want to check if they represent quantum channel.";


ExtendKraus::usage = "ExtendKraus[ch,n] - produces n-fold tensor products of Kraus operators from the list ch.";


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTransposeA::usage = "PartialTransposeA[\[Rho],m,n] performs partial transposition on the m-dimensional (first) subsystem of the m\[Cross]n-state.";


PartialTransposeB::usage = "PartialTransposeB[\[Rho],m,n] performs partial transposition on the n-dimensional (second) subsystem of the m\[Cross]n-state.";


PartialTraceA::usage = "PartialTraceA[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the m-demensional (first) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";


PartialTraceB::usage = "PartialTraceB[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the n-dimensional (second) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";


PartialTraceGeneral::usage = "PartialTraceGeneral[\[Rho],dim,sys] - Returns the partial trace, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}. See also: PartialTraceA, PartialTraceB.";


PartialTransposeGeneral::usage = "PartialTransposeGeneral[\[Rho],dim,sys] - Returns the partial transpose, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA,dimB}. ";


(* ::Subsection::Closed:: *)
(*Entanglement*)


Concurrence4::usage = "Concurrence4[\[Rho]] returns quantum concurrence of a density matrix \[Rho] representing a state of two-qubit system. This function uses Chop to provide numerical results.";


Negativity::usage = "Negativity[\[Rho],m,n] returns the sum of negative eigenvalues of the density matrix \[Rho]\[Element]\!\(\*SubscriptBox[\"\[DoubleStruckCapitalM]\", 
RowBox[{\"m\", \"\[Cross]\", \"n\"}]]\) after their partial transposition with respect to the first subsystem.";


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


\[Eta]::usage = "\[Eta][x] = -x Log[2,x].";


\[Eta]2::usage = "\[Eta]2[x] = \[Eta][x]+\[Eta][1-x].";


QuantumEntropy::usage = "QuantumEntropy[m] - von Neuman entropy for the matrix m.";


QuantumChannelEntropy::usage = "QuantumChannelEntropy[ch] - von Neuman entropy of the quantum channel calculated as a von Neuman entropy for the image of this channel in Jamiolkowski isomorphism. See also: Jamiolkowski, Superoperator.";


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


\[Delta]::usage = "\[Delta][a] is equivalent to \[Delta][x,''Dirac''] and it represents Dirac delta at x. If the second argument is ''Indicator'', \[Delta][x,''Indicator''] is equivalent to DiscreteDelta[x].";


VandermondeMatrix::usage = "VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}] - Vandermonde matrix for variables (\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)).";


ProdSum::usage = "ProdSum[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] gives \!\(\*SubsuperscriptBox[\"\[Product]\", 
RowBox[{\"i\", \"<\", \"j\"}], \"n\"]\)\!\(\*SubscriptBox[\"x\", \"i\"]\)+\!\(\*SubscriptBox[\"x\", \"j\"]\).";


ProdDiff2::usage = "ProdDiff2[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] is equivalent to Det[VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}]\!\(\*SuperscriptBox[\"]\", \"2\"]\) and gives a discriminant of the polynomial with roots {\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}.";


ProbBuresNorm::usage = "ProbBNorm[n] - Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Bures distance.";


ProbBures::usage = "ProbBures[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},\[Delta]] - Joint probability distribution of eigenvalues \[Lambda] = {\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}of a matrix according to Bures distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: ''Indicator''";


ProbHSNorm::usage = "Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Hilbert-Schmidt distribution.";


ProbHS::usage = "ProbHS[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},\[Delta]] Probability distribution of eigenvalues \[Lambda] = {\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)} of a matrix according to Hilbert-Schmidt distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: ''Indicator''";


(* ::Subsection::Closed:: *)
(*Random states and operations*)


RandomSimplex::usage = "RandomSimplex[d,\[Alpha]] generates a point on a d-dimensional simplex according to the Dirichlet distibution with parameter \[Alpha].\n RandomSimplex[d] uses the algorithm from the book 'Luc Devroye, Non-Uniform Random Variate Generation, Chapter 11, p. 568' and gives the flat distribution.";


RandomKet::usage = "RandomKet[d] - random ket in d-dimensional space. See: T. Radtke, S. Fritzsche, Comp. Phys. Comm., Vol. 179, No. 9, p. 647-664.";


RandomProductKet::usage = "RandomProductKet[{dim1,dim2,...,dimN}] - random pure state (ket vector) of the tensor product form with dimensions of subspaces specified dim1, dim2,...,dimN.";


RandomNormalMatrix::usage = "RandomNormalMatrix[d] - random normal matrix of dimension d.";


RandomDynamicalMatrix::usage = "RandomDynamicalMatrix[d,k] returns dynamical matrix of operation acting on d-dimensional states with k eigenvalues equal to 0. Thanks to Wojtek Bruzda.";


GinibreMatrix::usage = "GinibreMatrix[m,n] returns complex matrix of dimension m\[Cross]n with normal distribution of real and imaginary parts.";


RandomProductNumericalRange::usage = "RandomLocalNumericalRange[M,{dim1,dim2,...,dimN},n] returns n points from the product numerical range of the matrix M with respect to division specified as {dim1,dim2,...,dimN}. Note that dim1\[Cross]dim2\[Cross]...\[Cross]dimN must be equal to the dimension of matrix M.";


RandomMaximallyEntangledNumericalRange::usage = "RandomMaximallyEntangledNumericalRange[M,n] returns n points from the maximally entangled numerical range of the matrix M with respect to division Sqrt[dim[M]]\[Cross]Sqrt[dim[M]].";


RandomSpecialUnitary::usage = "Random special unitary matrix. Thanks to Rafal Demkowicz-Dobrzanski.";


RandomUnitary::usage = "Random unitary matrix. Thanks to Rafal Demkowicz-Dobrzanski.";


RandomState::usage = "RandomState[d,dist] - random density matrix of dimension d. Argument dist can be ''HS'' (default value) or ''Bures''. ''HS'' gives uniform distribution with respect to the Hilbert-Schmidt measure. ''Bures'' gives random state distributed according to Bures measure.";


(* ::Subsection::Closed:: *)
(*Random vectors*)


RandomComplexUnitVector::usage = "RandomComplexUnitVector[n] returns a normalized, n-dimensional vector of complex numbers.";


RandomRealUnitVector::usage = "RandomRealUnitVector[n] returns a normalized, n-dimensional vector of real numbers";


RandomUnitVector::usage = "RandomUnitVector[n] returns a normalized, n-dimensional vector of complex numbers. If the second argument is set to 'Real', thef unction will output a vector over \!\(\*SuperscriptBox[\"\[DoubleStruckCapitalR]\", \"n\"]\). See also: RandomKet.";


RandomEntangledUnitVector::usage = "RandomEntangledUnitVector[n] returns a maximally entangled unit vector on the n-dimensional vector space.";


RandomUnitVectorSchmidt::usage = "RandomUnitVectorSchmidt[n,r] returns a unit vector on n-dimensional space with a Schmidt rank r. Note that r has to smaller or equal \[Sqrt]n and n has to be a perfect square.";


(* ::Subsection::Closed:: *)
(*Numerical range*)


NumericalRangeBound::usage = "NumericalRangeBound[A,dx] - bound of numerical range of matrix A calculated with given step dx. Default value of dx is 0.01. Ref: Carl C. Cowen, Elad Harel, An Effective Algorithm for Computing the Numerical Range. Technical report, Dep. of Math. Purdue University, 1995.";


(* ::Subsection::Closed:: *)
(*Bloch Representation*)


BlochVector::usage = "BlochVector[A] - for a square matrix A returns a vector of coefficients obtained from expansion on normed generalized Pauli matrices. See also: GeneralizedPauliMatrices.";


StateFromBlochVector::usage = "StateFromBlochVector[v] - returns a matrix of appropriate dimension from Bloch vector, i.e. coefficients treated as coefficients from expansion on normalized generalized Pauli matrices. See also: GeneralizedPauliMatrices.";


(* ::Section::Closed:: *)
(*Private definitions*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Miscellaneous functions*)


qiGenDoc[docFile_:"qi_functions_list_alpha.tex",dir_:"~/zksi-repo/qi/doc"]:=Block[{latexHeader,latexFooter,f,txt,usage,name,lista},
	ExportString["","TeX"];
	SetDirectory[dir];

	lista=Table[{Names["QI`*"][[i]],ToExpression[Evaluate[Names["QI`*"][[i]]<>"::usage"]]},{i,1,Length[Names["QI`*"]]}];
	latexHeader="\\documentclass[a4paper,10pt]{scrartcl}
	\\usepackage{amsmath,amssymb,graphicx}
	\\usepackage{fullpage}
	\\parindent=0pt
	\\begin{document}
	\\title{QI Package for \\emph{Mathematica} 7.0 \\\\(version " <> QI`Private`qiVersion <> ")}" <>
	"\\author{Jaros{\\l}aw Adam Miszczak \\quad Piotr Gawron \\quad Zbigniew Pucha{\\l}a\\\\
	{The Institute of Theoretical and Applied Informatics}\\\\
	{Polish Academy of Sciences},\\\\
	{Ba{\\l}tycka 5, 44-100 Gliwice, Poland}}
	\\maketitle
	\\begin{abstract}"
	<> QI`Private`qiAbout <>
	"\\end{abstract}
	";

	latexFooter = "\\end{document}";
	f=OpenWrite[docFile];
	WriteString[f,latexHeader];
	For[i=1,i<= Length[lista],i++,
		name=ToString[TeXForm[lista[[i,1]]]];
		usage=ToString[TeXForm[DisplayForm[lista[[i,2]]]]];
		txt = QI`Private`qiFormatUsageMsg[name, usage];
		WriteString[f,txt];
		WriteString[f,"\\\\\n\n"];
	];
	WriteString[f,latexFooter];
	Close[f];
];

qiFormatUsageMsg[inName_,inMsg_] := Block[{name = "$ " <> inName <> " $ ",usage=inMsg,txt},
	usage=StringReplace[usage,"\{"->"LEWY"];
	usage=StringReplace[usage,"\}"->"PRAWY"];
	txt="\\textbf{"<>name<>"}"<>"-- "<>StringDrop[StringReplace[usage,RegularExpression["\\\\text{([^\}]{5,1000})}"]-> " $$1$ "],2] <> " $";
	txt=StringReplace[txt,"LEWY"-> "\{"];
	txt=StringReplace[txt,"PRAWY"->"\}"];
	txt
];


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


x_?MatrixQ \[CircleTimes] y_?MatrixQ := KroneckerProduct[x,y];
x_?VectorQ \[CircleTimes] y_?VectorQ := Flatten[KroneckerProduct[x,y]];
x_ \[CircleTimes] y_ \[CircleTimes] z_ := (x \[CircleTimes] y) \[CircleTimes] z;


KroneckerSum[A_?SquareMatrixQ,B_?SquareMatrixQ]:=Block[{dim=Length[A]},
	KroneckerProduct[A,IdentityMatrix[dim]]+KroneckerProduct[IdentityMatrix[dim],B]
];


x_ \[CirclePlus] y_ := KroneckerSum[x,y];
x_ \[CirclePlus] y_ \[CirclePlus] z_ := (x \[CirclePlus] y) \[CirclePlus] z;


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


ComplexToPoint[z_]:={Re[z],Im[z]};
SetAttributes[ComplexToPoint,Listable];


MatrixRe[A_?SquareMatrixQ]:=(A+A\[ConjugateTranspose])/2;


MatrixIm[A_?SquareMatrixQ]:=(A-A\[ConjugateTranspose])/2;


ExpectationValue[\[Rho]_,A_?SquareMatrixQ]:= If[SquareMatrixQ[\[Rho]], 
	Tr[\[Rho].A], 
	If[VectorQ[\[Rho]],
		Tr[Proj[\[Rho]].A],
		Message[ExpectationValue::argerr]
	]
];
ExpectationValue::argerr = "First argument should be a vector of a square matrix.";


Commutator[A_?SquareMatrixQ,B_?SquareMatrixQ] := If[ And@Dimensions[A]==Dimensions[B], A.B-B.A, Null];


















(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


MatrixSqrt[m_?SquareMatrixQ]:=MatrixPower[m,1/2];


MatrixAbs[a_?SquareMatrixQ]:=MatrixSqrt[a.(a\[ConjugateTranspose])];


Fidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=(Plus@@(Sqrt[Eigenvalues[a.b]]))^2;


Superfidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=Tr[a.b]+Sqrt[(1-Tr[a.a])]*Sqrt[(1-Tr[b.b])];


Subfidelity[a_?SquareMatrixQ,b_?SquareMatrixQ]:=Block[{prod = a.b}, Tr[prod] +Sqrt[2]Sqrt[(Tr[prod]*Tr[prod]-Tr[prod.prod])]];


TraceNorm[a_?SquareMatrixQ]:=Plus@@SingularValueList[a];


TraceDistance[a_?SquareMatrixQ,b_?SquareMatrixQ]:=1/2*Tr[MatrixAbs[a-b]];


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


sx = {{0,1},{1,0}};
sy = {{0,-I},{I,0}};
sz = {{1,0},{0,-1}};
\[Sigma]x = sx;
\[Sigma]y = sy;
\[Sigma]z = sz;
id = {{1,0},{0,1}};
wh = {{1,1},{1,-1}};


\[Lambda]1={{0,1,0},{1,0,0},{0,0,0}};
\[Lambda]2={{0,-I,0},{I,0,0},{0,0,0}};
\[Lambda]3={{1,0,0},{0,-1,0},{0,0,0}};
\[Lambda]4={{0,0,1},{0,0,0},{1,0,0}};
\[Lambda]5={{0,0,-I},{0,0,0},{I,0,0}};
\[Lambda]6={{0,0,0},{0,0,1},{0,1,0}};
\[Lambda]7={{0,0,0},{0,0,-I},{0,I,0}};
\[Lambda]8=Sqrt[1/3]{{1,0,0},{0,1,0},{0,0,-2}};


Proj[v_]:=Table[v[[i]]Conjugate[v[[j]]],{i,1,Length[v]},{j,1,Length[v]}];


BaseVectors[n_Integer]:=Table[UnitVector[n,k],{k,1,n}]


BaseMatrices[n_Integer]:=Table[Unres[UnitVector[n^2,k]],{k,1,n^2}];


KroneckerDeltaMatrix[m_,n_,dim_]:=Block[{mtx},
	mtx=Table[0,{dim},{dim}];
	mtx[[m,n]]=1;
	mtx
];


Lambda1[i_ ,j_,n_]:=Table[KroneckerDelta[j,\[Mu]]KroneckerDelta[i,\[Nu]] + KroneckerDelta[j,\[Nu]]KroneckerDelta[i,\[Mu]] ,{\[Mu],1,n},{\[Nu],1,n}];


Lambda2[i_ ,j_,n_]:=Table[-I(KroneckerDelta[i,\[Mu]]KroneckerDelta[j,\[Nu]] - KroneckerDelta[i,\[Nu]]KroneckerDelta[j,\[Mu]]) ,{\[Mu],1,n},{\[Nu],1,n}];


Lambda3[i_,n_]:=Sqrt[2/(i^2-i)]DiagonalMatrix[Join[Append[Table[1,{i-1}],-(i-1)],Table[0,{n-i}]]];


GeneralizedPauliMatrices[n_]:=Block[{l1,l2,l3,i,j},
	l1=Flatten[Table[Lambda1[i,j,n],{i,1,n},{j,i+1,n}],1];
	l2=Flatten[Table[Lambda2[i,j,n],{i,1,n},{j,i+1,n}],1];
	l3=Table[Lambda3[i,n],{i,2,n}];
	Join[l1,l2,l3]
];


Clear[\[Lambda]];
\[Lambda][i_,n_]:=GeneralizedPauliMatrices[n][[i]];


PauliMatrices = {sx,sy,sz};


GellMannMatrices = {\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5,\[Lambda]6,\[Lambda]7,\[Lambda]8};


(* ::Subsection::Closed:: *)
(*Quantum gates*)


Swap[dim_]:=Plus@@Flatten[Table[KroneckerProduct[Ketbra[i,j,Sqrt[dim]],Ketbra[j,i,Sqrt[dim]]],{i,0,Sqrt[dim]-1},{j,0,Sqrt[dim]-1}],1];


QFT[n_,method_:"Symbolic"]:=Block[{\[Omega]},
	If [method=="Numerical",\[Omega]=N[Exp[2 \[Pi] I/n]],\[Omega]=Exp[2 \[Pi] I/n]];
	Table[\[Omega]^(i k) ,{i,1,n},{k,1,n}]
];


cnot = {{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};


GeneralizedPauliX[d_]:=Sum[Ketbra[Mod[j-1,d],j,d],{j,0,d-1}];


GeneralizedPauliZ[d_]:=DiagonalMatrix[Table[Exp[2\[Pi] I j/d],{j,0,d-1}]];


(* ::Subsection::Closed:: *)
(*Special states*)


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


KetFromDigits[l_,b_:2]:=Ket[FromDigits[l,b],b^Length[l]];


MaxMix[n_Integer]:=1/n IdentityMatrix[n];


MaxEnt[dim_]:=Block[{subDim=Sqrt[dim]},
  If[IntegerQ[subDim],1/Sqrt[subDim] Plus@@Table[Ket[i,subDim]\[CircleTimes]Ket[i,subDim],{i,0,subDim-1}]]
];


WernerState[dim_,p_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],
	p Proj[MaxEnt[dim]] + (1-p) MaxMix[dim],
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


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


VectorSchmidtDecomposition[vec_,d1_,d2_]:=Block[{mtx, svd, vals, snum=Min[d1,d2]},
	If[VectorQ[vec],
		mtx = Partition[vec,d2];
		svd=SingularValueDecomposition[mtx];
		vals=Select[Diagonal[svd[[2]]],#!=0&];
		snum=Length[vals];
		Table[{vals[[i]],svd[[1]].BaseVectors[d1][[i]],svd[[3]]\[Conjugate].BaseVectors[d2][[i]]},{i,1,snum}],
		(*else*)
		Message[VectorSchmidtDecomposition::argerr,vec];
	]
];
VectorSchmidtDecomposition::argerr = "First argument should be a vector.";


OperatorSchmidtDecomposition[op_,d1_,d2_]:=Block[{mtx, svd, vals, snum=Min[d1*d1,d2*d2]},
	If[MatrixQ[op],
		mtx=ReshuffleGeneral[op,d1,d1,d2,d2];
		svd=SingularValueDecomposition[mtx];
		If[NumericQ[svd[[2,1]]],
			vals=Select[Diagonal[svd[[2]]],#!=0&],
			vals=Diagonal[svd[[2]]]
		];
		
		snum=Length[vals];
		Table[{vals[[i]],Unres[svd[[1]].Res[BaseMatrices[d1][[i]]]],Unres[svd[[3]]\[Conjugate].Res[BaseMatrices[d2][[i]]]]},{i,1,snum}],
		(*else*)
		Message[OperatorSchmidtDecomposition::argerr]
	]
];
OperatorSchmidtDecomposition::argerr = "First argument should be a matrix.";


SchmidtDecomposition[e_,d1_,d2_]:=Which[
	MatrixQ[e], OperatorSchmidtDecomposition[e,d1,d2],
	VectorQ[e], VectorSchmidtDecomposition[e,d1,d2],
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


Reshuffle[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
	If[dim1==0||dim2==0,
		dim=Length[\[Rho]];
		base1=BaseMatrices[Sqrt[dim]];
		Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])].Res[\[Rho]],{k,1,dim},{l,1,dim}],
		(* else *)
		base1=BaseMatrices[dim1];
		base2=BaseMatrices[dim2];
		Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])].Res[\[Rho]],{k,1,dim1 dim1},{l,1,dim2 dim2}]
	]
];


Reshuffle2[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
	If[dim1==0||dim2==0,
		dim=Length[\[Rho]];base1=BaseMatrices[Sqrt[dim]];Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])\[Transpose]].Res[\[Rho]],{l,1,dim},{k,1,dim}],
		base1=BaseMatrices[dim1];base2=BaseMatrices[dim2];Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])].Res[\[Rho]],{k,1,dim1 dim1},{l,1,dim2 dim2}]
	]
];


ReshuffleGeneral[A_,n1_,m1_,n2_,m2_]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]],{i1,0,n1 n2-1,n2},{i2,0,m1 m2-1,m2}]
,1];


ReshuffleGeneral2[A_,n1_,m1_,n2_,m2_]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]\[Transpose]],{i2,0,m1 m2-1,m2},{i1,0,n1 n2-1,n2}]
,1]\[Transpose];


MatrixElement[n_,\[Nu]_,m_,\[Mu]_,dim_,mtx_]:=mtx[[(n-1)*dim[[2]]+\[Nu],(m-1)*dim[[2]]+\[Mu]]];


ReshufflePermutation[dim1_,dim2_]:=Block[{initPos},
	initPos=Flatten[ReshuffleGeneral[Partition[Range[dim1*dim1*dim2*dim2],dim1*dim2],dim1,dim1,dim2,dim2]];
	Table[UnitVector[dim1*dim1*dim2*dim2,Position[initPos,i][[1,1]]],{i,1,dim1*dim1*dim2*dim2}]
];


ProductSuperoperator[m1_,m2_]:=Block[{dim1=Length[m1],dim2=Length[m2],perm},
	perm=ReshufflePermutation[Sqrt[dim1],Sqrt[dim2]];
	perm.(m1\[CircleTimes]m2).perm\[Transpose]
];


(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_]:=
Exp[I \[Alpha]] DiagonalMatrix[{Exp[-I \[Beta]/2],Exp[I \[Beta]/2]}].{{Cos[\[Gamma]/2],-Sin[\[Gamma]/2]},{Sin[\[Gamma]/2],Cos[\[Gamma]/2]}}.DiagonalMatrix[{Exp[-I \[Delta]/2],Exp[I \[Delta]/2]}];


SpecialUnitary2[\[Beta]_,\[Gamma]_,\[Delta]_]:=Unitary2[0,\[Beta],\[Gamma],\[Delta]];


Unitary3[al_,be_,ga_,th_,a_,b_,c_,ph_]:=MatrixExp[I *\[Lambda]3*al].MatrixExp[I*\[Lambda]2*be].
	MatrixExp[I*\[Lambda]3*ga].MatrixExp[I*\[Lambda]5*th].MatrixExp[I*\[Lambda]3*a].
	MatrixExp[I*\[Lambda]2*b].MatrixExp[I*\[Lambda]3*c].MatrixExp[I*\[Lambda]8*ph];


Unitary4Canonical[a1_,a2_,a3_]:=MatrixExp[I a1 KroneckerProduct[\[Sigma]x,\[Sigma]x]+a2 I KroneckerProduct[\[Sigma]y,\[Sigma]y]+a3 I KroneckerProduct[\[Sigma]z,\[Sigma]z]];


ProbablityVector[l_]:=Block[{ll,N},
	N=Length[l]+2;
	ll=Prepend[l,\[Pi]/2];
	Table[Sin[ll[[i-1]]]^2*Product[Cos[ll[[j-1]]]^2,{j,i+1,N}],{i,2,N}]
];


StateVector[l_]:=Block[{pr,ph,N},
	N=Length[l]/2;
	pr=ProbablityVector[l[[1;;N]]];
	ph=Prepend[Exp[I*l[[N+1;;2*N]]],1];
	FullSimplify[Sqrt[pr]*ph, Assumptions -> Table[0<l[[i]]<\[Pi]/2,{i,1,N}]]
];


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet[\[Alpha]_,\[Beta]_]:={Cos[\[Alpha]], Exp[I \[Beta]] Sin[\[Alpha]]};


QubitPureState[\[Alpha]_,\[Beta]_]:=Proj[QubitKet[\[Alpha],\[Beta]]];


QubitBlochState[a_,b_,c_]:=1/2id + a sx + b sy + c sz;


QubitGeneralState[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_,\[Lambda]_]:=Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]].DiagonalMatrix[{\[Lambda],1-\[Lambda]}].Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]]\[ConjugateTranspose];


(* ::Subsection::Closed:: *)
(*Quantum channels*)


IdentityChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]];


TransposeChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]\[Transpose]];


DepolarizingChannel=Function[{dim,p,\[Rho]},(1-p) IdentityMatrix[dim].\[Rho]+(p) Tr[\[Rho]]MaxMix[dim]];


QubitBitflipChannel=Function[{p,\[Rho]},(1-p) id.\[Rho]+p sx.\[Rho].sx\[ConjugateTranspose]];


QubitPhaseflipChannel=Function[{p,\[Rho]},(1-p) id.\[Rho]+p sz.\[Rho].sz\[ConjugateTranspose]];


QubitBitphaseflipChannel=Function[{dim,p,\[Rho]},(1-p) id.\[Rho]+p sy.\[Rho].sy\[ConjugateTranspose]];


HolevoWernerChannel=Function[{dim,p,\[Rho]},( p \[Rho]\[Transpose]+ (1-p)Tr[\[Rho]]MaxMix[dim])];


GeneralizedPauliKraus[d_,p_]:= Flatten[Table[Sqrt[ p[[i+1]][[j+1]]] (MatrixPower[GeneralizedPauliX[d],i].MatrixPower[GeneralizedPauliZ[d],j])\[ConjugateTranspose],{i,0,d-1},{j,0,d-1}],1];


ApplyKraus[ch_,\[Rho]_]:=Sum[ch[[k]].\[Rho].(ch[[k]]\[ConjugateTranspose]),{k,1,Length[ch]}];


ApplyUnitary[U_,\[Rho]_]:=U.\[Rho].U\[ConjugateTranspose];


ApplyChannel[f_,\[Rho]_] := Map[f,\[Rho],{0}];


ChannelToMatrix[fun_,dim_] := Table[
	Tr[(fun[BaseMatrices[dim][[i]]])\[ConjugateTranspose].BaseMatrices[dim][[j]]],{i,1,dim^2},{j,1,dim^2}
];


Superoperator[ch_List] := Sum[ch[[i]]\[CircleTimes]ch[[i]]\[Conjugate],{i,1,Length[ch]}];
Superoperator[fun_,dim_] := ChannelToMatrix[fun,dim];


DynamicalMatrix[ch_List] := Reshuffle[Superoperator[ch]];
DynamicalMatrix[fun_Function,dim_Integer] := Reshuffle[Superoperator[fun,dim]];


Jamiolkowski[ch_List] := 1/Length[ch[[1]]] DynamicalMatrix[ch];
Jamiolkowski[fun_Function,dim_Integer] := 1/dim DynamicalMatrix[fun,dim];


TPChannelQ[operators_] := Sum[operators[[i]]\[ConjugateTranspose].operators[[i]],{i,Length[operators]}] == IdentityMatrix[Length[operators[[1]]]];


ExtendKraus[operators_,n_] := Module[{tpl},tpl=Tuples[operators,n];Table[KroneckerProduct@@tpl[[i]],{i,1,Length[tpl]}]];


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTransposeA[\[Rho]_,m_,n_] := Reshuffle[Unres[(Swap[m m]\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho]]]]];


PartialTransposeB[\[Rho]_,m_,n_] := Reshuffle[Unres[(IdentityMatrix[m m]\[CircleTimes]Swap[n n]).Res[Reshuffle[\[Rho]]]]];


PartialTraceA[\[Rho]_,m_,n_]:=Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[m]Tr[#]&,m];
	Unres[((trMtx\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho],m,n]])[[1;;n n]]]
];


PartialTraceB[\[Rho]_,m_,n_] := Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[n]Tr[#]&,n];
	Unres[Unres[(IdentityMatrix[m m]\[CircleTimes]trMtx).Res[Reshuffle[\[Rho],m,n]]]\[Transpose][[1]]]
];


PartialTraceGeneral[\[Rho]_,dim_,sys_] := Block[{n,m},
	If[sys==1,
		Table[Sum[MatrixElement[m,\[Mu],m,\[Nu],dim,\[Rho]],{m,dim[[1]]}],{\[Mu],dim[[2]]},{\[Nu],dim[[2]]}],
		(* else *)
		Table[Sum[MatrixElement[n,\[Mu],m,\[Mu],dim,\[Rho]],{\[Mu],dim[[2]]}],{n,dim[[1]]},{m,dim[[1]]}]]
];


PartialTransposeGeneral[\[Rho]_,dim_,sys_]:=
If[sys==1,
	ArrayFlatten[Table[
		MatrixElement[n,\[Mu],m,\[Nu],dim,\[Rho]],{n,dim[[1]]},{m,dim[[1]]},{\[Nu],dim[[2]]},{\[Mu],dim[[2]]}
	]]
	,(*else*)
	ArrayFlatten[Table[
		MatrixElement[m,\[Nu],n,\[Mu],dim,\[Rho]],{n,dim[[1]]},{m,dim[[1]]},{\[Nu],dim[[2]]},{\[Mu],dim[[2]]}
	]]
];(*endif*)


(* ::Subsection::Closed:: *)
(*Entanglement*)


Concurrence4[m_]:=Block[{sqrtM=MatrixSqrt[m],evl},
	evl=Eigenvalues[MatrixSqrt[sqrtM.(sy\[CircleTimes]sy).m\[Conjugate].(sy\[CircleTimes]sy).sqrtM]];
	Max[0,Chop[Sqrt[evl[[1]]]-Sqrt[evl[[2]]]-Sqrt[evl[[3]]]-Sqrt[evl[[4]]]]]
];


Negativity[\[Rho]_,m_,n_]:=Plus@@Select[Eigenvalues[PartialTransposeA[\[Rho],m,n]],#>0&];


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


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


QubitDaviesSuperoperator[a_,c_,p_]:={{1 - a,0,0,(a p)/(1-p) },{0,c,0,0},{0,0,c,0},{a,0,0,1-(a p)/(1-p)}};


(* ::Subsection::Closed:: *)
(*One-qutrit channels*)


QutritSpontaneousEmissionKraus[A1_,A2_,t_]:={{{1,0,0},{0,Exp[-(A1 t/2)],0},{0,0,Exp[-(A2 t/2)]}},{{0,Sqrt[1-Exp[-(A1 t)]],0},{0,0,0},{0,0,0}},{{0,0,Sqrt[1-Exp[-(A2 t)]]},{0,0,0},{0,0,0}}};


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


ProdSum[l_]:=Times@@Table[Table[l[[i]] + l[[j]], {i, 1, j - 1}], {j, 2, Length[l]}];


ProdDiff2[l_]:=Block[{x},Discriminant[Times@@Table[(x-l[[i]]),{i,1,Length[l]}],x]];


ProbBuresNorm[N_]:=2^(N^2-N) Gamma[N^2/2]/(\[Pi]^(N/2) Product[Gamma[j+1],{j,1,N}]);


ProbBures[l_,delta_:"Dirac"]:=FullSimplify[ProbBuresNorm[Length[l]] \[Delta][\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)-1,delta]/\[Sqrt](\!\(
\*SubsuperscriptBox[\(\[Product]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)) Det[VandermondeMatrix[l]]^2/ProdSum[l]];


ProbHSNorm[N_]:=Gamma[N^2]/Product[Gamma[N-j] Gamma[N-j+1],{j,0,N-1}];


ProbHS[l_,delta_:"Dirac"]:=ProbHSNorm[Length[l]]\[Delta][1-(Plus@@l),delta] Det[VandermondeMatrix[l]]^2;


(* ::Subsection::Closed:: *)
(*Random states and operations*)


RandomSimplex[d_]:=Block[{r,r1,r2},
	r=Sort[Table[RandomReal[{0,1}],{i,1,d-1}]];
	r1=Append[r,1];r2=Prepend[r,0];r1-r2
];


RandomSimplex[d_,\[Alpha]_]:=Block[{gammaSample},
	gammaSample=RandomReal[GammaDistribution[\[Alpha],1],d];
	gammaSample/Plus@@gammaSample
];


RandomKet[n_]:=Block[{p,ph},
	p=Sqrt[RandomSimplex[n]];
	ph=Exp[I RandomReal[{0,2\[Pi]},n-1]];
	ph=Prepend[ph,1];
	p*ph
];


RandomProductKet[l_]:=Flatten[Apply[KroneckerProduct,Table[RandomKet[l[[i]]],{i,1,Length[l]}]]];


(* TODO: Find reference *)
RandomNormalMatrix[n_]:=Block[{DD,AA,QQ,RR},
	DD=DiagonalMatrix[RandomComplex [{-1-I,1+I},{n}]];
	AA=RandomComplex[{-1-I,1+I},{n,n}];
	{QQ,RR}=QRDecomposition[AA];
	QQ\[ConjugateTranspose].DD.QQ
];


RandomDynamicalMatrix[n_,m_:0]:=Block[{X,Y,sY},	
	X=GinibreMatrix[n^2,n^2-m];
	Y=PartialTraceGeneral[X.X\[ConjugateTranspose],{n,n},1];
	sY=MatrixPower[Y,-1/2];
	KroneckerProduct[IdentityMatrix[n],sY].X.X\[ConjugateTranspose].KroneckerProduct[IdentityMatrix[n],sY]
];


GinibreMatrix[m_,n_]:=RandomReal[NormalDistribution[0,1],{m,n}] + I RandomReal[NormalDistribution[0,1],{m,n}];


RandomProductNumericalRange[A_,sys_,noPoints_:1]:=Block[{prod},
	Table[prod=RandomProductKet[sys];Tr[Proj[prod].A],{noPoints}]
];


RandomMaximallyEntangledNumericalRange[A_,noPoints_]:=Block[{ent,dim},
	dim=Dimensions[A][[1]];
	Table[ent=RandomEntangledUnitVector[dim];ent\[Conjugate].A.ent,{noPoints}]
];


RandomSpecialUnitary[d_]:=Module[{psi,chi,r,s,phi,i,j,k,u,e,phi0,psi0,chi0},
    Do[psi[r,s]=2*Pi*Random[];,{r,1,d-1},{s,r+1,d}];
	Do[chi[r,s]=0;,{r,2,d-1},{s,r+1,d}];
	Do[chi[1,s]=2*Pi*Random[];,{s,2,d}];
	Do[phi[r,s]=ArcSin[(Random[])^(1/(2r))];,{r,1,d-1},{s,r+1,d}];
	e=Table[0,{r,1,d},{s,1,d},{i,1,d},{j,1,d}];
	Do[e[[r,s]]=IdentityMatrix[d];
    e[[r,s,r,r]]=Cos[phi0]*Exp[I*psi0];
    e[[r,s,s,s]]=Cos[phi0]*Exp[-I*psi0];
    e[[r,s,r,s]]=Sin[phi0]*Exp[I*chi0];
    e[[r,s,s,r]]=-Sin[phi0]*Exp[-I*chi0];,{r,1,d-1},{s,r+1,d}];
    u=IdentityMatrix[d];
    Do[u=(e[[r,r+1]] /. {phi0->phi[d-r,s+1],psi0->psi[d-r,s+1],chi0->chi[d-r,s+1]}).u;,{s,d-1,1,-1},{r,d-1,d-s,-1}];
    u
];


RandomUnitary[d_]:=Exp[I*RandomReal[2*\[Pi]]]*RandomSpecialUnitary[d];


RandomState[d_,dist_:"HS"]:=Block[{v,A,U},
	Switch[dist,
		"HS",
			v=Flatten[GinibreMatrix[d^2,1]];
			v=v/Norm[v,2];
			PartialTraceGeneral[Proj[v],{d,d},2]//Chop,
		"Bures",
			A=GinibreMatrix[d,d];
			U=RandomUnitary[d];
			(IdentityMatrix[d]+U).A.A\[ConjugateTranspose].(IdentityMatrix[d]+U)\[ConjugateTranspose]/Tr[(IdentityMatrix[d]+U).A.A\[ConjugateTranspose].(IdentityMatrix[d]+U)\[ConjugateTranspose]] //Chop,
		_, 
			Message[RandomState::argerr,dist]
	]

];
RandomState::argerr = "The second argument should be \"HS\" or \"Bures\", mesure \"`1`\" not implemented yet.";


(* ::Subsection::Closed:: *)
(*Random vectors*)


RandomComplexUnitVector[n_Integer]:=Block[{rv},
	rv=RandomReal[NormalDistribution[0,1],{n}]+I RandomReal[NormalDistribution[0,1],{n}];
	rv/Norm[rv]
];


RandomRealUnitVector[n_Integer]:=Block[{rv},
	rv = RandomReal[NormalDistribution[0,1],{n}];
	rv/Norm[rv]
];


RandomUnitVector[n_Integer,type_:"Complex"] :=Switch[type,
	"Complex", RandomComplexUnitVector[n],
	"Real", RandomRealUnitVector[n],
	_, RandomComplexUnitVector[n]
];


RandomEntangledUnitVector[n_Integer]:=Block[{u1,u2,d=Sqrt[n]},
	u1=RandomSpecialUnitary[d];
	u2=RandomSpecialUnitary[d];
	Plus@@Table[(u1[[i]])\[CircleTimes](u2[[i]]),{i,1,d}]/Sqrt[d]
];


RandomUnitVectorSchmidt[n_,r_]:=Block[{u1,u2,d,v},
	d=Sqrt[n];
	If[IntegerQ[d],
		u1=RandomSpecialUnitary[d];
		u2=RandomSpecialUnitary[d];
		v=Plus@@Table[(u1[[i]])\[CircleTimes](u2[[i]]),{i,1,r}];
		v/Norm[v],
		(* else *)
		Message[RandomUnitVectorSchmidt::argerr,n]
	]
];
RandomUnitVectorSchmidt::argerr = "`1` is not a perfect square!";


(* ::Subsection::Closed:: *)
(*Numerical range*)


NumericalRangeBound[A_?MatrixQ,step_:0.01]:=Block[{w,Ath,Hth,m,s,Kth,pKp,ee,rr,mm,sm,mM,sM,e,r},
	w={};
	Table[
	Ath=Exp[I*(-\[Theta])]*A;
	Hth=MatrixRe[Ath];
	{e,r}=Eigensystem[Hth];
	e=Re[e];
	m=Max[e];
	s=Position[e,m];
	If[
	Length[s]==1,(*then*)AppendTo[w,ArrayFlatten[Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose]]],
	(*else*)
	Kth=I*(Hth-Ath); pKp=Extract[r,s]\[Conjugate].Kth.Extract[r,s]\[Transpose]; {ee,rr}=Eigensystem[pKp]; ee=Re[ee]; mm=Min[ee]; sm=Position[ee,mm];
	AppendTo[w,ArrayFlatten[Extract[rr,sm]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sm]\[Transpose]]];
	mM=Max[ee];sM=Position[ee,mM];AppendTo[w,ArrayFlatten[Extract[rr,sM]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sM]\[Transpose]]]
	(*end if*)
	]
	,{\[Theta],0,2\[Pi],step}];
Flatten[w,2]
];


(* ::Subsection::Closed:: *)
(*Bloch Representation*)


BlochVector[A_]:=Block[{dim},
  dim=Length[A]; 1/Sqrt[2](Tr[A\[ConjugateTranspose].#]&/@GeneralizedPauliMatrices[dim])
];


StateFromBlochVector[vec_]:=Block[{dim},
	If[IntegerQ[Sqrt[Length[vec]+1]],
		dim= Sqrt[Length[vec]+1];
		1/dim IdentityMatrix[dim] + vec.GeneralizedPauliMatrices[dim]/Sqrt[2],
		Message[StateFromBlochVector::argerr, vec];
		Beep[];
	]
];
StateFromBlochVector::argerr= "Given vector (`1`) is not a Bloch vector of any dimension.";


End[];


(* ::Section::Closed:: *)
(*Package footer*)


Protect@@QI`Private`qiNames;


EndPackage[];
