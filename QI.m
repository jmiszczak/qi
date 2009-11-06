(* ::Package:: *)

(* ::Section:: *)
(*Package header*)


BeginPackage["QI`"];


(* File: QI.m *)
(* Description: Mathematica package for the analysis of quantum states and operations *)
(* Authors: Jaroslaw Miszczak <miszczak@iitis.pl>, Piotr Gawron <gawron@iiti.pl>, Zbigniew Puchala <z.puchala@iitis.pl> *)
(* License: GPLv3 *)
qiVersion = "0.3.0";
qiLastModification = "November 6, 2009";
qiHistory = {
	{"0.1.0", "Initial version"}, 
	{"0.1.1", "Fixed \[Eta] and \[Eta]2 functions, fixed problem with protected symbols."}, 
	{"0.1.2", "Added quantum channel parametrization for one qubit"}, 
	{"0.1.3", "Added alternative reshuffling."},
	{"0.1.4", "Changed default print output."}, 
	{"0.2.0", "Documentation generator added."},
	{"0.2.1", "Changed QubitGeneralState function."},
	{"0.2.2", "Added reshuffling permutation and product of superoperators."},
	{"0.2.3", "Minor update in documentation."},
	{"0.2.4", "Fixed \[Eta] function."},
	{"0.2.5", "Fixed Werner state and added IsotropicState."},
	{"0.2.6", "Spelling improvements"},
	{"0.2.7", "Improved \[CircleTimes] usage , added applyUnitary"},
	{"0.2.8", "Fixed small problem with MaxEnt"},
	{"0.2.9", "Some code cleanups, fixed SchmidtDecomposition"},
	{"0.3.0", "Added OperatorSchmidtDecomposition"}
};
qiAbout ="QI is a package of functions for Mathematica computer algebra system, which implements 
number of functions used in the analysis of quantum states. In contrast to many available 
packages for symbolic and numerical simulation of quantum computation presented package is focused 
on geometrical aspects of quantum information theory.";
qiGenDoc::usage = "Generate documentation for QI package. This function accepts up to two arguments - file name and output directory. By default the first one is set to 'qi_usage.tex' and the second one to '~/zksi-repo/qi/doc'.";
qiAbout::usage = "Display short information about QI package.";
qiVersion::usage = "Display version of QI package.";
qiHistory::usage = "Display history of modifications for QI package.";
Print["Package QI version ", qiVersion, " (last modification: ", qiLastModification, ")."];


$PrePrint = If[MatrixQ[#], MatrixForm[#], #]&;


(* ::Section::Closed:: *)
(*Help messages*)


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


KroneckerSum::usage = "KroneckerSum[A,B] returns the Kronecker sum of matrices A and B defined as A\[CircleTimes]1+1\[CircleTimes]B. Alternative syntax A\[CirclePlus]B for KroneckerSum[A,B] is provided.";


SymbolicMatrix::usage = "SymbolicMatrix[a,m,n] returns m\[Cross]n-matrix with elements a[i,j], i=1,...,m, j=1,...,n. If the second argument is ommited this function returns square n\[Cross]n matrix. This functions can save you some keystrokes and, thanks to TeXForm function, its results can be easily incorporated in LaTeX documents.";


SymbolicVector::usage = "SymbolicVector[a,m] is equivalent to Matrix[a,m,1] and it returns a vector with m elements a[i],i=1,...,m, j=1,...,n. This function is useful, for example, for generating lists of parameters.";


SymbolicHermitianMatrix::usage ="SymbolicHermitianMatrix[sym,d] produces d\[Cross]d hermitian matrix.";


ComplexToPoint::usage = "ComplexToPoint[z] returns real and imaginary parts of a complex number z as a pair of real numbers (point on real plane)).";


(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


MatrixSqrt::usage= "MatrixSqrt[M] returns square root for the matrix M.";


MatrixAbs::usage= "MatrixAbs[M] returns absolute value for matrix M defined as MatrixSqrt[M.M\!\(\*SuperscriptBox[\" \", \"\[Dagger]\"]\)].";


Fidelity::usage = "Fidelity[\[Rho],\[Sigma]] returns quantum fidelity between states \[Rho] and \[Sigma] calculated as tr\[Sqrt](\[Sqrt]\[Rho] \[Sigma]\[Sqrt]\[Rho]).";


Superfidelity::usage = "Superfidelity[A,B] calculates superfidelity between A and B defined as Tr[A.B] + Sqrt[1-Tr[A.A]]Sqrt[1-Tr[B.B]]."; 


Subfidelity::usage = "Subfidelity[A,B] returns subfidelity between states A and B calculated as \!\(\*SubscriptBox[\"tr\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)+\[Sqrt]2\[Sqrt]((\!\(\*SubscriptBox[\"tr\[Rho]\", \"2\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\))-\!\(\*SubscriptBox[\"tr\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)).";


TraceDistance::usage = "TraceDistance[A,B] returns trace distance between matrices A and B which is defined as \!\(\*FractionBox[\"1\", \"2\"]\)tr|A-B|.";


MatrixRe::usage = "Hermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A+ConjugateTranspose[A]).";


MatrixIm::usage = "Antyhermitian part of the matrix A i.e. \!\(\*FractionBox[\"1\", \"2\"]\)(A-ConjugateTranspose[A]).";


ExpectationValue::usage = "ExpectationValue[\[Rho],A] is just a alias for Tr[\[Rho].A].";


Commutator::usage = "Commutator[A,B] returns the commutator of matrices A and B i.e. [A,B]=AB - BA.";


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


sx::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\).";
sy::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\).";
sz::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\).";
\[Sigma]x::usage = sx::usage;
\[Sigma]y::usage = sy::usage;
\[Sigma]z::usage = sz::usage;
id::usage = "Identity matrix for one qubit. See also: IdentityMatrix";
wh::usage = "Hadamard gate for one qubit. See also: QFT.";


\[Lambda]1::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\).";
\[Lambda]2::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\).";
\[Lambda]3::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\).";
\[Lambda]4::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"4\"]\).";
\[Lambda]5::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"5\"]\).";
\[Lambda]6::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"6\"]\).";
\[Lambda]7::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"7\"]\).";
\[Lambda]8::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"8\"]\).";


Proj::usage = "Proj[{\!\(\*SubscriptBox[\"v\", \"1\"]\),\!\(\*SubscriptBox[\"v\", \"2\"]\),...,\!\(\*SubscriptBox[\"v\", \"n\"]\)}] returns projectors for the vectors in the input list.";


BaseVectors::usage = "BaseVectors[n] returns list with canonical basis in n-dimensional Hilbert space.";


BaseMatrices::usage = "BaseMatrices[n] returns list with canonical basis in n\[Cross]n-dimensional Hilbert-Schmidt space.";


KroneckerDeltaMatrix::usage = "KroneckerDeltaMatrix[i,j,d] returns d\[Cross]d matrix with 1 at position (i,j) and zeros elsewhere.";


Lambda1::usage = "Lambda1[i,j,n] generalized Pauli matrix. For example Lambda1[1,2,2] is equal to Pauli \[Sigma]x.";


Lambda2::usage = "Lambda2[i,j,n] generalized Pauli matrix. For example Lambda2[1,2,2] is equal to \[Sigma]y.";


Lambda3::usage ="Lambda3[i,n] generalized Pauli matrix. For example Lambda3[2,2] is equal to \[Sigma]z.";


GeneralizedPauliMatrices::usage = "GeneralizedPauliMatrices[n] returns list of generalized Pauli matrices for SU(n). For n=2 these are just Pauli matrices and for n=3 - Gell-Mann matrices. Note that identity matrix is not included in the list. See also: PauliMatrices, GellMannMatrices, \[Lambda], Lambda1, Lambda2, Lambda3.";


\[Lambda]::usage = "\[Lambda][i,n] is defined as GeneralizedPauliMatrices[n][[i]].";


PauliMatrices::usage = "Predefined list of Pauli matrices {\!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\)}. Use Map[MatrixForm[#]&,PauliMatrices] to get this list in more readible form.";


GellMannMatrices::usage = "List of Gell-Mann matrices. Use Map[MatrixForm[#]&,GellMannMatrices] to get this list in more readable form.";


(* ::Subsection::Closed:: *)
(*Quantum gates*)


Swap::usage="Swap[n] returns permutation operator \!\(\*UnderoverscriptBox[\"\[Sum]\", 
RowBox[{\"i\", \"=\", \"0\"}], 
RowBox[{\"n\", \"-\", \"1\"}]]\)\!\(\*UnderoverscriptBox[
RowBox[{\" \", \"\[Sum]\"}], 
RowBox[{\"j\", \"=\", \"0\"}], 
RowBox[{\"n\", \"-\", \"1\"}]]\)|i\[RightAngleBracket]\[LeftAngleBracket]j|\[CircleTimes]|j\[RightAngleBracket]\[LeftAngleBracket]i\[VerticalSeparator] acting on \!\(\*SuperscriptBox[\"n\", \"2\"]\)-dimensional space and exchanging two n-dimensional subsystems.";


QFT::usage = "QFT[n,method] - quantum Fourier transform of dimension n. This function accepts second optional argument, which specifies method used in calculation. Parameter method can be equal to 'Symbloic', which is default, or 'Numerical'. The second option makes this function much faster.";


CNot::usage = "Controlled not matrix for two qubits.";


H::usage = "Hadamard gate for one qubit.";


X::usage = "Generalized Pauli matrix X.";


Z::usage = "Generalized Pauli matrix Z.";


Circuit::usage = "Constructs quantum circuit from unitary gates.";


(* ::Subsection::Closed:: *)
(*Special states*)


Ket::usage = "Ket[i,d] returns |i\[RightAngleBracket] in d-dimensional Hilbert space. See also: StateVector for a different parametrization.";


Ketbra::usage = "Ketbra[i,j,d] returns \[VerticalSeparator]i\[RightAngleBracket]\[LeftAngleBracket]j\[VerticalSeparator] acting on d-dimensional space.";


KetFromDigits::usage = "KetFromDigits[list,base] - returns ket vector labeled by string of digits represented in given base.";


MaxMix::usage = "MaxMix[n] gies maximally mixed state in n-dimensional space of density matrices.";


MaxEnt::usage = "MaxEnt[N] - maximally entangled state in N dimensional vector space. Note that N must be perfect square.";


WernerState::"usage"="WernerState[d,p] - generalized Werner state d\[Cross]d-dimensional system with the mixing parameter p\[Element][0,1]. This state is defined as p Proj[MaxEnt[d]] + (1-p) MaxMix[d],. See also: MaxEnt, MaxMix..";


IsotropicState::usage = "IsotropicState[d,p] - isotropic state of dimensions d\[Cross]d with parameter p\[Element][0,1] defined as p Proj[MaxEnt[d]] + (1-p)/(\!\(\*SuperscriptBox[\"d\", \"2\"]\)-1)(I-Proj[MaxEnt[d]]]). This family of states is invariant for the operation of the form U\[CircleTimes]\!\(\*SuperscriptBox[\"U\", \"\[Star]\"]\).";


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


SchmidtDecomposition::usage = "SchmidtDecomposition[vec,d1,d2] - Schmidt decomposition of the vector vec in d1\[Cross]d2-dimensional Hilbert space.";


OperatorSchmidtDecomposition::usage = "OperatorSchmidtDecomposition[mtx,d1,d2] - Schmidt decomposition of mtx in the Hilbert-Schmidt space of matrices of dimension d1\[Cross]d2.";


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec::usage = "Vec[m] - vectorization of the matrix m column by column. See also: Res.";


Unvec::usage = "Unvec[v,c] - de-vectorization of the vector into the matrix with c columns. If the second parameter is omitted then it is assumed that v can be mapped into square matrix. See also: Unres, Vec.";


Res::usage = "Res[m] is equivalent to Vec[Transpose[m]]. Reshaping maps matrix m into vector row by row.";


Unres::usage = "de-reshaping of the vector into the matrix with c columns. If the second parameter is omitted then it is assumed that v can be mapped into square matrix. See also: Unvec, Res.";


Reshuffle::usage = "Reshuffle[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices. If  the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See also: ReshuffleGeneral, Reshuffle2.";


Reshuffle2::usage = "Alternative definition of the reshuffling operation. Reshuffle2[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices which are transposed versions of standard base matrices. If the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See: See also: ReshuffleGeneral, Reshuffle, BaseMatrices";


ReshuffleGeneral::usage = "ReshuffleGeneral[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix.";


ReshuffleGeneral2::usage = "ReshuffleGeneral2[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix - given by alternative definition of the reshuffling operation.";


MatrixElement::usage = "MatrixElement[n,\[Nu],m,\[Mu],div,M] - returns the matrix element of density matrix M indexed by two double indices n, \[Nu] and m, \[Mu] of the composite sytem of dimensions dim=dimA dimB.";


ReshufflePermutation::usage = "ReshufflePermutation[dim1,dim2] produces permutation matrix equivalent to the reshuffling operation on dim1\[Cross]dim2-dimensional system.";


ProductSuperoperator::usage = "ProductSuperoperator[m1,m2] computes product superoperator of superoperatos m1 and m2.";


(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2::usage = "Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]] returns the Euler parametrization of U(2).";


SpecialUnitary2::usage ="SpecialUnitary2[\[Beta],\[Gamma],\[Delta]] returns the Euler parametrization of SU(2). This is equivalent to Unitary2[0,\[Beta],\[Gamma],\[Delta]].";


Unitary3::usage = "Unitary3[\[Alpha],\[Beta],\[Gamma],\[Tau],a,b,c,ph] returns the Euler parametrization of U(3).";


Unitary4Canonical::usage = "Parametrization of non-local unitary matrices for two qubits. See: arXiv:quant-ph/0011050v1.";


ProbablityDistribution::usage = "ProbablityDistribution[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}] returns probability vectors of dimension n+1 parametrize with {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}. See also: StateVector.";


StateVector::usage = "StateVector[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\),\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"n\", \"+\", \"1\"}]]\),...,\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"2\", \" \", \"n\"}]]\)}] returns pure n+1-dimensional pure state (ket vector) constructed form probability distribution parametrize by numbers {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)} and phases {\!\(\*SubscriptBox[\"\[Phi]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \"n\"]\)}. See also: ProbablityDistribution, SymbolicVector.";


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet::usage = "QubitKet[\[Alpha],\[Beta]] parametrisation of the pure state (as a state vector) for one qubit as (Cos[\[Alpha]] Exp[i\[Beta]], Sin[\[Alpha]]). This is equivalent to StateVector[{\[Alpha],\[Beta]}]. See also: QubitPureState, StateVector.";


QubitPureState::usage = "QubitPureState[\[Alpha],\[Beta]] - parametrisation of the pure state as a density matrix for one qubit. This is just a alias for Proj[QubitKet[\[Alpha],\[Beta]]]. See also: QubitKet.";


QubitBlochState::usage = "Parametrization of the one-qubit mixed state on the Bloch sphere.";


QubitGeneralState::usage = "QubitGeneralState[\[Alpha],\[Beta],\[Gamma],\[Delta],\[Lambda]] - Parametrization of the one-qubit mixed state using rotations and eigenvalues. Returns one-qubits density matrix with eigenvalues \[Lambda] and 1-\[Lambda] rotated as U.diag(\[Lambda],1-\[Lambda]).\!\(\*SuperscriptBox[\"U\", \"\[Dagger]\"]\) with U defined by parameters \[Alpha],\[Beta],\[Gamma] and \[Delta].";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


IdentityChannel::usage="IdentityChannel[n, \[Rho]] - apply identity operation on n-dimensional density matrix \[Rho].";


TransposeChannel::usage="TransposeChannel[n, \[Rho]] - apply transposition operation on n-dimensional density matrix \[Rho]. This operations is not completely positive.";


DepolarizingChannel::usage = "DepolarizingChannel[n,p,\[Rho]] performs an action of the completely depolarizing channel with parameter p acting on n-dimensional input state \[Rho]. See also: QubitDepolarizingKraus, HolevoWernerChannel.";


BitflipChannel::usage  = "BitflipChannel[2,p,\[Rho]].";


PhaseflipChannel::usage  = "PhaseflipChannel[2,p,\[Rho]]";


BitphaseflipChannel::usage  = "BitphaseflipChannel[2,p,\[Rho]].";


HolevoWernerChannel::usage = "HolevoWernerChannel[n,p,\[Rho]] performs an action of the Holeve-Werner channel (also known as transpose depolarizing channel) with parameter p acting on n-dimensional input state \[Rho]. See also: DepolarizingChannel.";


ChannelToMatrix::usage = "ChannelToMatrix[E,d] returns matrix representation of a channel E acting on d-dimensional state space. First argument should be a pure function E such that E[\[Rho]] transforms input state according to the channel definition. For example for the Holevo-Werner channel one ca use ChannelToMatrix[HolevoWernerChannel[3,p,#]&,3] to obtain matrix representation of this channel acting on qutrits. See also: Superoperator.";


GeneralizedPauliKraus::usage = "GeneralizedPauliKraus[d,P] - list of Kraus operators for d-dimensional generalized Pauli channel with the d-dimesnional matrix of parameters P. See: M. Hayashi, Quantum Information An Introduction, Springer 2006, Example 5.8, p. 126.";


ApplyKraus::usage = "ApplyKraus[ch,\[Rho]] - apply channel ch, given as a list of Kraus operators, to the input state \[Rho].";


ApplyUnitary::usage="ApplyUnitary[U,\[Rho]] - apply unitary u to the input state \[Rho].";


Superoperator::usage = "Superoperator[kl] returns matrix representation of quantum channel given as a list of Kraus operators. Superoperator[fun,dim] is just am alternative name for ChannelToMatrix[fun,dim] and returns matrix representation of quantum channel, given as a pure function, acting on dim-dimensional space. So Superoperator[DepolarizingChannel[2,p,#]&,2] and Superoperator[QubitDepolarizingKraus[p]] returns the same matrix. See also: ChannelToMatrix.";


DynamicalMatrix::usage = "Dynamical matrix of quantum channel given as a list of Kraus operators (DynamicalMatrix[ch]) or as a function fun action on dim-dimensional space (DynamicalMatrix[fun,dim]). See alos: Superoperator, ChannelToMatrix.";


Jamiolkowski::usage = "Jamiolkowski[K] gives the image of the Jamiolkowski isomorphism for the channel given as the list of Karus operators K. Jamiolkowski[fun,dim] gives the image of the Jamiolkowski isomorphism for the channel given as a function fun action on dim-dimensional space. See alos: Superoperator, ChannelToMatrix, DynamicalMatrix.";


ChannelCondition::usage = "Performs some checks on Kraus operators. Use this if you want to check if they represent quantum channel.";


ExtendKraus::usage = "ExtendKraus[ch,n] - produces n-fold tensor products of Kraus operators from the list ch.";


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTransposeA::usage = "PartialTransposeA[\[Rho],m,n] performs partial transposition on the m-dimensional (first) subsystem of the m\[Cross]n-state.";


PartialTransposeB::usage = "PartialTransposeB[\[Rho],m,n] performs partial transposition on the n-dimensional (second) subsystem of the m\[Cross]n-state.";


PartialTraceA::usage = "PartialTraceA[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the m-demensional (first) subsystem.";


PartialTraceB::usage = "PartialTraceB[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the n-demensional (second) subsystem.";


PartialTraceGeneral::usage = "PartialTraceGeneral[\[Rho],dim,sys] - Returns the partial trace, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}. See alos: PartialTraceA, PartialTraceB.";


PartialTransposeGeneral::usage = "PartialTransposeGeneral[\[Rho],dim,sys] - Returns the partial transpose, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA,dimB}. ";


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


QubitDepolarizingKraus::usage = "Kraus operators of the depolarizing channel for one qubit. Note that it gives maximally mixed state for p=0.";


QubitDecayKraus::usage = "Kraus operators of the decay channel, also know as amplitude damping, for one qubit.";


QubitPhaseKraus::usage = "Kraus operators for one qubit phase damping channel.";


QubitBitflipKraus::usage = "Kraus operators for one qubit bit-flip channel.";


QubitPhaseflipKraus::usage = "Kraus operators for one qubit phase-flip channel.";


QubitBitphaseflipKraus::usage = "Kraus operators for one qubit bit-phase-flip channel.";


QubitDynamicalMatrix::usage = "QubitDynamicalMatrix[\!\(\*SubscriptBox[\"\[Kappa]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Kappa]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Kappa]\", \"z\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"x\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"y\"]\),\!\(\*SubscriptBox[\"\[Eta]\", \"z\"]\)] returns parametrization of one-qubit dynamical matrix. See: BZ Chapter 10, formula 10.81";


QubitDaviesDynamicalMatrix::usage = "Returns dynamical matrix for Davies channel with b = \!\(\*FractionBox[
RowBox[{\"a\", \" \", \"p\"}], 
RowBox[{\"1\", \"-\", \"p\"}]]\).";


(* ::Subsection::Closed:: *)
(*One-qutrit channels*)


QutritSpontaneousEmissionKraus::usage="QutritSpontaneousEmissionKraus[A1,A2,t] Kraus operators for qutrit epontaneous emission channel with parameters A1, A2, t >= 0, see 
A. Checinska, K. Wodkiewicz, Noisy Qutrit Channels, arXiv:quant-ph/0610127v2.";


(* ::Subsection::Closed:: *)
(*Entropies*)


Log0::usage = "Log0[x] is equal to Log[2,x] for x>0 and 1 for x=0.";


\[Eta]::usage = "\[Eta][x] = -x Log[2,x].";


\[Eta]2::usage = "\[Eta]2[x] = \[Eta][x]+\[Eta][1-x].";


QuantumEntropy::usage = "QuantumEntropy[m] - von Neuman entropy for the matrix m.";


S::usage = "S[m] = QuantumEntropy[m].";


QuantumChannelEntropy::usage = "QuantumChannelEntropy[ch] - von Neuman entropy of the quantum channel calculated as a von Neuman entropy for the image of this channel in Jamiolkowski isomorphism. See also: Jamiolkowski, Superoperator.";


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


\[Delta]::usage = "\[Delta][x] represents Dirac delta at x.";


VandermondeMatrix::usage = "VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}] - Vandermonde matrix for variables (\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)).";


ProdSum::usage = "ProdSum[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] gives \!\(\*SubsuperscriptBox[\"\[Product]\", 
RowBox[{\"i\", \"<\", \"j\"}], \"n\"]\)\!\(\*SubscriptBox[\"x\", \"i\"]\)+\!\(\*SubscriptBox[\"x\", \"j\"]\).";


ProdDiff2::usage = "ProdDiff2[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] is equivalent to Det[VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}]\!\(\*SuperscriptBox[\"]\", \"2\"]\) and gives a discriminant of the polynomial with roots {\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}.";


ProbBuresNorm::usage = "ProbBNorm[n] - Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Bures distance.";


ProbBures::usage = "ProbBures[\[Lambda]] - Joint probability distribution of eigenvalues \[Lambda] of a matrix according to Bures distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: \"Indicator\"";


ProbHSNorm::usage = "Normalization factor used for calculating probability distribution of eigenvalues of matrix of dimension N according to Hilbert-Schmidt distribution.";


ProbHS::usage = "ProbHS[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},] Probability distribution of eigenvalues of matrix according to Hilbert-Schmidt distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: \"Indicator\"";


(* ::Subsection::Closed:: *)
(*Random states and operations*)


RandomSimplex::usage = "RandomSimplex[d] - d-dimensional random simplex.";


RandomKet::usage = "RandomKet[d] - random ket in d-dimensional space. See: T. Radtke, S. Fritzsche / Computer Physics Communications 179 (2008) 647-664.";


RandomProductKet::usage = "RandomProductKet[{dim1,dim2,...,dimN}] - random pure state (ket vector) of the tensor product form with dimensions of subspaces specified dim1, dim2,...,dimN.";


RandomNormalMatrix::usage = "RandomNormalMatrix[d] - random normal matrix of dimension d.";


RandomDynamicalMatrix::usage = "RandomDynamicalMatrix[d,k] returns dynamical matrix of operation acting on d-dimensional states with k eigenvalues equalt to 0.";


GinibreMatrix::usage = "GinibreMatrix[m,n] returns complex matrix of dimension m\[Cross]n with normal distribution of real and imaginary parts.";


RandomProductNumericalRange::usage = "RandomLocalNumericalRange[M,{dim1,dim2,...,dimN},n] returns n points from the product numerical range of the matrix M with respect to division specified as {dim1,dim2,...,dimN}. Note that dim1\[Cross]dim2\[Cross]...\[Cross]dimN must be equal to the dimension of matrix M.";


RandomSpecialUnitary::usage = "Random spacial unitary matrix. Thanks to R.D-D.";


RandomUnitary::usage = "Random unitary matrix. Thanks to R.D-D.";


RandomState::usage = "RandomState[d] - random density matrix of dimension d. This gives uniform distribution with respect to Hilbert-Schmidt measure.";


(* ::Subsection::Closed:: *)
(*Numerical range*)


NumericalRangeBound::usage = "NumericalRangeBound[A_?MatrixQ,step_:0.01] - bound of numerical range of matrix A calculated with given step. Ref: Carl C. Cowen, Elad Harel, An Effective Algorithm for Computing the Numerical Range. Technical report, Dep. of Math. Purdue University, 1995.";


(* ::Subsection::Closed:: *)
(*Bloch Representation*)


BlochVector::usage = "BlochVector[A_MatrixQ] - for square matrix - vector of coefficients obtained from expansion on normed generalized pauli matrices, see function GeneralizedPauliMatrices"


StateFromBlochVector::usage = "StateFromBlochVector[vec_] - returns a matrix of appropriate dimension from Bloch vector (coefficients treated as coefficients from expansion on normed generalized Pauli matrices, see function GeneralizedPauliMatrices)"


(* ::Section::Closed:: *)
(*Private definitions*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Misc*)


qiGenDoc[docFile_:"qi_usage.tex",dir_:"~/zksi-repo/qi/doc"]:=Block[{latexHeader,latexFooter,f,txt,usage,name,lista},
ExportString["","TeX"];
SetDirectory[dir];
lista=Table[{Names["QI`*"][[i]],ToExpression[Evaluate[Names["QI`*"][[i]]<>"::usage"]]},{i,1,Length[Names["QI`*"]]}];
latexHeader="\\documentclass[a4paper,10pt]{scrartcl}
\\usepackage{amsmath,amssymb,graphicx}
\\usepackage{fullpage}
\\parindent=0pt
\\begin{document}
\\title{QI Package for \\emph{Mathematica} 7.0 (version " <> qiVersion <> ")}" <>
"\\author{Jaros{\\l}aw Miszczak, Piotr Gawron, Zbigniew Pucha{\\l}a\\\\ \\small{The Institute of Theoretical and Applied Informatics}\\\\
\\small{Polish Academy of Sciences},\\\\ \\small{Ba{\\l}tycka 5, 44-100 Gliwice, Poland}}
\\maketitle
\\begin{abstract}"
<> qiAbout <>
"\\end{abstract}
";
latexFooter = "\\end{document}";
f=OpenWrite[docFile];
WriteString[f,latexHeader];
For[i=1,i<= Length[lista],i++,
name=ToString[TeXForm[lista[[i,1]]]];
(*
name=StringReplace[name,RegularExpression["\\\\text{(.*?)}"]->"$1"];
name=StringReplace[name,RegularExpression["^\\\\delta"]->"$\\delta$"];
name=StringReplace[name,RegularExpression["^\\\\lambda"]->"$\\lambda$"];
name=StringReplace[name,RegularExpression["^\\\\eta"]->"$\\eta$"];
*)
name = "$ " <> name <> " $ ";
usage=ToString[TeXForm[DisplayForm[lista[[i,2]]]]];
usage=StringReplace[usage,"\{"->"LEWY"];
usage=StringReplace[usage,"\}"->"PRAWY"];
txt="\\textbf{"<>name<>"}"<>"-- "<>StringDrop[StringReplace[usage,RegularExpression["\\\\text{([^\}]{5,1000})}"]-> " $$1$ "],2] <> " $";
txt=StringReplace[txt,"LEWY"-> "\{"];
txt=StringReplace[txt,"PRAWY"->"\}"];
WriteString[f,txt];
WriteString[f,"\\\\\n\n"];
];
WriteString[f,latexFooter];
Close[f];
];


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


x_?MatrixQ \[CircleTimes] y_?MatrixQ := KroneckerProduct[x,y];
x_?VectorQ \[CircleTimes] y_?VectorQ := Flatten[KroneckerProduct[x,y]];
x_ \[CircleTimes] y_ \[CircleTimes] z_ := (x \[CircleTimes] y) \[CircleTimes] z;


Clear[KroneckerSum];
KroneckerSum[A_,B_]:=Block[{dim=Length[A]},
	KroneckerProduct[A,IdentityMatrix[dim]]+KroneckerProduct[IdentityMatrix[dim],B]
];


x_ \[CirclePlus] y_ := KroneckerSum[x,y];
x_ \[CirclePlus] y_ \[CirclePlus] z_ := (x \[CirclePlus] y) \[CirclePlus] z;


Clear[SymbolicMatrix];
SymbolicMatrix[sym_,d1_,d2_:0]:= Which[
	d2==0,Table[Subscript[sym, i,j],{i,1,d1},{j,1,d1}],
	d2==1,Table[Subscript[sym, i],{i,1,d1}],
	True,Table[Subscript[sym, i,j],{i,1,d1},{j,1,d2}] 
];


Clear[SymbolicVector];
SymbolicVector[sym_,d1_]:= SymbolicMatrix[sym,d1,1];


Clear[SymbolicHermitianMatrix];
SymbolicHermitianMatrix[sym_,d_]:=Block[{mtx},
	mtx=Table[0,{d},{d}];
	Table[mtx[[i,j]]=Subscript[sym, i,j],{i,1,d},{j,1,i}];
	mtx=mtx+mtx\[ConjugateTranspose];
	Table[mtx[[i,i]]=Subscript[sym, i,i],{i,1,d}];
	mtx
];


Clear[ComplexToPoint];
ComplexToPoint[z_]:={Re[z],Im[z]};
SetAttributes[ComplexToPoint,Listable];


(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


Clear[MatrixSqrt];
MatrixSqrt[m_]:=MatrixPower[m,1/2];


Clear[MatrixAbs];
MatrixAbs[a_]:=MatrixSqrt[a.(a\[ConjugateTranspose])];


Clear[Fidelity];
Fidelity[a_,b_]:=(Plus@@(Sqrt[Eigenvalues[a.b]]))^2;


Clear[Superfidelity];
Superfidelity[a_,b_]:=Tr[a.b]+Sqrt[(1-Tr[a.a])]*Sqrt[(1-Tr[b.b])];


Clear[Subfidelity];
Subfidelity[a_,b_]:=Block[{prod = a.b},Tr[prod] +\[Sqrt]2\[Sqrt](Tr[prod]^2-Tr[prod.prod])];


Clear[TraceDistance];
TraceDistance[a_,b_]:=1/2*Tr[MatrixAbs[a-b]];


Clear[MatrixRe];
MatrixRe[A_?MatrixQ]:=(A+A\[ConjugateTranspose])/2;


Clear[MatrixIm];
MatrixIm[A_?MatrixQ]:=(A-A\[ConjugateTranspose])/2;


Clear[ExpectationValue];
ExpectationValue[\[Rho]_,A_]:=Tr[\[Rho].A];


Clear[Commutator];
Commutator[A_,B_]:=A.B-B.A;


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


Unprotect[sx,sy,sz,\[Sigma]x,\[Sigma]y,\[Sigma]z,id,wh];
Clear[sx,sy,sz,\[Sigma]x,\[Sigma]y,\[Sigma]z,id,wh];
sx = {{0,1},{1,0}};
sy = {{0,-I},{I,0}};
sz = {{1,0},{0,-1}};
\[Sigma]x = sx;
\[Sigma]y = sy;
\[Sigma]z = sz;
id = {{1,0},{0,1}};
wh = {{1,1},{1,-1}};
SetAttributes[sx,Protected];
SetAttributes[sy,Protected];
SetAttributes[sz,Protected];
SetAttributes[\[Sigma]x,Protected];
SetAttributes[\[Sigma]y,Protected];
SetAttributes[\[Sigma]z,Protected];
SetAttributes[id,Protected];
SetAttributes[wh,Protected];


Unprotect[\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5,\[Lambda]6,\[Lambda]7,\[Lambda]8];
Clear[\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5,\[Lambda]6,\[Lambda]7,\[Lambda]8];
\[Lambda]1={{0,1,0},{1,0,0},{0,0,0}};
\[Lambda]2={{0,-I,0},{I,0,0},{0,0,0}};
\[Lambda]3={{1,0,0},{0,-1,0},{0,0,0}};
\[Lambda]4={{0,0,1},{0,0,0},{1,0,0}};
\[Lambda]5={{0,0,-I},{0,0,0},{I,0,0}};
\[Lambda]6={{0,0,0},{0,0,1},{0,1,0}};
\[Lambda]7={{0,0,0},{0,0,-I},{0,I,0}};
\[Lambda]8=Sqrt[1/3]{{1,0,0},{0,1,0},{0,0,-2}};
SetAttributes[\[Lambda]1,Protected];
SetAttributes[\[Lambda]2,Protected];
SetAttributes[\[Lambda]3,Protected];
SetAttributes[\[Lambda]4,Protected];
SetAttributes[\[Lambda]5,Protected];
SetAttributes[\[Lambda]6,Protected];
SetAttributes[\[Lambda]7,Protected];
SetAttributes[\[Lambda]8,Protected];


Clear[Proj];
Proj[v_]:=Table[v[[i]]Conjugate[v[[j]]],{i,1,Length[v]},{j,1,Length[v]}];


Clear[BaseVectors];
BaseVectors[n_Integer]:=Table[UnitVector[n,k],{k,1,n}]


Clear[BaseMatrices];
BaseMatrices[n_Integer]:=Table[Unres[UnitVector[n^2,k]],{k,1,n^2}];


Clear[KroneckerDeltaMatrix];
KroneckerDeltaMatrix[m_,n_,dim_]:=Block[{mtx},
	mtx=Table[0,{dim},{dim}];
	mtx[[m,n]]=1;
	mtx
];


Clear[Lambda1];
Lambda1[i_ ,j_,n_]:=Table[KroneckerDelta[j,\[Mu]]KroneckerDelta[i,\[Nu]] + KroneckerDelta[j,\[Nu]]KroneckerDelta[i,\[Mu]] ,{\[Mu],1,n},{\[Nu],1,n}];


Clear[Lambda2];
Lambda2[i_ ,j_,n_]:=Table[-I(KroneckerDelta[i,\[Mu]]KroneckerDelta[j,\[Nu]] - KroneckerDelta[i,\[Nu]]KroneckerDelta[j,\[Mu]]) ,{\[Mu],1,n},{\[Nu],1,n}];


Clear[Lambda3];
Lambda3[i_,n_]:=Sqrt[2/(i^2-i)]DiagonalMatrix[Join[Append[Table[1,{i-1}],-(i-1)],Table[0,{n-i}]]];


Clear[GeneralizedPauliMatrices];
GeneralizedPauliMatrices[n_]:=Block[{l1,l2,l3,i,j},
	l1=Flatten[Table[Lambda1[i,j,n],{i,1,n},{j,i+1,n}],1];
	l2=Flatten[Table[Lambda2[i,j,n],{i,1,n},{j,i+1,n}],1];
	l3=Table[Lambda3[i,n],{i,2,n}];
	Join[l1,l2,l3]
];


Unprotect[\[Lambda]];
Clear[\[Lambda]];
\[Lambda][i_,n_]:=GeneralizedPauliMatrices[n][[i]];
SetAttributes[\[Lambda],Protected];


Clear[PauliMatrices];
PauliMatrices = {sx,sy,sz};


PauliMatrices[GellMannMatrices];
GellMannMatrices = {\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5,\[Lambda]6,\[Lambda]7,\[Lambda]8};


(* ::Subsection::Closed:: *)
(*Quantum gates*)


Clear[Swap];
Swap[dim_]:=Plus@@Flatten[Table[KroneckerProduct[Ketbra[i,j,dim],Ketbra[j,i,dim]],{i,0,dim-1},{j,0,dim-1}],1];


Clear[QFT];
QFT[n_,method_:"Symbolic"]:=Block[{\[Omega]},
	If [method=="Numerical",\[Omega]=N[Exp[2 \[Pi] I/n]],\[Omega]=Exp[2 \[Pi] I/n]];
	Table[\[Omega]^(i k) ,{i,1,n},{k,1,n}]
];


Clear[CNot];
CNot=Ketbra[1,1,2]\[CircleTimes]sx+(IdentityMatrix[2]-Ketbra[1,1,2])\[CircleTimes]IdentityMatrix[2];


Unprotect[H];
Clear[H];
H = Sqrt[1/2]{{1,1},{1,-1}};
SetAttributes[H,Protected];


Unprotect[X];
Clear[X];
X[d_]:=Sum[Ketbra[Mod[j-1,d],j,d],{j,0,d-1}];
SetAttributes[X,Protected];


Unprotect[Z];
Clear[Z];
Z[d_]:=DiagonalMatrix[Table[Exp[2\[Pi] I j/d],{j,0,d-1}]];
SetAttributes[Z,Protected];


Clear[Circuit];
Circuit[s__]:=Block[{l},l=List[s];Fold[Dot,IdentityMatrix[Length[l[[1]]]],Reverse[l]]];


(* ::Subsection::Closed:: *)
(*Special states*)


Clear[Ket];
Ket[label_,dim_]:=Block[{vec},
	If[label<dim,
	vec =Table[0,{dim}];
	vec[[label+1]]=1;
	vec
	]
];


Clear[Ketbra];
Ketbra[i_,j_,dim_]:=KroneckerProduct[Ket[i,dim],Ket[j,dim]];


Clear[KetFromDigits];
KetFromDigits[l_,b_:2]:=Ket[FromDigits[l,b],b^Length[l]];


Clear[MaxMix];
MaxMix[n_Integer]:=1/n IdentityMatrix[n];


Clear[MaxEnt];
MaxEnt[dim_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],1/Sqrt[subDim] Plus@@Table[Ket[i,subDim]\[CircleTimes]Ket[i,subDim],{i,0,subDim-1}]]
];


Clear[WernerState];
WernerState[dim_,p_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],
	p Proj[MaxEnt[dim]] + (1-p) MaxMix[dim],
	Message[WernerState::argerr,dim];
]
];
WernerState::argerr = "The first `1` argument is not a perfect square.";


Clear[IsotropicState];
IsotropicState[dim_,p_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ[subDim],
	p Proj[MaxEnt[dim]] + (1-p)/(subDim^2-1)(IdentityMatrix[dim]-Proj[MaxEnt[dim]]),
	Message[WernerState::argerr,dim];
]
];
IsotropicState::argerr = "The first `1` argument is not a perfect square.";


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


Unprotect[SchmidtDecomposition];
Clear[SchmidtDecomposition];
SchmidtDecomposition[vec_,d1_,d2_]:=Block[{mtx, svd, vals, snum=Min[d1,d2]},
	mtx = Table[(BaseVectors[d1][[i]]\[CircleTimes]BaseVectors[d2][[j]]).vec,{i,1,d1},{j,1,d2}];
	svd=SingularValueDecomposition[mtx];
	vals=Select[Diagonal[svd[[2]]],#!=0&];
	snum=Length[vals];
	Table[{vals[[i]],svd[[1]].BaseVectors[d1][[i]],svd[[3]]\[Conjugate].BaseVectors[d2][[i]]},{i,1,snum}]
];
SetAttributes[SchmidtDecomposition,Protected];


Clear[OperatorSchmidtDecomposition];
OperatorSchmidtDecomposition[op_,d1_,d2_]:=Block[{mtx, svd, vals, snum=Min[d1*d1,d2*d2]},
	mtx = Table[Tr[(BaseMatrices[d1][[i]]\[CircleTimes]BaseMatrices[d2][[j]]).op\[Transpose]],{i,1,d1*d1},{j,1,d2*d2}];
	svd=SingularValueDecomposition[mtx];
	vals=Select[Diagonal[svd[[2]]],#!=0&];Table[{vals[[i]],Unres[svd[[1]].Res[BaseMatrices[d1][[i]]]],Unres[svd[[3]]\[Conjugate].Res[BaseMatrices[d2][[i]]]]},{i,1,snum}]
];
SetAttributes[OperatorSchmidtDecomposition,Protected];


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Clear[Vec];
Vec[m_]:=Flatten[Transpose[m]]; 


Clear[Unvec];
Unvec[v_List,cols_:0]:=Which[
 (cols== 0)&&IntegerQ[\[Sqrt]Length[v]],Transpose[Partition[v,\[Sqrt]Length[v]]],
 Mod[Length[v],cols]==0,Transpose[Partition[v,cols]]
];


Clear[Res];
Res[m_List]:=Flatten[m]; 


Clear[Unres];
Unres[v_List,cols_:0]:=Which[
 (cols== 0)&&IntegerQ[\[Sqrt]Length[v]],Partition[v,\[Sqrt]Length[v]],
 Mod[Length[v],cols]==0,Partition[v,cols]
];


Clear[Reshuffle];
Reshuffle[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
	If[dim1==0||dim2==0,
		dim=Length[\[Rho]];base1=BaseMatrices[Sqrt[dim]];Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])].Res[\[Rho]],{k,1,dim},{l,1,dim}],
		base1=BaseMatrices[dim1];base2=BaseMatrices[dim2];Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])].Res[\[Rho]],{k,1,dim1 dim1},{l,1,dim2 dim2}]
	]
];


Clear[Reshuffle2];
Reshuffle2[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
	If[dim1==0||dim2==0,
		dim=Length[\[Rho]];base1=BaseMatrices[Sqrt[dim]];Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])\[Transpose]].Res[\[Rho]],{l,1,dim},{k,1,dim}],
		base1=BaseMatrices[dim1];base2=BaseMatrices[dim2];Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])].Res[\[Rho]],{k,1,dim1 dim1},{l,1,dim2 dim2}]
	]
];


Clear[ReshuffleGeneral];
ReshuffleGeneral[A_,n1_,m1_,n2_,m2_]:=
Flatten[
 Table[
  Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]]
  ,{i1,0,n1 n2-1,n2},{i2,0,m1 m2-1,m2}]
,1]


Clear[ReshuffleGeneral2];
ReshuffleGeneral2[A_,n1_,m1_,n2_,m2_]:=
Flatten[
Table[
Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]\[Transpose]]
,{i2,0,m1 m2-1,m2},{i1,0,n1 n2-1,n2}]
,1]\[Transpose]


Clear[MatrixElement];
MatrixElement[n_,\[Nu]_,m_,\[Mu]_,dim_,M_]:=M[[(n-1)*dim[[2]]+\[Nu],(m-1)*dim[[2]]+\[Mu]]];


Clear[ReshufflePermutation];
ReshufflePermutation[dim1_,dim2_]:=Block[{initPos},
	initPos=Flatten[ReshuffleGeneral[Partition[Range[dim1*dim1*dim2*dim2],dim1*dim2],dim1,dim1,dim2,dim2]];
	Table[UnitVector[dim1*dim1*dim2*dim2,Position[initPos,i][[1,1]]]
,{i,1,dim1*dim1*dim2*dim2}]
];


Clear[ProductSuperoperator];
ProductSuperoperator[m1_,m2_]:=Block[{dim1=Length[m1],dim2=Length[m2],perm},
	perm=ReshufflePermutation[Sqrt[dim1],Sqrt[dim2]];
	perm.(m1\[CircleTimes]m2).perm\[Transpose]
];





(* ::Subsection::Closed:: *)
(*Parametrizations*)


Clear[Unitary2];
Unitary2[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_]:=
Exp[I \[Alpha]] DiagonalMatrix[{Exp[-I \[Beta]/2],Exp[I \[Beta]/2]}].{{Cos[\[Gamma]/2],-Sin[\[Gamma]/2]},{Sin[\[Gamma]/2],Cos[\[Gamma]/2]}}.DiagonalMatrix[{Exp[-I \[Delta]/2],Exp[I \[Delta]/2]}];


Clear[SpecialUnitary2];
SpecialUnitary2[\[Beta]_,\[Gamma]_,\[Delta]_]:=Unitary2[0,\[Beta],\[Gamma],\[Delta]];


Clear[Unitary3];
Unitary3[al_,be_,ga_,th_,a_,b_,c_,ph_]:=MatrixExp[I *\[Lambda]3*al].MatrixExp[I*\[Lambda]2*be].
	MatrixExp[I*\[Lambda]3*ga].MatrixExp[I*\[Lambda]5*th].MatrixExp[I*\[Lambda]3*a].
	MatrixExp[I*\[Lambda]2*b].MatrixExp[I*\[Lambda]3*c].MatrixExp[I*\[Lambda]8*ph];


Clear[Unitary4Canonical];
Unitary4Canonical[a1_,a2_,a3_]:=MatrixExp[I a1 KroneckerProduct[\[Sigma]x,\[Sigma]x]+a2 I KroneckerProduct[\[Sigma]y,\[Sigma]y]+a3 I KroneckerProduct[\[Sigma]z,\[Sigma]z]];


Clear[ProbablityDistribution];
ProbablityDistribution[l_]:=Block[{ll,N},
	N=Length[l]+2;
	ll=Prepend[l,\[Pi]/2];
	Table[Sin[ll[[i-1]]]^2*Product[Cos[ll[[j-1]]]^2,{j,i+1,N}],{i,2,N}]
];


Clear[StateVector];
StateVector[l_]:=Block[{pr,ph,N},
	N=Length[l]/2;
	pr=ProbablityDistribution[l[[1;;N]]];
	ph=Prepend[Exp[I*l[[N+1;;2*N]]],1];
	FullSimplify[Sqrt[pr]*ph, Assumptions -> Table[0<l[[i]]<\[Pi]/2,{i,1,N}]]
];


(* ::Subsection::Closed:: *)
(*One-qubit states*)


Clear[QubitKet];
QubitKet[\[Alpha]_,\[Beta]_]:={Cos[\[Alpha]], Exp[I \[Beta]] Sin[\[Alpha]]};


Clear[QubitPureState];
QubitPureState[\[Alpha]_,\[Beta]_]:=Proj[QubitKet[\[Alpha],\[Beta]]];


Clear[QubitBlochState];
QubitBlochState[a_,b_,c_]:=1/2id + a sx + b sy + c sz;


Clear[QubitGeneralState];
QubitGeneralState[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_,\[Lambda]_]:=Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]].DiagonalMatrix[{\[Lambda],1-\[Lambda]}].Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]]\[ConjugateTranspose];


(* ::Subsection::Closed:: *)
(*Quantum channels*)


Clear[IdentityChannel];
IdentityChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]];


Clear[TransposeChannel];
TransposeChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]\[Transpose]];


Clear[DepolarizingChannel];
DepolarizingChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p) Tr[\[Rho]]MaxMix[dim]];


Clear[BitflipChannel];
BitflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sx.\[Rho].sx\[ConjugateTranspose]];


Clear[PhaseflipChannel];
PhaseflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sz.\[Rho].sz\[ConjugateTranspose]];


Clear[BitphaseflipChannel];
BitphaseflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sy.\[Rho].sy\[ConjugateTranspose]];


Clear[HolevoWernerChannel];
HolevoWernerChannel=Function[{dim,p,\[Rho]},( p \[Rho]\[Transpose]+ (1-p)Tr[\[Rho]]MaxMix[dim])];


Clear[GeneralizedPauliKraus];
GeneralizedPauliKraus[d_,p_]:= Flatten[Table[Sqrt[ p[[i+1]][[j+1]]] (MatrixPower[X[d],i].MatrixPower[Z[d],j])\[ConjugateTranspose],{i,0,d-1},{j,0,d-1}],1];


Clear[ApplyKraus];
ApplyKraus[ch_,\[Rho]_]:=Sum[ch[[k]].\[Rho].(ch[[k]]\[ConjugateTranspose]),{k,1,Length[ch]}];


Clear[ApplyUnitary];
ApplyUnitary[U_,\[Rho]_]:=U.\[Rho].U\[ConjugateTranspose];


Clear[ChannelToMatrix];
ChannelToMatrix[fun_,dim_]:=
Table[Tr[(fun[BaseMatrices[dim][[i]]])\[ConjugateTranspose].BaseMatrices[dim][[j]]],{i,1,dim^2},{j,1,dim^2}];


Clear[Superoperator];
Superoperator[ch_List]:=Sum[ch[[i]]\[CircleTimes]ch[[i]]\[Conjugate],{i,1,Length[ch]}];
Superoperator[fun_,dim_]:=ChannelToMatrix[fun,dim];


Clear[DynamicalMatrix];
DynamicalMatrix[ch_List]:=Reshuffle[Superoperator[ch]];
DynamicalMatrix[fun_Function,dim_Integer]:=Reshuffle[Superoperator[fun,dim]];


Clear[Jamiolkowski];
Jamiolkowski[ch_List]:=1/Length[ch[[1]]] DynamicalMatrix[ch];
Jamiolkowski[fun_Function,dim_Integer]:=1/dim DynamicalMatrix[fun,dim];


Clear[ChannelCondition];
ChannelCondition[operators_]:=Sum[operators[[i]]\[ConjugateTranspose].operators[[i]],{i,Length[operators]}];


Clear[ExtendKraus];
ExtendKraus[operators_,n_]:=Module[{tpl},tpl=Tuples[operators,n];Table[KroneckerProduct@@tpl[[i]],{i,1,Length[tpl]}]];


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


Clear[PartialTransposeA];
PartialTransposeA[\[Rho]_,m_,n_]:=Reshuffle[Unres[(Swap[m]\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho]]]]];


Clear[PartialTransposeB];
PartialTransposeB[\[Rho]_,m_,n_]:=Reshuffle[Unres[(IdentityMatrix[m m]\[CircleTimes]Swap[n]).Res[Reshuffle[\[Rho]]]]];


Clear[PartialTraceA];
PartialTraceA[\[Rho]_,m_,n_]:=Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[m]Tr[#]&,m];
	Unres[((trMtx\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho],m,n]])[[1;;n n]]]
];


Clear[PartialTraceB];
PartialTraceB[\[Rho]_,m_,n_]:=Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[n]Tr[#]&,n];
	Unres[Unres[(IdentityMatrix[m m]\[CircleTimes]trMtx).Res[Reshuffle[\[Rho],m,n]]]\[Transpose][[1]]]
];


Clear[PartialTraceGeneral];
PartialTraceGeneral[\[Rho]_,dim_,sys_]:=Block[{n,m},
If[sys==1,Table[Sum[MatrixElement[m,\[Mu],m,\[Nu],dim,\[Rho]],{m,dim[[1]]}],{\[Mu],dim[[2]]},{\[Nu],dim[[2]]}],
	Table[Sum[MatrixElement[n,\[Mu],m,\[Mu],dim,\[Rho]],{\[Mu],dim[[2]]}],{n,dim[[1]]},{m,dim[[1]]}]]
];


Clear[PartialTransposeGeneral];
PartialTransposeGeneral[\[Rho]_,dim_,sys_]:=
If[sys==1,
	ArrayFlatten[Table[
		MatrixElement[n,\[Mu],m,\[Nu],dim,\[Rho]],{n,dim[[1]]},{m,dim[[1]]},{\[Nu],dim[[2]]},{\[Mu],dim[[2]]}
	]]
	,(*else*)
	ArrayFlatten[Table[
		MatrixElement[m,\[Nu],n,\[Mu],dim,\[Rho]],{n,dim[[1]]},{m,dim[[1]]},{\[Nu],dim[[2]]},{\[Mu],dim[[2]]}
	]]
](*endif*)


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


Clear[QubitDepolarizingKraus];
QubitDepolarizingKraus[p_]:={\[Sqrt]((3 p +1)/4)id,\[Sqrt]((1-p)/4)sx,\[Sqrt]((1-p)/4)sy,\[Sqrt]((1-p)/4)sz};


Clear[QubitDecayKraus];
QubitDecayKraus[p_]:={ {{1,0},{0,\[Sqrt]p}} , {{0,\[Sqrt](1-p)},{0,0}} };


Clear[QubitPhaseKraus];
QubitPhaseKraus[p_]:={ {{1,0},{0,\[Sqrt](1-p)}} , {{0,0},{0,\[Sqrt]p}} };


Clear[QubitBitflipKraus];
QubitBitflipKraus[p_]:={ \[Sqrt]p id,\[Sqrt](1-p) sx};


Clear[QubitPhaseflipKraus];
QubitPhaseflipKraus[p_]:={\[Sqrt]p id,\[Sqrt](1-p) sz};


Clear[QubitBitphaseflipKraus];
QubitBitphaseflipKraus[p_]:={\[Sqrt]p id,\[Sqrt](1-p) sy};


Clear[QubitDynamicalMatrix];
QubitDynamicalMatrix[kx_,ky_,kz_,nx_,ny_,nz_]:= 1/2{
	{1 + nz + kz, 0, kx + I ky, nx + ny},
	{0, 1 - nz + kz, nx - ny, kx + I ky},
	{kx - I ky, nx - ny, 1 - nz - kz, 0},
	{nx + ny, kx - I ky, 0, 1 + nz - kz}
}


QubitDaviesDynamicalMatrix[a_,b_,c_]:={{a,0,0,c},{0,b,0,0},{0,0,a,0},{c,0,0,1-b}};


(* ::Subsection::Closed:: *)
(*One-qutrit channels*)


Clear[QutritSpontaneousEmissionKraus];
QutritSpontaneousEmissionKraus[A1_,A2_,t_]:={{{1,0,0},{0,Exp[-(A1 t/2)],0},{0,0,Exp[-(A2 t/2)]}},{{0,Sqrt[1-Exp[-(A1 t)]],0},{0,0,0},{0,0,0}},{{0,0,Sqrt[1-Exp[-(A2 t)]]},{0,0,0},{0,0,0}}};


(* ::Subsection::Closed:: *)
(*Entropies*)


Log0[x_]:=If[x==0,0,Log[2,x]];
SetAttributes[Log0,Listable];


Unprotect[\[Eta]];
Clear[\[Eta]];
\[Eta][x_]:= -x Log[2,x];
SetAttributes[\[Eta],Protected];


Unprotect[\[Eta]2];
Clear[\[Eta]2];
\[Eta]2[x_]:= -x Log[2,x] - (1-x)Log[2,1-x];
SetAttributes[\[Eta]2,Protected];


Unprotect[QuantumEntropy];
Clear[QuantumEntropy];
QuantumEntropy[m_]:=Block[{eigvals,qe},eigvals=Chop[Eigenvalues[m]];
	qe=Sum[eigvals[[i]] Log0[eigvals[[i]]],{i,1,Length[eigvals]}];
	- Chop[qe]
];
SetAttributes[QuantumEntropy,Protected];


Unprotect[S];
Clear[S];
S[m_]:=QuantumEntropy[m];
SetAttributes[S,Protected];


Unprotect[QuantumChannelEntropy];
Clear[QuantumChannelEntropy];
QuantumChannelEntropy[ch_List]:=QuantumEntropy[Jamiolkowski[ch]];
QuantumChannelEntropy[fun_Function,dim_Integer]:=QuantumEntropy[Jamiolkowski[fun,dim]];
SetAttributes[QuantumChannelEntropy,Protected];


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


Unprotect[\[Delta]];
Clear[\[Delta]];
\[Delta][a_,type_:"Dirac"]=Switch[type,
	"Dirac",DiracDelta[a],
	"Indicator",DiscreteDelta[a]
]; 
SetAttributes[\[Delta],Protected];


Unprotect[VandermondeMatrix];
Clear[VandermondeMatrix];
VandermondeMatrix[l_]:=Table[Table[l[[j]]^i,{i,0,Length[l]-1}],{j,1,Length[l]}];
SetAttributes[VandermondeMatrix,Protected];


Unprotect[ProdSum];
Clear[ProdSum];
ProdSum[l_]:=Times@@\[AliasDelimiter][Table[Table[l[[i]]+l[[j]],{i,1,j-1}],{j,2,Length[l]}]];
SetAttributes[ProdSum,Protected];


Unprotect[ProdDiff2];
Clear[ProdDiff2];
ProdDiff2[l_]:=Block[{x},Discriminant[Times@@Table[(x-l[[i]]),{i,1,Length[l]}],x]];
SetAttributes[ProdDiff2,Protected];


Unprotect[ProbBuresNorm];
Clear[ProbBuresNorm];
ProbBuresNorm[N_]:=2^(N^2-N) Gamma[N^2/2]/(\[Pi]^(N/2) Product[Gamma[j+1],{j,1,N}]);
SetAttributes[ProbBuresNorm,Protected];


Unprotect[ProbBures];
Clear[ProbBures];
ProbBures[l_,delta_:"Dirac"]:=FullSimplify[ProbBuresNorm[Length[l]] \[Delta][\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)-1,delta]/\[Sqrt](\!\(
\*SubsuperscriptBox[\(\[Product]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)) Det[VandermondeMatrix[l]]^2/ProdSum[l]];
SetAttributes[ProbBures,Protected];


Unprotect[ProbHSNorm];
Clear[ProbHSNorm];
ProbHSNorm[N_]:=Gamma[N^2]/Product[Gamma[N-j] Gamma[N-j+1],{j,0,N-1}];
SetAttributes[ProbHSNorm,Protected];


Unprotect[ProbHS];
Clear[ProbHS];
ProbHS[l_,delta_:"Dirac"]:=ProbHSNorm[Length[l]]\[Delta][1-(Plus@@l),delta] Det[VandermondeMatrix[l]]^2;
SetAttributes[ProbHS,Protected];


(* ::Subsection::Closed:: *)
(*Random states and operations*)


Unprotect[RandomSimplex];
Clear[RandomSimplex];
RandomSimplex[n_]:=Block[{r,r1,r2},
	r=Sort[Table[RandomReal[{0,1}],{i,1,n-1}]];
	r1=Append[r,1];r2=Prepend[r,0];r1-r2
];
SetAttributes[RandomSimplex,Protected];


Unprotect[RandomKet];
Clear[RandomKet];
RandomKet[n_]:=Block[{p,ph},
	p=Sqrt[RandomSimplex[n]];
	ph=Exp[I RandomReal[{0,2\[Pi]},n-1]];
	ph=Prepend[ph,1];
	p*ph
];
SetAttributes[RandomKet,Protected];


Unprotect[RandomProductKet];
Clear[RandomProductKet];
RandomProductKet[l_]:=Flatten[Apply[KroneckerProduct,Table[RandomKet[l[[i]]],{i,1,Length[l]}]]];
SetAttributes[RandomProductKet,Protected];


(* TODO: Find reference *)
Unprotect[RandomNormalMatrix];
Clear[RandomNormalMatrix];
RandomNormalMatrix[n_]:=Block[{DD,AA,QQ,RR},
	DD=DiagonalMatrix[RandomComplex [{-1-I,1+I},{n}]];
	AA=RandomComplex[{-1-I,1+I},{n,n}];
	{QQ,RR}=QRDecomposition[AA];
	QQ\[ConjugateTranspose].DD.QQ
];
SetAttributes[RandomNormalMatrix,Protected];


Unprotect[RandomDynamicalMatrix];
Clear[RandomDynamicalMatrix];
RandomDynamicalMatrix[n_,m_:0]:=Block[{X,Y,sY},	
	X=GinibreMatrix[n^2,n^2-m];
	Y=PartialTraceGeneral[X.X\[ConjugateTranspose],{n,n},1];
	sY=MatrixPower[Y,-1/2];
	KroneckerProduct[IdentityMatrix[n],sY].X.X\[ConjugateTranspose].KroneckerProduct[IdentityMatrix[n],sY]
];
SetAttributes[RandomDynamicalMatrix,Protected];


Unprotect[GinibreMatrix];
Clear[GinibreMatrix];
GinibreMatrix[m_,n_]:=RandomReal[NormalDistribution[0,1],{m,n}] + I RandomReal[NormalDistribution[0,1],{m,n}];
SetAttributes[GinibreMatrix,Protected];


Unprotect[RandomProductNumericalRange];
Clear[RandomProductNumericalRange];
RandomProductNumericalRange[A_,sys_,noPoints_:1]:=Block[{prod},
	Table[prod=RandomProductKet[sys];Tr[Proj[prod].A],{noPoints}]
];
SetAttributes[RandomProductNumericalRange,Protected];


Unprotect[RandomSpecialUnitary];
Clear[RandomSpecialUnitary];
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
SetAttributes[RandomSpecialUnitary,Protected];


Unprotect[RandomUnitary];
Clear[RandomUnitary];
RandomUnitary[d_]:=Exp[I*RandomReal[2*\[Pi]]]*RandomSpecialUnitary[d];
SetAttributes[RandomUnitary,Protected];


Unprotect[RandomState];
Clear[RandomState];
RandomState[d_]:=Block[{v},
	v=RandomReal[NormalDistribution[0,1],d d] + I RandomReal[NormalDistribution[0,1],d d];
	v=v/Norm[v,2];
	PartialTraceGeneral[Proj[v],{d,d},2]
];
SetAttributes[RandomState,Protected];


(* ::Subsection::Closed:: *)
(*Numerical range*)


Unprotect[NumericalRangeBound];
Clear[NumericalRangeBound];
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
SetAttributes[NumericalRangeBound,Protected];


(* ::Subsection::Closed:: *)
(*Bloch Representation*)


Unprotect[BlochVector];
Clear[BlochVector];
BlochVector[A_]:=Block[{dim},
dim=Length[A]; 1/Sqrt[2](Tr[A\[ConjugateTranspose].#]&/@GeneralizedPauliMatrices[dim])];
SetAttributes[BlochVector,Protected];


Unprotect[StateFromBlochVector];
Clear[StateFromBlochVector];
StateFromBlochVector[vec_]:=Block[{dim},
	If[IntegerQ[Sqrt[Length[vec]+1]],
		dim= Sqrt[Length[vec]+1];
		1/dim IdentityMatrix[dim] + vec.GeneralizedPauliMatrices[dim]/Sqrt[2],
		Print["StateFromBlochVector: given vector is not a Bloch vector of any dimension"];
	]
];
SetAttributes[StateFromBlochVector,Protected];


(* ::Section::Closed:: *)
(*Package footer*)


End[];


EndPackage[];
