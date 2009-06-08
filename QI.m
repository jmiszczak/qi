(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package header*)


BeginPackage["QI`"];


(* File: QI.m *)
(* Description: Mathematica package for analysis of quantum states *)
(* Author: Jaroslaw Miszczak <miszczak@iitis.pl> *)
(* License: GPLv3 *)
(* Last modification: 4 Jun 2009 *)


(* ::Section:: *)
(*Help messages*)


(* ::Subsection::Closed:: *)
(*Kronecker sum and product, symbolic matrix*)


KroneckerSum::usage = "KroneckerSum[A,B] returns the Kronecker sum of A and B defined as A\[CircleTimes]1+1\[CircleTimes]B. Alternative syntax A\[CirclePlus]B for KroneckerSum[A,B] is provided.";


SymbolicMatrix::usage = "SymbolicMatrix[a,m,n] returns m\[Cross]n-matrix with elements \!\(\*SubscriptBox[\"a\", 
RowBox[{\"i\", \",\", \"j\"}]]\), i=1,...,m, j=1,...,n. If the second argument is ommited this function returns square n\[Cross]n matrix. This functions can save you some keystrokes and, thanks to TeXForm function, its results can be easily incorporeted in LaTeX documents.";


SymbolicVector::usage = "SymbolicVector[a,m] is equavalent to Matrix[a,m,1] and it returns a vector with m elements \!\(\*SubscriptBox[\"a\", \"i\"]\),i=1,...,m, j=1,...,n. This function is usefoul, for example, for generating lists of parameters.";


ComplexToPoint::usage = "ComplexToPoint[z] returns real and imaginary parts of a complex number z as a pair of real numbers (point in \!\(\*SuperscriptBox[\"R\", \"2\"]\)).";


(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


MatrixSqrt::usage= "MatrixSqrt[m] - square root for the matrix m.";


MatrixAbs::usage= "MatrixAbs[m] - absolute value for matrix m calculated as MatrixSqrt[m.m\[ConjugateTranspose]].";


Fidelity::usage = "Fidelity[\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\),\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)] returns quantum fidelity between states \!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\) and \!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\).";


Superfidelity::usage = "Superfidelity[A,B] calculates fuperfidelity between A and B defined as Tr[A.B] + Sqrt[1-Tr[A.A]]Sqrt[1-Tr[B.B]]."; 


Subfidelity::usage = "Subfidelity[A,B] returns superfidelity between states A and B calculated as \!\(\*SubscriptBox[\"tr\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)+\[Sqrt]2\[Sqrt]((\!\(\*SubscriptBox[\"tr\[Rho]\", \"2\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\))-\!\(\*SubscriptBox[\"tr\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"1\"]\)\!\(\*SubscriptBox[\"\[Rho]\", \"2\"]\)).";


TraceDistance::usage = "TraceDistance[A,B] = 1/2tr|A-B|.";


MatrixRe::usage = "Hermitian part of the matrix A - 1/2(A+A\[ConjugateTranspose]).";


MatrixIm::usage = "Antyhermitian part of the matrix A - 1/2(A-A\[ConjugateTranspose]).";


ExpectationValue::usage = "ExpectationValue[\[Rho],A] = Tr[\[Rho].A].";


Commutator::usage = "Comutator of matrices A and B.";


(* ::Subsection::Closed:: *)
(*Commonly used matrices*)


sx::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\).";
sy::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\).";
sz::usage = "Pauli matrix \!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\).";
\[Sigma]x::usage=sx::usage;
\[Sigma]y::usage=sy::usage;
\[Sigma]z::usage=sz::usage;
id::usage  = "Identity matrix for one qubit.";
wh::usage = "Hadamard gate for one qubit.";


\[Lambda]1::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\).";
\[Lambda]2::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\).";
\[Lambda]3::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\).";
\[Lambda]4::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"4\"]\).";
\[Lambda]5::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"5\"]\).";
\[Lambda]6::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"6\"]\).";
\[Lambda]7::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"7\"]\).";
\[Lambda]8::usage = "Gell-Mann matrix \!\(\*SubscriptBox[\"\[Lambda]\", \"8\"]\).";


Proj::usage = "Proj[{v1,...,v2}] returns projectors for the vectors in the input list.";


BaseVectors::usage = "BaseVectors[n] - canonical basis in n-dimensional Hilbert space";


BaseMatrices::usage = "BaseMatrices[n] - canonical basis in n\[Cross]n-dimensional Hilbert-Schmidt space.";


KroneckerDeltaMatrix::usage = "KroneckerDeltaMatrix[i,j,d] returns d\[Cross]d matrix with 1 at position (i,j) and zeroes elsewhere.";


Lambda1::usage = "Lambda1[i,j,n] generalized Pauli matrix {\!\(\*SubscriptBox[\"\[Lambda]\", \"1\"]\)(i,j)\!\(\*SubscriptBox[\"}\", 
RowBox[{\"k\", \",\", \"l\"}]]\)=\!\(\*SubscriptBox[\"\[Delta]\", \"jk\"]\)\!\(\*SubscriptBox[\"\[Delta]\", \"il\"]\)+\!\(\*SubscriptBox[\"\[Delta]\", \"jl\"]\)\!\(\*SubscriptBox[\"\[Delta]\", \"ik\"]\) for i<j. For example Lambda1[1,2,2] is equal to Pauli \!\(\*SubscriptBox[\"\[Sigma]\", \"x\"]\).";


Lambda2::usage = "Lambda2[i,j,n] generalized Pauli matrix {\!\(\*SubscriptBox[\"\[Lambda]\", \"2\"]\)(i,j)\!\(\*SubscriptBox[\"}\", 
RowBox[{\"k\", \",\", \"l\"}]]\)=-\[ImaginaryI](\!\(\*SubscriptBox[\"\[Delta]\", \"ik\"]\)\!\(\*SubscriptBox[\"\[Delta]\", \"jl\"]\)-\!\(\*SubscriptBox[\"\[Delta]\", \"il\"]\)\!\(\*SubscriptBox[\"\[Delta]\", \"jk\"]\)) for i<j. For example Lambda2[1,2,2] is equal to \!\(\*SubscriptBox[\"\[Sigma]\", \"y\"]\).";


Lambda3::usage ="Lambda3[i,n] generalized Pauli matrix {\!\(\*SubscriptBox[\"\[Lambda]\", \"3\"]\)}= \[Sqrt](\!\(\*FractionBox[\"2\", 
RowBox[{SuperscriptBox[\"n\", \"2\"], \"-\", \"n\"}]]\))diag(1,1,...,-1,0,...0), with -1 at n-th position, for i=2,...,n. For example Lambda3[2,2] is equal to \!\(\*SubscriptBox[\"\[Sigma]\", \"z\"]\).";


GeneralizedPauliMatrices::usage = "GeneralizedPauliMatrices[n] - list of generalized Pauli matrices for SU(n). For n=2 these are just Pauli matrices and for n=3 - Gell-Mann matrices. Note that identity matrix is not included in the list. See also: PauliMatrices, GellMannMatrices, \[Lambda], Lambda1, Lambda2, Lambda3.";


\[Lambda]::usage = "\[Lambda][i,n] is defined as GeneralizedPauliMatrices[n][[i]].";


PauliMatrices::usage = "List of Pauli matrices. Use Map[MatrixForm[#]&,PauliMatrices] to get this list in more readible form.";


GellMannMatrices::usage = "List of Gell-Mann matrices. Use Map[MatrixForm[#]&,GellMannMatrices] to get this list in more readible form.";


(* ::Subsection::Closed:: *)
(*Quantum gates*)


Swap::usage="Swap[n] returns permutation operator \!\(\*UnderoverscriptBox[\"\[Sum]\", 
RowBox[{\"i\", \"=\", \"0\"}], 
RowBox[{\"n\", \"-\", \"1\"}]]\)\!\(\*UnderoverscriptBox[
RowBox[{\" \", \"\[Sum]\"}], 
RowBox[{\"j\", \"=\", \"0\"}], 
RowBox[{\"n\", \"-\", \"1\"}]]\)|i\[RightAngleBracket]\[LeftAngleBracket]j|\[CircleTimes]|j\[RightAngleBracket]\[LeftAngleBracket]i\[VerticalSeparator] acting on \!\(\*SuperscriptBox[\"n\", \"2\"]\)-dimensional space and exchanging two n-dimensional subsystems.";


QFT::usage = "QFT[n,method] - quantum Fourier transform of dimension n. This function accepts second optional argument, which specifices method used in calculation. Parameter method can be equal to \"Symbloic\", which is default, or \"Numerical\". The second option makes this function much faster.";


CNot::usage = "Controlled not matrix for two qubits.";


H::usage = "Hadamard gates for one qubit.";


X::usage = "Generalized Pauli matrix X.";


Z::usage = "Generalized Pauli matrix Z.";


Circuit::usage = "Constructs quantum circuit from unitary gates.";


(* ::Subsection:: *)
(*Spacial states*)


Ket::usage = "Ket[i,d] returns |i\[RightAngleBracket] in d-dimensinal Hilbert space.";


Ketbra::usage = "Ketbra[i,j,d] returns \[VerticalSeparator]i\[RightAngleBracket]\[LeftAngleBracket]j\[VerticalSeparator] ascting on d-dimensional space.";


KetFromDigits::usage = "KetFromDigits[list,base] - returns ket vector labeled by string of digits represented in given base.";


MaxMix::usage = "MaxMix[n] gies maximally mixed state in n-dimensional space of density matrices.";


MaxEnt::usage = "MaxEnt[N] - maximally entangled state \!\(\*FractionBox[\"1\", 
RowBox[{\"\[Sqrt]\", \"N\"}]]\) \!\(\*UnderoverscriptBox[\"\[Sum]\", 
RowBox[{\"i\", \"=\", \"0\"}], 
RowBox[{\"n\", \"-\", \"1\"}]]\)|i\[RightAngleBracket]\[CircleTimes]|i\[RightAngleBracket] with N = \!\(\*SuperscriptBox[\"n\", \"2\"]\).";


WernerState::usage= "WernerState[p,n] - Werner state with parameter p\[Element][0,1] for n\[Cross]n-dimensional system. This state is defined as p \!\(\*FractionBox[\"2\", 
RowBox[{\"n\", 
RowBox[{\"(\", 
RowBox[{\"n\", \"+\", \"1\"}], \")\"}]}]]\)\!\(\*SubscriptBox[\"P\", \"sym\"]\) + (1-p) \!\(\*FractionBox[\"2\", 
RowBox[{\"n\", 
RowBox[{\"(\", 
RowBox[{\"n\", \"-\", \"1\"}], \")\"}]}]]\)\!\(\*SubscriptBox[\"P\", \"snty\"]\), where \!\(\*SubscriptBox[\"P\", \"sym\"]\) and \!\(\*SubscriptBox[\"P\", \"snty\"]\) are projectors for symmetric and anty-symmetric subspace.";


WernerState4::usage = "Werner state for two qubits.";


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


SchmidtDecomposition::usage = "SchmidtDecomposition[vec,d1,d2] - Schmidt decomposition of the vector vec in d1\[Cross]d2-dimensional Hilbert space.";


(* ::Subsection::Closed:: *)
(*Reshaping, vectorization and reshuffling*)


Vec::usage = "Vec[m] - vectorization of the matrix m column by column. See also: Res.";


Unvec::usage = "Unvec[v,c] - de-vectorization of the vector into the matrix with c colummns. If the second parameter is ommited then it is assumed that v can be mapped into square matrix. See also: Unres, Vec.";


Res::usage = "Res[m] is equvalent to Vec[m\[Transpose]]. Reshaping maps matrx m into vector row by row.";


Unres::usage = "de-reshaping of the vector into the matrix with c colummns. If the second parameter is ommited then it is assumed that v can be mapped into square matrix. See also: Unvec, Res.";


Reshuffle::usage = "Reshuffle[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices. If the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be ommited. In this case one obtains a reshuffle in the basis contrtucted by using two bases of d-dimensional Hilbert-Schmidt matrix spaces.";


MatrixElement::usage = "MatrixElement[n,\[Nu],m,\[Mu],div,M] - returns the matrix element of density matrix M indexed by two double indices n, \[Nu] and m, \[Mu] of the composite sytem of dimensions dim={dimA, dimB}";


(* ::Subsection:: *)
(*Parametrizations*)


Unitary2::usage = "Euler parametrization of U(2).";


SpecialUnitary2::usage ="Euler parametrization of SU(2).";


Unitary3::usage = "Euler parametrization of U(3)";


Unitary4Canonical::usage = "Parametrization of non-local unitary marices for two qubits. See: arXiv:quant-ph/0011050v1.";


ProbablityDistribution::usage = "ProbablityDistribution[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}] returns probability vectors of dimension n+1 parametrize with {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)}. See also: SymbolicVector.";


StateVector::usage = "StateVector[{\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\),\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"n\", \"+\", \"1\"}]]\),...,\!\(\*SubscriptBox[\"\[Phi]\", 
RowBox[{\"2\", \" \", \"n\"}]]\)}] returns pure n+1-dimensional pure state (ket vector) constructed form probability distribution parametrize by numbers {\!\(\*SubscriptBox[\"\[Theta]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Theta]\", \"n\"]\)} and phases {\!\(\*SubscriptBox[\"\[Phi]\", \"1\"]\),...,\!\(\*SubscriptBox[\"\[Phi]\", \"n\"]\)}. See also: ProbablityDistribution, SymbolicVector.";


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet::usage = "QubitKet[\[Alpha],\[Beta]] parametriation of the pure state (as a state vector) for one qubit as (cos(\[Alpha]) \!\(\*SuperscriptBox[\"\[ExponentialE]\", \"\[ImaginaryI]\[Beta]\"]\), sin(\[Alpha])). This is equivalent to StateVector[{\[Alpha],\[Beta]}]. See also: QubitPureState, StateVector.";


QubitPureState::usage = "QubitPureState[\[Alpha],\[Beta]] - parametriation of the pure state as a density matrix for one qubit. This is just a alias for Proj[QubitKet[\[Alpha],\[Beta]]]. See also: QubitKet.";


QubitBlochState::usage = "Parametrization of the one-qubit mixed state on the Bloch sphere.";


QubitState::usage = "QubitState[\[Alpha],\[Beta],\[Gamma],\[Delta],\[Lambda]] - Parametrization of the one-qubit mixed state using rotations and eigenvalues. Returns one-qubits density matrix with eigenvalues \[Lambda] and 1-\[Lambda] rotated as U.diag(\[Lambda],1-\[Lambda]).\!\(\*SuperscriptBox[\"U\", \"\[Dagger]\"]\) with U defined by parameters \[Alpha],\[Beta],\[Gamma] and \[Delta].";


QubitGeneralState::usage = "QubitGeneralState[a,b,c] - parametrization of the one-qubits mixed state using only normalization and self-adjointness.";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


IdentityChannel::usage="IdentityChannel[n, \[Rho]] - applay identity operation on n-dimensional density matrix \[Rho].";


TransposeChannel::usage="TransposeChannel[n, \[Rho]] - applay transposition operation on n-dimensional density matrix \[Rho]. This operations is not completely positive.";


DepolarizingChannel::usage = "DepolarizingChannel[n,p,\[Rho]] performs an action of the completely depolarizing channel with paramaeter p acting on n-dimensional input state \[Rho]. See also: QubitDepolarizingKraus, HolevoWernerChannel.";


BitflipChannel::usage  = "BitflipChannel[2,p,\[Rho]].";


PhaseflipChannel::usage  = "PhaseflipChannel[2,p,\[Rho]]";


BitphaseflipChannel::usage  = "BitphaseflipChannel[2,p,\[Rho]].";


HolevoWernerChannel::usage = "HolevoWernerChannel[n,p,\[Rho]] performs an action of the Holeve-Werner channel (also known as transpose depolarizing channel) with paramaeter p acting on n-dimensional input state \[Rho]. See also: DepolarizingChannel.";


ChannelToMatrix::usage = "ChannelToMatrix[E,d] returns matrix representation of a channel E acting on d-dimensional state space \!\(\*
StyleBox[\"i\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"e\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontSlant->\"Italic\"]\) matrix \!\(\*SubscriptBox[\"M\", \"E\"]\) such that \!\(\*SubscriptBox[\"M\", \"E\"]\) res(\[Rho])=\!\(\*SuperscriptBox[\"res\", 
RowBox[{\"-\", \"1\"}]]\)(E(\[Rho])). Here res(M) is a reshaping of matrix M implemented by Res[M] and \!\(\*SuperscriptBox[\"res\", 
RowBox[{\"-\", \"1\"}]]\)(M) is implemented as Unres[M]. First argument should be a pure function E such that E[\[Rho]] transforms input state according to the channel definition. For example for the Holevo-Werner channel one ca use ChannelToMatrix[HolevoWernerChannel[3,p,#]&,3] to obtain matrix representation of this channel acting on qutrits. See also: Superoperator.";


GeneralizedPauliKraus::usage = "GeneralizedPauliKraus[d,P] - list of Kraus operators for d-dimensional generalized Pauli channel with the d-dimesnional matrix of parameters P. See: M. Hayashi, Quantum Imformation An Introduction, Springer 2006, Exampl 5.8, p .126.";


ApplyKraus::usage = "ApplyKraus[ch,\[Rho]] - applay channel ch, given as a list of Kraus operators, to the input state \[Rho].";


Superoperator::usage = "Superoperator[kl] returns matrix representation of quantum channel given as a list of Kraus operators. Superoperator[fun,dim] is just am alternative name for ChannelToMatrix[fun,dim] and returns matrix representation of quantum channel, given as a pure function, acting on dim-deimensionla space. So Superoperator[DepolarizingChannel[2,p,#]&,2] and Superoperator[QubitDepolarizingKraus[p]] returns the same matrix. See also: ChannelToMatrix.";


DynamicalMatrix::usage = "Dynamical matrix of quantum channel given as a list of Kraus operators (DynamicalMatrix[ch]) or as a function fun action on dim-dimensional space (DynamicalMatrix[fun,dim]). See alos: Superoperator, ChannelToMatrix.";


Jamiolkowski::usage = "Jamiolkowski[K] gives the image of the Jamiolkowski isomorfizm for the channel given as the list of Karus operators K. Jamiolkowski[fun,dim] gives the image of the Jamiolkowski isomorfizm for the channel given as a function fun action on dim-dimensional space. See alos: Superoperator, ChannelToMatrix, DynamicalMatrix.";


ChannelCondition::usage = "Performs some checks on Kraus operators. Use this if you want to check if they represent quantum channel.";


ExtendKraus::usage = "ExtendKraus[ch,n] - produces n-fold tensor products of Kraus operators from the list ch.";


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTransposeA::usage = "PartialTransposeA[\[Rho],m,n] performs partial transposition on the m-dimensional (first) subsystem of the m\[Cross]n-state.";


PartialTransposeB::usage = "PartialTransposeB[\[Rho],m,n] performs partial transposition on the n-dimensional (second) subsystem of the m\[Cross]n-state.";


PartialTraceA::usage = "PartialTraceA[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the m-demensional (first) subsystem.";


PartialTraceB::usage = "PartialTraceB[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the n-demensional (second) subsystem.";


PartialTraceGeneral::usage = "PartialTraceGeneral[\[Rho],dim,sys] - Returns the partial trace, acording to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}. See alos: PartialTraceA, PartialTraceB.";


PartialTransposeGeneral::usage = "PartialTransposeGeneral[\[Rho],dim,sys] - Returns the partial transpose, acording to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}";


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


QubitDepolarizingKraus::usage = "Kraus operators of the depolarizing channel for one qubit. Note that it gives maximally mixed state for p=0.";


QubitDecayKraus::usage = "Kraus operators of the decay channel, also know as amplitude damping, for one qubit.";


QubitPhaseKraus::usage = "Kraus operators for one qubit phase damping channel.";


QubitBitflipKraus::usage = "Kraus operators for one qubit bit-flip channel.";


QubitPhaseflipKraus::usage = "Kraus operators for one qubit phase-flip channel.";


QubitBitphaseflipKraus::usage = "Kraus operators for one qubit bit-phase-flip channel.";


(* ::Subsection::Closed:: *)
(*Entropies*)


Log0::usage = "Log0[x] is equal to Log[2,x] for x>0 and 1 for x=0.";


\[Eta]::usage = "\[Eta][x] = -x Log[2,x].";


\[Eta]2::usage = "\[Eta]2[x] = \[Eta][x]+\[Eta][1-x].";


QuantumEntropy::usage = "QuantumEntropy[m] - von Neuman entropy for the matrix m.";


S::usage = "S[m] = QuantumEntropy[m].";


QuantumChannelEntropy::usage = "QuantumChannelEntropy[ch] - von Neuman entropy of the quantum channel calculated as a von Neuman entropy for the image of this channe in Jamiolkowski isomorphism. See also: Jamiolkowski, Superoperator.";


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


\[Delta]::usage = "\[Delta][x] represents Dirac delta at x. \[Delta][x,\"Indicator\"] is equal to 1 at x and zero elsewhere.";


VandermondeMatrix::usage = "VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)}] - Vandermonde matrix for variables (\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)).";


ProdSum::usage = "ProdSum[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] gives \!\(\*SubsuperscriptBox[\"\[Product]\", 
RowBox[{\"i\", \"<\", \"j\"}], \"n\"]\)(\!\(\*SubscriptBox[\"x\", \"i\"]\)+\!\(\*SubscriptBox[\"x\", \"j\"]\)).";


ProdDiff2::usage = "ProdDiff2[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}] is equivalent to Det[VandermondeMatrix[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}]\!\(\*SuperscriptBox[\"]\", \"2\"]\) and gives a discriminant of the polynomial with roots {\!\(\*SubscriptBox[\"x\", \"1\"]\),...,\!\(\*SubscriptBox[\"x\", \"n\"]\)}.";


ProbBuresNorm::usage = "ProbBNorm[n] - Normalization factor used for calculating probablity distribution of eigenvalues of matrix of dimension N according to Bures distance.";


ProbBures::usage = "ProbBures[\[Lambda]] - Joint probablity distribution of eigenvalues \[Lambda] of a matrix according to Bures distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: \"Indicator\"";


ProbHSNorm::usage = "Normalization factor used for calculating probablity distribution of eigenvalues of matrix of dimension N according to Hilbert-Schmidt distribution.";


ProbHS::usage = "ProbHS[{\!\(\*SubscriptBox[\"x\", \"1\"]\),...\!\(\*SubscriptBox[\"x\", \"n\"]\)},] Probablity distribution of eigenvalues of matrix according to Hilbert-Schmidt distance. By default \[Delta] is assumed to be Dirac delta. Other possible values: \"Indicator\"";


(* ::Subsection:: *)
(*Random states and operations*)


RandomSimplex::usage = "RandomSimplex[d] - d-dimesnional random simplex.";


RandomKet::usage = "RandomKet[d] - random ket in d-dimensional space. See: T. Radtke, S. Fritzsche / Computer Physics Communications 179 (2008) 647\[Dash]664.";


RandomProductKet::usage = "RandomProductKet[{dim1,dim2,...,dimN}] - random pure state (ket verctor) of the tensor product form with dimensions of subspaces specified dim1, dim2,...,dimN.";


RandomNormalMatrix::usage = "RandomNormalMatrix[d] - random normal matrix of dimension d.";


RandomDynamicalMatrix::usage = "RandomDynamicalMatrix[d,k] returns dynamical matrix of operation acting on d-dimensional states with k eigenvalues equalt to 0.";


GinibreMatrix::usage = "GinibreMatrix[m,n] returns complex matrix of dimension m\[Cross]n with normal distribution of real and imaginary parts.";


RandomProductNumericalRange::usage = "RandomLocalNumericalRange[M,{dim1,dim2,...,dimN},n] returns n points from the product numerical range of the matrix M with respect to division specified as {dim1,dim2,...,dimN}. Note that dim1\[Cross]dim2\[Cross]...\[Cross]dimN must be equal to the dimension of matrix M.";


RandomSpecialUnitry::usage = "R.D-D.";


RandomUnitary::usage = "Random unitary matrix. Thanks to R.D-D.";


RandomState::usage = "RandomState[d] - random density matrix of dimension d. This gives uniform distribution with respect to Hilbert-Schmidt measure.";


(* ::Subsection::Closed:: *)
(*Numerical range*)


NumericalRangeBound::usage = "NumericalRangeBound[A_?MatrixQ,step_:0.01] - bound of numerical range of matrix A calculated with given step. Ref: Carl C. Cowen, Elad Harel, An Effective Algorithm for Computing the Numerical Range. Technical report, Dep. of Math. Purdue University, 1995.";


(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Kronecker sum and product, symbolic matrix*)


x_ \[CircleTimes] y_ := KroneckerProduct[x,y];
x_ \[CircleTimes] y_ \[CircleTimes] z_ := (x \[CircleTimes] y) \[CircleTimes] z;


KroneckerSum[A_,B_]:=Block[{dim=Length[A]},
KroneckerProduct[A,IdentityMatrix[dim]]+KroneckerProduct[IdentityMatrix[dim],B]
];


x_ \[CirclePlus] y_ := KroneckerSum[x,y];
x_ \[CirclePlus] y_ \[CirclePlus] z_ := (x \[CirclePlus] y) \[CirclePlus] z;


SymbolicMatrix[sym_,d1_,d2_:0]:= Which[
	d2==0,Table[Subscript[sym, i,j],{i,1,d1},{j,1,d1}],
	d2==1,Table[Subscript[sym, i],{i,1,d1}],
	True,Table[Subscript[sym, i,j],{i,1,d1},{j,1,d2}] 
];


SymbolicVector[sym_,d1_]:= SymbolicMatrix[sym,d1,1];


ComplexToPoint[z_]:={Re[z],Im[z]};
SetAttributes[ComplexToPoint,Listable];


(* ::Subsection::Closed:: *)
(*Fidelity, trace distance etc.*)


MatrixSqrt[m_]:=Module[{u,w,v},{u,w,v}=SingularValueDecomposition[m];u.Sqrt[w].v\[ConjugateTranspose]]


MatrixAbs[a_]:=MatrixSqrt[a.(a\[ConjugateTranspose])];


Fidelity[a_,b_]:=Block[{sqrt=MatrixSqrt[a]},Tr[MatrixSqrt[sqrt.b.sqrt]]^2];


Superfidelity[a_,b_]:=Tr[a.b]+\[Sqrt](1-Tr[a.a])\[Sqrt](1-Tr[b.b]) ;


Subfidelity[a_,b_]:=Block[{prod = a.b},Tr[prod] +\[Sqrt]2\[Sqrt](Tr[prod]^2-Tr[prod.prod])];


TraceDistance[a_,b_]:=1/2Tr[MatrixAbs[a-b]];


MatrixRe[A_?MatrixQ]:=(A+A\[ConjugateTranspose])/2;


MatrixIm[A_?MatrixQ]:=(A-A\[ConjugateTranspose])/2;


ExpectationValue[\[Rho]_,A_]:=Tr[\[Rho].A];


Commutator[A_,B_]:=A.B-B.A;


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
mtx]


Lambda1[i_ ,j_,n_]:=Table[KroneckerDelta[j,\[Mu]]KroneckerDelta[i,\[Nu]] + KroneckerDelta[j,\[Nu]]KroneckerDelta[i,\[Mu]] ,{\[Mu],1,n},{\[Nu],1,n}];


Lambda2[i_ ,j_,n_]:=Table[-I(KroneckerDelta[i,\[Mu]]KroneckerDelta[j,\[Nu]] - KroneckerDelta[i,\[Nu]]KroneckerDelta[j,\[Mu]]) ,{\[Mu],1,n},{\[Nu],1,n}];


Lambda3[i_,n_]:=Sqrt[2/(i^2-i)]DiagonalMatrix[Join[Append[Table[1,{i-1}],-(i-1)],Table[0,{n-i}]]];


GeneralizedPauliMatrices[n_]:=Block[{l1,l2,l3,i,j},
l1=Flatten[Table[Lambda1[i,j,n],{i,1,n},{j,i+1,n}],1];
l2=Flatten[Table[Lambda2[i,j,n],{i,1,n},{j,i+1,n}],1];
l3=Table[Lambda3[i,n],{i,2,n}];
Join[l1,l2,l3]
];


\[Lambda][i_,n_]:=GeneralizedPauliMatrices[n][[i]];


PauliMatrices={sx,sy,sz};


GellMannMatrices = {\[Lambda]1,\[Lambda]2,\[Lambda]3,\[Lambda]4,\[Lambda]5,\[Lambda]6,\[Lambda]7,\[Lambda]8};


(* ::Subsection:: *)
(*Quantum gates*)


Swap[dim_]:=Plus@@Flatten[Table[KroneckerProduct[Ketbra[i,j,dim],Ketbra[j,i,dim]],{i,0,dim-1},{j,0,dim-1}],1];


QFT[n_,method_:"Symbolic"]:=Block[{\[Omega]},
	If [method=="Numerical",\[Omega]=N[Exp[2 \[Pi] I/n]],\[Omega]=Exp[2 \[Pi] I/n]];
	Table[\[Omega]^(i k) ,{i,1,n},{k,1,n}]
];


CNot=Ketbra[1,1,2]\[CircleTimes]sx+(IdentityMatrix[2]-Ketbra[1,1,2])\[CircleTimes]IdentityMatrix[2];


H = Sqrt[1/2]{{1,1},{1,-1}};
SetAttributes[H,Protected];


X[d_]:=Sum[Ketbra[Mod[j-1,d],j,d],{j,0,d-1}];
SetAttributes[X,Protected];


Z[d_]:=DiagonalMatrix[Table[Exp[2\[Pi] I j/d],{j,0,d-1}]];
SetAttributes[Z,Protected];


Circuit[s__]:=Block[{l},l=List[s];Fold[Dot,IdentityMatrix[Length[l[[1]]]],Reverse[l]]];


(* ::Subsection::Closed:: *)
(*Spacial states*)


Ket[label_,dim_]:=Block[{vec},
If[label<dim,
vec =Table[0,{dim}];
vec[[label+1]]=1;
vec
]
];


Ketbra[i_,j_,dim_]:=KroneckerProduct[Ket[i,dim],Ket[j,dim]];


KetFromDigits[l_,b_:2]:=Ket[FromDigits[l,b],b^Length[l]];


MaxMix[n_Integer]:=1/n IdentityMatrix[n];


MaxEnt[dim_]:=Block[{subDim=Sqrt[dim]},
If[IntegerQ [subDim],1/Sqrt[subDim] Plus@@Table[Flatten[{Ket[i,subDim]}\[CircleTimes]Ket[i,subDim]],{i,0,subDim-1}]]
];


WernerState[p_,dim_]:=Block[{subDim=\[Sqrt]dim,permut},
If[IntegerQ[subDim],permut=Swap[subDim];
p/(subDim (subDim+1))(IdentityMatrix[dim]+permut)+ p/(subDim (subDim-1)) (IdentityMatrix[dim]-permut)]
];


WernerState4[x_]:=(x*Proj[(1/\[Sqrt]2)(Ket[0,4]+Ket[3,4])]+(1-x)*(1/4)*IdentityMatrix[4]);


(* ::Subsection::Closed:: *)
(*Schmidt decomposition*)


SchmidtDecomposition[vec_,d1_,d2_]:=Block[{mtx,svd,vals,snum=Min[d1,d2]},
mtx = Table[(Flatten[{BaseVectors[d1][[i]]}\[CircleTimes]BaseVectors[d2][[j]]]).vec,{i,1,d1},{j,1,d2}];
svd=SingularValueDecomposition[mtx];
vals=Select[Diagonal[svd[[2]]],#!=0&];
snum=Length[vals];
Table[{vals[[i]],svd[[1]].BaseVectors[d1][[i]],svd[[3]]\[Conjugate].BaseVectors[d2][[i]]},{i,1,snum}]
];


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


Reshuffle[\[Rho]_,m_:0,n_:0]:=Block[{base1,base2,dim},
If[m==0 || n==0,dim=Length[\[Rho]];base1=BaseMatrices[Sqrt[dim]];Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])].Res[\[Rho]],{k,1,dim},{l,1,dim}],
	base1=BaseMatrices[m];base2=BaseMatrices[n];Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])].Res[\[Rho]],{k,1,m m},{l,1,n n}]]
];


MatrixElement[n_,\[Nu]_,m_,\[Mu]_,dim_,M_]:=M[[(n-1)*dim[[2]]+\[Nu],(m-1)*dim[[2]]+\[Mu]]];


(* ::Subsection::Closed:: *)
(*Parametrizations*)


Unitary2[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_]:=Exp[I \[Alpha]] DiagonalMatrix[{Exp[- I \[Beta]/2 ],Exp[I \[Beta]/2 ]}].{{Cos[\[Gamma]/2],-Sin[\[Gamma]/2]},{Sin[\[Gamma]/2],Cos[\[Gamma]/2]}}.DiagonalMatrix[{Exp[- I \[Delta]/2 ],Exp[I \[Delta]/2 ]}];


SpecialUnitary2[\[Beta]_,\[Gamma]_,\[Delta]_]:=Unitary2[0,\[Beta],\[Gamma],\[Delta]];


Unitary3[al_,be_,ga_,th_,a_,b_,c_,ph_]:=MatrixExp[I *\[Lambda]3*al].MatrixExp[I*\[Lambda]2*be].
MatrixExp[I*\[Lambda]3*ga].MatrixExp[I*\[Lambda]5*th].MatrixExp[I*\[Lambda]3*a].
MatrixExp[I*\[Lambda]2*b].MatrixExp[I*\[Lambda]3*c].MatrixExp[I*\[Lambda]8*ph];


Unitary4Canonical[a1_,a2_,a3_]:=MatrixExp[I a1 KroneckerProduct[\[Sigma]x,\[Sigma]x]+a2 I KroneckerProduct[\[Sigma]y,\[Sigma]y]+a3 I KroneckerProduct[\[Sigma]z,\[Sigma]z]];


ProbablityDistribution[l_]:=Block[{ll,N},
N=Length[l]+2;
ll=Prepend[l,\[Pi]/2];
Table[Sin[ll[[i-1]]]^2*Product[Cos[ll[[j-1]]]^2,{j,i+1,N}],{i,2,N}]
];


StateVector[l_]:=Block[{pr,ph,N},
N=Length[l]/2;
pr=ProbablityDistribution[l[[1;;N]]];
ph=Prepend[Exp[I*l[[N+1;;2*N]]],1];
FullSimplify[Sqrt[pr]*ph, Assumptions -> Table[0<l[[i]]<\[Pi]/2,{i,1,N}]]
];


(* ::Subsection::Closed:: *)
(*One-qubit states*)


QubitKet[\[Alpha]_,\[Beta]_]:={Cos[\[Alpha]], Exp[I \[Beta]] Sin[\[Alpha]]};


QubitPureState[\[Alpha]_,\[Beta]_]:=Proj[QubitKet[\[Alpha],\[Beta]]];


QubitBlochState[a_,b_,c_]:=1/2id + a sx + b sy + c sz;


QubitState[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Delta]_,\[Lambda]_]:=Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]].DiagonalMatrix[{\[Lambda],1-\[Lambda]}].Unitary2[\[Alpha],\[Beta],\[Gamma],\[Delta]]\[ConjugateTranspose];


QubitGeneralState [a_,b_,c_]:= {{a, b + I c},{b - I c ,1-a}};


(* ::Subsection:: *)
(*Quantum channels*)


IdentityChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]];


TransposeChannel=Function[{dim,\[Rho]},IdentityMatrix[dim].\[Rho]\[Transpose]];


DepolarizingChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p) Tr[\[Rho]]MaxMix[dim]];


BitflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sx.\[Rho].sx\[ConjugateTranspose]];


PhaseflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sz.\[Rho].sz\[ConjugateTranspose]];


BitphaseflipChannel=Function[{dim,p,\[Rho]},p IdentityMatrix[dim].\[Rho]+(1-p)sy.\[Rho].sy\[ConjugateTranspose]];


HolevoWernerChannel=Function[{dim,p,\[Rho]},1/(dim-1) ( p \[Rho]\[Transpose]+ (1-p)Tr[\[Rho]]MaxMix[dim])];


GeneralizedPauliKraus[d_,p_]:= Flatten[Table[Sqrt[ p[[i+1]][[j+1]]] (MatrixPower[X[d],i].MatrixPower[Z[d],j])\[ConjugateTranspose],{i,0,d-1},{j,0,d-1}],1];


ApplyKraus[ch_,\[Rho]_]:=Sum[ch[[k]].\[Rho].(ch[[k]]\[ConjugateTranspose]),{k,1,Length[ch]}];


ChannelToMatrix[fun_,dim_]:=
Table[Tr[(fun[BaseMatrices[dim][[i]]])\[ConjugateTranspose].BaseMatrices[dim][[j]]],{i,1,dim^2},{j,1,dim^2}];


Superoperator[ch_List]:=Sum[ch[[i]]\[CircleTimes]ch[[i]]\[Conjugate],{i,1,Length[ch]}];
Superoperator[fun_,dim_]:=ChannelToMatrix[fun,dim];


DynamicalMatrix[ch_List]:=Reshuffle[Superoperator[ch]];
DynamicalMatrix[fun_Function,dim_Integer]:=Reshuffle[Superoperator[fun,dim]];


Jamiolkowski[ch_List]:=1/Length[ch[[1]]] DynamicalMatrix[ch];
Jamiolkowski[fun_Function,dim_Integer]:=1/dim DynamicalMatrix[fun,dim];


ChannelCondition[operators_]:=Sum[operators[[i]]\[ConjugateTranspose].operators[[i]],{i,Length[operators]}];


ExtendKraus[operators_,n_]:=Module[{tpl},tpl=Tuples[operators,n];Table[KroneckerProduct@@tpl[[i]],{i,1,Length[tpl]}]];


(* ::Subsection::Closed:: *)
(*Partial trace and transposition*)


PartialTransposeA[\[Rho]_,m_,n_]:=Reshuffle[Unres[(Swap[m]\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho]]]]];


PartialTransposeB[\[Rho]_,m_,n_]:=Reshuffle[Unres[(IdentityMatrix[m m]\[CircleTimes]Swap[n]).Res[Reshuffle[\[Rho]]]]];


PartialTraceA[\[Rho]_,m_,n_]:=Block[{trMtx},
trMtx=ChannelToMatrix[IdentityMatrix[m]Tr[#]&,m];
Unres[((trMtx\[CircleTimes]IdentityMatrix[n n]).Res[Reshuffle[\[Rho],m,n]])[[1;;n n]]]
];


PartialTraceB[\[Rho]_,m_,n_]:=Block[{trMtx},
trMtx=ChannelToMatrix[IdentityMatrix[n]Tr[#]&,n];
Unres[Unres[(IdentityMatrix[m m]\[CircleTimes]trMtx).Res[Reshuffle[\[Rho],m,n]]]\[Transpose][[1]]]
];


PartialTraceGeneral[\[Rho]_,dim_,sys_]:=Block[{n,m},
If[sys==1,Table[Sum[MatrixElement[m,\[Mu],m,\[Nu],dim,\[Rho]],{m,dim[[1]]}],{\[Mu],dim[[2]]},{\[Nu],dim[[2]]}],
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
](*endif*)


(* ::Subsection::Closed:: *)
(*One-qubit quantum channels*)


QubitDepolarizingKraus[p_]:={\[Sqrt]((3 p +1)/4)id,\[Sqrt]((1-p)/4)sx,\[Sqrt]((1-p)/4)sy,\[Sqrt]((1-p)/4)sz};


QubitDecayKraus[p_]:={ {{1,0},{0,\[Sqrt]p}} , {{0,\[Sqrt](1-p)},{0,0}} };


QubitPhaseKraus[p_]:={ {{1,0},{0,\[Sqrt](1-p)}} , {{0,0},{0,\[Sqrt]p}} };


QubitBitflipKraus[p_]:={ \[Sqrt]p id,\[Sqrt](1-p) sx};


QubitPhaseflipKraus[p_]:={\[Sqrt]p id,\[Sqrt](1-p) sz};


QubitBitphaseflipKraus[p_]:={\[Sqrt]p id,\[Sqrt](1-p) sy};


(* ::Subsection:: *)
(*Entropies*)


Log0[x_]:=If[x==0,0,Log[2,x]];
SetAttributes[Log0,Listable];


\[Eta][x_]:= -x Log[2,x];
SetAttributes[\[Eta],Protected];


\[Eta]2[x_]:= -x Log[2,x] - (1-x)Log[2,1-x];
\[Eta]2::usage = "\[Eta]2[x] = \[Eta][x]+\[Eta][1-x]";


QuantumEntropy[m_]:=Block[{eigvals,qe},eigvals=Chop[Eigenvalues[m]];
qe=Sum[eigvals[[i]] Log0[eigvals[[i]]],{i,1,Length[eigvals]}];
- Chop[qe]
];


S[m_]:=QuantumEntropy[m];
SetAttributes[S,Protected];


QuantumChannelEntropy[ch_List]:=QuantumEntropy[Jamiolkowski[ch]];
QuantumChannelEntropy[fun_Function,dim_Integer]:=QuantumEntropy[Jamiolkowski[fun,dim]];


(* ::Subsection::Closed:: *)
(*Distribution of eigenvalues*)


\[Delta][a_,type_:"Dirac"]=Switch[type,
"Dirac",DiracDelta[a],
"Indicator",DiscreteDelta[a]]; 


VandermondeMatrix[l_]:=Table[Table[l[[j]]^i,{i,0,Length[l]-1}],{j,1,Length[l]}];


ProdSum[l_]:=Times@@Flatten[Table[Table[l[[i]]+l[[j]],{i,1,j-1}],{j,2,Length[l]}]];


ProdDiff2[l_]:=Block[{x},Discriminant[Times@@Table[(x-l[[i]]),{i,1,Length[l]}],x]];


ProbBuresNorm[N_]:=2^(N^2-N) Gamma[N^2/2]/(\[Pi]^(N/2) Product[Gamma[j+1],{j,1,N}]);


ProbBures[l_,delta_:"Dirac"]:=FullSimplify[ProbBuresNorm[Length[l]] \[Delta][\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)-1,delta]/\[Sqrt](\!\(
\*SubsuperscriptBox[\(\[Product]\), \(i = 1\), \(Length[l]\)]\(l[\([i]\)]\)\)) Det[VandermondeMatrix[l]]^2/ProdSum[l]];


ProbHSNorm[N_]:=Gamma[N^2]/Product[Gamma[N-j] Gamma[N-j+1],{j,0,N-1}];


ProbHS[l_,delta_:"Dirac"]:=ProbHSNorm[Length[l]]\[Delta][1-(Plus@@l),delta] Det[VandermondeMatrix[l]]^2;


(* ::Subsection:: *)
(*Random states and operations*)


RandomSimplex[n_]:=Block[{r,r1,r2},
	r=Sort[Table[RandomReal[{0,1}],{i,1,n-1}]];
	r1=Append[r,1];r2=Prepend[r,0];r1-r2
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
Table[prod=RandomProductState[sys];Tr[Proj[prod].A],{noPoints}]
];


RandomSpecialUnitary[d_]:=Module[{psi,chi,r,s,phi,i,j,k,u},
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


RandomState[d_]:=Block[{v},
	v=RandomReal[NormalDistribution[0,1],d d] + I RandomReal[NormalDistribution[0,1],d d];
	v=v/Norm[v,2];
	PartialTraceGeneral[Proj[v],{d,d},2]
];


(* ::Subsection:: *)
(*Numerical range*)


NumericalRangeBound[A_?MatrixQ,step_:0.01]:=Block[{w,Ath,Hth,m,s,Kth,pKp,ee,rr,mm,sm,mM,sM},
w={};
Table[
Ath=Exp[I*(-\[Theta])]*A;
Hth=MatrixRe[Ath];
{e,r}=Eigensystem[Hth];
e=Re[e];
m=Max[e];
s=Position[e,m];
If[Length[s]==1,
(*then*)
	AppendTo[w,MatrixToScalar[Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose]]],
(*else*)
	Kth=I*(Hth-Ath);
	pKp=Extract[r,s]\[Conjugate].Kth.Extract[r,s]\[Transpose];
	{ee,rr}=Eigensystem[pKp];
	ee=Re[ee];
	mm=Min[ee];
	sm=Position[ee,mm];
	AppendTo[w,MatrixToScalar[Extract[rr,sm]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sm]\[Transpose]]];
	mM=Max[ee];
	sM=Position[ee,mM];
	AppendTo[w,MatrixToScalar[Extract[rr,sM]\[Conjugate].Extract[r,s]\[Conjugate].A.Extract[r,s]\[Transpose].Extract[rr,sM]\[Transpose]]]
(*end if*)
]
,{\[Theta],0,2\[Pi],step}];
w
]


(* ::Section::Closed:: *)
(*Package footer*)


End[];


EndPackage[];
