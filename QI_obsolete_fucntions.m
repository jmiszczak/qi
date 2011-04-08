(* Mathematica Package *)

BeginPackage["QI`"]
(* Exported symbols added here with SymbolName::usage *)  

RandomUnitaryEuler::usage = "Random unitary matrix. Thanks to Rafal Demkowicz-Dobrzanski.";

MatrixElement::usage = "MatrixElement[n,\[Nu],m,\[Mu],dim,A] - returns the matrix element of a matrix A indexed by two double indices n, \[Nu] and m, \[Mu] of the composite sytem of dimensions dim=dim1*dim2.";

PartialTraceA::usage = "PartialTraceA[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the m-demensional (first) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";

PartialTraceB::usage = "PartialTraceB[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the n-dimensional (second) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";

PartialTraceGeneral::usage = "PartialTraceGeneral[\[Rho],dim,sys] - Returns the partial trace, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}. See also: PartialTraceA, PartialTraceB.";

PartialTransposeA::usage = "PartialTransposeA[\[Rho],m,n] performs partial transposition on the m-dimensional (first) subsystem of the m\[Cross]n-state.";

PartialTransposeB::usage = "PartialTransposeB[\[Rho],m,n] performs partial transposition on the n-dimensional (second) subsystem of the m\[Cross]n-state.";

PartialTransposeGeneral::usage = "BUGGY! PartialTransposeGeneral[\[Rho],dim,sys] - Returns the partial transpose, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA,dimB}. ";

ReshuffleBase::usage = "ReshuffleBase[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices. If  the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See also: Reshuffle, ReshuffleGeneral, Reshuffle2.";

ReshuffleBase2::usage = "Alternative definition of the reshuffling operation. Reshuffle2[\[Rho],m,n] returns representation of the m\[Cross]n-dimensional square matrix \[Rho] in the basis consisting of product matrices which are transposed versions of standard base matrices. If the matrix \[Rho] has dimension \!\(\*SuperscriptBox[\"d\", \"2\"]\) then two last arguments can be omitted. In this case one obtains a reshuffle in the basis constructed by using two bases of d-dimensional Hilbert-Schmidt matrix spaces. See: See also: Reshuffle2, ReshuffleGeneral, Reshuffle, BaseMatrices";

ReshuffleGeneral::usage = "ReshuffleGeneral[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix.";

ReshuffleGeneral2::usage = "ReshuffleGeneral2[\[Rho],n1,m1,n2,m2] for matrix of size (n1 n2)\[Times](m1 m2) returns a reshuffled matrix - given by alternative definition of the reshuffling operation.";

VectorSchmidtDecomposition::usage = "VectorSchmidtDecomposition[vec,d1,d2] - Schmidt decomposition of the vector vec in d1\[Cross]d2-dimensional Hilbert space.";

OperatorSchmidtDecomposition::usage = "OperatorSchmidtDecomposition[mtx,d1,d2] - Schmidt decomposition of mtx in the Hilbert-Schmidt space of matrices of dimension d1\[Cross]d2.";

SchmidtDecomposition::usage = "SchmidtDecomposition[e,d1,d2] - accepts a vector or a matrix as a first argument and returns apropriate Schmidt decomposition. See also: VectorSchmidtDecomposition, OperatorSchmidtDecomposition.";


Begin["`Private`"] (* Begin Private Context *) 


MatrixElement[n_,\[Nu]_,m_,\[Mu]_,dim_,mtx_]:=mtx[[(n-1)*dim[[2]]+\[Nu],(m-1)*dim[[2]]+\[Mu]]];

PartialTraceA[\[Rho]_,m_,n_]:=Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[m]Tr[#]&,m];
	Unres[((trMtx\[CircleTimes]IdentityMatrix[n n]).Res[ReshuffleGeneral[\[Rho],m,m,n,n]])[[1;;n n]]]
];


PartialTraceB[\[Rho]_,m_,n_] := Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[n]Tr[#]&,n];
	Unres[Unres[(IdentityMatrix[m m]\[CircleTimes]trMtx).Res[ReshuffleGeneral[\[Rho],m,m,n,n]]]\[Transpose][[1]]]
];

PartialTransposeA[\[Rho]_,m_,n_] := Reshuffle[Unres[(Swap[m*m]\[CircleTimes]IdentityMatrix[n*n]).Res[Reshuffle[\[Rho]]]]];


PartialTransposeB[\[Rho]_,m_,n_] := Reshuffle[Unres[(IdentityMatrix[m*m]\[CircleTimes]Swap[n*n]).Res[Reshuffle[\[Rho]]]]];

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


PartialTraceGeneral[\[Rho]_,dim_,sys_] := Block[{n,m},
	If[sys==1,
		Table[Sum[MatrixElement[m,\[Mu],m,\[Nu],dim,\[Rho]],{m,dim[[1]]}],{\[Mu],dim[[2]]},{\[Nu],dim[[2]]}],
		(* else *)
		Table[Sum[MatrixElement[n,\[Mu],m,\[Mu],dim,\[Rho]],{\[Mu],dim[[2]]}],{n,dim[[1]]},{m,dim[[1]]}]]
];

RandomUnitaryEuler[d_]:=Exp[I*RandomReal[2*\[Pi]]]*RandomSpecialUnitary[d];

RandomSpecialUnitary[d_]:=Module[{psi,chi,r,s,phi,i,j,u,e,phi0,psi0,chi0},
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

ReshuffleBase[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
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


ReshuffleBase2[\[Rho]_,dim1_:0,dim2_:0]:=Block[{base1,base2,dim},
	If[dim1==0||dim2==0,
		dim=Length[\[Rho]];
		base1=BaseMatrices[Sqrt[dim]];
		Table[Res[(base1[[k]]\[CircleTimes]base1[[l]])\[Transpose]].Res[\[Rho]],{l,1,dim},{k,1,dim}],
		(* else *)
		base1=BaseMatrices[dim1];
		base2=BaseMatrices[dim2];
		Table[Res[(base1[[k]]\[CircleTimes]base2[[l]])\[Transpose]].Res[\[Rho]],{l,1,dim1 dim1},{k,1,dim2 dim2}]
	]
];

ReshuffleGeneral[A_,n1_,m1_,n2_,m2_]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]],{i1,0,n1 n2-1,n2},{i2,0,m1 m2-1,m2}]
,1];

ReshuffleGeneral2[A_,n1_,m1_,n2_,m2_]:=Flatten[
	Table[Flatten[Part[A,1+i1;;n2+i1,1+i2;;m2+i2]\[Transpose]],{i2,0,m1 m2-1,m2},{i1,0,n1 n2-1,n2}]
,1]\[Transpose];

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
		mtx=Reshuffle[op,{d1,d2},{d1,d2}];
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


End[] (* End Private Context *)


EndPackage[]