(* Mathematica Package *)

BeginPackage["QI`"]
(* Exported symbols added here with SymbolName::usage *)  

RandomUnitaryEuler::usage = "Random unitary matrix. Thanks to Rafal Demkowicz-Dobrzanski.";

PartialTraceA::usage = "PartialTraceA[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the m-demensional (first) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";

PartialTraceB::usage = "PartialTraceB[\[Rho],m,n] performs partial trace on m\[Cross]n-dimensional density matrix \[Rho] with respect to the n-dimensional (second) subsystem. This function is implemented using composition of channels. Use PartialTraceGeneral for better performance.";

PartialTraceGeneral::usage = "PartialTraceGeneral[\[Rho],dim,sys] - Returns the partial trace, according to system sys, of density matrix \[Rho] composed of subsystems of dimensions dim={dimA, dimB}. See also: PartialTraceA, PartialTraceB.";


Begin["`Private`"] (* Begin Private Context *) 

PartialTraceA[\[Rho]_,m_,n_]:=Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[m]Tr[#]&,m];
	Unres[((trMtx\[CircleTimes]IdentityMatrix[n n]).Res[ReshuffleGeneral[\[Rho],m,m,n,n]])[[1;;n n]]]
];


PartialTraceB[\[Rho]_,m_,n_] := Block[{trMtx},
	trMtx=ChannelToMatrix[IdentityMatrix[n]Tr[#]&,n];
	Unres[Unres[(IdentityMatrix[m m]\[CircleTimes]trMtx).Res[ReshuffleGeneral[\[Rho],m,m,n,n]]]\[Transpose][[1]]]
];


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

End[] (* End Private Context *)


EndPackage[]