(* ::Package:: *)

(* ::Section:: *)
(*Package header*)


(* File: QIWalks.m *)
(* Description: Mathematica package for the analysis of quantum walks. *)
(* Authors: Jaroslaw Miszczak <miszczak@iitis.pl> *)
(* License: GPLv3 *)


BeginPackage["QIWalks`", { "QI`", "QIExtras`"}]
(* Exported symbols added here with SymbolName::usage *) 
Unprotect@@Names["QIWalks`*"]
Clear@@Names["QIWalks`*" ]


(* ::Section:: *)
(*Public definitions*)


UnitaryQuantumWalk::usage="UnitaryQuantumWalk[U,v,n] performs a quantum walk with n steps using operator U on an initial state v. This function does not assume anything about the structure of the coin, state space, and process memory.";


ProbsStateVector::usage="ProbsStateVector[v] returns probebilites for the computational base vectors in the state v.";


Probs::usage="Probs[q] returns probabilities for the computational base for the system in quantum state q, which can be represented by a density matrix or by a state vector.";


(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"] (* Begin Private Context *) 

QIDocRep = {"<v>" -> "\!\(\*StyleBox[\"" , "</v>" -> "\", \"TI\"]\)", "<f>"->"\!\(\*StyleBox[\"", "</f>" -> "\", \"Input\"]\)", "<s>" -> "", "</s>" -> ""} 
(MessageName[Evaluate[ToExpression[#]], "usage"] = StringReplace[MessageName[Evaluate[ToExpression[#]], "usage"],QIDocRep])& /@ Names["QINisq`*"];


UnitaryQuantumWalk[W_,inState_,steps_]:=Block[{i,state},
	state=W . inState;
	For[i=1,i<steps,i++,
		state=W . state;
	];
	state
];


ProbsStateVector[v_]:=Map[Abs[#]^2&,v];


ProbsMixedState[m_]:=Block[{dim},
	dim=Dimensions[m][[1]];
	Table[Tr[m . Ketbra[bs,bs,dim]],{bs,0,dim-1}]
];


Probs[inObj_?MatrixQ]:=ProbsMixedState[inObj];
Probs[inObj_?VectorQ]:=ProbsStateVector[inObj];


UnitaryQuantumWalkMemory[dim_,memState_:((Ket[0,2]+Ket[1,2])/Sqrt[2]),opCoin_:wh,steps_:10]:=Block[
	{opIdim,opI,opShift,opWalk,stInit,stParticle,maxStep,stateAtStep,i},
	maxStep=steps;
	opIdim=IdentityMatrix[dim];
	opI=Id;
	opShift=Plus@@
	Table[
		Ketbra[0,0,2]\[CircleTimes](Ketbra[Ket[0,2]\[CircleTimes]Ket[Mod[v-1,dim],dim],Ket[0,2]\[CircleTimes]Ket[v,dim]])+
		Ketbra[0,0,2]\[CircleTimes](Ketbra[Ket[1,2]\[CircleTimes]Ket[Mod[v+1,dim],dim],Ket[1,2]\[CircleTimes]Ket[v,dim]])+
		Ketbra[1,1,2]\[CircleTimes](Ketbra[Ket[1,2]\[CircleTimes]Ket[Mod[v+1,dim],dim],Ket[0,2]\[CircleTimes]Ket[v,dim]])+
		Ketbra[1,1,2]\[CircleTimes](Ketbra[Ket[0,2]\[CircleTimes]Ket[Mod[v-1,dim],dim],Ket[1,2]\[CircleTimes]Ket[v,dim]]),{v,0,dim-1}
	];
	opWalk=opShift . (opCoin\[CircleTimes]opI\[CircleTimes]opIdim);
	stParticle=Ket[0,dim];
	(*stInit=1/Sqrt[2](Ket[0,2]+Ket[1,2])\[CircleTimes]memState\[CircleTimes]stParticle;*)
	stInit=Ket[0,2]\[CircleTimes]memState\[CircleTimes]stParticle;
	stateAtStep=Table[Null,{maxStep}];
	stateAtStep[[1]]=opWalk . stInit;
	For[i=2,i<=maxStep,i++,
		stateAtStep[[i]]=opWalk . stateAtStep[[i-1]];
	];
	stateAtStep
];


(* ::Subsection:: *)
(*Internal functions (without usage strings)*)


Begin["`Private`"];

qiWalksAuthors = "Jaroslaw Miszczak <miszczak[at]iitis[dot]pl>.";

qiWalksLicense = "GPLv3 <http://www.gnu.org/licenses/gpl.html>";

qiWalksHistory = {
	{"0.0.1", "12/02/2023", "Jarek", "Initial import"}
};  

qiWalksVersion = Last[qiWalksHistory][[1]];

qiWalksLastModification = Last[qiWalksHistory][[2]];

qiWalksAbout = "QIWalks is a package of functions for Mathematica computer algebra system, based on QI package, which implements \
some of functions used in the analysis of quantum walks. The package is ";



(* ::Section:: *)
(*Package footer*)


Print["Package QIWalks ", QIWalks`Private`qiWalksVersion, " (last modification: ", QIWalks`Private`qiWalksLastModification, ")."];


End[] (* End Private Context *)

Protect@@Names["QIWalks`*"]

EndPackage[]



