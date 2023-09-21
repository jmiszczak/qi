(* ::Package:: *)

(* ::Section:: *)
(*Package header*)


(* File: QINisq.m *)
(* Description: Mathematica package for developing code for NISQ computers. *)
(* Authors: Jaroslaw Miszczak <jarek@miszczak.eu> *)
(* License: GPLv3 *)


(* Mathematica Package *)

BeginPackage["QINisq`", { "QI`", "QIExtras`"}]
(* Exported symbols added here with SymbolName::usage *) 
Unprotect@@Names["QINisq`*"]
Clear@@Names["QINisq`*" ]


(* ::Section:: *)
(*Public definitions*)


(* ::Subsection:: *)
(*Quantum gates*)


RX::usage = "RX[theta] one-qubit rotation wrt X axis.";


RY::usage = "RY[theta] one-qubit rotation wrt Y axis.";


RZ::usage = "RZ[theta] one-qubit rotation wrt Z axis.";


RXP::usage = "PXR[phi, theta] one-qubit rotation of the Bloch vector by an angle theta, where phi is the angle between the rotation axis and the X axis.";


S::usage = "S[] = RZ[pi/2]";


T::usage = "T[] = RZ[pi/4]";


H::usage = "H[] = Hadamard gate. Defined using wh.";


Toff::usage = "Toff[] returns Toffoli, or controlled-controlled-not, gate.";


XX::usage = "XX[\[Theta]] Ising gate with parameter \[Theta].";


Gate::usage = "Gate[g,t,q] returns quantum gate applied on systems specified in list t, on the computer with q qubits.";


CGate::usage = "CGate[g,c,t,q] returns controlled version of one-qubit gate g with control qubits c and target qubits q, acting on q qubits. Please note that currently this function does not check if the dimensions are appropriate and qubits specification do not overlap.";


CX::usage = "CX[c,t,qdim] generalized controlled not operation with control qubits c and target qubits t, acting on qdim qubits. Please note that this function does not check of the dimensions are appropriate and qubits specification do not overlap.";


CY::usage = "CY[c,t,qdim] generalized controlled sy. See CX usage for more details.";


CZ::usage = "CZ[c,t,qdim] generalized controlled sz. See CX usage for more details.";


(* ::Subsection:: *)
(*Quantum computer*)


QC::usage = "QC is the default list representing the quantum computer used during the calculations. Its first element is the number of qubits. Next elements are quantum gates, added using Q function. QC is initialized and reset using InitQC. See also: Q, InitQC.";


InitQC::usage = "InitQC[qdim,name] initializes a (virtual) quantum computer with qdim qubits of quantum memory. The name ncan be specified to create a quantum coputer with a nono-default name. When necessary, this function will resete the initilized quantum computer.";


Q::usage = "Q[qg, ops, qc] is used to change the instructions for the quantum computer. Itappends gate qg to the computer qc using options in the ops list. It is expected that the first element of qc list is an integer representing the number of available qubits.";


RunGate::usage = "RunGate[g,t,q] executes one-qubit operation operation g on qubits listed in t on the register of q qubits.";


RunCGate::usage = "RunCGate[g,ops,q] runs controlled gate g with control and target qubits specified in ops on a system with q qubits.";


QRun::usage = "QRun[qc] runs a simulation of the intructions in the list of gates, using the information about the paramters and the size of the computer qc.";


(* ::Subsection:: *)
(* Fidelity and friends *)


TruncatedFidelity::usage = "<f>TruncatedFidelity</f>[<v>r, s ,m, d</v>] returns quantum fidelity between <v>r</v> and <v>s</v> calculated by projecting <v>s</s> onto the <v>m</v> largest eigenvalues of <v>r</v>. Parameter <v>d</v> controls the accuarys of calculations of the eigenvalues and has default value set to 10E-6.";


(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"] (* Begin Private Context *) 

QIDocRep = {"<v>" -> "\!\(\*StyleBox[\"" , "</v>" -> "\", \"TI\"]\)", "<f>"->"\!\(\*StyleBox[\"", "</f>" -> "\", \"Input\"]\)", "<s>" -> "", "</s>" -> ""} 
(MessageName[Evaluate[ToExpression[#]], "usage"] = StringReplace[MessageName[Evaluate[ToExpression[#]], "usage"],QIDocRep])& /@ Names["QINisq`*"];


(* ::Subsection:: *)
(*Internal functions (without usage strings)*)


qiNisqAuthors = "Jaroslaw Miszczak <jarek[at]miszczak[dot]eu>";


qiNisqLicense = "GPLv3 <http://www.gnu.org/licenses/gpl.html>";


qiNisqHistory = {
	{"0.0.1", "15/03/2021", "Jarek", "Basic rotation gates."},
	{"0.0.2", "12/04/2021", "Jarek", "TruncatedFidelity moved from QIExtras."},
	{"0.0.2", "07/05/2021", "Jarek", "Added XX gate."},
	{"0.0.3", "17/05/2021", "Jarek", "Added Bloch verctor rotation wrt rotatet X axis."},
	{"0.0.4", "21/05/2021", "Jarek", "Added toffoli, S and T gates"},
	{"0.0.5", "22/05/2021", "Jarek", "Added general controlled one qubit gate and templates for CX, CY and CZ gates."},
	{"0.0.6", "25/05/2021", "Jarek", "Added template for function managing virtual quantum devices."},
	{"0.0.7", "28/05/2021", "Jarek", "Added function for extending one-qubit-gates and functions for executing non-parametric one-qubit gates and controlled gates."},
	{"0.0.8", "22/01/2023", "Jarek", "Updated description. Minor in usage messages."},
	{"0.0.9", "12/02/2023", "Jarek", "Minor update: description, usage messages."},
	{"0.0.10", "12/08/2023", "Jarek", "Better quantum computer initilization and management, fixed H gate, minor in usage messages."},
	{"0.0.11", "16/08/2023", "Jarek", "Fixed some usage messages."},
	{"0.0.12", "21/09/2023", "Jarek", "Improved demo file. Improved usage messages."}
};  


qiNisqVersion = Last[qiNisqHistory][[1]];


qiNisqLastModification = Last[qiNisqHistory][[2]];


qiNisqAbout = "QINisq is a Mathematica package encapsulating functions commonly used to developing quantum algorithms for noisy intermediate-scale quantum (NISQ) computers. The packages is based QI and QIExtras packages. Its goal is to providee functionality enabling the programming of quantum computers on a level similar to the one offered by Qiskit.";


(* ::Subsection:: *)
(*Quantum gates*)


RX[theta_]:=MatrixExp[- I theta/2 sx];


RY[theta_]:=MatrixExp[- I theta/2 sy];


RZ[theta_]:=MatrixExp[- I theta/2 sz];


S[] = MatrixPower[sz,1/2];


T[] = MatrixPower[sz,1/4];


H[] = 1/2 wh;


RXP[phi_,theta_]:={ {Cos[theta/2], -I Sin[theta/2] Exp[-I phi]}, {-I Sin[theta/2] Exp[I phi], Cos[theta/2]} };


Toff = Proj[Ket[1,2]]\[CircleTimes]Proj[Ket[1,2]]\[CircleTimes]sx + (Id[4]-Proj[Ket[1,2]]\[CircleTimes]Proj[Ket[1,2]])\[CircleTimes]id;


XX[theta_]:=Cos[theta](Id[2]\[CircleTimes]Id[2])-I Sin[theta](sx\[CircleTimes]sx);


Gate[op_,targ_,qdim_]:=Block[{res},
	res =Table[id,{qdim}];
	res[[Flatten[{targ}]]] = Table[op,{Flatten[{targ}]}];
	KroneckerProduct@@res
];


CGate[gate_,ctrl_,targ_,qdim_]:=Block[{idx,currIdx,idle},
	idle = Complement[Range[1,qdim],ctrl~Join~targ];
	idx = DeleteDuplicates[
	Table[
			currIdx=Reverse@IntegerDigits[i,2,qdim];
			currIdx[[idle]]=Table[id,{idle}];
			If[DeleteDuplicates[currIdx[[ctrl]]]=={1},
				currIdx[[targ]]=Table[gate,{targ}],
				currIdx[[targ]]=Table[id,{targ}]
			];
			currIdx[[ctrl]]=Map[Proj[Ket[#,2]]&,currIdx[[ctrl]]];
			currIdx,
			{i,0,2^qdim-1}
		]
	];
	Plus@@Table[KroneckerProduct@@currIdx,{currIdx,idx}]
];


CX[c_,t_,qdim_]:=CGate[sx, Flatten[{c}], Flatten[{t}],qdim];


CY[c_,t_,qdim_]:=CGate[sy, Flatten[{c}], Flatten[{t}],qdim];


CZ[c_,t_,qdim_]:=CGate[sz, Flatten[{c}], Flatten[{t}],qdim];


(* ::Subsection:: *)
(*Quantum computer*)


InitQC[dim_,name_:QC]:=Module[{}, 
	Print["[Info] Initializing quantum computer with " <> ToString[dim] <> " qubits."];
	Unprotect[QC];
	QC={dim};
	Protect[QC];
];


Q[gate_,ops_,qc_:QC]:=
Module[{}, 
	Unprotect[QC];
	AppendTo[qc,{gate,ops}];
	Protect[QC];
];
SetAttributes[Q,HoldAll];


RunGate[g_,ops_,qdim_]:=ToExpression[
	StringTemplate["Gate[`OP`[],`targ`,`qdim`]"][
		<|"OP"->g,
		"targ"->ops ,
		"qdim"->qdim|>]
];


RunCGate[g_,ops_,qdim_]:=ToExpression[
	StringTemplate["`OP`[`ctrl`,`targ`,`qdim`]"][
		<|"OP"->g,
		"ctrl"->ops[[1]] ,
		"targ"->ops[[2]],
		"qdim"->qdim|>]
];


QRun[qc_]:=Block[{dim=qc[[1]]},
(
	Dot@@Table[
	Switch[
	Head[qc[[i]][[2]][[1]]],
		Integer,RunGate[qc[[i]][[1]],qc[[i]][[2]],qc[[1]]],
		List,RunCGate[qc[[i]][[1]],{qc[[i]][[2,1]],qc[[i]][[2,2]]},qc[[1]]]
		],
	{i,2,Length[qc]}
])
];


(* ::Subsection:: *)
(* Fidelity and friends *)


TruncatedFidelity[rho_,sigma_,m_,delta_:10^-6]:=Block[{vec,val,proj},
	{val,vec}=Chop[Eigensystem[rho]];
	proj=Plus@@(Proj/@vec[[1;;m]]);
	Chop[Fidelity[proj . rho . proj,proj . sigma . proj],delta]
];


(* ::Section:: *)
(*Package footer*)


Print["Package QINisq ", QINisq`Private`qiNisqVersion, " (last modification: ", QINisq`Private`qiNisqLastModification, ")."];


End[] (* End Private Context *)

Protect@@Names["QINisq`*"]

EndPackage[]
