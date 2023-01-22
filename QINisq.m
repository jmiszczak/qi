(* ::Package:: *)

(* ::Section:: *)
(*Package header*)


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


RY::usage = "RX[theta] one-qubit rotation wrt Y axis.";


RZ::usage = "RX[theta] one-qubit rotation wrt Z axis.";


RXP::usage = "PXR[phi, theta] one-qubit rotation of the Bloch vector by an angle theta, where phi is the angle between the rotation axis and the X axis.";


S::usage = "S[] = RZ[pi/2]";


T::usage = "T[] = RZ[pi/4]";


H::usage = "H[] = Hadamard gate.";


Toff::usage = "Toff[] = Toffoli/controlled-controlled-not gate.";


XX::usage = "Ising gate with parameter \[Theta]."


Gate::usage = "Gate[g,t,q] returns quantum gate applied on systems specified in list t, on the computer with q qubits.";


CGate::usage = "CGate[g,c,t,q] returns controlled version of one-qubit gate g with control qubits c and target qubits q, acring on q qubits. Please note that this function does not check of the dimensions are appropriate and qubits specification do not overlap.";


CX::usage = "CX[c,t,qdim] generalized controlled not operation with control qubits c and target qubits t, acting on qdim qubits. Please note that this function does not check of the dimensions are appropriate and qubits specification do not overlap.";


CY::usage = "CX[c,t,qdim] generalized controlled sy. See CX usage for more details."


CZ::usage = "CX[c,t,qdim] generalized controlled sz. See CX usage for more details."


(* ::Subsection:: *)
(*Quantum computer*)


Q::usage = "Q[qg, ops, qc] appends gate qg to the computer qc using options in the ops list. It is expected that the first element of qc list is an integer representing the number of available qubits.";


RunGate::usage = "RunGate[g,t,q] executes one-qubit operation operation g on qubits listed in t on the register of q qubits.";


RunCGate::usage = "RunCGate[g,ops,q] runs controlled gate c with control and target qubits specified in ops on a system with q qubits.";


QRun::usage = "QRun[qc] runs a simulation of the intructions in the list of gates, using the information about the paramters and the size of the computer qc.";


(* ::Subsection:: *)
(* Fidelity and friends *)


TruncatedFidelity::usage = "<f>TruncatedFidelity</f>[<v>r, s ,m, d</v>] returns quantum fidelity between <v>r</v> and <v>s</v> calculated by projecting <v>s</s> onto the <v>m</v> largest eigenvalues of <v>r</v>. Parameter <v>d</v> controls the accuarys of calculations of the eigenvalues and has default value 10e-6.";


(* ::Section:: *)
(*Private definitions*)


Begin["`Private`"] (* Begin Private Context *) 

QIDocRep = {"<v>" -> "\!\(\*StyleBox[\"" , "</v>" -> "\", \"TI\"]\)", "<f>"->"\!\(\*StyleBox[\"", "</f>" -> "\", \"Input\"]\)", "<s>" -> "", "</s>" -> ""} 
(MessageName[Evaluate[ToExpression[#]], "usage"] = StringReplace[MessageName[Evaluate[ToExpression[#]], "usage"],QIDocRep])& /@ Names["QINisq`*"];


(* ::Subsection:: *)
(*Internal functions (without usage strings)*)


qiNisqAuthors = "Jaroslaw Miszczak <miszczak[at]iitis[dot]pl>";


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
	{"0.0.8", "22/01/2023", "Jarek", "Updated description. Minor in usage messages."}
};  


qiNisqVersion = Last[qiNisqHistory][[1]];


qiNisqLastModification = Last[qiNisqHistory][[2]];


qiNisqAbout = "QINisq is a Mathematica package encapsulating functions commonly used to developing quantum algorithms for noisy intermediate-scale quantum (NISQ) computers. The packages is based QI and QIExtras packages and provides functionality enabling the programming of quantum computers on a level similar to the one offered by Qiskit.";


(* ::Subsection:: *)
(*Quantum gates*)


RX[theta_]:=MatrixExp[- I theta/2 sx];


RY[theta_]:=MatrixExp[- I theta/2 sy];


RZ[theta_]:=MatrixExp[- I theta/2 sz];


S[] = MatrixPower[sz,1/2];


T[] = MatrixPower[sz,1/4];


H[] = 1/Sqrt[2]{{1,1},{1,-1}};


RXP[phi_,theta_]:={ {Cos[theta/2], -I Sin[theta/2] Exp[-I phi]}, {-I Sin[theta/2] Exp[I phi], Cos[theta/2]} };


toffoli = Proj[Ket[1,2]]\[CircleTimes]Proj[Ket[1,2]]\[CircleTimes]sx + (Id[4]-Proj[Ket[1,2]]\[CircleTimes]Proj[Ket[1,2]])\[CircleTimes]id;


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


Q[gate_,ops_,qc_]:=AppendTo[qc,{gate,ops}];
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
