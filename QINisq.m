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
	{"0.0.1", "15/03/2021", "Jarek", "Basic rotation gates."}
};  


qiNisqVersion = Last[qiNisqHistory][[1]];


qiNisqLastModification = Last[qiNisqHistory][[2]];


qiNisqAbout = "QINisq is a Mathematica package encapsulating functions useful for NISQ quantum computing. It is based QI and QIExtras packages.";


(* ::Subsection:: *)
(*Quantum gates*)


RX[theta_]:=MatrixExp[-I theta/2 sx];


RY[theta_]:=MatrixExp[-I theta/2 sy];


RZ[theta_]:=MatrixExp[-I theta/2 sz];


(* ::Section:: *)
(*Package footer*)


Print["Package QINisq ", QINisq`Private`qiNisqVersion, " (last modification: ", QINisq`Private`qiNisqLastModification, ")."];



End[] (* End Private Context *)

Protect@@Names["QINisq`*"]

EndPackage[]



