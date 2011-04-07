#!/usr/local/bin/MathematicaScript -script
Print[$ScriptCommandLine[[2]]]
Needs["QI`"];

qiFormatUsageMsg[inName_,inMsg_] := Block[{name = "$ " <> inName <> " $ ",usage=inMsg,txt},
	usage=StringReplace[usage,"\{"->"LEFTCURLY"];
	usage=StringReplace[usage,"\}"->"RIGHTCURLY"];
	txt="\\textbf{"<>name<>"}"<>"-- "<>StringDrop[StringReplace[usage,RegularExpression["\\\\text{([^\}]{5,1000})}"]-> " $$1$ "],2] <> " $";
	txt=StringReplace[txt,"LEFTCURLY"-> "\{"];
	txt=StringReplace[txt,"RIGHTCURLY"->"\}"];
	txt
];


qiGenDoc[docFile_]:=Block[{latexHeader,latexFooter,f,txt,usage,name,functionsList},
	ExportString["","TeX"];

	functionsList=Table[{Names["QI`*"][[i]],ToExpression[Evaluate[Names["QI`*"][[i]]<>"::usage"]]},{i,1,Length[Names["QI`*"]]}];
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
	For[i=1,i<= Length[functionsList],i++,
		name=ToString[TeXForm[functionsList[[i,1]]]];
		usage=ToString[TeXForm[DisplayForm[functionsList[[i,2]]]]];
		txt = qiFormatUsageMsg[name, usage];
		WriteString[f,txt];
		WriteString[f,"\\\\\n\n"];
	];
	WriteString[f,latexFooter];
	Close[f];
];


qiGenDoc[$ScriptCommandLine[[2]]];