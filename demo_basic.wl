#!/usr/local/bin/wolframscript -f
Needs["QI`"]

NewPar[] := Print@StringRepeat["-", 80]

NewPar[]
Print["Basic qubit states, gates, and channel action using QI.m."]
Print["Exmaple 0: A computational-basis qubit and its density matrix."]

ket0 = {1, 0};
rho0 = Proj[ket0];
Print["Base state as a vector:", ket0]

Print["Base state as a matrix:", MatrixForm[rho0]]


NewPar[]
Print["Example 1: Hadamard creates |+>; CNOT entangles |+> tensor |0>."]

ketPlus = 1/Sqrt[2] wh . ket0;
bellKet = ketPlus \[CircleTimes] ket0;
bellState = Proj[bellKet];

Print["Bell state as a vector:"]
Print[bellKet]

Print["Bell state as a density matrix:"]
Print[bellState]


NewPar[]
Print["Example 2: Apply the bit-flip channel with probability p."]
Print["           The channel is represented by two Kraus operators."]

p = 0.2;
bitFlipKraus = {Sqrt[1 - p] id, Sqrt[p] sx};
flippedMixture = ApplyKraus[bitFlipKraus, rho0];
statePlus = Proj[ketPlus]

Print["|+><+| = ", MatrixForm[statePlus]];
Print["Bell-state density matrix = ", MatrixForm[bellState]];
Print["Bit-flip channel on |0><0| = ", MatrixForm[flippedMixture]];

