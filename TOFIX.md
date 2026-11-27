# Outstanding code and API fixes

## Executable defects

~~-`QI.m`: `BlochToState` sends `StateFromBlochVector::argerr`, but only `BlochToState::argerr` is defined~~.
~~- `QI.m`: `RandomKet[{d1, d2}, "MaxEnt"]` applies an extra square root to its Schmidt coefficients and returns an unnormalized ket.~~
~~- `QI.m`: `Unvec[v]` for a non-square vector evaluates `Mod[Length[v], 0]` rather than reporting a controlled argument error.~~
~~- `QIExtras.m`: `QFT` omits the `1/Sqrt[n]` normalization, so it is not unitary.~~
- `QIExtras.m`: `QubitBitphaseflipChannel` requires an unused `dim` argument despite its documented two-argument API.
- `QIExtras.m`: `RandomMaximallyEntangledNumericalRange` calls undefined `RandomEntangledUnitVector`.


## Documentation/API inaccuracies

~~- `QI.m`: `RandomState` documents induced-measure parameter `K > 2`, while the implementation accepts `K >= d`.~~
~~- `QI.m`: the `ApplyChannel` usage string misspells its own exported name as `ApplayChannel`.~~
~~- `QI.m`: `ApplyKraus` and `ApplyChannel` usage strings reference nonexistent `ApplyUnitary`.~~
~~- `QIExtras.m`: the `Log0` usage text says `Log0[0]` is `1`; its implementation returns `0`.~~
~~- `QIExtras.m`: the `RandomProductNumericalRange` usage string documents `RandomLocalNumericalRange` instead.~~

## Future work

- `QIWalks.m`: usage-message conversion targets `QINisq` symbols instead of `QIWalks` symbols.
- `QIWalks.m`: two nested private contexts are opened but only one is closed, leaving the context stack unbalanced.
- `QIWalks.m`: a walk with zero steps still applies the walk operator once instead of returning its input state.

- `QINisq.m`: `QTransGate` and `QTransCGate` usage headers incorrectly name `QRunCGate`.
- `QINisq.m`: `Gate` localizes `id`, shadowing QI's identity matrix; untargeted qubits are therefore not represented by identity operators.
- `QINisq.m`: `H[]` uses `1/2 wh` rather than `wh/Sqrt[2]`, so it is not unitary.
- `QINisq.m`: `Toff` is a matrix value while the public API documents `Toff[]`; the documented call does not evaluate.
- `QINisq.m`: `InitQC[dim, name]` ignores `name` and always resets `QC`.
- `QINisq.m`: `QTrans` combines multiple strings with `Dot`, producing an unevaluated dot expression instead of transpiled output.
- `QINisq.m`: `TruncatedFidelity` documents `delta` as controlling eigensystem accuracy, but does not pass it to `Chop[Eigensystem[rho]]`.
