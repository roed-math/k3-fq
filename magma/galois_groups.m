// Galois groups of polynomials, as transitive group labels.
// Usage: magma -b infile:=in.txt outfile:=out.txt magma/galois_groups.m
// Each input line is a comma-separated coefficient list (constant term
// first) of an irreducible polynomial; each output line is
// "<coefficients> nTk".
SetColumns(0);
R<x> := PolynomialRing(Integers());
lines := [l : l in Split(Read(infile), "\n") | #l gt 0];
out := Open(outfile, "w");
for line in lines do
    f := R![StringToInteger(c) : c in Split(line, ",")];
    G := GaloisGroup(f);
    i, n := TransitiveGroupIdentification(G);
    fprintf out, "%o %oT%o\n", line, n, i;
end for;
delete out;
quit;
