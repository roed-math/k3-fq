// For each transitive group H = kTj with k <= 11, identify the imprimitive
// wreath product C2 wr H (degree 2k) as a transitive group.
// Output lines: k j n i, meaning C2 wr kTj = nTi (n = 2k).
// Run: magma -b magma/wreath_labels.m > k3fq/data/wreath_labels.txt
SetColumns(0);
for k in [1..11] do
    for j in [1..NumberOfTransitiveGroups(k)] do
        H := TransitiveGroup(k, j);
        W := WreathProduct(CyclicGroup(2), H);
        assert Degree(W) eq 2*k and IsTransitive(W);
        i, n := TransitiveGroupIdentification(W);
        printf "%o %o %o %o\n", k, j, n, i;
    end for;
end for;
quit;
