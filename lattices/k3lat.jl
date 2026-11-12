# Candidate pairs (N, phi) = (NS(X over Fqbar), Frobenius) for classes of K3 type over F_q.
#
#     julia --project=lattices lattices/k3lat.jl INPUT.jsonl OUTPUT.jsonl
#
# INPUT is written by prep.py (one class per line).  For a class with Weil polynomial
# P = (1-t)^rho C(t) S(t), C = prod Phi_d (d > 1), E = lcm of the d:
#
#  1. Coinvariant part (C_N, phi_C) = ((N^phi)^perp, phi): negative definite, even, rootless
#     (Lemma 3.1 of Auel-Shaikh, arXiv:2611.13468), char poly prod Phi_d, isometry of order E.
#     By their Theorem 4.1 (m = E), |disc C_N| divides gcd(E^rho, C(1)) * Res((t^E-1)/(t-1), S),
#     and it is at least delta_r (r = rank, Section 3).  Enumerated genus by genus with
#     Oscar's enumerate_classes_of_lattices_with_isometry.
#  2. Invariant part I = N^phi = NS(X): even, signature (1, rho-1), |disc I| = AT_1 / b^2.
#  3. (N, phi): equivariant primitive extensions of (I, id) and (C_N, phi_C), up to the
#     actions of O(I) and O(C_N, phi_C) (Oscar's equivariant_primitive_extensions).
#  4. Necessary conditions on (N, phi):
#     artin_tate   for every m | E, |Br X_m| = AT_m / |disc N^{phi^m}| is a square integer;
#     brauer_maps  for m | m' | E, the kernel of Br X_m -> Br X_m' is H^1(<phi^m>, N^{phi^m'})
#                  (Lemma 2.1), so |H^1| divides |Br X_m|, and |Br X_m| / |H^1| divides |Br X_m'|;
#     ell_length   for primes l != p, the l-part of D(N) has length <= 22 - rhobar;
#     ell_glue     for l != p, D(N)_l is glued Frobenius-equivariantly to the transcendental
#                  lattice in H^2(Z_l(1)), so S(phi) kills D(N) away from p.
#     The coinvariant lattices of all phi^m lie in C_N, so they are rootless too.
using Oscar
using JSON

# Lower bounds for the least discriminant of an even positive definite rootless lattice of
# rank k (Auel-Shaikh, Section 3; exact for k <= 8).
const DELTA = [4, 12, 32, 64, 128, 192, 256, 256, 298, 317, 312, 281, 230, 177, 128, 89, 58, 37, 24, 13, 8]

const Zx, x = polynomial_ring(ZZ, "x")

is_square_int(n::ZZRingElem) = n >= 0 && is_square(n)

"""Divisors d of n with n/d a square."""
square_cofactor_divisors(n::ZZRingElem) = [n ÷ b^2 for b in divisors(isqrt(n)) if n % b^2 == 0]

flat(M) = [Int(numerator(M[i, j])) for i in 1:nrows(M) for j in 1:ncols(M)]

has_roots(L::ZZLat) = rank(L) > 0 && !isempty(short_vectors(rescale(L, -1), 2))

# ---------------------------------------------------------------------------------------
# caches (shared by all classes of a run)

const COINV = Dict{Tuple{Vector{Int}, ZZRingElem}, Vector{ZZLatWithIsom}}()
const INV = Dict{Tuple{Int, ZZRingElem}, Vector{ZZLatWithIsom}}()
const GLUE = Dict{Tuple{UInt, UInt}, Vector{ZZLatWithIsom}}()

"""Rootless negative definite (C, phi) with char poly prod Phi_d (d in orders) and |disc| D."""
function coinvariant_candidates(orders::Vector{Int}, D::ZZRingElem)
    key = (sort(orders), D)
    haskey(COINV, key) && return COINV[key]
    r = sum(euler_phi(d) for d in orders)
    E = lcm(orders)
    chi = prod(cyclotomic(d, x) for d in orders)
    out = ZZLatWithIsom[]
    for G in integer_genera((0, r), (-1)^r * D; even=true)
        # root_test=true wrongly skips e.g. A2(-2) with Phi_3, so filter ourselves
        for Cf in enumerate_classes_of_lattices_with_isometry(G, E; char_poly=chi)
            has_roots(lattice(Cf)) || push!(out, Cf)
        end
    end
    COINV[key] = out
    return out
end

"""(I, id) for even I of signature (1, rho - 1) and |disc I| = dI."""
function invariant_candidates(rho::Int, dI::ZZRingElem)
    key = (rho, dI)
    haskey(INV, key) && return INV[key]
    out = ZZLatWithIsom[]
    for G in integer_genera((1, rho - 1), (-1)^(rho - 1) * dI; even=true)
        for I in representatives(G)
            push!(out, integer_lattice_with_isometry(I))
        end
    end
    INV[key] = out
    return out
end

function glue(If::ZZLatWithIsom, Cf::ZZLatWithIsom, orders::Vector{Int})
    key = (objectid(If), objectid(Cf))
    haskey(GLUE, key) && return GLUE[key]
    ok, V = equivariant_primitive_extensions(If, Cf; even=true)
    out = ZZLatWithIsom[L for (L, _, _) in V]
    GLUE[key] = out
    return out
end

# ---------------------------------------------------------------------------------------
# conditions

"""|H^1(<s>, L)| for an isometry s of finite order of L (Auel-Shaikh, Lemma 1.3)."""
function h1_order(Lf::ZZLatWithIsom)
    rank(Lf) == 0 && return ZZ(1)
    Linv = invariant_lattice(Lf)
    Cf = coinvariant_lattice(Lf)
    rank(Cf) == 0 && return ZZ(1)
    A2 = abs(det(lattice(Linv)) * det(lattice(Cf)) // det(lattice(Lf)))
    A = isqrt(numerator(A2))
    @assert A^2 == A2
    f = isometry(Cf)
    d = abs(det(identity_matrix(QQ, rank(Cf)) - f))
    h = d // A
    @assert denominator(h) == 1
    return numerator(h)
end

"""N^{phi^m} with the restriction of phi."""
fixed_by_power(Nf::ZZLatWithIsom, m::Int) = kernel_lattice(Nf, x^m - 1)

function check_conditions(Nf::ZZLatWithIsom, c::AbstractDict)
    E = c["E"]; p = c["p"]; AT = c["AT"]; rhobar = c["rhobar"]
    failed = String[]
    ms = divisors(E)
    sub = Dict(m => fixed_by_power(Nf, m) for m in ms)
    B = Dict{Int, ZZRingElem}()
    for m in ms
        dm = abs(det(lattice(sub[m])))
        q, r = divrem(AT[m], numerator(dm))
        if r != 0 || !is_square_int(q)
            push!(failed, "artin_tate")
            break
        end
        B[m] = q
    end
    if isempty(failed)
        for m in ms, m2 in ms
            (m < m2 && m2 % m == 0) || continue
            h = h1_order(sub[m2]^m)
            if B[m] % h != 0 || B[m2] % (B[m] ÷ h) != 0
                push!(failed, "brauer_maps")
                break
            end
        end
    end
    N = lattice(Nf)
    n = rank(N)
    Gz = change_base_ring(ZZ, gram_matrix(N))
    S_ = snf(Gz)
    for l in prime_divisors(abs(det(Gz)))
        l == p && continue
        if count(i -> S_[i, i] % l == 0, 1:n) > 22 - rhobar
            push!(failed, "ell_length")
            break
        end
    end
    # S(phi) kills D(N) away from p: S(F) maps the dual lattice into N[1/p]
    F = isometry(Nf)
    SF = zero_matrix(QQ, rank(N), rank(N))
    Fi = identity_matrix(QQ, rank(N))
    for s in c["S"]
        SF += s * Fi
        Fi *= F
    end
    W = inv(gram_matrix(N)) * SF
    for i in 1:n, j in 1:n
        den = denominator(W[i, j])
        while den % p == 0
            den = den ÷ p
        end
        if den != 1
            push!(failed, "ell_glue")
            break
        end
    end
    return failed
end

# ---------------------------------------------------------------------------------------

function process_class(c::AbstractDict)
    rho = c["rho"]; rhobar = c["rhobar"]; orders = Vector{Int}(c["orders"]); E = c["E"]
    AT = Dict(parse(Int, k) => ZZ(parse(BigInt, v)) for (k, v) in c["AT"])
    Res = Dict(parse(Int, k) => ZZ(parse(BigInt, v)) for (k, v) in c["Res"])
    q = c["q"]; p = Int(prime_divisors(ZZ(q))[1])
    cc = Dict("E" => E, "p" => p, "AT" => AT, "rhobar" => rhobar, "S" => c["S"])
    stats = Dict{String, Int}()
    rows = Dict{String, Any}[]
    t0 = time()
    dIs = square_cofactor_divisors(AT[1])
    pairs = ZZLatWithIsom[]
    if E == 1
        for dI in dIs
            append!(pairs, invariant_candidates(rho, dI))
        end
        stats["coinvariant_discs"] = 0
    else
        r = rhobar - rho
        chiC1 = prod(cyclotomic(d, x)(1) for d in orders)
        g = gcd(ZZ(E)^rho, chiC1)
        Ds = [D for D in divisors(g * Res[E]) if D >= DELTA[r]]
        stats["coinvariant_discs"] = length(Ds)
        ncoinv = 0
        for D in Ds
            Cs = coinvariant_candidates(orders, D)
            ncoinv += length(Cs)
            isempty(Cs) && continue
            for dI in dIs
                # |A|^2 = dI D / |disc N| with |disc N| = AT_E / c^2 and |A| | g
                ok = any(a -> (dI * D) % (a^2) == 0 && (AT[E] % ((dI * D) ÷ a^2) == 0) &&
                              is_square_int(AT[E] ÷ ((dI * D) ÷ a^2)), divisors(g))
                ok || continue
                for If in invariant_candidates(rho, dI), Cf in Cs
                    append!(pairs, glue(If, Cf, orders))
                end
            end
        end
        stats["coinvariant_lattices"] = ncoinv
    end
    stats["candidates"] = length(pairs)
    for Nf in pairs
        failed = check_conditions(Nf, cc)
        for f in failed
            stats[f] = get(stats, f, 0) + 1
        end
        isempty(failed) || continue
        N = lattice(Nf)
        I = lattice(invariant_lattice(Nf))
        push!(rows, Dict("gram" => flat(gram_matrix(N)), "frobenius" => flat(isometry(Nf)),
                         "disc" => string(det(N)),
                         "invariant_gram" => flat(gram_matrix(I)), "invariant_disc" => string(det(I)),
                         "coinvariant_disc" => string(det(lattice(coinvariant_lattice(Nf))))))
    end
    stats["lattices"] = length(rows)
    return Dict("label" => c["label"], "census" => c["census"], "realizable_k3" => c["realizable_k3"],
                "rho" => rho, "orders" => orders, "E" => E, "stats" => stats, "lattices" => rows,
                "secs" => round(time() - t0, digits=3))
end

function main(infile, outfile)
    done = Set{String}()
    if isfile(outfile)
        for line in eachline(outfile)
            push!(done, JSON.parse(line)["label"])
        end
    end
    open(outfile, "a") do OUT
        for line in eachline(infile)
            c = JSON.parse(line)
            c["label"] in done && continue
            res = process_class(c)
            println(OUT, JSON.json(res))
            flush(OUT)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2])
end
