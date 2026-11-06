# Monotonicity lemma

Source: the user's notes (`additional_notes.md` §3), written up here with the checks made
while adopting it. It generalises the n → n+1 argument of OPEN_QUESTIONS OQ1, which covers
only a norm-1 summand.

## Statement

Let K be a totally real field and L an admissible hyperbolic O_K-lattice. Let c ∈ O_K be
totally positive, and put L′ = L ⊥ O_K·u with (u,u) = c. Here W(·) is the group generated
by reflections in *all* roots of the lattice, and "reflective" means [O′(L) : W(L)] < ∞.

**Lemma.** If L′ is reflective, then L is reflective. Moreover, a fundamental polytope of
W(L) is the facet u^⊥ of a fundamental polytope of W(L′).

**Corollary.** For f_d^n = −d·x₀² + x₁² + … + xₙ², take c = 1 and u = v_{n+1}. Then:
- the set of n for which f_d^n is reflective is an interval {2, …, N(d)};
- non-reflective at n₀ implies non-reflective for every n ≥ n₀;
- reflective for some n implies reflective at n = 2.

## Proof

Write x ∈ L′ as x = x_L + a(x)·u with x_L ∈ L and a(x) ∈ O_K. Then

  (∗) (x,u)/c = a(x) ∈ O_K for every x ∈ L′.

1. **The roots of L are exactly the roots of L′ orthogonal to u.**
   - Let r be a root of L with norm m. Then (r,x) = (r,x_L) for every x ∈ L′, so
     2(r,x)/m ∈ O_K. Also r stays primitive in L′, so r is a root of L′.
   - Conversely, let r′ be a root of L′ with r′ ⊥ u. By (∗), a(r′) = 0, so r′ ∈ L. It is
     primitive there, and is a root of L.
   - The vector u is itself a root of L′, since 2(u,x)/c = 2a(x).
2. **A chamber with a facet on u^⊥.** Suppose L′ is reflective. Then W′ = W(L′) has
   finite-volume chambers with finitely many faces, and r_u ∈ W′. Choose a chamber P′
   with a facet F = P′ ∩ u^⊥.
3. **No odd angles at F.**
   - Suppose a facet with outward root e meets F at angle π/k with k odd. In the dihedral
     group ⟨r_e, r_u⟩ of order 2k, r_u and r_e are conjugate: w·r_u·w⁻¹ = r_e for some w.
   - Then wu ∈ ℝe and |wu| = |u|. The angle between e and u is π − π/k, so
     (wu,u) = ±c·cos(π/k).
   - Since wu ∈ L′, (∗) gives cos(π/k) ∈ O_K, so cos(π/k) is an algebraic integer.
   - But 2cos(π/k) = ζ⁻¹(1+ζ²) with ζ = e^{iπ/k}. For odd k, 1 + ζ² = (1−ζ⁴)/(1−ζ²) is
     a cyclotomic unit. So 2cos(π/k) is a unit, and cos(π/k) is not integral.
     Contradiction.
   - *Checked while adopting:* the conjugating w gives wu parallel to e itself (not to
     some other root in the dihedral group), because it carries the mirror u^⊥ onto e^⊥.
4. **Every ridge of F lies in a mirror of W(L).**
   - Take a ridge R = F ∩ F_j. Its angle is π/k with k even.
   - The dihedral group at R contains the reflection in the hyperplane through R
     perpendicular to u^⊥. That reflection lies in W′, so its root is a root of L′
     orthogonal to u, hence a root of L by step 1. (When k = 2, take the root of F_j.)
5. **F is a chamber of W(L).**
   - Every mirror of W(L) is r^⊥ ∩ u^⊥ with r^⊥ a mirror of W′ other than u^⊥ (step 1).
     So no mirror of W(L) meets the interior of the facet F.
   - Hence F lies in a single chamber C of W(L).
   - Each face hyperplane of F is a mirror of W(L) (step 4), so C ⊆ F, and therefore
     C = F.
   - So W(L) has a finite-volume chamber, and lies with finite index in O′(L). ∎

## Computational checks built on the lemma

- `tests/test_extension.py`:
  - no odd angle at the norm-1 root −v_{n+1} in any terminated run;
  - certificates are consistent with monotonicity in n.
- `tests/test_facet.py`, calibration (e): whenever (d, n+1) terminates, the facet
  −v_{n+1}^⊥ reproduces the polytope for (d, n). Its walls come from e itself when
  k = 2, and from the primitive projection of e onto u^⊥ when k is even.
