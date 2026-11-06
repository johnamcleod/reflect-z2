# Open questions and unverified steps

Each item lists what we would need to see before we stop treating it as an assumption.

## OQ1: Bugaenko's face theorem over ℤ[√2] (the n → n+1 step)

- **Claim.** If f = f′ ⊕ u with u unimodular of rank 1, then the fundamental polytope of the maximal reflection subgroup of O(f′) is a face of the one for O(f).
  - Source: Bugaenko 1990, Thm 2, as quoted in the thesis (Thm 3.0.8).
  - Consequence used: non-reflectivity of f_d^{(n)} implies non-reflectivity of f_d^{(m)} for all m ≥ n.
- **Status.** We have **not read the primary source**: the PDF is unavailable, and the user said to proceed without it.
  - The thesis applies the theorem over ℤ.
  - Bugaenko's own setting was rings of integers of totally real fields, so ℤ[√2] is probably covered, but the hypotheses are unchecked.
- **Partial check.** We will re-derive Mark's Prop. 3.1 argument over ℤ[√2] (see below) and record here whether it goes through.

## OQ2: elliptic subdiagrams ↔ faces for the infinitely-faceted chamber

- **Use.** The certificate argument (methods.md, Phase 3, step 4) uses the following: if a set S of n facets of the W-chamber P has an elliptic Gram matrix, then the mirrors of S meet in a vertex of P, and the cone of P at that vertex is cut out by S alone.
- **Status.** Vinberg proved this for acute-angled polyhedra; the thesis states it as Prop. 2.2.9 for Coxeter polytopes. For the chamber of a discrete reflection group with possibly infinitely many facets, local finiteness of the mirror arrangement should make the argument go through. That still has to be confirmed.
