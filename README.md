# Prime polynomials in Legendre intervals over F_{2^k}[t], II

Exact variances via quadratic towers and supersingular Artin–Schreier curves.
Sequel to [legendre-intervals-explicit-formula](https://github.com/Ruqing1963/legendre-intervals-explicit-formula)
(R. Chen, DOI 10.5281/zenodo.23219609), whose Open Problem 8.2 asked why the variance of the prime-power count
over Legendre intervals `{f^2 + s : deg s <= d}` in even characteristic departs from the Keating–Rudnick prediction.

Preprint DOI: 10.5281/zenodo.23224897.

## Results in one paragraph

* **d = 4, all q = 2^k (proved).** Closed formula for Ψ on every class of the Hayes group G_3 via the tower
  F_q ⊂ F_{q^2} ⊂ F_{q^4} ⊂ F_{q^8}. Consequences: Var_{G^2}(Ψ)/q^5 = q(q−1), Var_G(Ψ)/q^5 = 3 − 3/q,
  ratio q^2/3; the deviation is supported on {c_1 = 0}; ρ(ε) = 1 on the q−1 characters of conductor u^2 and
  (q−3)/(3(q−1)) otherwise for odd k; an elementary proof of Carlitz's cubic evaluation over F_{2^{8k}}.
* **d = 6, structure (proved).** The Legendre profile takes exactly three values, governed by [c_2 = 0] and
  Tr(c_4/c_2^2); it is determined by three sums A_6, X, Y over U = {Tr α = Tr α^3 = Tr α^5 = 0} ⊂ F_{q^12}.
* **d = 6, constants.** A_6 − q^9 and X are sums of Frobenius traces of the supersingular curves
  y^2 + y = λ_5 x^5 + λ_3 x^3 + λ_1 x. The pure cubic and pure quintic parts are evaluated for all k.
  For odd k the mixed part is evaluated completely (Theorem 6.12): W(λ_3, λ_5) = 2q^7 or 0 according as the
  Frobenius permutes the roots of the Gold-type quintic v^5 + v + c, c = (λ_5^3/λ_3^5)^{1/2}, as a 4-cycle or not,
  and the exact count of 4-cycles, (q−2)/2, is proved; hence **A_6 = q^9 + (q−1)q^5 for all odd k** (Corollary 6.13).
  The conjectured closed forms A_6 = q^9 + (−1)^{k+1}(q−1)q^5 (even k) and X = (−1)^{k+1} q^6 + q^5 (all k)
  are verified exactly for q ≤ 128. For X (odd k) the sums W' are order-4 character sums whose Frobenius classes
  on the 4-dimensional kernel of L' = L + x^4 run through S_6 ≅ Sp_4(F_2) with exact polynomial frequencies
  (Table 3); the structure, W'(0,0) and the pure parts are proved, the sign pattern is Conjecture 6.19.
  Y = q^3 m_k with m_k = 13, −1, −35, 47, 733 (k = 1..5) is a twisted Frobenius trace on a fixed nine-dimensional
  complete intersection X = V(e_1, e_3, e_5) ⊂ A^12 (no elementary closed form); Y = N_0 − N_1 is the difference of
  the twisted point counts of the 8-folds X ∩ {e_4 = 0}, X ∩ {e_4 = 1}, and C_6 = 1 is equivalent to
  #X_1^(σ)(F_q) = q^8 + q^5 − q^4 + o(q^{9/2}) (Proposition 7.9).
* **Scaling.** Var_{G^2}/q^{d+1} ≍ q^2 at d = 6 (not q^4); C_6 = lim Var_{G^2}/q^9 = 1 holds if m_k = o(q^{5/2}).
* **Spectral folding and general d.** The amplification is the coherence of the Fourier coefficients along the fibres
  of G^ → (G^2)^ (cosets of the 2-torsion G[2], of order q^{d/2}); at d = 6, Var_{G^2} = [(A_6−q^9)^2+(q−1)X^2]/q^4 + q^3(q−1)m_k^2
  exactly. The translation orbit K_d of the identity class has order 1 or q by Lucas' theorem. Exact variances for
  d ≤ 10 show the q^2 amplification is special to d = 4, 6: at d = 8 the Legendre/full-group ratio is ≈ 2.3 for q = 4, 8.

See `paper/paper_part2.pdf` (and the status table in its last section) for exactly what is proved, verified, or conjectural.

## Layout

```
paper/     paper_part2.tex / .pdf and eight figures (fig_profile, fig_scaling from make_figures.py; fig_gold_cycle,
           fig_singular_planes, fig_S6_classes, fig_spectral_bounds, fig_spectral_folding, fig_even_d_bifurcation
           from code/generate_figures.py)
code/      verification scripts (numpy; the Hayes-FFT driver also needs numba + scipy via the submodule)
data/      exact outputs quoted in the paper (text and CSV)
reports/   working notes (Chinese) that accompany the paper, parts 1–7
legendre-intervals-explicit-formula/   git submodule: the paper-I repository (its code/ is imported)
```

## Reproducing

```bash
git clone --recurse-submodules https://github.com/Ruqing1963/legendre-intervals-part2.git
cd legendre-intervals-part2/code
python legendre_var_check.py 2,4 4,4 8,4          # brute force over F_{q^8}, Table-4 values (q <= 8)
python legendre_d4_formula.py 7                   # d=4 closed formula, all checks, k = 1..7; Var_G identity to q = 2^22
python legendre_q16_hayes.py 16,4 32,4 64,4 128,4 16,6 32,6   # paper-I Hayes FFT, class-by-class comparison
python legendre_d6_tower.py                       # d=6 tower identities, three-value theorem, constants A6, X, Y
python legendre_d6_qf.py 1 2 3 4 5 6 7            # exact A6, X via F_2-quadratic forms (radical + Arf), q <= 128
python legendre_d6_qf_mixed.py 3 5 7              # mixed terms W per class j, radicals, all S(lambda_1) (odd k)
python quintic_check.py                           # factorisation types of v^5+v+c, k = 3..11 (Proposition 6.8)
python kasami_sign.py 3 5 7 9                     # Theorem 6.12 via L-polynomials of y^2+y = x^5+t x^3+b x
python legendre_d6_Wp_pattern.py 3 5 7            # W' (the sums behind X) by Frobenius class on ker L' (Table 3)
python L_function_spectral.py                     # m_k = 13,-1,-35,47,733: recurrences, Weil moduli, 2-adic constraints, N_0, N_1
python fourier_spectrum_folding.py                # spectral folding G^ -> (G^2)^: fibre sums, energies, coherence (Section 8)
python general_d_scan.py 2,8 4,8 8,8 2,10 4,10    # Legendre variance, identity class and its orbit K_d for d = 8, 10 (Section 9)
python generate_figures.py                        # six vector figures: gold cycle, singular planes, S_6 classes, spectral bounds, folding, even-d bifurcation
cd ../paper && python make_figures.py && pdflatex paper_part2.tex
```

`legendre_d6_qf.py` for k = 7 takes about four minutes and `kasami_sign.py 9` about a quarter of an hour; everything else runs in seconds to a minute.

## Citation

R. Chen, *Prime polynomials in Legendre intervals over F_{2^k}[t], II: exact variances via quadratic towers and
supersingular Artin–Schreier curves*, preprint (2026), DOI 10.5281/zenodo.23224897.

## License

Paper, data, figures and reports: CC BY 4.0. Code: MIT. See `LICENSE`.
