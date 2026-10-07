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
  y^2 + y = λ_5 x^5 + λ_3 x^3 + λ_1 x; the conjectured closed forms A_6 = q^9 + (−1)^{k+1}(q−1)q^5,
  X = (−1)^{k+1} q^6 + q^5 are verified exactly for q ≤ 128 (pure cubic part proved for all k).
  Y = q^3 m_k with m_k = 13, −1, −35, 47, 733 (k = 1..5) is a genuine exponential sum (no elementary closed form).
* **Scaling.** Var_{G^2}/q^{d+1} ≍ q^2 at d = 6 (not q^4); C_6 = lim Var_{G^2}/q^9 = 1 holds if m_k = o(q^{5/2}).

See `paper/paper_part2.pdf` (and the status table in its last section) for exactly what is proved, verified, or conjectural.

## Layout

```
paper/     paper_part2.tex / .pdf, figures (fig_profile.pdf, fig_scaling.pdf) and make_figures.py
code/      verification scripts (numpy; the Hayes-FFT driver also needs numba + scipy via the submodule)
data/      exact outputs quoted in the paper (text and CSV)
reports/   working notes (Chinese) that preceded the paper, parts 1–4
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
cd ../paper && python make_figures.py && pdflatex paper_part2.tex
```

`legendre_d6_qf.py` for k = 7 takes about four minutes; everything else runs in seconds to a minute.

## Citation

R. Chen, *Prime polynomials in Legendre intervals over F_{2^k}[t], II: exact variances via quadratic towers and
supersingular Artin–Schreier curves*, preprint (2026), DOI 10.5281/zenodo.23224897.

## License

Paper, data, figures and reports: CC BY 4.0. Code: MIT. See `LICENSE`.
