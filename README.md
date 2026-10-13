# HC-ts-cyclic

Code and data supporting the note *Cyclic symmetry and the Hodge conjecture for the fourfolds F(x)+F(y)=0*
(Kevin Saller, 2026), DOI: to be assigned.

**Subject.** For F very general in a cyclic family C[x0,x1,x2^n]_d (forms of degree d in which x2 occurs only
through x2^n), the note studies the fourfolds X_F = {F(x) + F(y) = 0} in P^5 and the surfaces
Y_F = {F(x) + u^d = 0} in P^3. Theorems A and B are proved by hand. Theorem C, the survey of 86 families of degree
6 to 10, is computer-assisted.

**Role of the computations.** The note's section "Computations" lists seven items. Items 1-5 and 7 are checks of
statements that the note proves by hand. Item 6, the survey, is part of the proof of Theorem C.
[`README.txt`](README.txt) gives the requirements, the commands, what each script checks, and the result of the
rerun made when this repository was assembled.

**Status:** not yet peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). The survey computation was done twice, by separately written implementations; the theory checks were
partly repeated by separately written programs. These were also produced with AI assistance.

## Licences

- **Code** (`code/`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs in `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
