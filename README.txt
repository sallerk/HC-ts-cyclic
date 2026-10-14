Supporting code and data for the note
"Cyclic symmetry and the Hodge conjecture for the fourfolds F(x)+F(y)=0"
=========================================================================

Status: accompanies the note (Kevin Saller, 2026), DOI to be assigned; not peer reviewed.

Setting. F is a very general member of a cyclic family C[x0,x1,x2^n]_d (forms of degree d in x0, x1, x2 in which
x2 occurs only through x2^n; d >= 6, 2 <= n <= d, n | d or n | d-1), X_F = {F(x) + F(y) = 0} in P^5 and
Y_F = {F(x) + u^d = 0} in P^3. The cohomology splits into pieces E_(a,chi) (a = character of u, chi = character of
x2). The "special point" is the Fermat form if n | d and the chain form x0^d + x1^d + x0 x2^(d-1) if n | d-1. The
survey (Theorem C) concerns 86 families of degree 6 to 10 with a cyclic group of diagonal symmetries.

Role of the computations. The note's section "Computations" lists seven items. Items 1-5 and 7 are CHECKS of
statements proved by hand in the note (Theorems A and B and the lemmas behind them). Item 6, the survey of the 86
families, is PART OF THE PROOF of Theorem C. It was computed twice, by separately written implementations
(components survey_enumeration_first and survey_second below). The theory checks were partly repeated by
separately written programs (component independent_checks).

Not checked by any script: the step from a computation modulo a prime at a random member (or at the special
point) to a statement about the very general complex member (semicontinuity, reduction modulo p), the lemmas that
turn computed bounds into Hodge classes, and the algebraicity of the classes are arguments of the note.

Labels. Script headers and some output lines use working labels of the author's notes (names of working
documents, of their sections, lemmas and rules, and short names of programs and of the people or runs that
wrote them). These documents are not part of this repository, and the labels are not numbers of the note.

Licence: the code in code/ is under the MIT licence; the outputs in outputs/ and this text are under CC BY 4.0
(see the LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3, tested with Python 3.12.9 on Windows 10. Third-party packages (tested versions):
  python-flint 0.9.0: refU_ring.py; enumeration/survey.py, ord3.py; verification/v1lib.py, probe.py;
                      residual/v1lib.py, ord3.py; survey_second/z3lib.py
  numpy 2.1.3:        independent_checks/cross_pairs/ref_xu.py; verification/probe.py;
                      survey_second/z3lib.py, z3run.py, z3tests.py
  sympy 1.13.1:       imported inside one function of survey_second/z3lib.py
All other scripts use the standard library only. Several scripts start worker processes: u1check.py (W, an
argument), refU_comm.py and refU_proofmoves.py (6), z3run.py (nproc, an argument), v1run.py and z1tri.py (8 by
default; environment variables V1_WORKERS and Z1_WORKERS), ord3drive.py (up to 8).


2. Layout and how to run
------------------------
code/<component>/...      the scripts (original file names)
outputs/<component>/...   the reference outputs; files named *.log or *_stdout.txt are saved screen output

All commands below are given from the top folder of this repository; OUT is the output file of your choice.

Where the scripts write. The scripts in code/theory_checks, code/independent_checks, code/appendix_B and
code/table take their output file as an argument or print to the screen. The scripts in
code/survey_enumeration_first/* and code/survey_second write into the folder "../outputs" relative to the script
file, that is code/survey_enumeration_first/outputs and code/outputs (not into outputs/). v1tests.py and
z3tests.py do not create that folder; create it first. Scripts that process stored results (summarize.py,
orbit_diag.py, secondorder.py, z1final.py, z1ord3plan.py, cond.py, curvecheck.py, z1curve.py, jcontrol.py,
z1show.py, z3tally.py, z3spot.py) read them from the same folder, so copy the stored outputs there first.
Not runnable as shipped (they read files under the folder names of the original working tree, or files that are
not included; their outputs are kept as records):
v1compare.py, z3compare.py, z3compare2.py, ref_cmp_table.py; compareB.py needs, for its comparison, data of an
earlier computation that is not included.


2.1 theory_checks/hodge_pieces (items 1 and 2)
    python -I -B code/theory_checks/hodge_pieces/u1check.py formulas DMAX OUT
    python -I -B code/theory_checks/hodge_pieces/u1check.py conn DMIN DMAX W OUT
    python -I -B code/theory_checks/hodge_pieces/u1check.py comm DMIN DMAX W OUT
    python -I -B code/theory_checks/hodge_pieces/u1check.py gen  DMIN DMAX W OUT
  formulas: for every pair (d,n) with 6 <= d <= DMAX: the closed form of c_n(D,g) (number of monomials of degree
    D with x2-exponent = g mod n) against a direct count, Hodge numbers, zero / pure / one-sided pieces, pure
    orbits, separation by Hodge numbers, Y(d,n), C(d,n) and dim Hdg^2(X_F) by the block rules against the closed
    formulas, and the cross-u closure count.
  conn: at the special point, the monomial basis of the Jacobian ring, the Hodge numbers from this basis against
    the formulas, and connectivity of the monomial graph of every non-pure weight-2 piece and every two-sided
    weight-1 piece. comm: the full first-order commutant of these pieces (exact rational arithmetic).
    gen: the generation step (each basis monomial of level q+1 is a basis monomial of level q times a family
    monomial).
  Stored: formulas_d6_100.txt (747 pairs, failures 0; screen output formulas_run.log), formulas_d6_60.txt (a run
    to d = 60, 382 pairs, failures 0); conn_d6_36.jsonl + conn_d37_60.jsonl (382 pairs), comm_d6_13.jsonl +
    comm_d14_26.jsonl (118 pairs), gen_d6_40.jsonl + gen_d41_60.jsonl (382 pairs); the *_run.log and *_run2.log
    files report failures=0 for each run.

2.2 theory_checks/cross_pairs (items 1, 4 and 5)
    python -I -B code/theory_checks/cross_pairs/w2xu.py DMAX OUT
    python -I -B code/theory_checks/cross_pairs/w2dim.py DMAX OUT
    python -I -B code/theory_checks/cross_pairs/w2fam.py DMAX_FAM DMAX_M345 OUT
  w2xu.py: for every pair with d <= DMAX, the cross pairs (pieces with different u-characters in non-pure orbits)
    that are not separated by the Hodge numbers of their Galois conjugates, against the predicted classification
    and against a criterion computed independently from fractional parts; for n = d the set KR(d) (ordered
    pairs), its swap part and the rest. Stored: xu_d6_100.txt (= xu_d6_100.log), xu_d6_30.txt. The last line
    "# DONE bad=1" of both comes from d = 6 only: the line "6 6 ... SWAPCOUNT 8 16" compares the swap count 8 with
    the closed form (d-1)(d-3) + [2|d] = 16, which the script tests at every d but which holds only for d >= 7.
    kr_table.txt lists d:|KR(d)|(non-swap part) for d = 6..100; it is a reformatting of the KR_ordered and
    nonswap values of the n = d lines of xu_d6_100.txt (checked at assembly: identical). |KR(d)| is the number
    k(d) of item 5.
  w2dim.py: the closed form of dim Hdg^2 (without the KR term) and dim != 2k; dim_d6_100.txt: 747 pairs, bad=0.
  w2fam.py: the explicit families of non-separated cross pairs (n = d/2 and the sporadic pairs) up to DMAX_FAM and
    the classes m in {3,4,5} up to DMAX_M345; fam_d400_m345_d180.txt (= fam_run.log): bad=0.

2.3 theory_checks/planes (item 3)
    python -I -B code/theory_checks/planes/planes_gram.py > OUT
  No arguments (6 <= d <= 20 and all n are fixed in the script). Gram matrix of h^2 and the n*d planes
  {y = eta*gamma(x)} on X_F from the intersection numbers, rank modulo 2^61-1, against the expected rank and
  dim Hdg^2(X_F). Stored: planes_gram_d6_20.txt (78 pairs, last line "bad 0"); this output was produced when the
  repository was assembled (about 40 s).

2.4 independent_checks/pieces_commutant (items 1 and 2; separately written check programs)
    python -I -B code/independent_checks/pieces_commutant/refU_formulas.py DMAX_HS DMAX_ALL OUT
    python -I -B code/independent_checks/pieces_commutant/refU_ring.py DMIN DMAX OUT
    python -I -B code/independent_checks/pieces_commutant/refU_comm.py DMIN DMAX OUT [control|sub2]
    python -I -B code/independent_checks/pieces_commutant/refU_proofmoves.py DMIN DMAX OUT
  refU_formulas.py: c_n against a direct count, the Hilbert series of the Jacobian ring at the special point
    (d <= DMAX_HS), Hodge numbers, pure pieces, Y, C, dim and the non-separated cross pairs (d <= DMAX_ALL).
    out_formulas.txt: DMAX_HS = 40, DMAX_ALL = 100; 747 pairs, failures 0; 152 pairs with non-separated cross
    pairs (all n = d, all n = d/2, and (18,6), (24,6), (36,12), (48,12), (54,18), (72,18), (72,24), (90,30),
    (96,24)).
  refU_ring.py: the Jacobian ring at the chain point by exact linear algebra over Z; out_ring.txt: d = 5..10 OK.
  refU_comm.py: the full first-order commutant of every non-pure weight-2 and two-sided weight-1 piece at the
    special point, with all monomials as unknowns, modulo 2^61-1; out_comm_6_11_control.jsonl (with positive
    controls), out_comm_12_22.jsonl, out_comm_23_32.jsonl; out_comm_sub2_control.jsonl (positive control: the
    subfamily C[x0,x1,x2^(2n)]); the two out_comm_*.log files are empty.
  refU_proofmoves.py: the explicit moves of a proof step for every ordered pair of level-0 basis monomials;
    out_proofmoves_6_40.txt: 41,247,998 pairs, failures 0.

2.5 independent_checks/cross_pairs (items 1, 4 and 5; separately written check programs)
    python -I -B code/independent_checks/cross_pairs/ref_xu.py full DMIN DMAX OUT [mlist]
    python -I -B code/independent_checks/cross_pairs/ref_xu.py fam DMAX OUT
    python -I -B code/independent_checks/cross_pairs/ref_summ.py FILE [FILE ...]
    python -I -B code/independent_checks/cross_pairs/ref_lemmas.py
  ref_xu.py full: from the definitions only (direct count c_n, Hodge numbers, Galois action), for every pair with
    DMIN <= d <= DMAX the non-separated cross pairs, the dimension check and KR(d); stored for d <= 150 in
    out_full_6_110.txt, out_full_111_130.txt, out_full_131_142.txt, out_full_143_150.txt, and for 151 <= d <= 360
    (restricted list of classes, see the file name) in out_ext_m3456_151_360.txt; out_test.txt is a short
    earlier run in the same format. ref_xu.py fam: explicit families; out_fam_1000.txt (d <= 1000, bad 0).
  ref_summ.py: summaries summ_full_6_150.txt (1240 pairs, 0 rows with a failed check; non-separated cross pairs
    exactly for n = d, n = d/2 and the pairs (18k,6k), (24k,6k) listed there; the table d:|KR(d)|(kappa) for
    d <= 150; kappa = 0 iff gcd(d,6) = 1; kappa = 48(p-1) for d = 2p and 24(p-1) for d = 3p, p >= 11 prime) and
    summ_ext.txt (d = 152..360).
  ref_lemmas.py: the closed form N* (all pieces, d <= 150) and further lemma checks; out_lemmas.txt.
  ref_cmp_table.py compares a table of a working manuscript (not included) with a summary file.

2.6 survey_enumeration_first (item 6: enumeration and first implementation)
  enumeration/ (primes near 3*10^4, e.g. 30097 and 30169):
    python -I -B code/survey_enumeration_first/enumeration/survey.py d [only=NAME] [allub1=J1,J2] [max=K] [tag=TAG]
    python -I -B code/survey_enumeration_first/enumeration/summarize.py TAG d1 d2 ...
    python -I -B code/survey_enumeration_first/enumeration/ord3.py FAMILY d i p seed narcs order [small]
    python -I -B code/survey_enumeration_first/enumeration/orbit_diag.py TAG d1 d2 ...
    survey.py enumerates the families of degree d (deduplicated by monomial set; Hilbert function checked at
    random members modulo two primes) and computes lower and upper bounds for the Hodge classes block by block.
    Stored (final run, TAG = _v4): survey_d<d>_v4.json/.txt and run_d<d>_v4.log for d = 6..10; summary_v4.txt
    (screen output of summarize.py _v4 6 7 8 9 10), summary_v4_d678.txt, summary_v4_d9.txt; ord3_*.json (orders
    1-3 for single pieces, written by ord3.py when called from secondorder.py) and secondorder_d678.txt,
    secondorder_d8_d10.txt (screen output of secondorder.py); orbit_diag_v4*.txt; survey_d7_chk_v4_d7.json/.txt (a check run at d = 7).
  verification/ (a separate implementation of the same tests; primes above 5*10^4, a second run above 7*10^4):
    python -I -B code/survey_enumeration_first/verification/v1run.py [all|...]    (V1_START, V1_WORKERS)
    python -I -B code/survey_enumeration_first/verification/v1tests.py
    python -I -B code/survey_enumeration_first/verification/v1robust.py
    v1lib.py: Jacobian rings mod p, Griffiths reduction, flat Taylor expansion, hull ranks, Higgs fields,
    first-order block bounds. Stored: v1_<family>_d<d>.json, v1_summary_all.txt, v1_summary_sel.txt; the same
    with suffix _s70000 (V1_START=70000), compared in robustness_s70000.txt; unit_tests.txt (Hilbert function
    against the Koszul formula, Taylor coefficients along trivial arcs); compare_S1.txt (v1compare.py).
    probe.py, s1probe.py, s1probe2.py, v1inspect.py, v1orbits.py are small inspection helpers.
  residual/ (the families left open above; primes near 5*10^4):
    python -I -B code/survey_enumeration_first/residual/z1tri.py [all|...]      (Z1_WORKERS)
    python -I -B code/survey_enumeration_first/residual/ord3drive.py fam:d:i [fam:d:i ...]
    python -I -B code/survey_enumeration_first/residual/z1final.py
    z1tri.py imports the local copies v1lib.py (byte-identical to verification/v1lib.py) and v1run.py (a modified
    copy); ord3drive.py runs the local ord3.py (a modified copy of enumeration/ord3.py) at two primes; refU_comm.py
    is a byte-identical copy of independent_checks/pieces_commutant/refU_comm.py. Stored, named after the scripts:
    z1_<family>_d<d>.json, z1_summary*.txt, z1_stdout.txt (z1tri.py); ord3_*.json, ord3log_*.txt (ord3drive.py);
    z1o3_*.json (z1ord3plan.py); jt_*.json, jtest_stdout.txt (jtest.py); jcontrol.txt, cond.txt, curvecheck.txt,
    z1curve.txt; z1_final_table.txt (z1final.py, the merged result); v1_*.json, v1_summary_sel.txt and
    v1rerun_stdout.txt (v1run.py of this folder).

2.7 survey_second (item 6: second implementation, primes between 10^5 and 1.02*10^5)
    python -I -B code/survey_second/z3tests.py
    python -I -B code/survey_second/z3run.py [--p0|--p1] [nproc] [family:d ...]
    python -I -B code/survey_second/z3tally.py
  z3lib.py: Jacobian ring mod p by degree and character, Griffiths reduction, pieces, Higgs and multiplication
    maps, first-order (Hom-)commutants, flat Taylor coefficients up to order 4 along random arcs, hull ranks,
    certificates. z3tests.py: unit tests (Jacobian ring against the Hilbert series, Taylor coefficients along
    trivial arcs); unit_tests.txt. z3run.py: per family the lower and upper bound of Y for every block orbit and
    whether the family is closed (z3_<family>_d<d>.json/.txt; with --p0 or --p1 only one of the two primes,
    files *_only_p0 / *_only_p1; summary z3_summary_<tag>.txt). Without family arguments it runs its built-in
    list (31 families and the two controls cyc12_017 d6 and cyc14_018 d7, which stay open here and are treated
    in Appendix B of the note).
  outputs/survey_second/main: that list (both primes and each prime alone; run_all.log, run_p0.log, run_p1.log),
    unit_tests.txt, z3_table.txt (z3tally.py), z3_arcs.txt (z3arcs.py), and the comparisons with the first
    implementation z3_compare.txt, z3_compare_allblocks.txt, z3_spot.txt.
  outputs/survey_second/further_families: the other survey families, computed with byte-identical copies of
    z3lib.py and z3run.py (included once, in code/survey_second); run.log is the screen output of the run of the
    15 families that the first implementation alone had closed, run2.log that of the run of the remaining 38
    families (all closed).

2.8 table (item 6)
    python -I -B code/table/survey_table2.py [OUT]
    python -I -B code/table/survey_crosscheck.py
    python -I -B code/table/check_primes.py
  survey_table2.py reads outputs/survey_enumeration_first/enumeration/summary_v4.txt and the z3_<family>_d<d>.json
  files in outputs/survey_second/{main,further_families}; it prints the number of rows per route, missing runs,
  inconsistent runs and the primes used, and the LaTeX rows of the survey table (also written to OUT).
  survey_crosscheck.py compares Y of the two implementations on every family closed by both, and Y of the
  families cyc<n>_001 with the formula of Theorem A. check_primes.py checks that every run of the second
  implementation used two distinct primes.
    python -I -B code/table/smooth_criterion_enum.py ROWS
  smooth_criterion_enum.py (added after a referee report): completeness of the list. It enumerates the families
  S_{n,w}, 6 <= d <= 10, n <= 2d, up to permutations of the coordinates, with the exact criterion of Section 9.1
  of the note for smooth general members ((i) no coordinate divides every monomial; (ii) for each k, x_k^d or some
  x_k^(d-1) x_j is in S), and compares with the rows of the survey table read from ROWS, the file written by
  survey_table2.py ROWS. Stored: outputs/table/survey_table_rows.tex (that file) and
  outputs/table/smooth_criterion_enum.txt (86 families with smooth general member, 86 table families, no
  difference, all with moduli; "RESULT: complete").

2.9 appendix_B (item 7)
    python -I -B code/appendix_B/z2a_struct.py
    python -I -B code/appendix_B/z2a_enum.py DMIN DMAX
  z2a_struct.py (screen output = z2a_struct.txt): the dictionary between the pieces of the two decomposable
  families and tensor products, the Hodge numbers of all pieces against the Jacobian ring, the Hodge numbers of
  the eigenspaces V_chi and the types of the lines used, by exact counts. z2a_enum.py: enumeration of
  decomposable families; z2a_enum.txt = screen output of z2a_enum.py 6 12.

2.10 theory_checks/kd_fermat (item 5; added after a referee report)
    python -I -B code/theory_checks/kd_fermat/kd_fermat_check.py > OUT
  No arguments (6 <= d <= 100 fixed). Counts KR(d) directly from its definition (n = d), and the Hodge characters
  B_d of the Fermat surface of degree d, rho = 1 + #B_d; checks k(d) = rho - 2d^2 + 5d - 4 for d >= 7 (Remark 10.3
  of the note), the listed values of k(d), and that B_d equals the set of decomposable characters exactly when
  gcd(d,6) = 1 (d >= 7). Stored: outputs/theory_checks/kd_fermat/kd_fermat_d6_100.txt (mismatches 0; the d = 6
  line shows k(6) = 16 and 40 for the formula, which is claimed only for d >= 7).


3. Rerun when assembling this repository
----------------------------------------
Every script below was run from this repository's copy (the scripts that write into "../outputs" from a
byte-identical copy of their code folder in a scratch folder, so that nothing was written into the repository),
with Python 3.12.9, python-flint 0.9.0, numpy 2.1.3, on small ranges, and compared with the stored outputs.
"Identical" means identical lines or records; fields with run times were ignored.

  u1check.py formulas 30: 146 pairs; all 147 output lines (146 pairs + c_n line) identical to formulas_d6_100.txt.
  u1check.py conn 6 14 2 / comm 6 10 2 / gen 6 14 2: 42 / 22 / 42 records, all identical to the stored jsonl.
  w2xu.py 30: output byte-identical to xu_d6_30.txt (including "# DONE bad=1", see 2.2).
  w2dim.py 30: 146 pairs, bad=0. w2fam.py 60 60: the 5 m345 lines identical to fam_d400_m345_d180.txt; bad=0.
  planes_gram.py: full run d = 6..20 (stored as planes_gram_d6_20.txt), bad 0, no MISMATCH.
  refU_formulas.py 12 20: all 78 pair lines identical to out_formulas.txt; failures 0.
  refU_ring.py 5 7: identical to the first 3 lines of out_ring.txt.
  refU_comm.py 12 12: 6 records identical to out_comm_12_22.jsonl.
  refU_proofmoves.py 6 10: 6984 pairs, failures 0 (the stored file has only the total for d = 6..40).
  ref_xu.py full 6 30: 146 records identical to out_full_6_110.txt; ref_summ.py on it: 0 failed checks, KR table
    for d <= 30 identical to summ_full_6_150.txt. ref_xu.py fam 60: 64 lines identical to out_fam_1000.txt.
  ref_lemmas.py: full run (25 s), output identical to out_lemmas.txt.
  summarize.py _v4 6 7 8 9 10 (on the stored survey_d<d>_v4.json): screen output identical to summary_v4.txt.
  survey.py 6 only=cyc2_001 tag=_smoke: the cyc2_001 block identical to survey_d6_v4.txt (time field ignored).
  v1tests.py: unit_tests.txt identical to the stored one (41 lines).
  z1final.py (on the stored residual outputs): z1_final_table.txt byte-identical.
  z3tests.py: unit_tests.txt identical to the stored one (22 lines), 0 failures.
  z3run.py 1 cyc3_012:6 cyc12_017:6: JSON identical up to the run-time field, Y lb/ub/closed 157/157/True and
    1009/1017/False as stored, same primes 100057, 100129; the .txt outputs identical.
  z3tally.py (on the stored outputs in main): z3_table.txt byte-identical.
  z2a_struct.py, z2a_enum.py 6 12: screen output identical to z2a_struct.txt, z2a_enum.txt.
  survey_table2.py: 86 rows (routes P 5, A 22, C 57, W 2), no missing or inconsistent run, primes
    100049..101281; the LaTeX rows are identical to the rows of the survey table of the note.
    survey_crosscheck.py and check_primes.py: same screen output as before the move into this repository (51
    families closed by both implementations, all with the same Y; every one of the 86 runs used 2 primes).
  smooth_criterion_enum.py and kd_fermat_check.py (added later): full runs; the stored outputs are their output.
Not rerun: the long runs (u1check.py beyond the ranges above, refU_comm.py beyond d = 12, the full survey runs,
v1run.py, z1tri.py, ord3.py/ord3drive.py, jtest.py, z3run.py on the other families) and the scripts listed in
section 2 as not runnable as shipped.
