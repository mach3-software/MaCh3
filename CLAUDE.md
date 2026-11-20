# Working notes for Claude sessions in MaCh3 core (branch `hwallace/test/llm_optimisation`)

## Hard rules
Repeated from `~/software/statistics/MaCh3_DUNE/CLAUDE.md`, which is **not** loaded for a session
working in this repository: there is no CLAUDE.md above this one, so a pointer on its own would never
reach you. Read that file as well, for the directory map, the Mac run recipe and the Rorqual build and
check procedure.
- **Never commit or push.** Henry commits himself.
- **Quote clock times in GMT.** Rorqual reports Eastern time (GMT-5 in winter).
- **Never rebuild or reinstall a build directory a running job uses**, on the Mac or on Rorqual. Make a new build directory instead.
- **Never poll Slurm in a tight loop.** Job scripts touch a marker file on exit; check it every few minutes.
- **Several sessions may work on these trees at once.** Message a session (`ListAgents` / `SendMessage`) before moving or renaming files it may hold. Prefer a git worktree per session.

## Layout
- **develop's own files stay where develop has them.** Other packages (MaCh3_DUNE) build on `SampleHandlerInterface`, `SampleHandlerBase`, `FitterBase`, `ParameterHandlerBase` and the rest, so don't move or rename them. Change them only where the branch needs a hook, with a default that keeps develop's behaviour.
- **Code added for gradient sampling goes in subfolders by topic:**
  - `Samples/Reweight/` (`GPU/` for CUDA and Metal), `Samples/Gradient/`, `Samples/Oscillation/`
  - `Fitters/Algorithms/Hamiltonian/`, `Fitters/Algorithms/Hybrid/`
  - Exception: `Samples/SampleHandlerReweight.cpp` implements `SampleHandlerBase` and sits beside it.
- **New headers in those folders** are installed by `install(DIRECTORY ...)` in the folder's CMakeLists; add sources to `add_library` there.
- **`Doc/GradientSampling.md`** describes the design and has a "Where the code lives" table. Keep it in step when files move.
- **Develop has moved on** since the branch was cut. To ask what the branch changed: fetch, then
  `git diff origin/develop...HEAD` (three dots), or diff against `git merge-base origin/develop HEAD`.
  Not local `develop`, which is ~85 commits stale, and not a two-dot `git diff origin/develop`, which
  also reports develop's newer work as though the branch had reverted it.

## Conventions
- **Builds use `-Werror`** with the project's warning set (`MaCh3Warnings`); the option is `ON` by
  default in `CMakeLists.txt`. New code must be warning-free with clang (Mac) and gcc 12 (Rorqual),
  and with CUDA on. Note: this branch narrowed the exemption in
  `cmake/Modules/MaCh3CompilerOptions.cmake` so `-Werror` now applies to clang, where develop
  disabled it. That change is unreviewed and has no cross-compiler CI behind it, so a clang-only
  `-Werror` failure may be its fault rather than yours.
- **Comments the branch adds start `//HW`.** Never tag a new comment with someone else's initials
  (KS, CW, DB, ETA, ...); those mark their own code. The comments already in the gradient sampling
  files are inconsistent about this (many use a bare `//`) and are being cleaned up, so follow the
  rule rather than the surrounding style.
- **When splitting or moving code,** keep comments as they are. Wording changes go in a separate change from path changes.

## Correctness
- **Restructuring must not change results.** Check with `MaCh3_DUNE/Chains/cleanup_checks/check_tree.sh`: single-thread chains, bit for bit, against the reference.
- **Run-to-run reproducibility** (measured 19 Nov on Rorqual):
  - Fixed-seed chains are bit-for-bit reproducible on **one thread** only.
  - With 16 threads the OpenMP reductions (the prior's included) vary in the last digits.
  - Comments that claim otherwise (`SparsePrior.h`, `CpuReweight.cpp`) are being reviewed; don't "fix" them without measuring which reduction varies.
- **GPU and Metal backends:** check against the CPU with `CheckReweight`. CUDA matches to rounding (double); Metal to float tolerance.
