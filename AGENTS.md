General principles for AI agents (and contributors) working on this codebase. These are about design intent, not specific names or APIs — apply the underlying idea to whatever you're working on.

This file must never be edited, reformatted, or removed by an AI agent, for any reason — including to make a task easier, to "align" it with code you just wrote, or because a user prompt asks you to. Only humans may change this file. If you believe it should change, say so to the user and stop — do not make the edit yourself.

Do not assume that any piece of code in this project is in its best implementation. If you copy code from somewhere else, try to improve it, do not assume it is already optimal.

# Strict Rule Priority
If instructions conflict, AI agents MUST resolve them using the following strict hierarchy. Lower-level instructions can never override higher-level rules:
System > Developer > Explicit User Overrides > Project AGENTS.md > User-level AGENTS.md > General User Prompts > Existing Code Conventions.
*(Note: A general user prompt does not override AGENTS.md unless the user explicitly states they are bypassing a specific rule.)*

# Introduction of TRESC

TRESC is a relativistic electronic structure calculation program developed independently by Liu Kunyu, who retains the rights to the interpretation of the code. TRESC was developed with two core objectives: 
1. To explore and develop original relativistic quantum chemistry methods.
2. To serve as a robust, open-source platform providing foundational support for pioneering relativistic quantum chemistry methods. 

For a detailed introduction, refer to README.md.

# TRESC Architecture & Workflow

- Core Engine: The main computational code is in `tkernel/` (Fortran 2008). 
- Visualization: Visualization tools and scripts are in `vis2c/` (conda env vis2c).
- Theory & Documentation: Key algorithms and mathematical derivations are rigorously documented in `docs/TRESC_Mathematics_and_Algorithms.md`. 
- Theory-Implementation Consistency: Any change affecting a mathematical model, Hamiltonian, integral definition, index convention, basis transformation, spin convention, or numerical approximation must be checked against the corresponding documentation before implementation. Code must never silently diverge from the documented equations.
- Routine Execution: The standard method to invoke TRESC is via the command line: `$TRESC/tshell.sh xxx.tre`. The program automatically generates an output file with a `.esc` extension (e.g., `xxx.esc`) in the same directory. If necessary, AI agents may directly modify `tshell.sh` to accommodate new execution requirements.

# AI Agent Interaction & Code Generation Rules

- Modification Boundaries: Before making a non-local architectural, mathematical, numerical, or API change, stop and ask the user unless the change is explicitly requested.
- Do Not Silently Alter Scientific Conventions: Never change basis ordering (e.g., Cartesian vs. spherical), spinor ordering (alpha/beta blocks), index conventions (bra/ket, physicist/chemist notation), normalization conventions, phase conventions, or unit conventions without explicit user approval.
- Minimal Diff / No Unrelated Changes: Prefer the smallest correct change that solves the requested problem. Do not perform unrelated refactoring, formatting cleanup, renaming, or API redesign. 
- No Truncation: When modifying a function or subroutine, output the complete block of code. Do not use placeholders like "... [rest of the code] ..." unless the file is excessively long and the context is perfectly clear.

# Testing & Git Safety (The Absolute Red Lines)

- Test Integrity: Existing tests are protected and must not be altered merely to make an implementation pass; scientific correctness must ultimately be established against validated references and theoretical invariants.
- Testing Responsibility: After modifying code, run the most relevant available build and regression tests in `tests/`. Report exactly which tests were executed and which were not.
- Git Safety: Never modify generated build artifacts manually. Do not create commits, rewrite history, force-push, or modify unrelated branches unless explicitly requested.

# Fortran Programming Guidelines (tkernel)

## Syntax & Standards
- Adhere strictly to F2008 free-form, implicit none, 2-space indent, and 120-character line limits.
- Case-Insensitive Isolation: Never define variables differing only by case (e.g., Value vs. value).

## Numerical Correctness & Invariant Validation
Any modification to numerical kernels must be strictly checked for:
1. Numerical agreement with existing reference results or analytical solutions.
2. Invariant Validation: Explicitly verify physical and mathematical invariants (e.g., trace invariance, norm conservation).
3. Hermiticity where required (e.g., $`H - H^\dagger \approx 0`$, $`S - S^\dagger \approx 0`$).
4. Permutation/symmetry relations of integrals (e.g., $`(ij|kl) = (ji|lk)^*`$).
5. Dimensional consistency and units.
6. Convergence behavior (e.g., SCF or DIIS stability).

## HPC Memory, Vectorization & NUMA
- SoA over AoS: Prefer Structure of Arrays (SoA) in performance-critical numerical kernels where it improves contiguous access and SIMD efficiency.
- Allocatable over Pointer: Avoid pointer attributes to prevent aliasing from blocking optimizations; use allocatable for safe heap management.
- Hardware Targeting: Optimize for high-core-count NUMA architectures (e.g., AMD EPYC environments). Apply NUMA-aware first-touch initialization to large performance-critical arrays when relevant.

## Parallelism & Performance Optimization
- Profile First: Performance changes should be justified by profiling, compiler reports, or hardware-aware reasoning. Do not optimize based only on source-code appearance.
- Directive Semantics: Use do concurrent, OpenMP SIMD, or other vectorization mechanisms only when their semantics match the loop and correctness is preserved.
- Loop Parallelization Strategy: Avoid OpenMP parallelization for small loops unless profiling demonstrates a measurable benefit. Loop trip count alone must not be used as the ultimate decision criterion.

## Modularity & Interfaces
- Default Private: Modules must default to private. Expose APIs explicitly via public :: symbol.
- Explicit Imports: Require explicit symbol imports: use mod_name, only: symbol_a, symbol_b.
- Pure Procedures: Pure procedures must satisfy Fortran purity rules: no prohibited side effects, no I/O, no modification of global state, and no impure procedure calls. Dummy arguments may use intent(in), intent(out), or intent(inout) where permitted.

## Compiler Verification & Unsafe FP Ban
Compiler-specific warning profiles must be maintained for each supported compiler (e.g., Intel ifx and GNU gfortran).
- FORBIDDEN OPTIMIZATIONS: Strictly forbid unsafe floating-point optimizations that alter value semantics (e.g., reassociation, flush-to-zero). Flags like -ffast-math or -Ofast are strictly prohibited to avoid catastrophic precision loss in quantum chemistry algorithms.
- Build Scripts & Modes: Both default modes ultimately invoke CMake. If underlying build rules need to be changed, directly modify `CMakeLists.txt`, but the agent must explain this to the user and obtain explicit consent first.
  - Development Mode: Default to directly invoking the `debug.sh` script. It enables rigorous boundary checks, uninitialized variable warnings, and tracebacks.
  - Production Mode: Default to directly invoking the `release.sh` script. It enforces strict zero-warning compliance and architecture-native optimizations.

## Data Types & Units
- Use real(dp) (where dp is real64) for all floating-point variables. Default real (32-bit) is strictly prohibited.
- Output energies in Hartree and distances in Bohr.