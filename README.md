## Many 0-1 proofs

# Files

The repository contains the following files:
- `group_params.py`: contains the parameters of the group used in the
  computations, that is the order of G_p and G_q and a generator of G_p.
- `group.py`: contains a class PowRadix used to precompute tables for faster
  exponentiations. 
- `multiEG.py`: contains an implementation of a multi-ElGamal encryption
  scheme, as well as various 0-1 proof methods.
- `montgomery.py`: contains an implementation of a Montgomery representation of
  G_p.

These files require the module gmpy2 (tested with version 2.0.8).

# Benchmarks

The repository contains the following benchmarks:
- `bench_group.py`
- `bench_multiEG.py`
- `bench_proofs.py`
- `bench_log.py`
- `bench_log_batch_comparison_pre.py`
- `bench_log_batch_comparison_nopre.py`

Each one of these benchmarks can be launched with `python3 filename.py`.

`bench_multiEG.py`, `bench_log_batch_comparison_*.py` require the module
matplotlib to plot the graphs (tested with version 3.5.3).
