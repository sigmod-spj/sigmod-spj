# sigmod-spj

`spj-full.pdf` is the full version of the paper that contains all proofs and SQL's.
==It has been updated for revision.==

- Synthetic: Experiments for synthetic datasets (uniform and zipf)
- Benchmark: Experiments for benchmark datasets (TPCDS Q28 and Q54)
- Real: Experiments for real datasets (IMDb and Network)

To run the experiments:
1. Generate synthetic data by excecuting `data/main.cpp`. The other datasets are available online.
2. Compile the code (example makefile is given under `synthetic` folder). For synthetic dataset, uncomment the appropriate filter before compiling.
3. Run the binary. Output will be written to `./result/`

An execution outputs the following files:
- `time.txt`: Computation time
- `exact.txt`: Exact value of each query
- `ub.txt`: upperbound of each query
- `uds.txt`: uniform distinct sampling estimates
- `wds.txt`: weighted distinct sampling estimates
- `histo.txt`: histogram approximation for wds
- `rw.txt`: random walk algorithm estimates
