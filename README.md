## Running the autotuning and benchmarking scripts

Autotuning: 
```sh
python3 tiled_vectorized_autotuning.py
```  

Benchmarking of all implementations, also generates best unrolling code:
```sh
python3 benchmark_all.py ../autotuning/autotuned_parameters.csv
```

Generate plots with the results:
```sh
python3 plot_results.py benchmark_results.csv
```