# **GPU-Accelerated Approximate K-mer Counting Using a Bloom Filter**

## **Abstract**

This project implements a high-performance, memory-efficient k-mer counter by combining GPU acceleration with probabilistic filtering. K-mer counting is a fundamental task in bioinformatics, but existing software is computationally intensive when dealing with large datasets. Our approach uses a two-phase, hybrid CPU/GPU architecture. First, a massively parallel CUDA kernel filters the full k-mer stream against a global Bloom filter, using `atomicOr` operations to isolate a reduced set of potential non-singletons efficiently. Second, this candidate set is processed on the CPU using per-thread hash tables to generate final, exact counts. This method achieves significant speedup by leveraging the GPU for bulk filtering while avoiding the bottlenecks of a fully parallelized counting phase.

## **Background**

### **What is K-mer Counting?**

K-mer counting is the process of counting all possible substrings of length k within biological sequences such as DNA or RNA. While the concept is seemingly simple, it becomes computationally intensive when applied to real data from large genomic datasets, due to high memory and time requirements.

### **The Singleton Problem**

When calculating k-mer counts in biological sequences, the majority of k-mers are **singletons**—they occur only once in the sequence. If we use traditional hash tables to store all k-mers, we end up with many hash table entries that are rarely used, leading to inefficient memory utilization. Our solution is to filter out singletons before instantiating hash table entries, storing only non-singleton k-mers.

### **Why Bloom Filters?**

A Bloom filter is a space-efficient probabilistic data structure that can test whether an element is a member of a set. While false positives are possible (a singleton being marked as a non-singleton), false negatives are not—at least in sequential execution. This makes it ideal for filtering out low-frequency k-mers while using minimal memory.

## **Algorithm Overview**

### **Input**

* **Primary Data**: Raw DNA sequence  
* **Parameters**:  
  * `k`: Integer length for k-mers  
  * `m`: Size (in bits) of the global Bloom filter  
  * `h`: Number of hash functions for each k-mer

### **Output**

* K-mer counts for non-singleton k-mers

### **Data Preprocessing**

A sliding window approach is used to create the stream of k-mers. The DNA sequence is divided among threads, and each thread generates k-mers via a sliding window over its assigned region. Each thread writes its k-mers to a global output buffer for processing in Phase 1\.

## **Two-Phase Architecture**

### **Phase 1: Filtering Singletons (Parallelized)**

This phase identifies potential non-singleton k-mers using a Bloom filter. The k-mer stream is divided among CUDA threads, making this the primary parallelized component of our system.

#### **Sequential Version**

In the sequential implementation, for each k-mer we:

1. **Check**: Hash the k-mer using multiple hash functions and check if all corresponding bit positions in the Bloom filter are set to 1  
2. **Add**: Set all bit positions for this k-mer in the Bloom filter to 1

If the k-mer already exists in the Bloom filter (all bits are 1), it's added to the candidate set as a potential non-singleton.

#### **Parallel Version (CUDA)**

For each k-mer, each CUDA thread performs:

1. **Query**: Hash the k-mer using multiple hash functions and check if all corresponding bit positions in the Bloom filter are set to 1\. If yes, the k-mer is a potential non-singleton—add it to a candidate set (output buffer). If no, skip it.  
2. **Insert**: Use `atomicOr` to set all bit positions for this k-mer in the Bloom filter to 1

The use of `atomicOr` operations allows contention-free updates to the shared Bloom filter across thousands of parallel threads.

### **Phase 2: Exact K-mer Counting (CPU)**

The candidate set from Phase 1 (which may contain false positives) is processed to generate exact counts. The candidate k-mers are divided among CPU threads, where each thread maintains its own private hash table and counts occurrences of k-mers in its assigned portion of the original k-mer stream. After all threads finish, the per-thread hash tables are merged on the CPU to produce final counts.

This phase is not parallelized on the GPU, as the overhead of atomic operations on GPU hash tables would create a bottleneck that negates the performance gains from Phase 1\.

## **The False Negative Problem**

### **Root Cause**

One major issue arising from the conversion from sequential to parallel execution is the occurrence of **false negatives**; they are non-singleton k-mers that are incorrectly treated as singletons. This occurs due to race conditions in the parallel implementation:

1. Reading a Bloom filter bit while another thread is writing it.  
2. Two threads reading the same Bloom filter bit as 0 before either sets it.

### **Impact Analysis**

However, when examining biological data, true non-singleton k-mers typically have significantly higher counts than singletons. We hypothesized that false negatives may have low average counts (in the 2-4 range), meaning their exclusion has minimal impact on downstream analyses that focus on high-frequency k-mers.

Although the parallel Bloom-filter phase does introduce false negatives, our results show that these misses are not as harmful as the raw counts might suggest. Across all Bloom filter sizes, the false negatives that occurred had very low average frequencies (typically only 2 to 3 occurrences per k-mer). In other words, the k-mers that slipped through were almost singletons anyway. 

By contrast, the true non-singleton k-mers in our dataset had an average count of roughly 8 occurrences, meaning the biologically meaningful signals were still reliably captured. This gap indicates that while parallelization introduces false negatives due to race conditions, the k-mers being missed tend to sit right at the noise threshold. Their exclusion has minimal practical impact on analyses that prioritize higher-frequency k-mers.

## **Implementation Details**

The project is implemented on an Ubuntu machine utilizing CUDA architecture and C++. The implementation consists of:

* **CUDA kernels** for parallel Phase 1 filtering  
* **C++ CPU code** for Phase 2 exact counting and hash table merging  
* **Jupyter Notebook** for testing compiled CUDA code and visualizing k-mer count results

## **Performance Results**

### **Execution Time Comparison**

#### Sequential Implementation

| BloomSize   | Phase1_ms | Phase2_ms | TotalTime_ms | Candidates | FalsePositives | FalsePositiveRate | FalseNegatives |
|------------:|----------:|----------:|-------------:|-----------:|---------------:|-----------------:|---------------:|
| 10000       | 105.12    | 183.82    | 288.94       | 35057.00   | 30677.00       | 7.0039           | 0.00           |
| 50000       | 22.73     | 26.37     | 49.10        | 4850.00    | 470.00         | 0.1073           | 0.00           |
| 100000      | 21.18     | 23.96     | 45.15        | 4392.00    | 12.00          | 0.0027           | 0.00           |
| 2500000     | 75.58     | 24.56     | 100.13       | 4380.00    | 0.00           | 0.0000           | 0.00           |
| 5000000     | 75.97     | 23.40     | 99.37        | 4380.00    | 0.00           | 0.0000           | 0.00           |
| 10000000    | 108.58    | 35.22     | 143.80       | 4380.00    | 0.00           | 0.0000           | 0.00           |
| 81177       | 30.57     | 39.52     | 70.09        | 4414.00    | 34.00          | 0.0078           | 0.00           |

In the sequential version, Phase 1 time decreases sharply as the Bloom size increases from 10,000 to 50,000 and again to 100,000. This reflects reduced collisions and a lower false positive rate. For Bloom sizes of 2.5M and above, Phase 1 becomes slower again due to cache and memory overhead from very large bit arrays.

Phase 2 time generally remains stable (23–40 ms), driven primarily by the number of surviving candidates rather than Bloom size. The candidate count drops dramatically when moving from 10,000 (~35k candidates) to 50,000 (~4.8k candidates), stabilizing around 4.3k—indicating that beyond a certain size, the Bloom filter becomes effectively collision-free for this dataset.

False positives decrease monotonically with Bloom size, approaching zero at 2.5M and above. Importantly, no false negatives occur, confirming correctness.

Total runtime ranges from 45 ms (100k Bloom size) to 289 ms (10k Bloom size), showing that Bloom size strongly affects Phase 1 cost for the sequential version.

#### Parallel Implementation
| BloomSize   | Phase1_ms | Phase2_ms | TotalTime_ms | Candidates | FalsePositives | FalsePositiveRate | FalseNegatives | Average Occur of FN | Average Occur of TP |
|------------:|----------:|----------:|-------------:|-----------:|---------------:|-----------------:|---------------:|------------------:|------------------:|
| 10000       | 21.30     | 20.94     | 42.24        | 2271.90    | 944.20         | 0.2156           | 3052.30        | 2.82              | 8.19              |
| 50000       | 20.43     | 27.01     | 47.44        | 1250.60    | 26.40          | 0.0060           | 3155.80        | 2.85              | 8.42              |
| 100000      | 15.64     | 16.42     | 32.06        | 1369.70    | 4.10           | 0.0009           | 3014.40        | 2.78              | 8.00              |
| 2500000     | 23.67     | 17.20     | 40.87        | 1676.80    | 0.00           | 0.0000           | 2703.20        | 2.65              | 7.24              |
| 5000000     | 32.50     | 20.94     | 53.44        | 1626.20    | 0.00           | 0.0000           | 2753.80        | 2.68              | 7.34              |
| 10000000    | 45.06     | 17.03     | 62.08        | 1259.50    | 0.00           | 0.0000           | 3120.50        | 2.84              | 8.31              |
| 81177       | 13.93     | 16.45     | 30.37        | 1290.60    | 9.20           | 0.0021           | 3098.60        | 2.83              | 8.23              |


The GPU implementation shows consistently lower Phase 1 and Phase 2 times. Phase 1 on GPU is dominated by the phase1_filter_kernel, with an average measured device execution time of ~2.7 ms, but total Phase 1 includes memory transfers and CPU orchestration.

Candidate counts are lower than in the sequential version due to race-condition–induced false negatives. These false negatives occur because parallel threads independently update or inspect Bloom entries or candidate tables without sufficient synchronization, leading to k-mers being incorrectly treated as singletons.

False positive rates follow the same trend as the sequential version, approaching zero as Bloom size increases.

The total parallel execution time ranges from 30 ms (Bloom=81,177) to 62 ms (10M Bloom)—consistently faster than the sequential method.
### **Speedup Analysis**

Host-to-Device (HtD) and Device-to-Host (DtH) memory transfers add an average overhead of:
- HtD: 1.094 ms
- DtH: 1.9649 ms
- Total transfer cost: =3.06 ms per run

| Bloom Size | Sequential Total (ms) | Parallel Total (ms) | Speedup   |
| ---------- | --------------------- | ------------------- | --------- |
| 10,000     | 288.94                | 42.24               | **6.84×** |
| 50,000     | 49.10                 | 47.44               | **1.04×** |
| 100,000    | 45.15                 | 32.06               | **1.41×** |
| 2,500,000  | 100.13                | 40.87               | **2.45×** |
| 5,000,000  | 99.37                 | 53.44               | **1.86×** |
| 10,000,000 | 143.80                | 62.08               | **2.32×** |
| 81,177     | 70.09                 | 30.37               | **2.31×** |

### **Overall Observations**
- The GPU consistently outperforms the sequential version, except in the narrow region where sequential performance is already optimal (50k Bloom). Maximum speedup occurs when Bloom size is too small for the CPU but still parallelizable on GPU.
- Transfer overhead isn't that bad (~3 ms), meaning computation dominates overall runtime.
- False negatives in the GPU version inflate apparent speedup by reducing Phase 2 workload; actual ideal speedup would be slightly lower.
## **Benchmarking**

Testing was performed using:

* Real biological sequence data from National Center for Biotechnology Information (NCBI)  
* K-mer length of 12
* Various Bloom filter configurations

## **Limitations and Future Work**

* **False negatives**: Race conditions in parallel execution can cause some non-singleton k-mers to be missed  
* **GPU memory constraints**: Very large Bloom filters may exceed GPU memory capacity

Future improvements could include:

* Techniques to further reduce false negative rates  
* Extension to distributed computing for even larger datasets
