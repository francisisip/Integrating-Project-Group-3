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

\[**Analysis of false negative rates and their k-mer count distributions**\]

## **Implementation Details**

The project is implemented on an Ubuntu machine utilizing CUDA architecture and C++. The implementation consists of:

* **CUDA kernels** for parallel Phase 1 filtering  
* **C++ CPU code** for Phase 2 exact counting and hash table merging  
* **Jupyter Notebook** for testing compiled CUDA code and visualizing k-mer count results

## **Performance Results**

### **Execution Time Comparison**

#### Sequential Implementation

| BloomSize | Phase1\_ms | Phase2\_ms | TotalTime\_ms | Candidates | FalsePositives | FalsePositiveRate | FalseNegatives |
|-----------|-----------|-----------|--------------|------------|----------------|-------------------|----------------|
| 10000     | 142.33    | 303.66    | 445.99       | 43148      | 40949          | 18.6216           | 0              |
| 50000     | 35.96     | 38.69     | 74.66        | 3129       | 930            | 0.4229            | 0              |
| 100000    | 32.68     | 31.53     | 64.21        | 2244       | 45             | 0.0205            | 0              |
| 2500000   | 127.16    | 31.42     | 158.57       | 2199       | 0              | 0.0000            | 0              |
| 5000000   | 99.89     | 21.93     | 121.83       | 2199       | 0              | 0.0000            | 0              |
| 10000000  | 107.13    | 18.88     | 126.01       | 2199       | 0              | 0.0000            | 0              |
| 81168     | 19.96     | 16.61     | 36.58        | 2289       | 90             | 0.0409            | 0              |

#### Parallel Implementation
| BloomSize | Phase1\_ms | Phase2\_ms | TotalTime\_ms | Candidates | FalsePositives | FalsePositiveRate | FalseNegatives | Average Occurences of FN |  
|-----------|-----------|-----------|--------------|------------|----------------|-------------------|----------------|---------------------------|  
| 10000     | 38.68     | 118.46    | 157.13       | 15182      | 13934          | 6.3365            | 951            | 2.43                      |  
| 50000     | 13.27     | 27.02     | 40.28        | 2003       | 1023           | 0.4652            | 1219           | 2.51                      |  
| 100000    | 11.93     | 23.82     | 35.75        | 1229       | 252            | 0.1146            | 1222           | 2.49                      |  
| 2500000   | 15.10     | 23.52     | 38.62        | 1049       | 0              | 0.0000            | 1150           | 2.57                      |  
| 5000000   | 18.94     | 24.01     | 42.95        | 1254       | 0              | 0.0000            | 945            | 2.49                      |  
| 10000000  | 26.80     | 28.53     | 55.33        | 1192       | 0              | 0.0000            | 1007           | 2.57                      |  
| 81168     | 12.36     | 23.74     | 36.09        | 1205       | 326            | 0.1482            | 1320           | 2.60                      |

### **Speedup Analysis**

\[**Chart/Table: Speedup factor achieved by GPU acceleration for various values of k and dataset sizes**\]

Host \-\> Device and Device \-\> Host

## **Benchmarking**

Testing was performed using:

* Real biological sequence data from National Center for Biotechnology Information (NCBI)  
* Synthetic datasets of varying sizes  
* K-mer length of 12
* Various Bloom filter configurations

## **Limitations and Future Work**

* **False negatives**: Race conditions in parallel execution can cause some non-singleton k-mers to be missed  
* **GPU memory constraints**: Very large Bloom filters may exceed GPU memory capacity

Future improvements could include:

* Techniques to further reduce false negative rates  
* Extension to distributed computing for even larger datasets
