# **Analysis of AluY Pattern Matching in Human Genome Assemblies**

## **1. Implementation Description**
The goal of this assignment was to identify occurrences of the AluY sequence in two human genome assemblies: **GRCh38 (hg38)** and **T2T CHM13v2.0**. Two pattern-matching approaches were implemented:

- **Exact Matching:** Searches for occurrences of AluY with a strict character-by-character comparison.
- **Fuzzy Matching (1-Mismatch):** Allows for up to **one mismatch**, including substitutions, insertions, or deletions.

The implementations were written in **Python** without using built-in string-matching functions. The scripts (**exact_match.py** and **non_exact_match.py**) were executed on **Ibex HPC**, using SLURM job scheduling. Performance analysis was conducted based on execution time and peak memory usage.

## **2. Results: AluY Counts in Two Genome Assemblies**
The results for both **exact** and **fuzzy** matching approaches are as follows:

| Genome Assembly | Exact Matches | Fuzzy Matches (≤1 mismatch) |
|----------------|--------------|------------------------------|
| GRCh38 (hg38)  | **3** | **13** |
| T2T CHM13v2.0  | **2** | **10** |

The **fuzzy search** identified more matches, highlighting potential variations or sequencing errors in the genome. The output files with location of
each match are found in the files: **exact_Ch38.out**, **exact_T2T.out**, **fuzzy_Ch38.out** and **fuzzy_T2T.out**.

## **3. Performance Analysis**
The execution time and peak memory usage was analyzed for both pattern-matching algorithms, and are presented in the files: (**time_ch38_exact.txt**, **time_T2T_exact.txt**, **time_ch38_fuzzy.txt** and **time_T2T_fuzzy.txt**).

### **3.1 Execution Time**
| Genome Assembly | Exact Matching Time | Fuzzy Matching Time |
|----------------|--------------------|---------------------|
| GRCh38 (hg38)  | **8 min 56 sec**  | **56 min 16 sec**  |
| T2T CHM13v2.0  | **8 min 22 sec**  | **53 min 26 sec**  |

### **3.2 Peak Memory Usage**
| Genome Assembly | Exact Matching (MaxRSS) | Fuzzy Matching (MaxRSS) |
|----------------|------------------------|-------------------------|
| GRCh38 (hg38)  | **3.24 GB**             | **3.24 GB**             |
| T2T CHM13v2.0  | **3.08 GB**             | **3.08 GB**             |

### **3.3 Observations**
1. **Exact matching is significantly faster (~8 minutes) compared to fuzzy matching (~53-56 minutes).** This is expected as the fuzzy search involves additional comparisons and checks for mismatches.
2. **Memory usage is similar for both approaches (~3.0-3.2 GB)**, indicating that the computational cost is mainly in processing time rather than memory consumption.

