# Benchmarks

## Table of Contents

- [Benchmark Results](#benchmark-results)
  - [Polynomial Multiplication Benchmarks](#polynomial-multiplication-benchmarks)

## Benchmark Results

### Polynomial Multiplication Benchmarks

|             | `Old NTT`                  | `New NTT`                         | `Direct NTT`                      |
| :---------- | :------------------------- | :-------------------------------- | :-------------------------------- |
| **`4`**     | `1.84 us` (✅ **1.00x**)   | `933.95 us` (❌ _507.09x slower_) | `507.46 us` (❌ _275.53x slower_) |
| **`8`**     | `31.83 us` (✅ **1.00x**)  | `1.06 ms` (❌ _33.21x slower_)    | `555.49 us` (❌ _17.45x slower_)  |
| **`16`**    | `60.56 us` (✅ **1.00x**)  | `1.32 ms` (❌ _21.81x slower_)    | `611.75 us` (❌ _10.10x slower_)  |
| **`32`**    | `142.43 us` (✅ **1.00x**) | `1.95 ms` (❌ _13.72x slower_)    | `688.96 us` (❌ _4.84x slower_)   |
| **`64`**    | `304.66 us` (✅ **1.00x**) | `3.15 ms` (❌ _10.34x slower_)    | `801.15 us` (❌ _2.63x slower_)   |
| **`128`**   | `661.51 us` (✅ **1.00x**) | `5.58 ms` (❌ _8.44x slower_)     | `1.12 ms` (❌ _1.69x slower_)     |
| **`256`**   | `1.45 ms` (✅ **1.00x**)   | `10.54 ms` (❌ _7.28x slower_)    | `1.75 ms` (❌ _1.21x slower_)     |
| **`512`**   | `3.20 ms` (✅ **1.00x**)   | `20.24 ms` (❌ _6.31x slower_)    | `2.69 ms` (✅ **1.19x faster**)   |
| **`1024`**  | `6.97 ms` (✅ **1.00x**)   | `41.83 ms` (❌ _6.00x slower_)    | `4.44 ms` (✅ **1.57x faster**)   |
| **`2048`**  | `15.77 ms` (✅ **1.00x**)  | `81.65 ms` (❌ _5.18x slower_)    | `10.19 ms` (✅ **1.55x faster**)  |
| **`4096`**  | `34.78 ms` (✅ **1.00x**)  | `165.38 ms` (❌ _4.76x slower_)   | `17.26 ms` (🚀 **2.01x faster**)  |
| **`8192`**  | `75.87 ms` (✅ **1.00x**)  | `337.18 ms` (❌ _4.44x slower_)   | `33.25 ms` (🚀 **2.28x faster**)  |
| **`16384`** | `164.70 ms` (✅ **1.00x**) | `778.55 ms` (❌ _4.73x slower_)   | `68.76 ms` (🚀 **2.40x faster**)  |
| **`32768`** | `349.04 ms` (✅ **1.00x**) | `1.42 s` (❌ _4.08x slower_)      | `145.98 ms` (🚀 **2.39x faster**) |

---

Made with [criterion-table](https://github.com/nu11ptr/criterion-table)
