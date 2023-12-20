# Benchmarks

## Table of Contents

- [Benchmark Results](#benchmark-results)
    - [Polynomial Multiplication Benchmarks](#polynomial-multiplication-benchmarks)

## Benchmark Results

### Polynomial Multiplication Benchmarks

|             | `Old NTT`                 | `New NTT`                          | `Direct NTT`                        |
|:------------|:--------------------------|:-----------------------------------|:----------------------------------- |
| **`4`**     | `1.91 us` (✅ **1.00x**)   | `973.41 us` (❌ *508.39x slower*)   | `527.70 us` (❌ *275.61x slower*)    |
| **`8`**     | `30.60 us` (✅ **1.00x**)  | `943.25 us` (❌ *30.83x slower*)    | `556.09 us` (❌ *18.17x slower*)     |
| **`16`**    | `61.02 us` (✅ **1.00x**)  | `1.00 ms` (❌ *16.45x slower*)      | `611.82 us` (❌ *10.03x slower*)     |
| **`32`**    | `140.19 us` (✅ **1.00x**) | `1.22 ms` (❌ *8.69x slower*)       | `717.35 us` (❌ *5.12x slower*)      |
| **`64`**    | `321.59 us` (✅ **1.00x**) | `1.65 ms` (❌ *5.14x slower*)       | `828.52 us` (❌ *2.58x slower*)      |
| **`128`**   | `668.54 us` (✅ **1.00x**) | `2.62 ms` (❌ *3.92x slower*)       | `1.18 ms` (❌ *1.77x slower*)        |
| **`256`**   | `1.50 ms` (✅ **1.00x**)   | `4.54 ms` (❌ *3.02x slower*)       | `1.82 ms` (❌ *1.21x slower*)        |
| **`512`**   | `3.22 ms` (✅ **1.00x**)   | `8.23 ms` (❌ *2.56x slower*)       | `2.81 ms` (✅ **1.14x faster**)      |
| **`1024`**  | `7.10 ms` (✅ **1.00x**)   | `16.34 ms` (❌ *2.30x slower*)      | `4.85 ms` (✅ **1.46x faster**)      |
| **`2048`**  | `16.21 ms` (✅ **1.00x**)  | `32.75 ms` (❌ *2.02x slower*)      | `9.00 ms` (🚀 **1.80x faster**)      |
| **`4096`**  | `34.82 ms` (✅ **1.00x**)  | `66.27 ms` (❌ *1.90x slower*)      | `19.25 ms` (🚀 **1.81x faster**)     |
| **`8192`**  | `76.88 ms` (✅ **1.00x**)  | `138.69 ms` (❌ *1.80x slower*)     | `34.25 ms` (🚀 **2.24x faster**)     |
| **`16384`** | `167.22 ms` (✅ **1.00x**) | `293.78 ms` (❌ *1.76x slower*)     | `69.72 ms` (🚀 **2.40x faster**)     |
| **`32768`** | `362.20 ms` (✅ **1.00x**) | `656.02 ms` (❌ *1.81x slower*)     | `153.91 ms` (🚀 **2.35x faster**)    |

---
Made with [criterion-table](https://github.com/nu11ptr/criterion-table)

