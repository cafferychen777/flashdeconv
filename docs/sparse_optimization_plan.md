# FlashDeconv Sparse Matrix Optimization Plan

## Executive Summary

当前 FlashDeconv 代码在处理大规模数据时存在内存瓶颈，主要原因是在多个位置将稀疏矩阵转换为稠密矩阵。本文档详细分析问题并提出修复方案。

---

## 1. 问题诊断

### 1.1 内存瓶颈位置

| 位置 | 文件:行号 | 问题代码 | 严重性 |
|------|----------|----------|--------|
| **P1** | `genes.py:46-47` | `Y_dense = Y.toarray()` | 🔴 致命 |
| **P2** | `deconv.py:240-241` | `Y_subset = Y[:, gene_idx].toarray()` | 🟡 严重 |
| **P3** | `deconv.py:144` | `Y / Y.sum(axis=1, keepdims=True)` | 🟡 严重 |

### 1.2 内存消耗估算

假设 `N = 1,000,000` spots, `G = 30,000` genes, `G' = 4,000` selected genes:

| 问题 | 矩阵维度 | 内存 (float64) |
|------|----------|----------------|
| P1: 全矩阵转 dense | N × G | **240 GB** |
| P2: 子集转 dense | N × G' | **32 GB** |
| P3: numpy 除法输出 | N × G' | **32 GB** |
| 理想: Sketch 输出 | N × 512 | **4 GB** |

### 1.3 当前数据流

```
Y_sparse (N × 30K, sparse)
    │
    ├──▶ select_hvg()
    │       └── Y.toarray()           ← P1: 240GB 内存爆炸！
    │
    ├──▶ gene_idx (2000-4000 genes)
    │
    ▼
Y[:, gene_idx].toarray()              ← P2: 32GB
    │
    ▼
_preprocess_data()
    │   Y / Y.sum(axis=1) * 1e4       ← P3: numpy 广播产生 dense
    │   np.log1p(Y_cpm)
    ▼
Y_tilde (dense, N × 4K)
    │
    ▼
sketch_data()
    │   Y_tilde @ Omega
    ▼
Y_sketch (dense, N × 512)             ← 这一步是合理的
```

---

## 2. 修复方案

### 2.1 目标数据流

```
Y_sparse (N × 30K, sparse CSR)
    │
    ├──▶ select_hvg_sparse()          ← 全程 sparse 操作
    │       └── mean, var 直接从 sparse 计算
    │
    ├──▶ gene_idx
    │
    ▼
Y_sparse[:, gene_idx]                 ← 保持 sparse (N × 4K)
    │
    ▼
_preprocess_data_sparse()
    │   diags(1/lib_size) @ Y         ← sparse @ sparse = sparse
    │   Y.data = np.log1p(Y.data)     ← 原位操作，保持稀疏
    ▼
Y_tilde (sparse, N × 4K)
    │
    ▼
sketch_data()
    │   Y_sparse @ Omega_dense        ← sparse @ dense = dense (高效)
    ▼
Y_sketch (dense, N × 512)             ← 最终输出，可接受
```

---

## 3. 具体修改

### 3.1 修改 P1: `genes.py` - HVG 选择

**当前代码 (lines 45-49):**
```python
# Convert to dense if sparse
if sparse.issparse(Y):
    Y_dense = Y.toarray()  # ❌ 240GB for 1M spots
else:
    Y_dense = np.asarray(Y)
```

**修改方案 (数值完全等价):**

**关键洞察: `log1p(0) = 0`，所以 log 变换后稀疏结构保持不变！**

```python
def select_hvg_sparse(Y, n_top=2000):
    """
    Sparse-friendly HVG selection - NUMERICALLY EQUIVALENT to dense version.

    Key insight: log1p(0) = 0, so sparsity is preserved after log transform!
    """
    N, G = Y.shape

    if sparse.issparse(Y):
        # Step 1: Row normalize (CPM-like) using diagonal matrix
        lib_size = np.array(Y.sum(axis=1)).flatten()
        lib_size = np.maximum(lib_size, 1.0)
        D = diags(10000.0 / lib_size)
        Y_norm = D @ Y  # Still sparse!

        # Step 2: Log1p transform - zeros stay zeros!
        Y_log = Y_norm.copy()
        Y_log.data = np.log1p(Y_log.data)

        # Step 3: Compute mean and variance per gene
        gene_means = np.array(Y_log.mean(axis=0)).flatten()
        mean_sq = np.array(Y_log.power(2).mean(axis=0)).flatten()
        # Sample variance (ddof=1): Var = N/(N-1) * (E[X^2] - E[X]^2)
        gene_vars = N / (N - 1) * (mean_sq - gene_means ** 2)
        gene_vars = np.maximum(gene_vars, 0)

        # Step 4: Bin-based normalization (same as original, on small arrays)
        n_bins = 20
        bins = np.percentile(gene_means[gene_means > 0], np.linspace(0, 100, n_bins + 1))
        bins = np.unique(bins)

        gene_bins = np.digitize(gene_means, bins) - 1
        gene_bins = np.clip(gene_bins, 0, len(bins) - 2)

        normalized_dispersion = np.zeros(G)
        for i in range(len(bins) - 1):
            mask = gene_bins == i
            if np.sum(mask) > 1:
                bin_vars = gene_vars[mask]
                bin_mean = np.mean(bin_vars)
                bin_std = np.std(bin_vars) + 1e-10
                normalized_dispersion[mask] = (bin_vars - bin_mean) / bin_std

        # Step 5: Select top genes
        hvg_idx = np.argsort(normalized_dispersion)[::-1][:n_top]
        return np.sort(hvg_idx)
    else:
        # Dense path (original logic)
        ...
```

**验证结果 (见 `validation/sparse_hvg_correctness_test.py`):**
- HVG 重叠率: **100%** (2000/2000)
- 均值差异: 9.99e-16
- 方差差异: 6.11e-14
- 内存减少: **4.2x**
- 速度提升: **4.0x**

### 3.2 修改 P2: `deconv.py` - 删除 toarray()

**当前代码 (lines 239-244):**
```python
# Subset to selected genes
if sparse.issparse(Y):
    Y_subset = Y[:, gene_idx].toarray()  # ❌ 32GB
else:
    Y_subset = Y[:, gene_idx]
X_subset = X[:, gene_idx]
```

**修改方案:**
```python
# Subset to selected genes (keep sparse!)
if sparse.issparse(Y):
    Y_subset = Y[:, gene_idx]  # ✓ Still CSR, N × 4000
    # Ensure CSR format for efficient row operations
    if not sparse.isspmatrix_csr(Y_subset):
        Y_subset = Y_subset.tocsr()
else:
    Y_subset = Y[:, gene_idx]
X_subset = X[:, gene_idx]
```

### 3.3 修改 P3: `deconv.py` - Sparse Log-CPM

**当前代码 (lines 142-146):**
```python
if method == "log_cpm":
    # CPM normalization + log1p (recommended)
    Y_cpm = Y / (Y.sum(axis=1, keepdims=True) + 1e-10) * 1e4  # ❌ dense output
    X_cpm = X / (X.sum(axis=1, keepdims=True) + 1e-10) * 1e4
    return np.log1p(Y_cpm), np.log1p(X_cpm)
```

**修改方案:**
```python
from scipy.sparse import diags, issparse

def _preprocess_data(self, Y, X, method):
    if method == "log_cpm":
        # ========== Y (可能是 sparse) ==========
        if issparse(Y):
            # 1. Compute library size (row sums)
            lib_size = np.array(Y.sum(axis=1)).flatten()
            lib_size[lib_size == 0] = 1.0

            # 2. Build diagonal scaling matrix D = diag(1e4 / lib_size)
            scale_factors = 1e4 / lib_size
            D = diags(scale_factors)

            # 3. CPM = D @ Y (sparse @ sparse = sparse)
            Y_cpm = D @ Y

            # 4. Log1p: operate on .data directly (preserves sparsity!)
            # Because log1p(0) = 0, zeros stay zeros
            Y_norm = Y_cpm.copy()
            Y_norm.data = np.log1p(Y_norm.data)
        else:
            # Dense path
            lib_size = Y.sum(axis=1, keepdims=True)
            lib_size[lib_size == 0] = 1.0
            Y_cpm = Y / lib_size * 1e4
            Y_norm = np.log1p(Y_cpm)

        # ========== X (通常很小，dense 没问题) ==========
        if issparse(X):
            X = X.toarray()  # K × G' is small (~800KB)
        X_lib = X.sum(axis=1, keepdims=True)
        X_lib[X_lib == 0] = 1.0
        X_cpm = X / X_lib * 1e4
        X_norm = np.log1p(X_cpm)

        return Y_norm, X_norm
```

### 3.4 确认 Sketching 兼容性

**当前代码 (`sketching.py:204`):**
```python
Y_sketch = Y_tilde @ Omega
```

**分析:**
- 如果 `Y_tilde` 是 sparse CSR, `Omega` 是 dense ndarray
- `scipy.sparse` 的 `@` 运算符会自动调用优化的 sparse-dense 乘法
- 输出是 dense ndarray (N × d)，这是我们期望的

**建议:** 确保 `Omega` 是 dense ndarray（而不是 sparse）：
```python
def sketch_data(Y_tilde, X_tilde, Omega, ...):
    # Ensure Omega is dense for optimal sparse @ dense multiplication
    if sparse.issparse(Omega):
        Omega = Omega.toarray()

    # This handles both sparse and dense Y_tilde
    Y_sketch = Y_tilde @ Omega  # ✓ Works for both
    X_sketch = X_tilde @ Omega

    return Y_sketch, X_sketch
```

---

## 4. 实现检查清单

- [ ] **P1: `genes.py:select_hvg()`**
  - [ ] 使用 `Y.mean(axis=0)` 代替 `np.mean(Y_dense, axis=0)`
  - [ ] 使用 `Y.power(2).mean(axis=0)` 计算 E[Y²]
  - [ ] 移除 `Y.toarray()` 调用
  - [ ] 验证: 1M spots 时内存不超过输入矩阵大小

- [ ] **P2: `deconv.py:fit()` line 240**
  - [ ] 移除 `.toarray()` 调用
  - [ ] 确保返回 CSR 格式

- [ ] **P3: `deconv.py:_preprocess_data()`**
  - [ ] 对 sparse Y 使用 `diags()` 进行行归一化
  - [ ] 使用 `Y.data = np.log1p(Y.data)` 原位操作
  - [ ] 保持 X 的 dense 处理（因为 K 很小）

- [ ] **P4: `sketching.py:sketch_data()`**
  - [ ] 确保 Omega 是 dense
  - [ ] 验证 sparse @ dense 工作正常

---

## 5. 验证方案

### 5.1 单元测试

每个修改点都需要验证：
1. **正确性**: 输出与原 dense 实现一致（允许浮点误差 < 1e-6）
2. **稀疏性**: 中间结果保持 sparse
3. **内存**: 峰值内存符合预期

### 5.2 集成测试

```python
# 测试场景
test_cases = [
    {"N": 5000, "G": 20000, "expected_peak_mb": 500},    # Visium
    {"N": 100000, "G": 20000, "expected_peak_mb": 5000}, # Xenium
    {"N": 1000000, "G": 20000, "expected_peak_mb": 12000}, # Atlas (target!)
]
```

### 5.3 Benchmark 脚本

见 `validation/sparse_benchmark.py`

---

## 6. 验证结果

### 6.1 Demo 脚本验证

**Log-CPM 验证** (`validation/sparse_optimization_demo.py --demo logcpm`):
```
Sparse Log-CPM: 178.6 MB, 0.065s
Dense Log-CPM:  640.1 MB, 0.308s
Memory reduction: 3.6x
Max difference: 8.88e-16 (numerically identical)
```

**HVG 选择验证** (`validation/sparse_hvg_correctness_test.py`):
```
N=5000, G=10000, density=0.15, n_top=2000

Numerical Comparison:
  Gene means max diff: 9.99e-16
  Gene vars max diff:  6.11e-14
  Norm dispersion max diff: 9.00e-12

HVG Selection:
  Overlap: 2000/2000 (100.0%)  ← 完全一致！
  Identical HVG sets: True

Memory/Speed (N=20000, G=20000):
  Sparse: 2285 MB, 0.92s
  Dense:  9601 MB, 3.69s
  Reduction: 4.2x memory, 4.0x faster
```

### 6.2 完整 Pipeline Benchmark

运行 `validation/sparse_memory_benchmark.py`:

| Scale | N | Current Peak | Optimized Peak | Reduction |
|-------|---|--------------|----------------|-----------|
| Small | 5,000 | 1,837 MB | 368 MB | **5.0x** |
| Medium | 50,000 | 19,332 MB | 3,779 MB | **5.1x** |
| Large | 200,000 | ~32,000 MB (OOM) | 9,847 MB | **3.3x** |

### 6.3 关键发现

1. **优化版本在 200K spots 仅需 ~10GB 内存**（当前版本需要 32GB，会 OOM）
2. **所有阶段都有显著改善**，HVG 阶段改善最大
3. **数值精度完全保持**：
   - Log-CPM 差异 < 1e-15
   - HVG 选择 **100% 一致**（关键洞察：`log1p(0) = 0`）
4. **速度提升 4x**：sparse 操作比 dense 更快

---

## 7. 风险评估

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| 数值精度差异 | 结果略有不同 | 设置合理的 tolerance |
| scipy 版本兼容 | 旧版本可能缺少某些 API | 文档说明最低版本要求 |
| Pearson residuals 不兼容 | 该方法需要全矩阵操作 | 仅对 log_cpm 做优化 |

---

## 7. 参考资料

- [scipy.sparse documentation](https://docs.scipy.org/doc/scipy/reference/sparse.html)
- [Efficient sparse matrix operations (Stack Overflow)](https://stackoverflow.com/questions/12305021/efficient-way-to-normalize-a-scipy-sparse-matrix)
- [scanpy sparse handling](https://scanpy.readthedocs.io/)
