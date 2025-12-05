# FlashDeconv Sparse 优化调研报告总结

## 问题诊断

当前代码存在 **3 个内存瓶颈**，导致无法处理大规模数据：

| 位置 | 问题代码 | N=1M 时内存 |
|------|----------|-------------|
| `genes.py:47` | `Y.toarray()` | **240 GB** |
| `deconv.py:241` | `Y[:, idx].toarray()` | 32 GB |
| `deconv.py:144` | `Y / Y.sum(axis=1)` | 32 GB |

## 验证结果

运行 `validation/sparse_memory_benchmark.py` 和 `validation/sparse_optimization_demo.py`:

### 内存对比

| 规模 | N | 当前实现 | 优化后 | 改善 |
|------|---|----------|--------|------|
| Small | 5K | 1.8 GB | 0.4 GB | 5.0x |
| Medium | 50K | 19.3 GB | 3.8 GB | 5.1x |
| Large | 200K | OOM (~32GB) | 9.8 GB | ✓ |

### 正确性验证

- **Log-CPM**: 数值差异 < 1e-15 (完全一致)
- **Sketching**: 数值差异 = 0 (完全一致)
- **HVG 选择**: **100% 一致** (关键洞察: `log1p(0) = 0`)

## 核心修改方案

### 1. HVG 选择 (`genes.py`)

**关键洞察: `log1p(0) = 0`，所以 log 变换后稀疏结构保持！**

```python
# Before (OOM for N=1M)
Y_dense = Y.toarray()

# After (sparse friendly, 100% equivalent)
# Step 1: Row normalize using diagonal matrix
D = diags(10000.0 / np.array(Y.sum(axis=1)).flatten())
Y_norm = D @ Y  # Still sparse!

# Step 2: Log1p on non-zero values only
Y_log = Y_norm.copy()
Y_log.data = np.log1p(Y_log.data)  # zeros stay zeros!

# Step 3: Compute mean/var
gene_means = np.array(Y_log.mean(axis=0)).flatten()
mean_sq = np.array(Y_log.power(2).mean(axis=0)).flatten()
gene_vars = N / (N-1) * (mean_sq - gene_means ** 2)
```

### 2. Log-CPM 预处理 (`deconv.py`)

```python
# Before (produces dense output)
Y_cpm = Y / Y.sum(axis=1, keepdims=True) * 1e4

# After (stays sparse)
from scipy.sparse import diags
D = diags(1e4 / np.array(Y.sum(axis=1)).flatten())
Y_cpm = D @ Y
Y_cpm.data = np.log1p(Y_cpm.data)  # in-place, preserves sparsity
```

### 3. 基因子集 (`deconv.py`)

```python
# Before
Y_subset = Y[:, gene_idx].toarray()

# After
Y_subset = Y[:, gene_idx]  # keep sparse
```

## scipy.sparse 关键操作参考

| 操作 | 输出类型 | 保持稀疏? |
|------|----------|-----------|
| `Y.mean(axis=0)` | matrix | 本身 dense，但很小 |
| `Y.power(2)` | sparse | ✓ |
| `diags(d) @ Y` | sparse | ✓ |
| `Y @ dense_matrix` | dense | (预期行为) |
| `Y / Y.sum(axis=1)` | **dense** | ❌ 危险！ |
| `Y.data = np.log1p(Y.data)` | sparse | ✓ in-place |

## 文件清单

- `docs/sparse_optimization_plan.md` - 完整修改方案
- `validation/sparse_optimization_demo.py` - 各优化点验证脚本
- `validation/sparse_memory_benchmark.py` - 内存基准测试脚本
- `validation/sparse_hvg_correctness_test.py` - **HVG 正确性验证** (100% 一致)

## 下一步

1. 确认修改方案无误
2. 修改 `flashdeconv/utils/genes.py`
3. 修改 `flashdeconv/core/deconv.py`
4. 运行完整测试验证
5. 更新论文中的 scalability 声明
