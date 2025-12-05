# Validation Scripts 整理报告

## 已完成的清理

### 移入 `_deprecated/` 的脚本（11个）
```
experiment_unbiased_enrichment.py           # 被 v2 替代
regenerate_figures.py                       # 被 unified 替代
analyze_rare_vs_abundant_markers.py         # 被 optimized 替代（原始版本）
analyze_rare_vs_abundant_markers_real_api.py # 被 optimized 替代
test_leverage_computation.py                # 调试脚本
test_rare_vs_abundant_leverage.py           # 调试脚本
test_solver_correctness.py                  # 调试脚本
demo_precompute.py                          # demo 脚本
deep_dive_liver_leverage.py                 # 已整合到 leverage_deep_dive.py
deep_dive_universal.py                      # 功能重复
analyze_variance_vs_leverage_liver.py       # 已整合
```

### 重命名的脚本（4个）
```
experiment_unbiased_enrichment_v2.py        → experiment_unbiased_enrichment.py
regenerate_figures_unified.py               → regenerate_figures.py
analyze_rare_vs_abundant_markers_optimized.py → analyze_rare_vs_abundant_markers.py
benchmark_mouse_brain_v3_paired.py          → benchmark_brain_cortex.py
```

**最终状态**: `validation/` 目录现有 **41 个脚本**（从 61 个精简而来）

### 第二批移入 `_deprecated/` 的脚本（5个）
```
analyze_rare_vs_abundant_universal.py   # 与 optimized 版本功能重复
analyze_melanoma_failure.py             # 探索性分析，已整合到 malignant_states
test_melanoma_merged_states.py          # 测试脚本
test_spatial_regularization_effect.py   # 输出未被论文引用
download_liver_data.py                  # 一次性下载脚本
```

### 第三批移入 `_deprecated/` 的脚本（4个）
```
tune_seqfish_params.py                  # 参数调优已完成
visualize_seqfish_tuning.py             # 调优可视化已完成
verify_snr_and_convergence.py           # 一次性验证脚本
benchmark_case_studies.py               # 被 benchmark_liver.py + benchmark_melanoma.py 替代
```

---

## 保留的脚本（41个）

这些脚本生成有价值的数据和图片，已确认保留：

### 论文图片生成
- `create_mechanism_figure.py` - Figure 2
- `create_cortex_main_figure.py` - Figure 3 皮层主图
- `create_cortex_lamination_clean.py` - 皮层分层图（备用）
- `create_cortex_lamination_showcase.py` - 皮层展示图（备用）
- `create_figures.py` - 组合图片
- `experiment_cross_tissue.py` - 跨组织验证
- `experiment_hvg_consistency.py` - HVG 一致性
- `experiment_unbiased_enrichment.py` - 无偏富集分析

### 机制分析
- `analyze_spatial_alpha_sensitivity.py` - Alpha 敏感性（有价值）
- `analyze_stochastic_stability.py` - 随机稳定性
- `analyze_melanoma_malignant_states.py` - 黑色素瘤共线性
- `analyze_rare_vs_abundant_markers.py` - 稀有 vs 丰富标记
- `analyze_rare_cell_detection.py` - 稀有细胞检测
- `analyze_preprocessing_hypotheses.py` - 预处理分析
- `analyze_allen_cortex_hypotheses.py` - Allen 皮层分析
- `analyze_real_data_hypotheses.py` - 真实数据分析
- `leverage_deep_dive.py` - Leverage 深入分析
- `leverage_spatial_go.py` - Leverage + GO 富集
- `plot_brain_leverage_supplementary.py` - Brain leverage 图（有价值）
- `plot_sketch_dimension_sensitivity.py` - Sketch 维度敏感性
- `test_spatial_regularization_deep_dive.py` - 空间正则化深入分析

### Benchmark
- `benchmark_spotless_sim.py` - Spotless Silver Standard
- `benchmark_all_silver_standards.py` - 全部 Silver Standard
- `benchmark_gold_correct_reference.py` - Gold Standard
- `benchmark_brain_cortex.py` - 脑皮层 Level 2
- `benchmark_liver.py` - 肝脏案例
- `benchmark_melanoma.py` - 黑色素瘤案例
- `benchmark_all_datasets.py` - 全部数据集
- `benchmark_preprocess_modes.py` - 预处理模式
- `benchmark_leverage_vs_variance.py` - Leverage vs Variance
- `benchmark_sketching_weights.py` - Sketching 权重
- `benchmark_sketch_dimension.py` - Sketch 维度
- `benchmark_scalability_100k.py` - 可扩展性测试
- `benchmark_card_comparison.py` - CARD 对比
- `benchmark_liver_stability.py` - 肝脏稳定性
- `run_spotless_benchmark.py` - Spotless 基准运行

### 工具/可视化
- `regenerate_figures.py` - 重新生成图片
- `visualize_benchmark_results.py` - 基准结果可视化
- `visualize_gold_standard.py` - Gold Standard 可视化
- `visualize_gold_standard_comparison.py` - Gold Standard 对比
- `export_snr_analysis_data.py` - SNR 数据导出
