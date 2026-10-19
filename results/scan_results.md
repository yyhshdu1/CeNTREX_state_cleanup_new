# Scan results

本文件夹存放 CeNTREX state cleanup 的全部模拟数据。本页由 `scripts/write_results_index.py` 从数据文件自动生成；新增扫描后重新运行即可更新。

## 文件

| 文件 | 内容 | 结果个数 | 大小 |
|---|---|---|---|
| `SPB_J1_cleanup_scans.h5` | CU1：X J=1，R(1) F1'=5/2 F'=2，150 V/cm，80 mW | 17 | 12.1 MB |
| `SPB_J2_cleanup_scans.h5` | CU2：X J=2，R(2) F1'=7/2 F'=3 + F'=4，80 + 24 mW，17.5 MHz EOM β=1.5 | 22 | 2.4 MB |

## 约定

- 每个结果存为 `runs/<名字>`：`results` 数组、每个坐标轴一个数据集（带单位），以及 JSON 元数据（`fixed_parameters_json`、`model_json`、`solver_json`、`notes`）。schema 为 `centrex-cleanup-scan-v1`。
- 除特别说明外，`results` 是清除后留在该 J 的布居，即 1 − 效率。表中“最好”是该结果里效率最高的点；对轨迹、多列的逐态结果和表格类结果不适用，记为 —。
- “4.8 cm / 2.4 cm / 3 cm”指激光光斑高度（4σ 全宽）。实验中 CU1、CU2 现在都是约 3 cm；4.8 cm（CU1）和 2.4 cm（CU2）是早期基线。
- 横向图只算了 y0 ≥ 0，问题关于 y0 对称。“接受范围内平均”的权重：25 kV 束流模拟（σ_vx = 0.65 m/s；σ_y = 10.95 mm @ CU1，10.35 mm @ CU2），截断在 |vx| ≤ 1 m/s、|y0| ≤ 8.5 mm（探测光高 1.7 cm）。
- 报告图表编号指 `docs/state_cleanup/main.pdf`。

读取示例：

```python
from scripts.scan_hdf import load_scan_hdf
run = load_scan_hdf("results/SPB_J2_cleanup_scans.h5", "doppler_cu2_R1_uniform_3cm_coarse")
run["results"], run["axes"], run["fixed_parameters"]
```

## `SPB_J1_cleanup_scans.h5`

### CU1 轴上扫描（2026-08）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `detuning_board` | 失谐扫描 −20…50 MHz | 4.8 cm，26 能级，均匀 12 态 | 1×71 | laser_power_mW = 80；laser_detuning_MHz -20…50 (71) | 99.9269 % |  | 2026-08-24 |
| `detuning_board_-20to_100MHz` | 失谐扫描 −20…100 MHz，步长 1 MHz | 4.8 cm，26 能级（含 Λ 双重态伙伴），均匀 12 态 | 1×121 | laser_power_mW = 80；laser_detuning_MHz -20…100 (121) | 99.9269 % | 图 5.2a | 2026-08-25 |
| `detuning_withoutP-_-20to_100MHz` | 同上，不含 Λ 双重态伙伴 | 21 能级；元数据 retain_opposite_parity_levels 记录有误 | 1×121 | laser_power_mW = 80；laser_detuning_MHz -20…100 (121) | 99.9275 % | 图 5.2a | 2026-08-26 |
| `detuning_fine` | 失谐 19…22 MHz，步长 0.1 MHz | 4.8 cm，26 能级 | 1×31 | laser_power_mW = 80；laser_detuning_MHz 19…22 (31) | 99.9628 % | 表 5.1，图 5.2b | 2026-08-24 |
| `power_detuning_no_phase_eom` | 功率 10–110 mW × 失谐 −16…30 MHz | 4.8 cm，不加 EOM | 6×46 | laser_power_mW 10…110 (6)；laser_detuning_MHz -16…30 (46) | 99.9546 % | 表 5.2 | 2026-08-23 |
| `power_detuning_eom` | 功率 × 失谐 × EOM 频率（0–30 MHz） | 早期 EOM 测试 | 5×6×3 | laser_power_mW 10…110 (5)；laser_detuning_MHz -15…30 (6)；eom_frequency_MHz 0…30 (3) | 99.9474 % |  | 2026-08-23 |
| `vertical_position_beam_diameter` | y0 −12…12 mm × 光斑高度 10–50 mm | 旧配置（1 MHz 偏振切换），恒定峰值光强扫描，不可直接使用 | 50×5 | vertical_position_mm -12…12 (50)；beam_diameter_4sigma_mm 10…50 (5) | 99.8773 % | 图 5.4 | 2026-08-23 |

### CU1 横向速度/位置（4.8 cm）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `doppler_cu1_spb` | vx −1…1 m/s × y0 0…15 mm，细网格 | +20.2 MHz，SPB 相关 8 态 | 21×11 | vx_m_per_s -1…1 (21)；y0_mm 0…15 (11) | 99.9698 % | 表 5.3，图 5.5a | 2026-10-06 |
| `doppler_cu1_uniform` | 同上 | +20.2 MHz，均匀 12 态 | 21×11 | vx_m_per_s -1…1 (21)；y0_mm 0…15 (11) | 99.9682 % | 图 5.5b | 2026-10-07 |
| `doppler_cu1_spb_eom1p6MHz_beta1p5_det20MHz` | 同上，加 1.6 MHz EOM | β=1.5，+20.0 MHz，SPB 8 态 | 21×11 | vx_m_per_s -1…1 (21)；y0_mm 0…15 (11) | 99.4854 % | 图 5.6d | 2026-10-07 |
| `eom_beta_detuning_vx_cu1_spb_4p8cm` | 1.6 MHz EOM：β × 失谐 × vx（y0=0） | SPB 8 态；NaN = 未计算的组合 | 11×17×11 | eom_beta 0…4 (11)；laser_detuning_MHz 16…23.5 (17)；vx_m_per_s -1…1 (11) | 99.9787 % | 图 5.6a | 2026-10-09 |
| `doppler_cu1_spb_4p8cm_coarse_settings` | 9 个设置的 vx × y0 粗网格 | 设置 (β, 失谐) 见元数据 settings；SPB 8 态 | 9×11×6 | setting_index (9)；vx_m_per_s -1…1 (11)；y0_mm 0…15 (6) | 99.9647 % | 表 5.4 | 2026-10-09 |

### CU1 当前条件（3 cm）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `detuning_vx_y_cu1_spb_3cm_coarse` | 失谐 17–23 MHz × vx × y0（0–9 mm） | 不加 EOM，SPB 8 态 | 13×11×4 | laser_detuning_MHz 17…23 (13)；vx_m_per_s -1…1 (11)；y0_mm 0…9 (4) | 99.999 % | 表 5.5 | 2026-10-09 |
| `doppler_cu1_spb_3cm_det19MHz` | vx × y0 细网格（设定点） | +19.0 MHz，不加 EOM，SPB 8 态；接受范围内平均 98.7 % | 21×11 | vx_m_per_s -1…1 (21)；y0_mm 0…15 (11) | 99.9942 % | 图 5.8 | 2026-10-08 |

### CU1 逐态结果与轨迹

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `per_level_cu1_20p2MHz_4p8cm` | 从每个 J=1 子能级出发的残留（轴上） | +20.2 MHz，4.8 cm | 12 | start_level_index (12) | 99.9943 % | 能级图（slides） | 2026-10-09 |
| `per_state_cu1_spb_20p2MHz_4p8cm` | 8 个 SPB 相关态：残留和去向（J=0/2/3/4） | +20.2 MHz，4.8 cm，轴上 | 8×5 | state_index (8)；quantity_index (5) | — |  | 2026-10-09 |
| `trajectories_cu1_vx0_19p0_20p2MHz_4p8cm` | 各能级布居随时间变化，每个起始子能级一条 | +19.0 / +20.2 MHz，vx=0，4.8 cm；能级标签在元数据 | 2×12×2719×26 | laser_detuning_MHz 19…20.2 (2)；start_level_index (12)；time_us 0…271.739 (2719)；level_index (26) | — | 图 5.7 | 2026-10-09 |


## `SPB_J2_cleanup_scans.h5`

### CU2 250 V/cm（旧电场，仅供参考）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `Fprime3_Fprime4_detuning` | δ0 × δ1 小网格 | 2.4 cm，两束独立失谐 | 3×3 | Fprime3_detuning_MHz 17.5…19.5 (3)；Fprime4_detuning_MHz -7.5…-5.5 (3) | 99.9909 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset` | 两束同时偏移 −20…20 MHz | 中心 (18.5, −6.5) MHz | 41 | shared_detuning_offset_MHz -20…20 (41) | 99.9909 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset_board` | 两束同时偏移 −40…79 MHz | 中心 (−18.5, −43.5) MHz | 120 | shared_detuning_offset_MHz_board -40…79 (120) | 99.9909 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset_boarder` | 同上（重复一次） | 中心 (−18.5, −43.5) MHz | 120 | shared_detuning_offset_MHz_board -40…79 (120) | 99.9909 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset_boarder_extension` | 偏移 43…82 MHz | 中心 (18.5, −6.5) MHz | 40 | shared_detuning_offset_MHz_board 43…82 (40) | 99.646 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset_boarder_84-124MHz` | 偏移 84…123 MHz | 中心 (18.5, −6.5) MHz | 40 | shared_detuning_offset_MHz_board 84…123 (40) | 89.9137 % | 附录 B | 2026-08-24 |
| `correlated_detuning_shared_offset_boarder_combined` | 以上偏移扫描拼接，δ0 −58.5…141.5 MHz | δ1 − δ0 = −25 MHz | 200 | Fprime3_detuning_MHz -58.5…141.5 (200) | 99.9909 % | 附录 B | 2026-08-25 |

### CU2 200 V/cm 轴上（2.4 cm）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `correlated_detuning_shared_offset_200V` | 固定 f_AOM=186.83 MHz 扫激光频率 | 偏移 −10…8 MHz，δ1 − δ0 = −25 MHz | 10 | shared_detuning_offset_MHz_board -10…8 (10) | 99.3267 % | 表 6.1 | 2026-09-02 |
| `laser_AOM_2D_200V_coarse_2MHz` | δ0 5–23 MHz × f_AOM 185–215 MHz，步长 2 MHz | 均匀 20 态 | 10×16 | Fprime3_detuning_MHz 5…23 (10)；AOM_frequency_MHz 185…215 (16) | 99.7914 % | 图 6.3 | 2026-10-05 |
| `laser_AOM_2D_200V_fine_0p5MHz_region1` | R1 附近 0.5 MHz 细网格 | f_AOM 207–211 MHz | 9×9 | Fprime3_detuning_MHz 9…13 (9)；AOM_frequency_MHz 207…211 (9) | 99.8034 % | 图 6.4 | 2026-10-06 |
| `laser_AOM_2D_200V_fine_0p5MHz_region2` | R2 附近 0.5 MHz 细网格 | f_AOM 189–193 MHz | 9×9 | Fprime3_detuning_MHz 9…13 (9)；AOM_frequency_MHz 189…193 (9) | 99.8535 % | 图 6.4 | 2026-10-06 |
| `laser_AOM_2D_200V_fine_0p5MHz_region3` | 第三个极大附近 0.5 MHz 细网格 | f_AOM 195–199 MHz | 9×9 | Fprime3_detuning_MHz 21…25 (9)；AOM_frequency_MHz 195…199 (9) | 99.5559 % | 图 6.4 | 2026-10-06 |
| `laser_AOM_2D_200V_confirm_rtol1e-7` | 3 个最佳点用 rtol 1e-7 复算 | f_AOM 见元数据 | 3 | Fprime3_detuning_MHz 11…12 (3) | 99.8508 % | 表 6.3 | 2026-10-06 |
| `trajectories_cu2_R1_R2_2p4cm` | R1、R2 下各能级布居随时间变化 | 均匀 20 态；能级标签、光子数在元数据 | 2×2719×56 | f_AOM_MHz 191…208.5 (2)；time_us 0…271.739 (2719)；level_index (56) | — | 图 6.5，6.6 | 2026-10-09 |
| `per_level_cu2_R1_R2_2p4cm` | 从每个 J=2 子能级出发的残留和泄漏到 J=1 | R1、R2 | 2×20×2 | f_AOM_MHz 191…208.5 (2)；start_level_index (20)；quantity_index (2) | — | 能级图（slides） | 2026-10-09 |
| `per_state_cu2_spb_R1_R2_2p4cm` | 4 个 SPB 源态：残留和去向（J=1/3/4/5） | R1、R2 | 2×4×5 | f_AOM_MHz 191…208.5 (2)；state_index (4)；quantity_index (5) | — |  | 2026-10-09 |
| `level_scheme_cu2_200V_sublevels` | 激发态子能级表：位置、耦合强度、展宽 | 每束光各自的失谐轴 | 32×6 | row_index (32)；column_index (6) | — | 图 6.7 | 2026-10-09 |
| `level_scheme_cu2_200V_pairs` | 耦合对（激发态, 基态）表：共振位置、\|M\|² | X、Z 偏振合计 | 228×6 | pair_index (228)；column_index (6) | — | 图 6.7 | 2026-10-09 |

### CU2 当前条件（3 cm）

| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |
|---|---|---|---|---|---|---|---|
| `doppler_cu2_R1_uniform_3cm_coarse` | vx × y0 粗网格 | R1（f_AOM 208.5 MHz），均匀 20 态；接受范围内平均 96.7 % | 11×6 | vx_m_per_s -1…1 (11)；y0_mm 0…15 (6) | 99.5198 % | 表 6.7，图 6.8a | 2026-10-07 |
| `doppler_cu2_R1_spb_3cm_coarse` | vx × y0 粗网格 | R1，SPB 源态（lens 权重）；接受范围内平均 99.6 % | 11×6 | vx_m_per_s -1…1 (11)；y0_mm 0…15 (6) | 99.9736 % | 表 6.7，图 6.8b | 2026-10-07 |
| `doppler_cu2_R2_uniform_3cm_coarse` | vx × y0 粗网格 | R2（f_AOM 191.0 MHz），均匀 20 态；接受范围内平均 96.1 % | 11×6 | vx_m_per_s -1…1 (11)；y0_mm 0…15 (6) | 99.5856 % | 表 6.7 | 2026-10-07 |
| `doppler_cu2_R2_spb_3cm_coarse` | vx × y0 粗网格 | R2，SPB 源态；接受范围内平均 99.35 % | 11×6 | vx_m_per_s -1…1 (11)；y0_mm 0…15 (6) | 99.9747 % | 表 6.7 | 2026-10-07 |

## 不在本文件夹的中间文件

- `tmp/` 下的 `*.npz` 是扫描的断点文件，内容都已写进上面的结果；作图和续算脚本仍直接读取它们，请保留。
- `tmp/h5_backup_20261009/`：2026-10-09 写入新结果之前两个数据文件的备份。
