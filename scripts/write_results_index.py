"""Write results/scan_results.md: an index of every stored run in the two scan files.

    python scripts/write_results_index.py

Shape, axes, best value and creation date are read from the HDF5 files; the short description, the conditions
and the report reference of each run are kept in DESCRIPTIONS below. A run without a description is still
listed (with its stored notes), so the index never silently misses a run. Re-run after adding scans.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FILES = {"SPB_J1_cleanup_scans.h5": "CU1：X J=1，R(1) F1'=5/2 F'=2，150 V/cm，80 mW",
         "SPB_J2_cleanup_scans.h5": "CU2：X J=2，R(2) F1'=7/2 F'=3 + F'=4，80 + 24 mW，17.5 MHz EOM β=1.5"}

# name -> (group, description, conditions, report reference)
DESCRIPTIONS = {
    # ---------------- CU1 ----------------
    "detuning_board": ("CU1 轴上扫描（2026-08）", "失谐扫描 −20…50 MHz", "4.8 cm，26 能级，均匀 12 态", ""),
    "detuning_board_-20to_100MHz": ("CU1 轴上扫描（2026-08）", "失谐扫描 −20…100 MHz，步长 1 MHz",
                                    "4.8 cm，26 能级（含 Λ 双重态伙伴），均匀 12 态", "图 5.2a"),
    "detuning_withoutP-_-20to_100MHz": ("CU1 轴上扫描（2026-08）", "同上，不含 Λ 双重态伙伴",
                                        "21 能级；元数据 retain_opposite_parity_levels 记录有误", "图 5.2a"),
    "detuning_fine": ("CU1 轴上扫描（2026-08）", "失谐 19…22 MHz，步长 0.1 MHz", "4.8 cm，26 能级", "表 5.1，图 5.2b"),
    "power_detuning_no_phase_eom": ("CU1 轴上扫描（2026-08）", "功率 10–110 mW × 失谐 −16…30 MHz",
                                    "4.8 cm，不加 EOM", "表 5.2"),
    "power_detuning_eom": ("CU1 轴上扫描（2026-08）", "功率 × 失谐 × EOM 频率（0–30 MHz）", "早期 EOM 测试", ""),
    "vertical_position_beam_diameter": ("CU1 轴上扫描（2026-08）", "y0 −12…12 mm × 光斑高度 10–50 mm",
                                        "旧配置（1 MHz 偏振切换），恒定峰值光强扫描，不可直接使用", "图 5.4"),
    "doppler_cu1_spb": ("CU1 横向速度/位置（4.8 cm）", "vx −1…1 m/s × y0 0…15 mm，细网格",
                        "+20.2 MHz，SPB 相关 8 态", "表 5.3，图 5.5a"),
    "doppler_cu1_uniform": ("CU1 横向速度/位置（4.8 cm）", "同上", "+20.2 MHz，均匀 12 态", "图 5.5b"),
    "doppler_cu1_spb_eom1p6MHz_beta1p5_det20MHz": ("CU1 横向速度/位置（4.8 cm）", "同上，加 1.6 MHz EOM",
                                                   "β=1.5，+20.0 MHz，SPB 8 态", "图 5.6d"),
    "eom_beta_detuning_vx_cu1_spb_4p8cm": ("CU1 横向速度/位置（4.8 cm）", "1.6 MHz EOM：β × 失谐 × vx（y0=0）",
                                           "SPB 8 态；NaN = 未计算的组合", "图 5.6a"),
    "doppler_cu1_spb_4p8cm_coarse_settings": ("CU1 横向速度/位置（4.8 cm）", "9 个设置的 vx × y0 粗网格",
                                              "设置 (β, 失谐) 见元数据 settings；SPB 8 态", "表 5.4"),
    "detuning_vx_y_cu1_spb_3cm_coarse": ("CU1 当前条件（3 cm）", "失谐 17–23 MHz × vx × y0（0–9 mm）",
                                         "不加 EOM，SPB 8 态", "表 5.5"),
    "doppler_cu1_spb_3cm_det19MHz": ("CU1 当前条件（3 cm）", "vx × y0 细网格（设定点）",
                                     "+19.0 MHz，不加 EOM，SPB 8 态；接受范围内平均 98.7 %", "图 5.8"),
    "per_level_cu1_20p2MHz_4p8cm": ("CU1 逐态结果与轨迹", "从每个 J=1 子能级出发的残留（轴上）",
                                    "+20.2 MHz，4.8 cm", "能级图（slides）"),
    "per_state_cu1_spb_20p2MHz_4p8cm": ("CU1 逐态结果与轨迹", "8 个 SPB 相关态：残留和去向（J=0/2/3/4）",
                                        "+20.2 MHz，4.8 cm，轴上", ""),
    "trajectories_cu1_vx0_19p0_20p2MHz_4p8cm": ("CU1 逐态结果与轨迹", "各能级布居随时间变化，每个起始子能级一条",
                                                "+19.0 / +20.2 MHz，vx=0，4.8 cm；能级标签在元数据", "图 5.7"),
    # ---------------- CU2 ----------------
    "Fprime3_Fprime4_detuning": ("CU2 250 V/cm（旧电场，仅供参考）", "δ0 × δ1 小网格",
                                 "2.4 cm，两束独立失谐", "附录 B"),
    "correlated_detuning_shared_offset": ("CU2 250 V/cm（旧电场，仅供参考）", "两束同时偏移 −20…20 MHz",
                                          "中心 (18.5, −6.5) MHz", "附录 B"),
    "correlated_detuning_shared_offset_board": ("CU2 250 V/cm（旧电场，仅供参考）", "两束同时偏移 −40…79 MHz",
                                                "中心 (−18.5, −43.5) MHz", "附录 B"),
    "correlated_detuning_shared_offset_boarder": ("CU2 250 V/cm（旧电场，仅供参考）", "同上（重复一次）",
                                                  "中心 (−18.5, −43.5) MHz", "附录 B"),
    "correlated_detuning_shared_offset_boarder_extension": ("CU2 250 V/cm（旧电场，仅供参考）", "偏移 43…82 MHz",
                                                            "中心 (18.5, −6.5) MHz", "附录 B"),
    "correlated_detuning_shared_offset_boarder_84-124MHz": ("CU2 250 V/cm（旧电场，仅供参考）", "偏移 84…123 MHz",
                                                            "中心 (18.5, −6.5) MHz", "附录 B"),
    "correlated_detuning_shared_offset_boarder_combined": ("CU2 250 V/cm（旧电场，仅供参考）",
                                                           "以上偏移扫描拼接，δ0 −58.5…141.5 MHz", "δ1 − δ0 = −25 MHz",
                                                           "附录 B"),
    "correlated_detuning_shared_offset_200V": ("CU2 200 V/cm 轴上（2.4 cm）", "固定 f_AOM=186.83 MHz 扫激光频率",
                                               "偏移 −10…8 MHz，δ1 − δ0 = −25 MHz", "表 6.1"),
    "laser_AOM_2D_200V_coarse_2MHz": ("CU2 200 V/cm 轴上（2.4 cm）", "δ0 5–23 MHz × f_AOM 185–215 MHz，步长 2 MHz",
                                      "均匀 20 态", "图 6.3"),
    "laser_AOM_2D_200V_fine_0p5MHz_region1": ("CU2 200 V/cm 轴上（2.4 cm）", "R1 附近 0.5 MHz 细网格",
                                              "f_AOM 207–211 MHz", "图 6.4"),
    "laser_AOM_2D_200V_fine_0p5MHz_region2": ("CU2 200 V/cm 轴上（2.4 cm）", "R2 附近 0.5 MHz 细网格",
                                              "f_AOM 189–193 MHz", "图 6.4"),
    "laser_AOM_2D_200V_fine_0p5MHz_region3": ("CU2 200 V/cm 轴上（2.4 cm）", "第三个极大附近 0.5 MHz 细网格",
                                              "f_AOM 195–199 MHz", "图 6.4"),
    "laser_AOM_2D_200V_confirm_rtol1e-7": ("CU2 200 V/cm 轴上（2.4 cm）", "3 个最佳点用 rtol 1e-7 复算",
                                           "f_AOM 见元数据", "表 6.3"),
    "trajectories_cu2_R1_R2_2p4cm": ("CU2 200 V/cm 轴上（2.4 cm）", "R1、R2 下各能级布居随时间变化",
                                     "均匀 20 态；能级标签、光子数在元数据", "图 6.5，6.6"),
    "per_level_cu2_R1_R2_2p4cm": ("CU2 200 V/cm 轴上（2.4 cm）", "从每个 J=2 子能级出发的残留和泄漏到 J=1",
                                  "R1、R2", "能级图（slides）"),
    "per_state_cu2_spb_R1_R2_2p4cm": ("CU2 200 V/cm 轴上（2.4 cm）", "4 个 SPB 源态：残留和去向（J=1/3/4/5）",
                                      "R1、R2", ""),
    "level_scheme_cu2_200V_sublevels": ("CU2 200 V/cm 轴上（2.4 cm）", "激发态子能级表：位置、耦合强度、展宽",
                                        "每束光各自的失谐轴", "图 6.7"),
    "level_scheme_cu2_200V_pairs": ("CU2 200 V/cm 轴上（2.4 cm）", "耦合对（激发态, 基态）表：共振位置、|M|²",
                                    "X、Z 偏振合计", "图 6.7"),
    "doppler_cu2_R1_uniform_3cm_coarse": ("CU2 当前条件（3 cm）", "vx × y0 粗网格",
                                          "R1（f_AOM 208.5 MHz），均匀 20 态；接受范围内平均 96.7 %", "表 6.7，图 6.8a"),
    "doppler_cu2_R1_spb_3cm_coarse": ("CU2 当前条件（3 cm）", "vx × y0 粗网格",
                                      "R1，SPB 源态（lens 权重）；接受范围内平均 99.6 %", "表 6.7，图 6.8b"),
    "doppler_cu2_R2_uniform_3cm_coarse": ("CU2 当前条件（3 cm）", "vx × y0 粗网格",
                                          "R2（f_AOM 191.0 MHz），均匀 20 态；接受范围内平均 96.1 %", "表 6.7"),
    "doppler_cu2_R2_spb_3cm_coarse": ("CU2 当前条件（3 cm）", "vx × y0 粗网格",
                                      "R2，SPB 源态；接受范围内平均 99.35 %", "表 6.7"),
}


ORDER = list(DESCRIPTIONS)                      # runs are listed in the order of DESCRIPTIONS


def cell(text):
    return str(text).replace("|", "\\|")


def axis_text(name, values, unit):
    v = np.asarray(values)
    if v.size == 1:
        return f"{name} = {v.ravel()[0]:g}"
    if unit == "index":
        return f"{name} ({v.size})"
    return f"{name} {v.min():g}…{v.max():g} ({v.size})"


def best_text(g):
    if "results" not in g:
        return "—"
    quantity = g.attrs.get("quantity", "")
    if not str(quantity).startswith("final_population_remaining"):
        return "—"
    r = np.asarray(g["results"][...]).real
    if not np.isfinite(r).any():
        return "—"
    return f"{100 * (1 - np.nanmin(r)):.4f}".rstrip("0").rstrip(".") + " %"


def main():
    lines = [
        "# Scan results",
        "",
        "本文件夹存放 CeNTREX state cleanup 的全部模拟数据。本页由 `scripts/write_results_index.py` 从数据文件自动生成；"
        "新增扫描后重新运行即可更新。",
        "",
        "## 文件",
        "",
        "| 文件 | 内容 | 结果个数 | 大小 |",
        "|---|---|---|---|",
    ]
    tables = []
    for fn, what in FILES.items():
        path = RESULTS / fn
        with h5py.File(path, "r") as h:
            runs = sorted(h["runs"], key=lambda n: (ORDER.index(n) if n in ORDER else len(ORDER), n))
            lines.append(f"| `{fn}` | {what} | {len(runs)} | {path.stat().st_size / 1e6:.1f} MB |")
            groups = {}
            for name in runs:
                g = h["runs"][name]
                order = json.loads(g.attrs["axis_order_json"])
                axes = "；".join(axis_text(a, g["axes"][a][...], g["axes"][a].attrs.get("unit", "")) for a in order)
                shape = "×".join(str(s) for s in g["results"].shape) if "results" in g else "—"
                created = str(g.attrs.get("created_utc", ""))[:10]
                grp, desc, cond, ref = DESCRIPTIONS.get(
                    name, ("未分类", str(g.attrs.get("notes", ""))[:80], "", ""))
                groups.setdefault(grp, []).append(
                    f"| `{name}` | {cell(desc)} | {cell(cond)} | {shape} | {cell(axes)} | {best_text(g)} | {ref} "
                    f"| {created} |")
            tables.append((fn, groups))
    lines += [
        "",
        "## 约定",
        "",
        "- 每个结果存为 `runs/<名字>`：`results` 数组、每个坐标轴一个数据集（带单位），以及 JSON 元数据"
        "（`fixed_parameters_json`、`model_json`、`solver_json`、`notes`）。schema 为 `centrex-cleanup-scan-v1`。",
        "- 除特别说明外，`results` 是清除后留在该 J 的布居，即 1 − 效率。表中“最好”是该结果里效率最高的点；"
        "对轨迹、多列的逐态结果和表格类结果不适用，记为 —。",
        "- “4.8 cm / 2.4 cm / 3 cm”指激光光斑高度（4σ 全宽）。实验中 CU1、CU2 现在都是约 3 cm；"
        "4.8 cm（CU1）和 2.4 cm（CU2）是早期基线。",
        "- 横向图只算了 y0 ≥ 0，问题关于 y0 对称。“接受范围内平均”的权重：25 kV 束流模拟（σ_vx = 0.65 m/s；"
        "σ_y = 10.95 mm @ CU1，10.35 mm @ CU2），截断在 |vx| ≤ 1 m/s、|y0| ≤ 8.5 mm（探测光高 1.7 cm）。",
        "- 报告图表编号指 `docs/state_cleanup/main.pdf`。",
        "",
        "读取示例：",
        "",
        "```python",
        "from scripts.scan_hdf import load_scan_hdf",
        'run = load_scan_hdf("results/SPB_J2_cleanup_scans.h5", "doppler_cu2_R1_uniform_3cm_coarse")',
        'run["results"], run["axes"], run["fixed_parameters"]',
        "```",
    ]
    for fn, groups in tables:
        lines += ["", f"## `{fn}`", ""]
        for grp, rows in groups.items():
            lines += [f"### {grp}", "",
                      "| 名字 | 内容 | 条件 | 形状 | 坐标轴 | 最好 | 报告 | 创建 |",
                      "|---|---|---|---|---|---|---|---|", *rows, ""]
    lines += [
        "## 不在本文件夹的中间文件",
        "",
        "- `tmp/` 下的 `*.npz` 是扫描的断点文件，内容都已写进上面的结果；作图和续算脚本仍直接读取它们，请保留。",
        "- `tmp/h5_backup_20261009/`：2026-10-09 写入新结果之前两个数据文件的备份。",
        "",
    ]
    out = RESULTS / "scan_results.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    missing = [n for fn, groups in tables for n in groups.get("未分类", [])]
    print(f"wrote {out}" + (f"; {len(missing)} runs without description" if missing else ""))


if __name__ == "__main__":
    sys.exit(main())
