# P1CRC Visium HD 错误分割修正复现

本目录提供 Proseg、ResolVI 和 SPLIT 在同一个 P1CRC Visium HD 样本上的运行脚本。原始输入来自 10x Genomics 的 `Visium_HD_Human_Colon_Cancer_P1` 2 µm binned outputs、H&E BTF 和空间配准文件；单细胞参考为已有的 `adata_sc_all_reanno.h5ad` 中 `Patient == P1CRC` 部分。

## 方法输入与比较口径

| 方法 | 输入 | 输出 | 备注 |
| --- | --- | --- | --- |
| Proseg | 2 µm bins、H&E 上 StarDist 的细胞核先验、Space Ranger 坐标 | 修正后的细胞分配、计数矩阵、细胞元数据 | 在整个组织区域运行 |
| ResolVI | 同一 StarDist 核先验聚合得到的原始细胞×基因 UMI、细胞坐标 | 潜变量、解码表达矩阵 | 无监督模式；不读取 Proseg 输出 |
| SPLIT | 同一原始细胞 UMI、坐标、P1CRC scRNA 参考 | RCTD 双细胞分解、SPLIT 纯化计数 | 正式运行取中心连续区域约 1 万细胞，随机种子 42，保存筛选名单和 RCTD 拒绝细胞统计 |

这三个方法解决的问题并不完全相同：Proseg 改变细胞归属和边界；ResolVI 修正环境 RNA 与细胞内表达；SPLIT 依据参考和 RCTD 权重去除错误归属的转录本。因此结果应分别报告细胞数、有效基因、运行时间和下游指标，不把三者的输出当成同一种数值解释。

压缩包中的 `P1CRC_HD.h5ad` 包含 **507,684 个 8 µm bin**，而非 507,684 个分割细胞；其 `obsm['spatial']` 是 BTF 全分辨率像素坐标。脚本可用 `--spatial-scale 0.27380817798463214` 把它换算为 µm，并用 `--spatial-sampling random` 对整张切片固定随机抽样。该输入可作为旧版 benchmark 的对照口径，但必须与本流程由 H&E 核先验得到的约 22.6 万个细胞分开命名和报告。ResolVI 需要先过滤低于 5 UMI 的 bin；RCTD 默认 `UMI_min=100`，直接运行旧 8 µm 输入会过滤较多 bin。

输入口径预检：H&E 细胞输入的 300 个空间观测中，SPLIT 输出 235 个；旧 8 µm 文件的全组织随机 300 个 bin 中，RCTD 因低 UMI 排除 120 个，SPLIT 最终输出 159 个。这只是流程和保留率检查，不是正式性能比较。聊天记录中约 1 万输入、7 千多最终用于指标的比例更接近前者；旧版具体预处理脚本未找到，不能据此断定历史使用了哪一种输入。

正式 SPLIT 输入的独立预检从中心连续区域选出 10,000 个 StarDist 细胞、5,443 个 P1CRC 参考细胞和 18,071 个共同基因；空间细胞的 UMI 中位数为 210，其中 8,876 个达到 RCTD 默认的 100 UMI 下限，1,124 个低于下限。实际 RCTD 和 SPLIT 还可能进一步拒绝细胞，最终保留数以正式输出为准。

## 当前 qz 路径

```bash
ROOT=/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark/misseg
PY="$ROOT/../envs/revise/bin/python"
RESOLVI_PY="$ROOT/envs/resolvi/bin/python"
R="$ROOT/envs/split/bin/Rscript"
```

原始数据位于 `$ROOT/raw_p1crc/`，参考在 `$ROOT/adata_sc_all_reanno.h5ad`。结果统一放在 `$ROOT/results/`。运行前确认 BTF 已下载完整、`binned_outputs/` 已解压。完整 H&E 的 StarDist 掩膜由 `segment_p1crc_he.py` 生成，裁剪覆盖组织坐标 `(x=0..26000, y=11000..38000)`，每 2 个 BTF 像素预测 1 个掩膜像素。`*.npy.json` 保存从掩膜像素到 µm 的仿射变换，不能丢弃。

全图分割使用 StarDist `prob_thresh=0.01`、`nms_thresh=0.001`，这是为了保留密集组织中的核候选，应在结果中作为预处理参数报告。同一 1024×1024 小区域里，该设定检出 3,284 个核；模型默认阈值检出 970 个、明显漏掉可见核；`prob_thresh=0.2,nms_thresh=0.3` 检出 2,790 个。宽松阈值耗时较长，性能指标可能受多检出的核影响。

## 运行命令

```bash
# 1. H&E 细胞核先验。大图需要较长时间；运行结束后检查 JSON 中 num_nuclei。
"$ROOT/envs/stardist/bin/python" "$ROOT/segment_p1crc_he.py" \
  --image "$ROOT/raw_p1crc/Visium_HD_Human_Colon_Cancer_P1_tissue_image.btf" \
  --output "$ROOT/results/proseg/p1crc_stardist_masks.npy.gz" \
  --microns-per-pixel 0.27380817798463214 \
  --x0 0 --y0 11000 --x1 26000 --y1 38000 --downsample 2 --block-size 2048

# 2. Proseg：运行前应在官方源码上应用 proseg-coordinate-fix.patch 并重新编译。
"$PY" "$ROOT/run_proseg.py" \
  --binary "$ROOT/methods/proseg/target/release/proseg" \
  --binned-outputs "$ROOT/raw_p1crc/binned_outputs" \
  --mask "$ROOT/results/proseg/p1crc_stardist_masks.npy.gz" \
  --output "$ROOT/results/proseg/full" --nthreads 48

# 3. 对同一个先验聚合未经过 Proseg 修正的原始计数，供另两种方法使用。
"$PY" "$ROOT/aggregate_stardist_cells.py" \
  --binned-outputs "$ROOT/raw_p1crc/binned_outputs" \
  --mask "$ROOT/results/proseg/p1crc_stardist_masks.npy.gz" \
  --output "$ROOT/results/common/star_dist_cells_5um.h5ad" \
  --max-expansion-um 5

# 4. ResolVI。
"$RESOLVI_PY" "$ROOT/run_resolvi.py" \
  --input "$ROOT/results/common/star_dist_cells_5um.h5ad" \
  --output "$ROOT/results/resolvi/full" --epochs 100 --batch-size 256

# 5. SPLIT：先用 P1CRC scRNA 构建参考，然后运行 RCTD 双细胞模型与 SPLIT。
"$PY" "$ROOT/prepare_split_inputs.py" \
  --spatial "$ROOT/results/common/star_dist_cells_5um.h5ad" \
  --reference "$ROOT/adata_sc_all_reanno.h5ad" \
  --output "$ROOT/results/split/full_input" \
  --max-spatial-cells 10000 --max-reference-per-type 500 --seed 42
"$R" "$ROOT/run_split.R" \
  "$ROOT/results/split/full_input" "$ROOT/results/split/full" 16
```

若需单独复现压缩包中旧 8 µm bin 的输入口径，可在独立结果目录运行：

```bash
"$RESOLVI_PY" "$ROOT/run_resolvi.py" \
  --input "$ROOT/P1CRC_HD.h5ad" \
  --output "$ROOT/results/resolvi/legacy_8um" \
  --spatial-scale 0.27380817798463214 --min-counts 5 \
  --epochs 100 --batch-size 256
"$PY" "$ROOT/prepare_split_inputs.py" \
  --spatial "$ROOT/P1CRC_HD.h5ad" \
  --reference "$ROOT/adata_sc_all_reanno.h5ad" \
  --output "$ROOT/results/split/legacy_8um_input" \
  --spatial-scale 0.27380817798463214 \
  --spatial-sampling random --max-spatial-cells 10000 \
  --max-reference-per-type 500 --seed 42
"$R" "$ROOT/run_split.R" \
  "$ROOT/results/split/legacy_8um_input" \
  "$ROOT/results/split/legacy_8um" 16
```

### 运行与结果检查

检查每个脚本的退出码。Proseg 应产出 `counts.mtx.gz`、`cell_metadata.parquet` 和 `proseg-output.zarr/`；ResolVI 应产出 `model/`、`resolvi_latent.h5ad`、`resolvi_expression.h5ad`；SPLIT 应产出 `rctd.rds`、`split_result.rds`、`purified_counts.mtx`、`cell_metadata.csv` 和记录输入、线程数、随机种子的 `run_parameters.txt`。每个方法需记录实际参与指标计算的细胞数。旧聊天中的 SPLIT 大约 1 万输入、最终 7 千多用于指标只是历史观察值，不能代替本次运行的真实数值。

流水线完整结束后，运行 `"$PY" "$ROOT/verify_full.py" --root "$ROOT"`。该命令检查各阶段完成标记、矩阵维度、空间坐标、Proseg 退出码和 ResolVI/SPLIT 输出，并写入 `$ROOT/results/full_verification.json`。它还将 StarDist 原始细胞标签与 Proseg 细胞矩阵行号对应，保存为 `$ROOT/results/common/proseg_cell_correspondence.parquet`，供后续在同一批细胞上比较。

## 已验证的小区域

H&E 上 `(10000,20000)..(12048,22048)` 全分辨率像素裁剪得到 3,284 个 StarDist 核。对应的 2 µm bin 小数据有 78,600 个 bin、18,085 个基因。Proseg 试跑产出 3,284 个细胞的计数矩阵，其质心 `x=2739..3297 µm, y=5478..6035 µm` 与该裁剪区域吻合。ResolVI 在 500 个细胞上完成训练并输出解码表达。这些只是流程验证，不作为正式性能结果。

## Proseg 坐标修补

当前克隆的 Proseg 3.2.0（commit `4caa6f3`）在 Visium HD 读取时把 `pxl_row_in_fullres` 用作 x，把 `pxl_col_in_fullres` 用作 y。H&E 掩膜的变换使用列作为 x、行作为 y；`proseg-coordinate-fix.patch` 将读取改为 `x = col_px × microns_per_pixel`、`y = row_px × microns_per_pixel`。未应用这个补丁时核先验与转录本坐标错位。正式运行保留该版本的默认迭代数：200 次 burn-in、200 次采样、50 次 hill-climb。SPLIT 使用 commit `e880e39`（0.3.0），ResolVI 使用 scvi-tools 1.3.3。

## 可比性限制

StarDist 掩膜源于 H&E，采用 2 倍降采样。聚合原始细胞计数时向细胞核外扩 5 µm，这是明确的预处理参数，不能事后与 10x 自带分割的计数混称为相同输入。SPLIT 的 RCTD 会拒绝部分细胞；统计应同时给出输入数和输出数。ResolVI 的解码表达是连续值，不能当作原始 UMI 计数。
