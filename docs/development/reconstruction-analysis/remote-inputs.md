# 两仓库远程输入

本入口覆盖 P2 Xenium、P1 HD、P2 Visium 全量重建与分析，以及 Xenium 四方法比较。代码和配置随 Git 同步；下面的六个本地资源单独传输。两个远程仓库保持并列，名称为 `REVISE`、`REVISE_Analysis_Agent`。

## 准备上传目录

在 REVISE 根目录运行：

```bash
python scripts/prepare_remote_inputs.py
```

默认生成 `upload/remote-inputs/`；目录已存在时拒绝覆盖，可用 `--output` 指定新目录。脚本不修改源文件，不建立符号链接。该目录包含两仓库的目标路径、`manifest.json`、`SHA256SUMS` 和放置说明。`upload/` 不进入 Git。

| 仓库 | 数据相对路径 | 用途 |
|---|---|---|
| REVISE | `raw_data/Real_application/P2CRC_Xenium.h5ad` | iST 输入 |
| REVISE | `raw_data/Real_application/P1CRC_HD.h5ad` | hST 输入 |
| REVISE | `raw_data/Real_application/P2CRC_Visium.h5ad` | sST 输入 |
| REVISE | `raw_data/Real_application/adata_sc_all_reanno.h5ad` | 三路线参考；P2 过滤沿用配置 |
| REVISE | `output/mini-acceptance/20260920/assembly/original_spatial.h5ad` | 49279×0 历史分群与坐标基准；不是全量重建结果 |
| REVISE_Analysis_Agent | `resources/h.all.v2025.1.Hs.symbols.gmt` | Hallmark 通路分析资源 |

2026-09-20 六文件总计 3,114,093,428 字节，约 2.90 GiB；每次实际大小与 SHA256 以新生成的清单为准。清单中的本地源路径用于追溯，服务器放置使用相对目标路径。

## 传输和放置

后续上传工作流负责将整个目录送到远端。`qz-runner` 用于执行与核验；大数据须通过可用的远端存储或文件传输通道到位，不嵌入提交脚本。

先在收到的目录内运行：

```bash
sha256sum -c SHA256SUMS
# macOS 可使用 shasum -a 256 -c SHA256SUMS
```

将其中 `REVISE/`、`REVISE_Analysis_Agent/` 的内容分别合入已由 Git 配置好的同名仓库；只安装缺失文件或校验和相同的文件。遇到同路径不同内容时保留远端文件并停止放置，不直接覆盖。在两个仓库的共同父目录，用该清单的绝对路径再次运行 `sha256sum -c /path/to/remote-inputs/SHA256SUMS`，核验最终位置。

无需上传 `.venv`、旧分析库 `data/P2CRC_Xenium`、历史结果、benchmark 数据或完整 mini 交付。环境按两库现有安装说明配置。

## 运行顺序

1. 在 REVISE 根目录按[服务器全量入口](acceptance.md#服务器全量入口)运行六份全量配置，生成三路线交付与四种 Xenium assembly。环境、线程与容量限制沿用该入口。
2. 在 Analysis Agent 根目录运行 `configs/p2_full_project.yaml`、`p1_hd_full_project.yaml`、`p2_visium_full_project.yaml`。它们直接消费上游生成的 `sample.yaml + raw.h5ad + SVC.h5ad`，无需再复制一份原始数据。
3. 在 REVISE 根目录执行全量比较：

```bash
python scripts/execute_acceptance_notebook.py \
  --repo . \
  --source reproduce/case/assembly_comparison.ipynb \
  --comparison-config examples/assembly-comparison-full.json \
  --output output/full-acceptance/assembly/notebook
```

[全量比较配置](../../../examples/assembly-comparison-full.json)读取 `results/acceptance-full/P2CRC_Xenium/{mean,random,within_cluster,outside_cluster}/SVC.h5ad`，保存至 `output/full-acceptance/assembly/comparison/`。这些 SVC 是服务器运行后产物，未包含在上传目录中。原 [mini 配置](../../../examples/assembly-comparison-real.json)继续指向原有证据。

全量运行与科学解释尚未验收。hST 表达来源仍为 unknown；sST 坐标尺度与真零支持限制沿用当前验收说明，文件整理不改变这些科学条件。

## 本地清理记录

逐文件核验和删除回执保存在本地 `upload/cleanup-audit/`；处理汇总见同目录下的 `summary.md`。原始数据、历史 carrier、当前结果和两库 `mini-acceptance/20260920`、`20260920-closeout` 保留。只清理本次明确列出的旧分析快照、普通缓存及有保留原件覆盖的内容一致副本。
