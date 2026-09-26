# 超分辨率 benchmark 复现交接文档

## 1. 交接范围与当前状态

本目录整理了 REVISE 的 Sim2Real-ST spot-size 数据上 **TESLA、iStar** 的输入适配、运行、预测对齐和统一评分代码，并提供 REVISE 自身 spot 路线的运行脚本。范围为 P1CRC 的 `part1`、`part2`、`part3`，每个区域有 50、100、150、200 µm 四种伪 spot，共 12 个数据案例，基因数均为 324。

截至 2026-09-26，qz 上 TESLA 和 iStar 共 **24/24 次运行已完成并通过结果校验**。REVISE 另有两个成功的单例试测，且 `part2` 四个尺寸的整条路线已成功运行；由于下载包缺少两份原始输入，REVISE 的试测使用了重建输入，**不能称为论文结果的严格复现**。仓库只包含代码、说明和小体积汇总表，不包含 Xenium/H&E 原始数据、模型权重或大量细胞级预测。

| 方法 | 已完成 | qz 计时口径与结果 |
| --- | --- | --- |
| TESLA | 12/12 | 推断中位数 85.5 秒，范围 26.4–1686 秒；不含预处理和统一评分 |
| iStar | 12/12 | 每例推断及评分中位数 278.5 秒，范围 215–1286 秒；另需每区域约 0.5–2 分钟、可复用的 HIPT 图像特征提取 |
| REVISE | 2 个单例试测 + `part2` 四尺寸整条路线 | 单例 `part1/50` 为 25 秒、`part2/200` 为 34 秒；`part2` 四尺寸合计 58 秒，均包含启动、输出和评分 |

这些数字不是控制了输入和硬件条件的算法速度排名。REVISE spot 路线读取细胞 ID、spot 归属和 Xenium 细胞坐标；TESLA/iStar 输入为伪 spot 表达、spot 坐标和 H&E 图像。REVISE 的真实细胞类型标签用于诊断匹配率；其 spot SR 分配步骤使用预测类型。iStar 运行中有部分 CPU 线程竞争异常值，见下文。

## 2. 推荐的 GitHub 交接方式

把本目录作为 REVISE 仓库中的 `benchmark/spot_sr/`，通过 **fork 分支 + Pull Request** 交给维护者审核；师兄可立即从 fork 下载源码或 ZIP，仓库维护者合并后即可从 REVISE 主仓库获取。不建议把数据、H&E、HIPT 权重、Python 虚拟环境或细胞级输出提交到 Git。当前 GitHub 账号对 `wuys13/REVISE` 只有 `READ` 权限；若需要直接向其分支 `git push`，须由仓库管理员邀请该账号为 collaborator。

仓库代码可追踪和审阅，压缩包适合一次性交接；大规模 benchmark 应以 Git commit 作为代码版本标识，而不是只依赖压缩包。结果文件应记录仓库 commit、各方法 commit、输入校验和、随机种子、CPU/GPU、线程数、是否命中图像特征缓存和计时边界。

## 3. qz 目录与数据

本次的工作根目录：

```bash
export BENCHMARK_ROOT=/inspire/hdd/global_user/zhangyinghao-240107100021/REVISE_benchmark
cd "$BENCHMARK_ROOT"
```

运行脚本假设如下文件存在：

```text
spot/part1, part2, part3/
  selected_xenium.h5ad
  real_sc_ref_part.h5ad
  spot_50, spot_100, spot_150, spot_200/xenium_spot.h5ad
Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_image.ome.tif
Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_he_imagealignment.csv
```

H&E 原图及配准 CSV 用于 TESLA/iStar。`prepare_inputs.py` 将 Xenium 微米坐标通过 10x 对齐点映射至 H&E 像素，按区域加 250 µm 边界裁切，再重采样至 0.5 µm/像素。已测文件的对齐点残差约 2.87 个原始 H&E 像素。扩展至新样本时，应检查 `prepared/<part>/he-raw.jpg` 与 `image_transform.json`，尤其是 spot 是否落在组织区域。

`selected_xenium.h5ad` 对 TESLA/iStar 只在**模型推断完成后**用于取细胞坐标、采样预测并计算指标；其表达矩阵不能进入模型推断。

## 4. 固定源码与环境

| 方法 | 已测 commit |
| --- | --- |
| REVISE | `c83dc97d25b6d513b59cc301255e5bdc7e9c7cd9` |
| iStar | `3cb0e5352a86df337f41c5841d0808da0003457b` |
| TESLA | `8c4dcf896497dd4a34a6e720f03ec68c1502d571` |

```bash
mkdir -p methods envs
git clone https://github.com/daviddaiweizhang/istar.git methods/istar
git clone https://github.com/jianhuupenn/TESLA.git methods/TESLA
git -C methods/istar checkout 3cb0e5352a86df337f41c5841d0808da0003457b
git -C methods/TESLA checkout 8c4dcf896497dd4a34a6e720f03ec68c1502d571

# 本目录放在 methods/REVISE/benchmark/spot_sr/ 后：
SCRIPT="$BENCHMARK_ROOT/methods/REVISE/benchmark/spot_sr"
git -C methods/istar apply "$SCRIPT/istar_thread_limit.patch"
```

本次 qz 使用基于 `/opt/conda/envs/spacec` 的 Python 3.10 venv。基底已有 PyTorch 2.5.1、NumPy 1.26.4、Scanpy 1.11.5、AnnData 0.11.4、pyvips 等。若在**同一个 qz 镜像**重建额外依赖，可参考：

```bash
/opt/conda/envs/spacec/bin/python -m venv --system-site-packages envs/baselines
envs/baselines/bin/pip install --no-deps \
  einops==0.6.1 pytorch-lightning==2.0.8 torchvision==0.20.1 \
  lightning-utilities==0.15.2 torchmetrics==0.11.4

/opt/conda/envs/spacec/bin/python -m venv --system-site-packages envs/revise
envs/revise/bin/pip install 'POT==0.9.5' 'squidpy==1.6.5' 'scikit-misc==0.5.2'
envs/revise/bin/pip install --no-deps 'numpy==1.26.4'
```

这是一份**已测 qz 镜像的环境说明**，不是对全新服务器通用的完整 lockfile。安装 Squidpy 时，pip 曾把 NumPy 升至 2.x，导致 qz 原有二进制包导入失败；因此要固定回 1.26.4，并重新检查 Scanpy/Squidpy/POT 导入。HIPT 的两个官方 checkpoint 要放在 `methods/istar/checkpoints/`；源仓库的 Box 链接本次不可用，实际使用 [HIPT 官方仓库](https://github.com/mahmoodlab/HIPT)的权重。权重 SHA-256 见 [README.md](README.md)。

## 5. TESLA 与 iStar 完整运行

```bash
cd "$BENCHMARK_ROOT"
SCRIPT="$BENCHMARK_ROOT/methods/REVISE/benchmark/spot_sr"
bash "$SCRIPT/run_all.sh" prepare
bash "$SCRIPT/run_all.sh" tesla
GPU=0 bash "$SCRIPT/run_all.sh" istar

# 如果某个案例失败，可以只运行它；有完整评分输出的案例自动跳过。
PARTS=part2 SIZES=200 bash "$SCRIPT/run_all.sh" tesla
PARTS=part2 SIZES=200 GPU=0 bash "$SCRIPT/run_all.sh" istar

envs/baselines/bin/python "$SCRIPT/summarize.py" \
  --prepared prepared --output results/baselines_summary.csv
envs/baselines/bin/python "$SCRIPT/validate_results.py" \
  --prepared prepared --output results/validation.json
```

TESLA 采用 spot 的 `log1p` 计数、10 µm 超像素格点及 10 个邻居，并把每个评估细胞映射到最近的预测超像素。iStar 每个区域只提取一次默认 shifted HIPT 图像特征，四种 spot 尺寸共用；逐例以 `--epochs=400 --n-states=5` 运行，再从包含细胞中心的 16 像素图像 patch 取预测。iStar 预设圆形捕获范围，但本数据的伪 spot 为方格，适配器以方格半宽作为半径代理；解释指标时应保留该限制。

预测、每基因指标和日志保存在 `$BENCHMARK_ROOT/prepared/` 与 `$BENCHMARK_ROOT/logs/`。统一评分按 REVISE 的规范：每细胞总量归一到 10,000、每基因分别做 min-max 归一化，然后计算 PCC、SSIM、MSE、NRMSE。`validate_results.py` 验证 24 个案例的细胞数、基因集、指标有限性和汇总一致性。`results/` 下的小表是本次完成结果的快照。

## 6. REVISE 自身路线：先确认输入来源

REVISE 上游 spot 配置要求每个案例的 `real_sc_ref_all.h5ad`、数据根目录的 `PM_on_cell.csv`，路径格式为 `spot/P1CRC/cut_partN/`。当前下载包只有三个区域各自的 `real_sc_ref_part.h5ad`，也没有 PM 文件。因此严格复现需要找到作者原始的完整参考和 PM；**不要把下面的替代输入称为原始 benchmark 输入。**

如果目的是先确认代码能运行并估算耗时，可以显式构造适配输入：

```bash
cd "$BENCHMARK_ROOT"
SCRIPT="$BENCHMARK_ROOT/methods/REVISE/benchmark/spot_sr"
envs/revise/bin/python "$SCRIPT/prepare_revise_reference.py" --spot-root spot
envs/revise/bin/python "$SCRIPT/prepare_revise_pm.py" --spot-root spot --seed 0
PARTS='part1 part2 part3' bash "$SCRIPT/run_revise.sh"
```

参考构造脚本对三个子集按 cell ID 去重，得到 15,280 个细胞、18,071 个基因，并生成上游所需的路径链接。PM 脚本仿照 REVISE Visium 示例，用参考数据中的细胞类型比例加极小随机扰动产生 47,606 × 11 的先验矩阵，并写入 `PM_on_cell.provenance.json`。已有 PM 文件时默认拒绝覆盖。REVISE 自身的 `run_revise.sh` 使用上游四个尺寸的配置逐区域运行，输出在 `results/revise/`、日志在 `logs/revise_partN.log`。

两个单例试测分别是 `part1/spot_50`（37,837 细胞，25 秒）和 `part2/spot_200`（6,292 细胞，34 秒）；`part2` 四尺寸整条路线运行成功，总耗时 58 秒。`part1/spot_50` 没有执行图聚合更新，所以单例时间与细胞数不单调；不要由此估算完整 12 例或新数据的总时间。

## 7. 扩展到大规模 benchmark 的建议

1. **先定输入协议。** 记录每个方法允许使用的 spot counts、H&E、参考单细胞、细胞坐标/分割、spot-cell 归属；若输入信息不同，应分别呈现结果并解释，不能只做一个混合排名。
2. **以案例为单位保存清单。** 建议每个新案例包含样本 ID、组织/平台、spot 数、细胞数、基因面板、H&E 分辨率、坐标变换和校验和。扩展新样本时修改 `prepare_inputs.py` 中的样本识别与配准参数，而不是仅改目录名。
3. **固定计时边界。** 同时报告冷启动总耗时、预处理/共享图像特征耗时、逐案例推断耗时、评分耗时，并记录 GPU 型号、线程数及并发任务。现有三个方法的计时口径不完全一致。
4. **保留可检查产物。** 每例保存基因顺序一致的预测 AnnData、每基因指标、日志和方法版本；汇总时按区域/样本聚合，避免把不同区域的细胞直接合并后重新计算 PCC。
5. **先小规模试跑再扩张。** 用 1 个小、1 个密集案例检验内存与时间；之后按区域切分并行。TESLA 的 `part1/spot_50` 本次耗时约 28 分钟，是现有数据中的明显重负载案例。

`annotate_execution.py` 仅用于标注本次历史 iStar 结果中哪些案例在加线程限制补丁之前启动，不应用它标注以后新增的运行。
