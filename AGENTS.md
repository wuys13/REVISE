# REVISE 工作入口

- 本仓库说明：[README](README.md)。
- 公共工作区背景：[根目录 README](../../README.md)。
- 本地与启智协同协议：[协同协议](../../docs/collaboration.md)。
- REVISE 负责参考准备、重建、Raw 保护与 `sample.yaml + raw.h5ad + SVC.h5ad` 交付；下游分析由并列的 `../REVISE_Analysis_Agent` 负责。
- 保留用户和并发修改；历史 provenance、receipt、冻结结果与 Notebook 输出不做路径字符串替换。

## 接手与交付

- 按[根级接手规则](../../AGENTS.md)读取本仓状态入口 `docs/development/reconstruction-analysis/acceptance.md`，核对 `revise-2.0`、输入身份和允许修改范围。
- 上游输入由工作区登记的 SC/空间数据交付提供；下游交给并列的 `../REVISE_Analysis_Agent`，边界仍是 `sample.yaml + raw.h5ad + SVC.h5ad`。
- 结束时保存运行 receipt 或 acceptance 事实，并记录 commit、执行位置、解释器、输入指纹、输出入口和验证范围；重建完成不等于下游分析或科学接受完成。
