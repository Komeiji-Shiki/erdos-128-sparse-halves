# 外部数学证书的来源

`sarid-clebsch-tangent.json` 由 **Amir Sarid** 发布，不是本项目原创。

- 来源：[aimir/erdos-128-sparse-halves](https://github.com/aimir/erdos-128-sparse-halves)
- 固定提交：`ba85f319ebb33ec2194ef4c2929a12916329e1e7`
- 原文件：`clebsch_tangent_certificate.json`
- SHA-256：`606b26d835d5b7fe70d077230f8091bbd9677283049a1d0b46a3e3ad30a200ce`
- 读取日期：2026-09-05

证书声明 `t(C4) >= rho³ - (3/64 + 1/50000000000) rho`。本地原验证程序运行通过；本项目另写的 `src/verify_four_cycle.py` 也用不同实现重建全部六顶点无三角形图和系数，验证通过。我们没有运行或借用原项目的 Lean 形式化作为本稿的证明。

该文件按原字节保存；`.gitattributes` 禁止对它及区间树作换行转换，以保持公布的哈希。论文明确引用原作者，未将其已有四边形不等式或证书据为新成果。
