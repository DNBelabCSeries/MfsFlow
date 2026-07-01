# MfsFlow 性能与健壮性优化计划

## Context

MfsFlow 是一条单细胞 RNA-seq 流水线。通读核心代码后定位到 3 处真实性能瓶颈与 3 处低风险健壮性问题。本次优化目标是在**不改变任何分析结果、不破坏下游契约、保持现有测试通过**的前提下，消除 Filtering 阶段的冗余 I/O 与 CPU 开销，并修掉两处潜在 bug。

### 已核实的契约（必须保留，来自下游消费者审计）

- BAM tag 按名查找：下游全部用 `has_tag`/`get_tag`/`set_tag`，**无位置依赖** → 改写 BAM 写入方式安全
- flag 值是硬契约：SE=`4`，PE read1=`77`、read2=`141`
- `query_qualities` 必须是 int list 且长度与 `query_sequence` 一致（否则 [barcode_corrector._adjust_read_sequence](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/barcode_corrector.py#L101) 会置 None）
- 文件名契约：`{project}{suffix}.raw.tagged.bam`（[mapping_analysis.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/mapping_analysis.py#L354) 靠 `.raw.tagged.bam` 后缀判断是否流式纠错）
- 文本格式契约：
  - `BCstats.txt`：`bc\tcount`，无表头（[barcode_detection.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/barcode_detection.py) `pd.read_csv(sep='\t', header=None, names=['XC','n'])`）
  - `Q30stats.txt`：首行 `metric\ttotal_bases\tq30_bases`，后续 `metric\ttotal\tq30`（[test_q30_stats.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/tests/test_q30_stats.py) 断言行格式 `["R1","40","30","0.750000"]`）
- BAM header：可仅 `@HD` + `@PG`，SQ 可空（下游全部 `check_sq=False` 打开）

### 已核实的测试覆盖

- `test_fqfilter_logic.py`：仅测 `extract_seq`、`hamming_distance` — **不涉及 BAM 输出/check_qual/fastq_iter**
- `test_q30_stats.py`：测 `merge_q30_stats` 行格式 — **Q30 输出格式必须保留**
- fqfilter 的 BAM 写入、check_qual、fastq_iter、split_fastq、merge_bam_stats — **无测试覆盖**
- pysam 已是项目依赖（pyproject.toml `pysam>=0.21`），`unset_tag` 在 0.14+ 可用

---

## Phase A — 性能优化（高价值）

### A1. fqfilter.py：用 pysam 直接写二进制 BAM，替换「拼文本 SAM → samtools view -Sb」管道

**文件**：[mfsflow/scripts/fqfilter.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/fqfilter.py)（L272-278 子进程管道 + L375-402 手工拼 SAM 行）

**当前问题**：
- 手工拼 `b"\t".join([rid, b"4", b"*", ..., seq1_out, qual1_out])` 文本 SAM，再 pipe 给 `samtools view -Sb -` 转二进制。文本序列化 + samtools 文本解析双重开销
- 手工拼字段对 rid 中含空格/特殊字符不健壮

**改法**：
- 删除 `samtools_proc = subprocess.Popen([samtools, 'view', '-@', ..., '-Sb', '-'], stdin=PIPE, stdout=out_bam_fh)` 与 `out_bam_fh`
- 改为 `import pysam`；`outfile = pysam.AlignmentFile(out_bam, "wb", header=header_dict)`
  - `header_dict = {"HD": {"VN": "1.6"}, "PG": [{"ID": "MfsFlow-fqfilter", "PN": "MfsFlow-fqfilter", "VN": "3.0", "CL": " ".join(sys.argv[1:])}]}`
- `process_records` 中，对每条记录构造 `pysam.AlignedSegment()`：
  - `a.query_name = rid.decode()`（rid 已 strip 过 `\n\r`；rid 头部 `@` 已在原代码剥离，保留此逻辑）
  - `a.flag = 4`（SE）或 `77`（PE R1）/ `141`（PE R2）— 与原代码完全一致
  - `a.query_sequence = final_cdna.decode() if final_cdna else None`（保留 `*` 等价语义：None → SAM 输出 `*`）
  - `a.query_qualities = [b - 33 for b in final_cdna_q] if final_cdna_q else None`（bytes 逐字节减 33 得 int list；长度自动与 seq 一致）
  - `a.set_tag("CR", final_bc.decode(), "Z")`、`("UR", final_umi.decode(), "Z")`、`("CY", final_bc_q.decode(), "Z")`、`("UY", final_umi_q.decode(), "Z")`
  - `outfile.write(a)`
- PG header 改为在 header_dict 里声明，不再手工写 `@PG` 行
- `samtools` 与 `samtools_threads` 参数从 fqfilter 命令行移除（不再需要）；同步更新 [stages/filtering.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/stages/filtering.py#L143-L148) 构造 fqfilter 命令处：去掉 `samtools` 位置参数与 `--samtools-threads`，pigz 仍保留（用于解压输入 chunk）

**契约保留检查**：flag 4/77/141 ✓、tag 名 CR/UR/CY/UY ✓、`.raw.tagged.bam` 文件名 ✓（`out_bam` 路径不变）、PG 行 ✓（在 header 里）、BCstats/Q30stats 输出不变 ✓（这两段在 BAM 写入之后，与改写无关）。

### A2. fqfilter.py：check_qual 用 bytes.translate 查表，替换逐字节 Python 循环

**文件**：[mfsflow/scripts/fqfilter.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/fqfilter.py#L289-L298)

**当前问题**：每条 read 对 BC（20bp）和 UMI（10bp）各跑一次 `for q in q_str: if q < limit` 逐字节 Python 循环。

**改法**：在 main() 开头按 phred 阈值预建查表（与现有 `Q30_TABLE` 同模式）：
```python
def _lowq_table(phred_threshold):
    limit = phred_threshold + 33
    return bytes(1 if i < limit else 0 for i in range(256))

BC_LOWQ_TABLE = _lowq_table(bc_filter[1])
UMI_LOWQ_TABLE = _lowq_table(umi_filter[1])

def check_qual(q_str, threshold_count, table):
    return sum(q_str.translate(table)) < threshold_count
```
调用处改传对应 table。**注意**：原实现是「累计到 threshold_count 即返回 False」早退；translate+sum 全量扫描。对 20bp/10bp 的短串差异在纳秒级，translate 在 C 层更快。语义等价（都是「低质量碱基数 ≥ threshold_count 则过滤」）。

### A3. pipeline_modules.split_fastq：SeqKit `-s` 单遍切，替换 `-p` 两遍

**文件**：[mfsflow/pipeline_modules.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/pipeline_modules.py#L181-L195)

**当前问题**：`seqkit split2 -p {split_parts}` 会先 count 再 split（两遍读），代码注释自己也承认。`lines_per_chunk` 已经算好但只用于 GNU split fallback。

**改法**：SeqKit 分支改用 `-s`（每份记录数，单遍）：
```python
records_per_chunk = max(1, lines_per_chunk // 4)  # lines_per_chunk 已 4 对齐
if mode == "PE":
    cmd = f"{seqkit_cmd} split2 -s {records_per_chunk} -1 {shlex.quote(fpath1)} -2 {shlex.quote(fpath2)} -O {shlex.quote(out_dir)} -f -j {seqkit_threads} {ext_flag}"
else:
    cmd = f"{seqkit_cmd} split2 -s {records_per_chunk} -O {shlex.quote(out_dir)} -f {shlex.quote(fpath1)} -j {seqkit_threads} {ext_flag}"
```
**契约保留**：输出文件名仍为 `.part_NNN` 模式，后续重命名逻辑不变；`split_parts` 变量仍用于估算 `lines_per_chunk`，保留。

---

## Phase B — 健壮性 / 清理（低风险）

### B1. barcode_corrector.py：`set_tag("CC", None)` → `unset_tag`

**文件**：[mfsflow/scripts/barcode_corrector.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/scripts/barcode_corrector.py#L90-L91)

**问题**：`read.set_tag("CC", None)` / `read.set_tag("CB", None)` 在 pysam 不同版本下行为不一致（可能写空值而非删除 tag）。下游 `stream_corrector.get_or_apply_correction` 用 `read.has_tag("CC")` 判断是否已纠错——若空值 tag 存在则 has_tag 返回 True，会跳过纠错。

**改法**：
```python
else:
    read.unset_tag("CC") if read.has_tag("CC") else None
    read.unset_tag("CB") if read.has_tag("CB") else None
```
或更简洁，封装一个 `safe_unset_tag` helper。pysam 0.21+ 支持 `unset_tag`。

### B2. 删除 process_fastq_inputs 死代码

**文件**：[mfsflow/bootstrap.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/bootstrap.py#L145-L154) + [mfsflow/cli.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/cli.py#L98-L99)

**问题**：`process_fastq_inputs` 是空占位函数，CLI 每次都调用但只 `return None`，误导读者以为有处理。

**改法**：删除函数定义 + 删除 cli.py 中 `process_fastq_inputs(config)` 调用与下一行 `log_info('Fastq processed.')`。

### B3. stages/filtering.py：移除无用的 multiprocessing.Pool 包装

**文件**：[mfsflow/stages/filtering.py](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/stages/filtering.py#L101-L121)

**问题**：`Pool(processes=min(2, num_threads))` 里只 `apply_async` 一个任务再 `get()`——纯进程 fork 开销，无并行收益。真正的并行在 `split_fastq` 内部用 Popen 实现。

**改法**：去掉 Pool，直接调用：
```python
chunk_suffixes = pipeline_modules.split_fastq(
    fq1_files, num_threads, lines_per_chunk, tmp_merge_path, project,
    pigz, seqkit, fq2_files, False, split_parts,
)
```

---

## 不在本次范围（已评估，暂缓）

- **split_fastq 文件重命名启发式**（[pipeline_modules.py L291-466](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/pipeline_modules.py#L291)）：脆但能用，重写收益不抵风险
- **merge_bam_stats 运行中改写 YAML + flag==4 推断 layout**（[pipeline_modules.py L579-600](file:///Users/lishuangshuang/Documents/project/Mhsflt_toolkit/mfsflow/pipeline_modules.py#L579)）：行为性变更，需更深重构，单独立项
- **scripts/ 双导入 fallback、两份 path_layout.py、which_Stage 命名不一致**：清理项，触及多文件，价值低
- **run_featurecounts.py / report.py 拆分**：单体文件大但功能正常，无性能问题

---

## 验证

1. **单元测试**：`python3 -m unittest discover -s tests` 必须 100% 通过（重点 `test_fqfilter_logic.py`、`test_q30_stats.py`、`test_barcode_corrector.py`、`test_stream_corrector.py`、`test_stage_state.py`）
2. **白盒核查**：
   - fqfilter 改写后，手工核查生成的 BAM：`pysam.AlignmentFile(out_bam, "rb", check_sq=False)` 能打开、首 read 的 `flag`/`query_sequence`/`get_tag("CR")` 与原文本 SAM 一致
   - check_qual：构造 qual bytes 使低质量碱基数恰为 threshold_count-1 / threshold_count / threshold_count+1，验证三种边界
   - split_fastq `-s`：确认输出文件数与原 `-p` 接近（records_per_chunk 估算准的话应等于 split_parts）
3. **whitespace**：`git diff --check` 无报错
4. **端到端冒烟**（受限）：本机为 macOS，bundled STAR/samtools/seqkit 为 Linux ELF 无法本地跑全流程。需在 Linux 环境用小 FASTQ 跑一次 `mfsflow --fastqs ... --stage Filtering`，比对改写前后 `*.raw.tagged.bam` 的 read 数与 tag 内容一致

## 修改文件清单

| 文件 | 改动 |
|------|------|
| mfsflow/scripts/fqfilter.py | A1 pysam 直写 BAM + A2 check_qual translate + 移除 samtools 参数 |
| mfsflow/stages/filtering.py | A1 同步去掉 fqfilter 命令的 samtools/samtools-threads 参数 + B3 去 Pool |
| mfsflow/pipeline_modules.py | A3 SeqKit `-s` 单遍切 |
| mfsflow/scripts/barcode_corrector.py | B1 unset_tag |
| mfsflow/bootstrap.py | B2 删 process_fastq_inputs |
| mfsflow/cli.py | B2 删调用 |
