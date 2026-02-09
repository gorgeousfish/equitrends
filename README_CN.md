# equitrends

双重差分估计中预趋势的等价性检验

[![Stata](https://img.shields.io/badge/Stata-16%2B-blue)](https://www.stata.com/)
![Version](https://img.shields.io/badge/version-0.1.0-informational)
![License](https://img.shields.io/badge/license-AGPL--3.0-blue)

![逆转举证责任](images/image.png)

EQUITRENDS 是一个 Stata 软件包，实现了双重差分（DiD）设计中预趋势的等价性检验，基于 [Dette & Schumann (2024)](https://doi.org/10.1080/07350015.2024.2308121)，发表于 *Journal of Business & Economic Statistics*。

## 目录

- [概述](#概述)
- [系统要求](#系统要求)
- [安装](#安装)
- [获取帮助](#获取帮助)
- [数据要求](#数据要求)
- [快速开始](#快速开始)
- [命令参考](#命令参考)
- [检验选择指南](#检验选择指南)
- [方法论](#方法论)
- [RMS 检验的 alpha 限制](#rms-检验的-alpha-限制)
- [存储结果](#存储结果)
- [可视化](#可视化)
- [示例](#示例)
- [引用](#引用)
- [作者](#作者)
- [许可证](#许可证)
- [未来路线图](#未来路线图)
- [相关软件包](#相关软件包)
- [推荐资源](#推荐资源)
- [支持](#支持)

## 概述

DiD 设计中的标准预趋势检验使用*精确*平行趋势的零假设（H₀: β = 0）。这种方法存在根本性的局限：

1. **未拒绝 ≠ 支持证据**：统计功效不足可能无法检测到实际的违反。未能拒绝并不支持平行趋势假设。
2. **条件偏差放大**：Roth (2022) 证明，当违反存在时，以通过传统预检验为条件反而会*放大* DiD 偏差。
3. **无明确阈值**：传统检验没有提供判断何为"可忽略"的平行趋势偏差的框架。

EQUITRENDS 实现了**逆转举证责任**的*等价性*检验：零假设为偏差*较大*（H₀: ‖β‖ ≥ 阈值），拒绝则提供偏差*较小*的统计证据。这种方法：

- 要求明确论证等价性阈值
- 控制第一类错误（错误地得出等价结论）
- 统计功效随样本量增加而提高
- 允许研究者量化等价性成立的最小阈值

### 主要特性

- **三种等价性假设**：最大值、均值和均方根（Dette & Schumann, 2024, 第 3.1 节）
- **最小等价性阈值**：计算能够得出等价结论的最小阈值
- **最大值检验的多种推断方法**：IU（解析法）、球形 Bootstrap 和 Wild Bootstrap
- **可视化**：带等价性界限的系数图（`equivtest_plot`）

## 系统要求

- **Stata 16.0** 或更高版本
- **Python（Stata Python 集成）** 以及 **numpy** 和 **scipy**，用于基于 Bootstrap 的最大值检验：
  - `equivtest ..., type(max) method(boot)`
  - `equivtest ..., type(max) method(wild)`

**注意**：`type(max)` 的 IU 方法以及 `type(mean)` / `type(rms)` 检验不需要 Python。Bootstrap 方法（`method(boot)` 和 `method(wild)`）需要安装 Python 及 **numpy** 和 **scipy**。

## 安装

### 从 GitHub 安装

```stata
* 安装软件包
net install equitrends, from("https://raw.githubusercontent.com/gorgeousfish/equitrends/main") replace
```

### 本地安装

```stata
* 安装软件包
net install equitrends, from("/path/to/equitrends-main") replace
```

### 加载示例数据

```stata
* 加载软件包自带的示例数据集
equitrends_data, clear
```

该数据集随软件包一起安装，从本地安装目录加载。

## 获取帮助

安装后，可访问 Stata 内置文档：

```stata
help equitrends       // 软件包概述和命令列表
help equivtest        // 主要统一检验接口
help equivtest_plot   // 可视化选项
help maxequivtest     // 最大值检验详情
help meanequivtest    // 均值检验详情
help rmsequivtest     // RMS 检验详情
help equivsim         // 蒙特卡洛模拟
help equitrends_data  // 加载示例数据集
```

## 数据要求

**使用本软件包前，请确保您的数据满足以下要求：**

| 要求                          | 说明                                                                                                 |
| :---------------------------- | :--------------------------------------------------------------------------------------------------- |
| **面板结构**            | 包含个体（`id`）和时间（`time`）标识符的面板数据。**支持平衡面板和非平衡面板。**           |
| **最少处理前期数**      | 至少一个处理前时期（*T* ≥ 1）。更多时期可提高统计功效。                                           |
| **处理组指示变量**      | 二值变量，编码为 0（对照组）或 1（处理组）。                                                         |
| **统一处理设计**        | 所有处理单元必须在同一时间接受处理。交错处理需要按队列分析（参见 Dette & Schumann, 2024, 第 5 节）。 |
| **完整的时间-组别单元** | 每个时间段必须同时包含处理组和对照组的观测值。                                                       |

### 非平衡面板

本软件包自动检测和处理非平衡面板（即个体具有不同数量的观测时期）。当检测到非平衡面板时：

- `e(is_balanced)` 返回 0
- `e(T_min)` 和 `e(T_max)` 存储每个个体的最小和最大时期数
- 输出显示 "Panel type: Unbalanced" 及时期范围

无需特殊语法——按常规方式运行命令即可：

```stata
xtset id time  // 如果面板不平衡，会显示 "unbalanced"
equivtest y, type(max) id(id) group(treat) time(time) pretreatment(1 2 3 4 5) baseperiod(5)
```

### 数据完整性

当某些时期缺少某一组别的观测值时，安慰剂回归无法正确估计。在这种情况下：

- 命令将报错或产生缺失值
- 请检查每个时间-组别单元至少包含一个观测值

**常见原因：**

- 处理变量、结果变量或标识变量中存在缺失值
- 样本限制导致某一组别在整个时期内没有观测值

**解决方案：** 确保每个时间-处理单元至少有一个观测值，或将分析限制在具有完整覆盖的时期内。在运行等价性检验前，使用 `xtset` 验证面板结构。

## 快速开始

### 实证示例：Di Tella & Schargrodsky (2004)

本示例复现了 Dette & Schumann (2024, 第 7 节) 的实证应用，使用 Di Tella & Schargrodsky (2004) 的布宜诺斯艾利斯犯罪数据。该数据集包含 876 个布宜诺斯艾利斯街区的月度汽车盗窃数量（1994 年 4 月至 12 月），其中 37 个街区在 7 月恐怖袭击后获得了警察保护。

```stata
* 加载并准备数据
equitrends_data, clear
drop if mes == 72 | mes == 73                        // 删除处理后月份
drop if mes > 7                                      // 仅保留处理前时期
rename (observ totrob mes) (ID Y period)
gen G = (distanci == 0)                              // 处理指示变量
xtset ID period

* 最大值检验（IU 方法，聚类稳健标准误）
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)

* 均值检验
equivtest Y, type(mean) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) vce(cluster) cluster(ID)

* RMS 检验
equivtest Y, type(rms) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) seed(2024)

* 可视化结果
equivtest_plot, ci
```

### Bootstrap 方法（仅限最大值检验）

```stata
* 球形 Bootstrap（假设球形误差；定理 1）
equivtest Y, type(max) method(boot) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(12345)

* Wild Bootstrap（推荐用于非球形误差；注释 1(c)）
equivtest Y, type(max) method(wild) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) nboot(1000) seed(12345)
```

### 控制变量（条件平行趋势）

`x()` 选项在 TWFE 安慰剂回归中加入额外的控制变量（Dette & Schumann, 2024, 第 5 节）。时间不变的协变量在双重去均值过程中被个体固定效应吸收。对于完整的条件 PTA 规范（公式 5.5），请在数据中构造协变量与时间的交互项，然后通过 `x()` 传入。

```stata
equivtest Y, type(max) method(iu) id(ID) group(G) time(period) ///
    pretreatment(4 5 6 7) baseperiod(7) x(edpub estserv banco) ///
    vce(cluster) cluster(ID)
```

## 命令参考

| 命令                | 说明                             |
| :------------------ | :------------------------------- |
| `equivtest`       | 三种检验的统一接口（推荐使用）   |
| `maxequivtest`    | 最大绝对系数检验（IU/Boot/Wild） |
| `meanequivtest`   | 均值系数检验                     |
| `rmsequivtest`    | 均方根（RMS）检验                |
| `equivtest_plot`  | 带等价性界限的系数图             |
| `equivsim`        | 蒙特卡洛模拟（功效分析）         |
| `equitrends_data` | 加载软件包自带的示例数据集       |

### 统一语法

```stata
equivtest depvar, type(max|mean|rms) id(varname) group(varname) time(varname) [options]
```

### 核心选项

| 选项                      | 说明                                       |
| :------------------------ | :----------------------------------------- |
| `type(max/mean/rms)`    | 检验类型（必需）                           |
| `id(varname)`           | 面板标识符（必需）                         |
| `group(varname)`        | 处理组指示变量 0/1（必需）；也接受 `g()` |
| `time(varname)`         | 时间变量（必需）；也接受 `period()`      |
| `threshold(#)`          | 等价性阈值；省略则计算最小阈值             |
| `alpha(#)`              | 显著性水平；默认 0.05                      |
| `pretreatment(numlist)` | 纳入分析的处理前时期                       |
| `baseperiod(#)`         | 安慰剂构造的基期                           |
| `x(varlist)`            | 控制变量                                   |

`type(max)` 特有选项：

| 选项                     | 说明                          |
| :----------------------- | :---------------------------- |
| `method(iu/boot/wild)` | 推断方法；默认 `iu`         |
| `nboot(#)`             | Bootstrap 重复次数；默认 1000 |
| `seed(#)`              | Bootstrap 随机种子            |
| `nodots`               | 抑制 Bootstrap 进度显示       |

`type(rms)` 特有选项：

| 选项            | 说明                             |
| :-------------- | :------------------------------- |
| `nolambda(#)` | 自归一化使用的子样本数量；默认 5 |
| `seed(#)`     | 子抽样随机种子                   |

稳健标准误选项（适用于 IU/均值检验；不适用于 Bootstrap 或 RMS）：

| 选项                 | 说明                                |
| :------------------- | :---------------------------------- |
| `vce(vcetype)`     | 方差估计量；见下表                  |
| `cluster(varname)` | 聚类变量（聚类稳健 VCE 类型时必需） |

**方差估计量类型（`vcetype`）：**

| 类型                | 说明                                                            |
| :------------------ | :-------------------------------------------------------------- |
| `ols`             | 同方差 OLS 方差（默认）                                         |
| `robust`/`hc1`  | HC1 异方差稳健（White, 1980）                                   |
| `hc2`             | HC2 杠杆调整（MacKinnon & White, 1985）；有限样本性质更优       |
| `hc3`             | HC3 更保守的杠杆调整（Davidson & MacKinnon, 1993）              |
| `hac`             | 面板数据的 Arellano (1987) HAC 估计量                           |
| `cluster`/`cr0` | CR0 聚类稳健，无小样本调整；需要 `cluster()`                  |
| `cr1`             | CR1 聚类稳健，含 G/(G-1) 调整（Stata 默认）；需要 `cluster()` |
| `hc1_cluster`     | HC1 聚类稳健，含小样本调整；需要 `cluster()`                  |

## 检验选择指南

EQUITRENDS 提供三种具有不同性质的等价性检验。请根据您的研究情境选择：

| 特征               | 最大值检验       | 均值检验             | RMS 检验     |
| :----------------- | :--------------- | :------------------- | :----------- |
| **假设**     | max\|βₜ\| < δ | \|β̄\| < τ        | β_RMS < ζ  |
| **度量**     | 最大单一违反     | 平均违反             | 均方根       |
| **抵消效应** | 无               | 有（相反符号会抵消） | 无           |
| **敏感性**   | 任何单一大偏差   | 系统性方向偏差       | 各偏差间均衡 |
| **保守程度** | 最保守           | 最不保守             | 适中         |

### 建议

1. **最大值检验（`type(max)`）**：作为默认的保守选择，建议首先使用。

   - 检测任何单一大违反
   - 使用 `method(iu)` 进行解析推断，或使用 `method(wild)` 处理非球形误差
   - 当您希望排除*任何*实质性预趋势违反时推荐使用
2. **均值检验（`type(mean)`）**：当预期违反为单调（同号）时使用。

   - 当偏差方向一致时功效更高
   - **注意**：相反方向的违反可能相互抵消，导致错误的等价结论
3. **RMS 检验（`type(rms)`）**：作为通用替代方案使用。

   - 对所有安慰剂系数的敏感性均衡
   - 无抵消问题
   - 自归一化（无需方差估计）

### 解读最小阈值

当省略 `threshold()` 时，EQUITRENDS 报告在指定显著性水平下能够得出等价结论的最小等价性阈值（δ*、τ* 或 ζ*）。将其与您估计的处理效应进行比较：

- **δ* << 估计的 ATT**：预趋势可忽略的强证据
- **δ* ≈ 估计的 ATT**：预趋势违反可能解释处理效应
- **δ* >> 估计的 ATT**：平行趋势的证据不足；考虑替代设计

## 方法论

设 $\beta = (\beta_1,\ldots,\beta_T)'$ 为 TWFE 安慰剂回归（Dette & Schumann, 2024, 公式 (2.5)）的安慰剂（处理前）系数向量。EQUITRENDS 实现了三种等价性假设（第 3.1 节）：

1. **最大偏差（公式 (3.1)）**：

$$
H_0: \|\beta\|_{\infty} \ge \delta \quad \text{vs.} \quad H_1: \|\beta\|_{\infty} < \delta, \qquad \|\beta\|_{\infty}=\max_{l\in\{1,\ldots,T\}}|\beta_l|
$$

2. **均值偏差（公式 (3.2)）**：

$$
\bar{\beta}=\frac{1}{T}\sum_{l=1}^{T}\beta_l, \qquad H_0: |\bar{\beta}| \ge \tau \quad \text{vs.} \quad H_1: |\bar{\beta}| < \tau
$$

3. **均方根偏差（公式 (3.3)）**：

$$
\beta_{\mathrm{RMS}}=\sqrt{\frac{1}{T}\sum_{l=1}^{T}\beta_l^2}, \qquad H_0: \beta_{\mathrm{RMS}} \ge \zeta \quad \text{vs.} \quad H_1: \beta_{\mathrm{RMS}} < \zeta
$$

### 最大值检验的推断方法

- **IU（交叉-联合，解析法）**：对每个安慰剂系数 *t* = 1, ..., *T*，当且仅当所有 |beta_t| < Q(alpha) 时拒绝 H0，其中 Q 为均值为 delta、方差为 sigma_tt/n 的折叠正态分布的 alpha 分位数（Dette & Schumann, 2024, 公式 (4.4)）。计算效率高但对大 *T* 较为保守。折叠正态 CDF/分位数在 Mata 中实现。
- **Bootstrap**（`method(boot)`）：在 beta 约束下使用约束 OLS 生成 Bootstrap 样本，然后计算经验 alpha 分位数作为临界值（Dette & Schumann, 2024, 定理 1）。假设球形误差。对 *T* > 1 比 IU 更有功效。需要 Python 及 numpy 和 scipy。
- **Wild Bootstrap**（`method(wild)`）：用 Rademacher 加权残差替代 i.i.d. Bootstrap 误差，使检验对异方差和序列相关稳健（Dette & Schumann, 2024, 注释 1(c)）。推荐用于非球形误差。需要 Python 及 numpy 和 scipy。

### 均值检验的推断

当安慰剂系数的绝对样本均值低于均值为 tau、方差为 1'Sigma1/(nT^2) 的折叠正态分布的 alpha 分位数时，均值检验拒绝 H0（Dette & Schumann, 2024, 公式 (4.12)）。

### RMS 检验的推断

RMS 检验使用基于子抽样的自归一化统计量（Dette & Schumann, 2024, 定理 2）。当 beta_RMS^2 < zeta^2 + Q_W(alpha) * V_n 时拒绝 H0，其中 Q_W(alpha) 是极限分布（布朗运动泛函）的 alpha 分位数，V_n 由子样本估计计算（公式 (4.18)）。该检验是枢轴的，不需要方差估计。

## RMS 检验的 alpha 限制

RMS 检验仅支持：

$$
\alpha \in \{0.01, 0.025, 0.05, 0.1, 0.2\}
$$

这反映了基于 Dette & Schumann (2024, 定理 2-3) 中极限分布临界值的实现。

## 存储结果

`equivtest` 将结果存储在 `e()` 中（另见 `src/sthlp/equivtest.sthlp`）。

### 标量

| 结果                       | 说明                                             |
| :------------------------- | :----------------------------------------------- |
| `e(N)`                   | 观测数                                           |
| `e(N_g)`                 | 个体（面板）数                                   |
| `e(no_placebos)`         | 安慰剂系数数量（*T*）                          |
| `e(alpha)`               | 使用的显著性水平                                 |
| `e(base_period)`         | 安慰剂构造的基期                                 |
| `e(is_balanced)`         | 平衡面板为 1，否则为 0                           |
| `e(threshold_specified)` | 指定了 `threshold()` 为 1，否则为 0            |
| `e(threshold)`           | 等价性阈值（如已指定）                           |
| `e(reject)`              | 拒绝 H₀ 为 1，否则为 0（如已指定阈值）          |
| `e(min_threshold)`       | 等价性成立的最小阈值 δ*/τ*/ζ*（如未指定阈值） |

`type(max)` 特有标量：

| 结果                 | 说明                                               |
| :------------------- | :------------------------------------------------- |
| `e(max_abs_coef)`  | 最大绝对安慰剂系数                                 |
| `e(nboot)`         | Bootstrap 重复次数（`method(boot)` 或 `wild`） |
| `e(boot_critical)` | Bootstrap 临界值（`method(boot)` 或 `wild`）   |

`type(mean)` 特有标量：

| 结果                       | 说明                           |
| :------------------------- | :----------------------------- |
| `e(abs_mean_placebo)`    | 安慰剂系数均值的绝对值         |
| `e(var_mean_placebo)`    | 安慰剂系数均值的方差           |
| `e(se_mean_placebo)`     | 安慰剂系数均值的标准误         |
| `e(p_value)`             | p 值（如已指定阈值）           |
| `e(mean_critical_value)` | 均值检验临界值（如已指定阈值） |

`type(rms)` 特有标量：

| 结果                      | 说明                       |
| :------------------------ | :------------------------- |
| `e(rms_placebo_coefs)`  | 安慰剂系数的均方根         |
| `e(nolambda)`           | Lambda 子样本数量          |
| `e(rms_critical_value)` | RMS 临界值（如已指定阈值） |

### 宏

| 结果                      | 说明                                                       |
| :------------------------ | :--------------------------------------------------------- |
| `e(cmd)`                | 命令名称（`equivtest`）                                  |
| `e(cmdline)`            | 完整的输入命令                                             |
| `e(type)`               | 检验类型：`max`、`mean` 或 `rms`                     |
| `e(method)`             | 推断方法：`iu`、`boot` 或 `wild`（仅 `type(max)`） |
| `e(depvar)`             | 因变量名称                                                 |
| `e(idvar)`              | 面板标识变量                                               |
| `e(groupvar)`           | 处理组变量                                                 |
| `e(timevar)`            | 时间变量                                                   |
| `e(vce)`                | 方差估计量类型（适用时）                                   |
| `e(clustvar)`           | 聚类变量（如使用聚类稳健 VCE）                             |
| `e(preperiods)`         | 指定的处理前时期（如使用了 `pretreatment()`）            |
| `e(placebo_coef_names)` | 安慰剂系数时期名称                                         |

### 矩阵

| 结果                        | 说明                                                               |
| :-------------------------- | :----------------------------------------------------------------- |
| `e(b_placebo)`            | 安慰剂系数向量（*T* × 1）；适用于 IU 最大值和均值检验         |
| `e(V_placebo)`            | 安慰剂方差-协方差矩阵（*T* × *T*）；适用于 IU 最大值/均值检验 |
| `e(se_placebo)`           | 安慰剂标准误（*T* × 1）；适用于 IU 最大值检验                   |
| `e(IU_critical_values)`   | 各系数的临界值（type=max, method=iu, 已指定阈值）                  |
| `e(min_equiv_thresholds)` | 各系数的最小阈值（type=max, method=iu, 未指定阈值）                |

## 可视化

`equivtest_plot` 创建带等价性界限的安慰剂系数图。在运行 `equivtest` 后使用它来可视化结果。

### 基本语法

```stata
equivtest_plot [, options]
```

### 主要选项

| 选项             | 说明                                                                |
| :--------------- | :------------------------------------------------------------------ |
| `threshold(#)` | 水平线的等价性阈值；默认为 `e(threshold)` 或 `e(min_threshold)` |
| `ci`           | 显示置信区间（仅 IU 方法）                                          |
| `level(#)`     | 置信水平；默认 95                                                   |
| `connect`      | 用线连接各点                                                        |
| `noline`       | 抑制零参考线                                                        |
| `nobase`       | 抑制基期参考点                                                      |
| `nothreshold`  | 抑制等价性阈值线                                                    |

### 样式选项

| 选项                                 | 说明                            |
| :----------------------------------- | :------------------------------ |
| `msymbol(symbolstyle)`             | 标记符号；默认 `O`（圆形）    |
| `msize(markersizestyle)`           | 标记大小；默认 `medium`       |
| `mcolor(colorstyle)`               | 标记颜色；默认 `navy`         |
| `basemsymbol(symbolstyle)`         | 基期标记符号；默认 `S`        |
| `basemcolor(colorstyle)`           | 基期标记颜色                    |
| `threshlcolor(colorstyle)`         | 阈值线颜色；默认 `red`        |
| `threshlwidth(linewidthstyle)`     | 阈值线宽度；默认 `medium`     |
| `threshlpattern(linepatternstyle)` | 阈值线样式；默认 `dash`       |
| `cilcolor(colorstyle)`             | 置信区间线颜色；默认 `navy`   |
| `cilwidth(linewidthstyle)`         | 置信区间线宽度；默认 `medium` |

### 输出选项

| 选项                       | 说明               |
| :------------------------- | :----------------- |
| `title(string)`          | 图形标题           |
| `subtitle(string)`       | 图形副标题         |
| `xtitle(string)`         | X 轴标题           |
| `ytitle(string)`         | Y 轴标题           |
| `xlabel(rule_or_values)` | X 轴标签           |
| `ylabel(rule_or_values)` | Y 轴标签           |
| `note(string)`           | 图形注释           |
| `scheme(schemename)`     | 图形方案           |
| `saving(filename)`       | 保存图形到文件     |
| `replace`                | 保存时替换已有文件 |
| `name(windowname)`       | 图形窗口名称       |

### 示例

```stata
* 运行 equivtest 后的基本绘图
equivtest_plot

* 带 95% 置信区间的绘图
equivtest_plot, ci

* 出版质量的绘图
equivtest_plot, ci level(95) title("预趋势分析") ///
    msymbol(O) mcolor(navy) threshlpattern(dash) scheme(s2color)

* 保存图形到文件
equivtest_plot, ci saving(pretrend_plot) replace
```

## 示例

### 模拟面板数据

本示例使用满足平行趋势的模拟数据演示软件包的使用。

```stata
* 生成模拟面板数据
clear
set seed 12345
set obs 1000
gen id = ceil(_n/10)
gen time = mod(_n-1, 10) + 1
gen treat = (id <= 50)

* 生成结果变量：处理前（time <= 5）满足平行趋势
* 处理后（time > 5）处理效应为 0.5
gen y = rnormal() + 0.1*time + treat*(time > 5)*0.5

* 使用 IU 方法运行最大值检验
equivtest y, type(max) method(iu) id(id) group(treat) time(time) ///
    pretreatment(1 2 3 4) baseperiod(5)

* 显示最小阈值
display "最小等价性阈值 (delta*): " %6.4f e(min_threshold)
```

### 解读结果

当运行 `equivtest` 而未指定 `threshold()` 时，命令报告在指定显著性水平下能够拒绝非等价零假设的最小等价性阈值：

```
Equivalence Test for Pre-Trends (Maximum)
──────────────────────────────────────────────────────────────────────────────
Type:             max                    Observations:       1000
Method:           iu                     Individuals:         100
VCE:              ols                    Placebo coefs:         4
──────────────────────────────────────────────────────────────────────────────

Hypothesis Test:
  H0: max|placebo effect| >= delta  (non-equivalence)
  H1: max|placebo effect| <  delta  (equivalence)

Minimum Equivalence Threshold Search
──────────────────────────────────────────────────────────────────────────────
        Period |  Abs. Estimate   Std. Error   Min. Threshold
───────────────+──────────────────────────────────────────────────────────────
    placebo_1  |       0.034521     0.100234       0.199543
    placebo_2  |       0.012345     0.098765       0.174832
    placebo_3  |       0.056789     0.101234       0.223456
    placebo_4  |       0.023456     0.099876       0.187654
──────────────────────────────────────────────────────────────────────────────

Minimum equivalence threshold delta* =   0.2235
```

**解读**：在 alpha = 0.05 下，我们可以得出最大绝对安慰剂系数小于 delta* 的结论。如果该阈值相对于您估计的处理效应较小，则您有证据支持预趋势违反可忽略。

### 使用预设阈值进行检验

```stata
* 检验最大违反是否 < 0.2
equivtest y, type(max) method(iu) id(id) group(treat) time(time) ///
    pretreatment(1 2 3 4) baseperiod(5) threshold(0.2)

* 检查拒绝决策
if e(reject) == 1 {
    display "等价结论成立：max|beta_t| < 0.2"
}
else {
    display "无法在阈值 0.2 下得出等价结论"
}
```

## 引用

如果您在研究中使用了本软件包，请同时引用方法论文和 Stata 实现：

**APA 格式：**

> Dette, H., & Schumann, M. (2024). Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation. *Journal of Business & Economic Statistics*, 42(4), 1289–1301. https://doi.org/10.1080/07350015.2024.2308121
>
> Cai, X., & Xu, W. (2025). *Equitrends: Stata module for equivalence tests for pre-trends in DiD* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/equitrends

**BibTeX：**

```bibtex
@article{dette2024testing,
  title={Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation},
  author={Dette, Holger and Schumann, Martin},
  journal={Journal of Business \& Economic Statistics},
  volume={42},
  number={4},
  pages={1289--1301},
  year={2024},
  publisher={Taylor \& Francis},
  doi={10.1080/07350015.2024.2308121}
}

@software{equitrends2025stata,
  title={Equitrends: Stata module for equivalence tests for pre-trends in DiD},
  author={Cai, Xuanyu and Xu, Wenli},
  year={2025},
  version={0.1.0},
  url={https://github.com/gorgeousfish/equitrends}
}
```

## 作者

**Stata 实现：**

- **蔡宣宇 (Xuanyu Cai)**，澳门城市大学
  邮箱：[xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **许文理 (Wenli Xu)**，澳门城市大学
  邮箱：[wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**方法论：**

- **Holger Dette**，波鸿鲁尔大学数学系
- **Martin Schumann**，马斯特里赫特大学商业与经济学院

## 许可证

AGPL-3.0 许可证。详见 [LICENSE](LICENSE)。

## 未来路线图

开发团队正在评估以下扩展方向：

- **交错处理支持**：通过队列堆叠将等价性检验框架扩展到交错处理设计（Dette & Schumann, 2024, 第 5 节）。
- **阈值敏感性分析**：可视化检验决策和最小阈值如何随等价性阈值的连续范围变化。
- **额外的 Bootstrap 方法**：实现乘数 Bootstrap 和子抽样替代方案。
- **Bootstrap 并行化**：通过 `parallel` 包实现多核 Bootstrap 计算，以减少大数据集的运行时间。

## 相关软件包

- **EquiTrends (R)**：[TiesBos/EquiTrends](https://github.com/TiesBos/EquiTrends) — 预趋势等价性检验的原始 R 实现（Dette & Schumann, 2024）
- **pretest (Stata)**：[gorgeousfish/pretest](https://github.com/gorgeousfish/pretest) — DiD 的条件外推预检验与偏差调整置信区间（Mikhaeil & Harshaw, 2025）

## 推荐资源

因果推断和计量经济学入门推荐：

- [Causal Inference for the Brave and True](https://matheusfacure.github.io/python-causality-handbook/landing-page.html) — Matheus Facure 撰写的优秀因果推断入门教程
- [Causal Inference for the Brave and True（中文版）](https://ci-book.huangwz.com/intro) — 黄文喆、许文立 中文翻译
- [What&#39;s Trending in Difference-in-Differences? A Synthesis of the Recent Econometrics Literature](https://doi.org/10.1016/j.jeconom.2023.03.008) — Roth, Sant'Anna, Bilinski & Poe (2023) 的综合综述

## 支持

如有问题或发现 Bug，请在 [GitHub](https://github.com/gorgeousfish/equitrends/issues) 上提交 Issue。
