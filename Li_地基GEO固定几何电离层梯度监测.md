# 地基 GEO 固定几何的电离层梯度监测

> **Li, Zhiyao; Zhong, Jiahao; Wang, Ningbo; et al.**  
> *Ionospheric gradient estimation using ground-based GEO observations for monitoring multi-scale ionospheric dynamics*  
> **Satellite Navigation**, 2025, 6:30 · DOI: [10.1186/s43020-025-00187-4](https://doi.org/10.1186/s43020-025-00187-4)

---

## 1. 研究背景

### 1.1 现有电离层观测的矛盾

| 观测方式 | 优点 | 局限 |
|---|---|---|
| MEO/LEO 卫星、GNSS 掩星 | 覆盖广，适合全球或大尺度结构 | 对某一地点的重复时间有限；IPP 持续移动 |
| 地基雷达、气辉成像 | 时间连续性好 | 区域覆盖有限 |
| GIM 等格网产品 | 全球、易获取 | 插值与平滑会削弱瞬态和细尺度结构 |

电离层扰动同时跨越秒到小时、几十到上千公里的尺度。传统方案往往只能在“空间覆盖、时间连续性、空间分辨率”之间取舍。

### 1.2 本文的核心问题

将 GEO–接收机形成的稳定视线几何充分利用起来，直接、连续地监测 TEC 的空间与时间变化，而不依赖卫星运动校正或大范围插值

### 1.3 最核心的物理与几何洞见

对于同一颗 GEO 卫星：

- 每个地面站对应一个近乎固定的 IPP；
- 两个 IPP 之间的基线长度 $S_{ij}$ 和方位角 $	\theta_{ij}$ 随时间近乎不变；
- 同一历元的两点 VTEC 之差除以基线长度，就是沿该方向的空间梯度；
- 单点相邻历元的 TEC 差除以时间间隔，就是时间梯度。

因此，**固定的 IPP 对不是“两个离散点”，而是一台具有固定尺度与方向的差分传感器。** 若有 $N$ 个 IPP，可形成最多 $N(N-1)/2$ 条 IPP 对基线，观测冗余度和空间采样密度都会显著提高。

---

## 2. 数据准备

| 项目 | 论文设置 |
|---|---|
| 接收机网络 | 207 个地基 GNSS 接收机，主要位于东亚中低纬；平均站距约 200 km |
| 卫星 | 5 颗北斗 GEO：C01、C02、C03、C04、C05；星下点经度覆盖 $59.0^\circ\mathrm{E}$–$160.0^\circ\mathrm{E}$ |
| TEC 数据 | 北斗 B1I–B3I 双频载波相位；载波相位观测码为 L2I、L6I |
| 质量控制 | 卫星高度角 $>20^\circ$ |
| 时间分辨率 | 30 s |
| 电离层薄壳 | 固定高度 $h=400\ \mathrm{km}$ |
| DCB | 使用中国科学院（CAS）日 DCB 产品改正卫星和接收机 DCB |

---

## 3. 核心公式与梯度参量

### 3.1 双频载波相位与斜向 TEC（STEC）

$$
T_S=
\frac{f_1^2 f_2^2}{40.3\left(f_1^2-f_2^2\right)}
\left[
\left(L_1-L_2\right)-\left(N_1-N_2\right)+\left(d_1-d_2\right)
\right].
$$

| 符号 | 含义 |
|---|---|
| $T_S$ | 沿卫星—接收机视线的斜向总电子含量（STEC） |
| $f_1,f_2$ | 两个载波频率；本文为 B1I 与 B3I |
| $L_1,L_2$ | 两个频点的载波相位观测值 |
| $N_1,N_2$ | 载波相位整周模糊度 |
| $d_1,d_2$ | 两频点相关硬件偏差 |

载波相位噪声小，但含未知的整周模糊度；作者随后用载波—伪距平滑（CCL）确定该常数，并用 DCB 产品给出绝对 TEC 标尺。

### 3.2 单层映射：由 STEC 得到 VTEC

$$
M(\delta)=
\frac{1}{\sqrt{1-\left(\frac{R_e\cos\delta}{R_e+h}\right)^2}},
$$

$$
T_S^i=M(\delta_i)\,T_V^i+b_{\mathrm{R,DCB}}+b_{\mathrm{S,DCB}}+\varepsilon,
$$

$$
T_V^i=
\frac{T_S^i-b_{\mathrm{R,DCB}}-b_{\mathrm{S,DCB}}}{M(\delta_i)}+\xi.
$$

| 符号 | 含义 |
|---|---|
| $T_V^i$、$T_S^i$ | 第 $i$ 个 IPP 的垂直 / 斜向 TEC |
| $M(\delta_i)$ | 与高度角 $delta_i$ 有关的映射函数 |
| $R_e$、$h$ | 地球半径与薄壳高度；本文 $h=400\ \mathrm{km}$ |
| $\xi$ | 多路径、测量噪声、映射误差等残余误差 |

作者指出北斗 VTEC 的残余误差通常 $|\xi|\leq1\ \mathrm{TECU}​$；本文主要研究梯度，做差分后这类残余影响会进一步减弱。

### 3.3 TEC时间梯度ROT：同一固定 IPP 在相邻历元的变化

$$
\frac{\partial T_V^i}{\partial t}
=\frac{1}{M(\delta_i)}\frac{\Delta T_S^i}{\Delta t}
=\frac{1}{M(\delta_i)}I_{\mathrm{ROT}}^i,
$$

$$
I_{\mathrm{ROT}}^i=
\frac{T_S^i(t_2)-T_S^i(t_1)}{t_2-t_1},
\qquad \Delta t=30\ \mathrm{s}.
$$

前后两个历元的差主要反映真实时间演化，而不是卫星运动导致的采样位置变化。它可显示 EIA 的增长/衰减、EPB 的快速波动或 TID 波前的通过。

### 3.4 TEC波动率ROTI：识别电离层不规则体

$$
I_{\mathrm{ROTI}}^i=
\sqrt{
\left\langle\left(I_{\mathrm{ROT}}^i\right)^2\right\rangle
-\left\langle I_{\mathrm{ROT}}^i\right\rangle^2
}.
$$

尖括号表示 5 min 时间窗内的平均。ROT 是 TEC 的变化率；ROTI 是 ROT 在短时间窗内的标准差。ROTI 大说明 TEC 变化很不规则，常用于标记 EPB 等小尺度扰动。它能说明“扰动强不强”，但不能直接给出**扰动的空间方向与边界**；这正是空间梯度场的补充价值。

### 3.5 空间梯度：体现扰动的空间方向

$$
\left.\frac{\partial T_V^{(i,j)}}{\partial S_{ij}}\right|_t
=
\left(
\frac{T_S^i}{M(\delta_i)}-
\frac{T_S^j}{M(\delta_j)}
\right)\frac{1}{S_{ij}}
+\frac{\Delta\xi}{S_{ij}}.
$$

$$
\nabla T_V^{(i,j)}(t)=
\left|
\left.\frac{\partial T_V^{(i,j)}}{\partial S_{ij}}\right|_t
\right|.
$$

| 符号 | 含义 |
|---|---|
| $S_{ij}$ | IPP $i$ 与 $j$ 在 400 km 薄壳上的大圆距离 |
| $\theta_{ij}$ | 从正北顺时针量起的 IPP 基线方位角 |
| $\Delta\xi$ | 两个 IPP 的残余误差之差 |
| $\nabla T_V^{(i,j)}(t)$ | 沿该向量方向的空间TEC梯度**大小** |

同一时刻两个固定 IPP 的 VTEC 差除以它们的距离，即为沿基线方向的空间变化率。残余误差项被 $S_{ij}$ 除，以较长基线会降低这部分影响；但过长基线又会平滑掉小尺度结构，实际应用需权衡基线长度与目标尺度。

 **将梯度分解为东西、南北分量**
$$
\left.\frac{\partial T_V^{(i,j)}}{\partial x}\right|_t
=\nabla T_V^{(i,j)}(t)\sin\theta_{ij},
\qquad
\left.\frac{\partial T_V^{(i,j)}}{\partial y}\right|_t
=\nabla T_V^{(i,j)}(t)\cos\theta_{ij}.
$$

其中 $x$ 为向东、$y$ 为向北。

梯度大小只说明大小，分量则说明方向。例如：EPB 的 TEC 空洞两侧常出现相反符号的东西向梯度；EIA 峰两侧则常呈明显的南北向梯度。这是仅用 ROTI 无法提供的信息。

### 3.6 混合时空导数：检查EPB引起的线性演化与非线性扰动

在 GEO 几何固定、时间间隔很小的条件下，作者构造：

$$
\frac{\partial^2 T_V^{(i,j)}}{\partial x\,\partial t}
\approx
\frac{
\left|
\frac{\partial T_V^i}{\partial t}
-\frac{\partial T_V^j}{\partial t}
\right|
\sin\theta_{ij}}
{S_{ij}}
\approx
\frac{\partial^2 T_V^{(i,j)}}{\partial t\,\partial x},
$$

$$
\frac{\partial^2 T_V^{(i,j)}}{\partial y\,\partial t}
\approx
\frac{
\left|
\frac{\partial T_V^i}{\partial t}
-\frac{\partial T_V^j}{\partial t}
\right|
\cos\theta_{ij}}
{S_{ij}}
\approx
\frac{\partial^2 T_V^{(i,j)}}{\partial t\,\partial y}.
$$

**直观理解**：不同位置的 TEC 变化速度相差多少，即时间演化在空间上的不均匀性。若 TEC 场在时空上平滑演化，混合偏导可以交换次序，两种计算应近似一致；文中将它们的差作为非线性或不连续过程的指标。

**物理解释**：平静背景或 EIA 早期的平滑发展更接近线性、可交换；EPB 的产生、漂移和耗散会造成强烈的时空耦合，混合导数不对称性增大，会破坏其本身的连续性。因此该量是一个**时空连续性/非线性诊断量**。

### 3.7 格网内标准差：分辨率是否足以解析结构

$$
\sigma\!\left(\nabla T_V\right)=
\sqrt{
\left\langle\left(\nabla T_V\right)^2\right\rangle
-\left\langle\nabla T_V\right\rangle^2
}.
$$

多条不同方向、不同长度的 IPP 间基线，都在该格网附近对 TEC 空间变化作出了估计，对同一格网内这些梯度样本求平均，得到该格网的空间梯度场。若同一格网内多个梯度样本的离散度很大，说明格网内仍有未被完全分辨的细结构；该统计量帮助判断是否足以描述局地变化，而不是把格网均值误认为真实均匀场。**如果主导尺度比较大这个标准差应该表现出0值附近，如果尺度比较小有很多细节则会表现出较大值**

---

## 4. 应用实例

### 案例 I：EIA 的日变化（2024-09-21，平静日，13:40）

**研究对象**：赤道电离层异常（EIA）。赤道上空等离子体上升后沿磁力线向两侧扩散，形成赤道附近低谷与约 $15^\circ$–$20^\circ$ 磁纬两侧峰区。

| 观察量 | 作者观察到的现象 | 具体表现 |
|---|---|---|
| VTEC 与空间梯度 | 在约 $25^\circ\mathrm{N}$ 的 EIA 峰两侧出现显著梯度带 | 梯度直接给出峰区两侧的陡变位置 |
| 南北梯度 | 显示峰区向高、低纬的扩散 | 可刻画喷泉效应相关的经向运输 |
| 东西梯度 | 显示电离分布的经向不对称 | 可观察局地电动力背景差异 |
| 时间梯度与混合导数 | 11:40–14:40 LT 的早期增长更平滑 | 早期 EIA 主要受输运/局地产生等近线性过程控制，后期就乱起来了 |
| 晚间 ROTI 与导数不对称 | 19:40–01:40 LT 出现 EPB 相关扰动，南侧梯度结构被侵蚀 | 判断何时被EPB破坏以及EPB的移动方向 |

该案例说明框架不仅能指出EIA 峰是位置，还能显示峰区如何形成、向何处扩散、何时被 EPB 破坏（会造成TEC空洞）。

### 案例 II：EPB 的生命周期（2024-09-21 22:40）

**研究对象**：赤道等离子体泡（EPB），即低纬电离层中电子密度亏损的不规则结构，会损害通信与导航信号。

| 关键发现 | 具体表现 |
|---|---|
| 出现与漂移 | EPB 约在 20:00 LT、$110^\circ\mathrm{E}$ 附近出现；至 23:40 LT 漂移到约 $115^\circ\mathrm{E}$，总体向东漂移 |
| 形态 | TEC 空洞范围约 $105^\circ$–$115^\circ\mathrm{E}$，向北可延伸至 $25^\circ\mathrm{N}$，呈西倾形态 |
| 边界定位 | 空洞两侧出现正、负相反的东西向梯度极值；梯度符号变化可定位空洞槽底与两侧泡壁 |
| 强度与演化 | 梯度强度约在 22:40 LT 最强，之后减弱；倾斜角在增长与衰减阶段呈不同转动趋势 |
| 非线性特征 | ROTI 升高且混合导数不对称显著增加 |

ROTI 可以提示出现了不规则体，而固定几何梯度场进一步说明不规则体的边界在哪里、朝哪边漂移、形态怎样变化。

### 案例 III：磁暴期 LSTID 的传播（2024-10-11）

**研究对象**：大尺度行进电离层扰动（LSTID）。事件日 Dst = $-308\ \mathrm{nT}$，属于强磁暴背景。

| 结果 | 含义 |
|---|---|
| 时间梯度中出现交替正负带 | 对应向南传播的 LSTID 波前 |
| 主要表现于南北梯度，东西分量较弱 | 与 LSTID 以子午向传播为主一致 |
| 混合导数描出波相位转换区 | $\partial T_V/\partial t$ 给出当前相位状态；$\partial T_V/\partial y$ 给出空间斜率；$\partial^2T_V/\partial y\partial t$ 给出相位演化在空间上的不均匀性 |

同一框架能够同时描述 EIA 的低纬背景、EPB 的细尺度不规则体和 LSTID 的大尺度波动，说明其具有跨尺度能力。

---

## 5. 创新点

1. **固定几何的传感器化思维**：不是把 GEO 当作普通 GNSS 卫星使用，而是把稳定的 IPP–IPP 基线定义为具有固定长度和方向的梯度传感单元。
2. **直接的空间差分，避免运动伪影**：GEO IPP 基本固定，因此空间梯度不需要对卫星运动、IPP 漂移或高度角快速变化作额外校正。
3. **空间、时间、混合导数一体化**：除 VTEC 外，同时输出空间梯度、时间梯度、ROTI 和混合导数不对称性，分别描述结构边界、演化速度、短时扰动和非线性。
4. **多 GEO 视角与大量 IPP 对带来的密采样**：五颗 GEO 卫星和 207 个接收机形成高冗余、连续的观测架构，而非少数 IPP 的个案比较。
5. **过程描述**：不仅发现 EIA/EPB/LSTID，还尝试刻画 EIA 的输运、EPB 的漂移与形变、LSTID 的相位和传播。

---

## 6. 我提出的建议

画图应该放一个正常情况下，无特殊情况的对比图，这样看着可能更清晰一些。