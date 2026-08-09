# GEO 观测电离层梯度监测

## 1. 主题：电离梯度监测

**如何利用地基 GNSS 观测，及时识别电离层中影响导航定位的快速变化、空间梯度与不规则结构。**

```math
\text{ROTI：识别不规则体}
\;\longrightarrow\;
\text{DIX：量化时空扰动}
\;\longrightarrow\;
\text{GIX：直接估计空间 TEC 梯度}.
```
对应文献：

1. Pi et al. (1997)：提出 ROT 与 ROTI，用于监测电离层不规则体；
2. Jakowski et al. (2012)：提出 DIX，用多链路 TEC 变化率构建时空扰动指数；
3. Jakowski & Hoque (2019)：提出 GIX，直接以 IPP 对之间的 VTEC 差分估计空间梯度。

---

## 2. 基本概念

对传统 MEO GNSS 卫星，电离层穿刺点（IPP）会随卫星运动。在一个移动 IPP 上观测到的 TEC 变化率为：

```math
\frac{d\mathrm{TEC}}{dt}
=
\frac{\partial\mathrm{TEC}}{\partial s}v_{\mathrm{IPP}}
+
\frac{\partial\mathrm{TEC}}{\partial t}.
```
其中：

- 第一项表示 IPP 扫过空间 TEC 梯度产生的变化；
- 第二项表示固定位置上的 TEC 真正随时间发生的变化；
- $v_{\mathrm{IPP}}$ 为 IPP 在电离层薄壳上的运动速度。

因此，传统单站相邻历元 TEC 差分通常是**空间变化和时间变化的混合量**。三篇文章分别采用 ROTI、DIX 和 GIX，从不同角度应对这一问题。

---

## 3. Pi et al. (1997)：ROT 与 ROTI

> Pi, X., Mannucci, A. J., Lindqwister, U. J., & Ho, C. M. (1997). *Monitoring of global ionospheric irregularities using the Worldwide GPS Network*. Geophysical Research Letters, 24(18), 2283–2286.

### 3.1 文章核心

本文首次系统展示：可以利用全球 GPS 接收机网，依据双频载波相位的快速起伏监测电离层不规则体。其关注对象主要是低纬 EPB、极区不规则体及磁暴期间的相位闪烁。

它回答的不是“哪里有多大的二维梯度”，而是：

> 某条卫星—接收机视线附近，是否出现了快速相位起伏和小尺度电离层不规则体？

### 3.2 ROT：TEC 的变化率

```math
\mathrm{ROT}(t)
=
\frac{\mathrm{STEC}(t_2)-\mathrm{STEC}(t_1)}{t_2-t_1}.
```
ROT 表示沿卫星—接收机视线的 TEC 变化速率，常用单位为 $\mathrm{TECU/min}$。

### 3.3 ROTI：ROT 的短时标准差

在一个短时间窗口内（经典做法为 5 分钟），定义：

```math
\mathrm{ROTI}
=
\sqrt{
\left\langle \mathrm{ROT}^{2}\right\rangle
-
\left\langle \mathrm{ROT}\right\rangle^{2}
}.
```
其中 $\langle\cdot\rangle$ 表示时间窗口内的平均。

### 3.4 物理含义

- ROTI 小：相位变化平稳，LOS 附近没有显著的小尺度不规则体；
- ROTI 大：TEC 变化率强烈起伏，常对应 EPB、闪烁、极区不规则体或磁暴扰动；
- ROTI 是**起伏强度指标**，不保留梯度的正负号和方向。

在“冻结不规则体”近似下，有：

```math
\mathrm{ROT}
\approx
\frac{\partial\mathrm{TEC}}{\partial s}v.
```
所以若已知不规则体相对视线的扫过速度 $v$，ROT 可与沿该方向的空间梯度联系起来；但实际中 $v$ 往往未知。因此，ROTI 不能直接替代严格的空间梯度估计。

---

## 4. Jakowski et al. (2012)：DIX 时空扰动指数

> Jakowski, N., Borries, C., & Wilken, V. (2012). *Introducing a disturbance ionosphere index*. Radio Science, 47, RS0L14.

### 4.1 文章核心

DIX（Disturbance Ionosphere Index）旨在用一个可近实时更新的数值，描述某个区域电离层的**时空扰动强度**，为导航、通信和空间天气服务提供预警信息。DIX **不是直接的空间梯度**，也不是某一个固定地点的纯时间梯度。它利用许多运动 GNSS 链路的 TEC 变化率，通过统计组合来判断区域内是否存在强烈扰动。

### 4.2 基础观测量

对双频载波相位几何无关组合 $\Delta\Phi$，可得到近似的垂直 TEC 变化率：

```math
\frac{\Delta\mathrm{TEC}}{\Delta t}
\approx
\frac{1}{aM}
\frac{\Delta\Phi}{\Delta t},
```
其中 $a$ 为频率相关系数，$M$ 为映射函数。将单条链路的观测变化率简记为：

```math
r_i
=
\frac{\Delta\mathrm{TEC}_i}{\Delta t}.
```
由于 MEO 卫星在运动，IPP 也在电离层薄壳上移动。因此 $r_i$ 同时包含空间项和时间项：

```math
r_i
=
\frac{\partial\mathrm{TEC}}{\partial s}v_{\mathrm{IPP},i}
+
\frac{\partial\mathrm{TEC}}{\partial t}.
```
其中第一项表示 IPP 移动时扫过空间 TEC 坡度产生的变化，第二项表示固定位置上的 TEC 真正随时间发生的变化。因此，单条链路的 TEC 变化率是**时空混合量**。

### 4.3 两条链路的作差：突出空间差异

在同一区域、同一时刻，选取两条链路 $k$ 与 $l$：

```math
r_k=\frac{\Delta\mathrm{TEC}_k}{\Delta t},
\qquad
r_l=\frac{\Delta\mathrm{TEC}_l}{\Delta t}.
```
计算两者之差：

```math
r_k-r_l.
```
如果发生太阳耀斑等大范围同步事件，两处的时间变化通常近似相同，公共时间项会在作差时被削弱，**其差值较小，两处变化很接近，空间差异不明显**。因此，$r_k-r_l$ 更敏感于两条链路经历的空间结构差异。空间 DIX 分量可概念化为：

```math
\mathrm{DIX}_{S}
\propto
\sqrt{
\frac{1}{N_P}
\sum_{k,l}
\left(r_k-r_l\right)^2
}.
```
用于描述整个空间的差异情况。

### 4.4 两条链路的平均：突出共同变化

计算两条链路变化率的平均：

```math
\frac{r_k+r_l}{2}.
```
若发生大范围 TEC 同步增强，两条链路通常同号且数值接近，因此该平均值会明显增大。若变化主要来自 IPP 扫过局地空间结构，不同链路的空间项可能正负不同；在大量链路平均后，这部分会部分抵消。故该平均量更偏向区域共同的时间变化。

### 4.5 总 DIX：时空扰动的综合指标

完整 DIX 将时间项与空间项结合：

```math
\mathrm{DIX}
\propto
\sqrt{
\frac{1}{N_P}
\sum_{k,l}
\left[
\left(\frac{r_k+r_l}{2}\right)^2
+
h\left(r_k-r_l\right)^2
\right]
},
```
其中 $h$ 为控制空间项权重的调节参数；$h$ 越大，指标越重视空间结构扰动。第一项偏向区域共同时间变化，第二项偏向区域空间结构差异。DIX值大只能说明此区域发生了强烈的TEC时空扰动，但是并不能说明方向问题。

| 情形 | $r_k$ 与 $r_l$ 的关系 | $r_k-r_l$ | $(r_k+r_l)/2$ | DIX 的主要响应 |
|:---|:---|:---|:---|:---|
| 大范围同步 TEC 增强 | 同号且接近 | 小 | 大 | 时间项增强 |
| 两条链路扫过强 TEC 锋面 | 差异大，甚至异号 | 大 | 小或中等 | 空间项增强 |

例如，若 $r_k=0.8\ \mathrm{TECU/min}$、$r_l=0.7\ \mathrm{TECU/min}$，两者差为 $0.1\ \mathrm{TECU/min}$，而平均为 $0.75\ \mathrm{TECU/min}$，更像区域同步时间变化。若 $r_k=0.8\ \mathrm{TECU/min}$、$r_l=-0.6\ \mathrm{TECU/min}$，差值为 $1.4\ \mathrm{TECU/min}$、平均仅为 $0.1\ \mathrm{TECU/min}$，则更像强空间结构或传播扰动。

### 4.6 DIX 的尺度与局限

DIX 通过设定 IPP 对最大距离 $L_c$ 来控制关注的空间尺度：

```math
L_c\ \text{较小}
\Rightarrow
\text{更关注局地、小尺度扰动},
```
```math
L_c\ \text{较大}
\Rightarrow
\text{更关注中、大尺度扰动}.
```
DIX 的优势是无需绝对 TEC 标定、实时性强、对 DCB 残差更稳健；局限是依赖多条运动链路的统计组合，不能直接给出局地梯度的准确位置、方向和幅度。

---

## 5. Jakowski & Hoque (2019)：GIX 直接空间梯度估计

> Jakowski, N., & Hoque, M. M. (2019). *Estimation of spatial gradients and temporal variations of the Total Electron Content using ground-based GNSS measurements*. Space Weather, 17(2), 339–356.

### 5.1 文章核心

利用同一时刻两个 IPP 的 VTEC 差除以其距离，直接估计沿基线方向的 TEC 空间梯度。**把一个区域内很多条 IPP 基线梯度，汇总成一个区域整体梯度强度的指标。**

### 5.2 单条 IPP 基线的梯度

```math
\nabla\mathrm{TEC}_{ij}
=
\frac{\mathrm{VTEC}_i-\mathrm{VTEC}_j}{\Delta s_{ij}}.
```
其中 $\Delta s_{ij}$ 为两个 IPP 在电离层薄壳上的距离。

### 5.3 梯度方向分解

若 IPP 对方位角为 $\delta$，则可分解为东西、南北分量：

```math
\nabla\mathrm{TEC}_{x,ij}
=
\nabla\mathrm{TEC}_{ij}\sin\delta,
```
```math
\nabla\mathrm{TEC}_{y,ij}
=
\nabla\mathrm{TEC}_{ij}\cos\delta.
```
方向信息很关键：

- EIA 峰两侧通常具有显著南北向梯度；
- EPB 空洞两侧通常具有强烈、符号相反的梯度；
- LSTID 则常呈现连续推进的正负梯度带。

### 5.4 区域梯度指标

对区域内多个 IPP 对进行统计，得到：

**平均梯度 GIX：**

```math
\mathrm{GIX}
=
\sqrt{
\left\langle\nabla\mathrm{TEC}_{x}\right\rangle^2
+
\left\langle\nabla\mathrm{TEC}_{y}\right\rangle^2
}.
```
**梯度离散度 GIXS：**

```math
\mathrm{GIXS}
=
\sqrt{
\left\langle(\nabla\mathrm{TEC})^2\right\rangle
-
\left\langle\nabla\mathrm{TEC}\right\rangle^2
}.
```
**极端梯度指标 GIXP95：**

```math
\mathrm{GIXP}_{95}
=
P_{95}
\left(
\left|\nabla\mathrm{TEC}_{ij}\right|
\right).
```
GIX 描述整体平均梯度；GIXS 描述区域内梯度是否高度不均匀；GIXP95 描述接近极端的梯度水平，适合完整性和安全风险评估。

**为何引入GIXS和GIXP95：**

空洞左边，TEC 可能向东减小：$g_x<0$,空洞右边，TEC 又向东增加：$g_x>0$,如果把两侧梯度平均，可能发生抵消使得$g_x$平均后约等于0，此时即使 EPB 边界很陡，区域平均的 GIX 仍可能不大。

因此，EPB 情况下常见的图景是：GIX不一定大，但是另外两个指标可能很大。

### 5.5 方法局限

该方法已经直接估计了空间梯度，但其 IPP 来自 MEO 卫星，仍然会移动：

- 相邻历元的梯度并不严格位于同一基线和同一位置；
- 时间梯度仍需要通过多链路平均来压低空间项影响；
- 绝对 VTEC 依赖 DCB 改正和映射函数，短基线时残余偏差可能被放大；
- 区域平均可能平滑局地极端结构。

---

## 6. 与 GEO 固定几何梯度框架的关系

GEO 卫星相对地面测站近似固定，因此：

```math
\frac{dS_{ij}}{dt}\approx 0,
\qquad
\frac{d\theta_{ij}}{dt}\approx 0,
\qquad
v_{\mathrm{IPP}}\approx 0.
```
这带来三个直接改进：

1. **时间梯度更纯粹**：固定 IPP 的相邻历元差分主要反映局地 TEC 的时间变化；
2. **空间梯度可重复比较**：同一 IPP 对在每个历元对应近似相同的位置、距离与方位；
3. **可进一步讨论混合梯度**：可以比较空间梯度的时间变化与时间梯度的空间变化。

因此，GEO 框架可同时给出：

```math
\frac{\partial T_V}{\partial x},
\qquad
\frac{\partial T_V}{\partial y},
\qquad
\frac{\partial T_V}{\partial t},
\qquad
\frac{\partial^2T_V}{\partial x\partial t},
\qquad
\frac{\partial^2T_V}{\partial y\partial t}.
```
---

## 7. 三类核心量区别表

| 核心量 | 它回答的问题 | 主要适用对象 | 主要局限 |
|:---|:---|:---|:---|
| ROTI | 是否存在快速相位起伏和小尺度不规则体？ | EPB、闪烁、极区不规则体 | 无方向、无正负号、不是严格空间梯度 |
| DIX | 区域整体时空扰动有多强？ | 磁暴、TID、空间天气预警 | 是统计扰动指标，无法精确定位梯度边界 |
| GIX / GIXS / GIXP95 | 区域梯度平均强度、离散性和极端程度如何？ | TEC 锋面、强梯度区、定位完整性 | MEO IPP 运动导致时空混叠 |
| GEO 空间梯度 | 梯度在哪、方向为何、幅度多大？ | EIA、EPB 空洞边界、LSTID 波前 | 依赖 GEO 覆盖、薄壳假设与 VTEC 改正 |

