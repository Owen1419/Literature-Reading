# 中国区域 VTEC 一维剖面 B 样条拟合：文献阅读与研究建议

## 1. B 样条函数

### 1.1 直观理解

B-spline（Basis spline，B 样条）不是一条预先指定的 S 型曲线，而是一组局部、平滑、彼此重叠的“基函数”。拟合曲线是这些基函数的加权和：

$$
f(x)=\sum_{j=1}^{p} c_j B_{j,q}(x).
$$

其中，$x$ 是位置（本研究中可为归一化纬度），$B_{j,q}$ 是第 $j$ 个 $q$ 次 B 样条基函数，$c_j$ 是待估系数。直觉上，某个 $c_j$ 只主要影响其附近的一小段区域；因此站点缺失或异常值不会像全局多项式/球谐那样轻易影响整个中国区域。

### 1.2 节点、次数与局部支撑

给定非递减节点序列

$$
\tau_0\le \tau_1\le\cdots\le\tau_{p+q},
$$

零次基函数为

$$
B_{j,0}(x)=
\begin{cases}
1,&\tau_j\le x<\tau_{j+1},\\
0,&\text{其他}.
\end{cases}
$$

高次基函数由 Cox–de Boor 递推关系给出：

$$
B_{j,q}(x)=
\frac{x-\tau_j}{\tau_{j+q}-\tau_j}B_{j,q-1}(x)
+
\frac{\tau_{j+q+1}-x}{\tau_{j+q+1}-\tau_{j+1}}B_{j+1,q-1}(x),
$$

分母为 0 时对应项定义为 0。工程上最常用三次 B 样条（degree $q=3$）：函数值、一次导数和二次导数在普通内部节点处连续，既足够平滑，也能表达区域结构。

关键性质：

1. **局部支撑**：$B_{j,q}$ 只在 $[\tau_j,\tau_{j+q+1}]$ 内非零；
2. **非负与分割统一性**：$B_{j,q}\ge0$，且在定义域内部 $\sum_jB_{j,q}(x)=1$；
3. **稀疏设计矩阵**：任意观测点只激活至多 $q+1$ 个基函数。三次样条每行至多 4 个非零量；
4. **不等于“不会过拟合”**：节点太密时它同样会追随噪声，仍需要正则化和交叉验证。

### 1.3 一维 VTEC 拟合的矩阵形式

将一列中位于 $x_i$ 的观测记为 $y_i=\mathrm{VTEC}_i$，构造

$$
\mathbf y=\mathbf B\mathbf c+\boldsymbol\varepsilon,
\qquad B_{ij}=B_{j,q}(x_i).
$$

基础加权最小二乘为

$$
\hat{\mathbf c}
=\arg\min_{\mathbf c}
(\mathbf y-\mathbf B\mathbf c)^T\mathbf W(\mathbf y-\mathbf B\mathbf c).
$$

仅靠这一式在站点少、站距不均或节点过密时会病态。推荐的平滑样条目标函数是

$$
\hat{\mathbf c}=
\arg\min_{\mathbf c}
\left[
(\mathbf y-\mathbf B\mathbf c)^T\mathbf W(\mathbf y-\mathbf B\mathbf c)
+\lambda\int\left(f''(x)\right)^2dx
\right].
$$

若把曲率惩罚写成 $\mathbf c^T\mathbf R\mathbf c$，则正规方程为

$$
(\mathbf B^T\mathbf W\mathbf B+\lambda\mathbf R)\hat{\mathbf c}
=\mathbf B^T\mathbf W\mathbf y.
$$

$\lambda$ 控制“贴数据”和“平滑”的平衡：$\lambda=0$ 容易在稀疏列中振荡；过大则会抹平 EIA 峰和真实扰动。

### 1.4 梯度如何从 B 样条直接得到

你的目标含“线性梯度”，而 B 样条的一个优势是可解析求导：

$$
\frac{\partial\widehat{\mathrm{VTEC}}}{\partial\phi}
=\sum_j\hat c_j\frac{dB_{j,q}(\phi)}{d\phi}.
$$

若要报告物理梯度而非 TECU/deg，沿南北方向可近似转换为

$$
g_N\,[\mathrm{TECU/km}]
=\frac{1}{R_E\pi/180}
\frac{\partial\widehat{\mathrm{VTEC}}}{\partial\phi[\mathrm{deg}]}
\approx\frac{1}{111.2}
\frac{\partial\widehat{\mathrm{VTEC}}}{\partial\phi[\mathrm{deg}]}.
$$

这比对 0.25° 栅格做生硬差分更稳定；同时仍应报告原始/未平滑差分作为敏感性分析。

---

## 2. 三篇 B 样条论文

### 2.1 Schmidt et al. (2011)：局部 B 样条优于全局球谐处理不均匀数据

**模型**：以 VTEC 相对于背景模型的改正数 $\Delta\mathrm{VTEC}$ 为对象，构造三维时空张量积模型：

$$
\Delta\mathrm{VTEC}(\lambda,\phi,t)
=\sum_{k_\lambda}\sum_{k_\phi}\sum_{k_t}
d_{k_\lambda k_\phi k_t}
\,B_{k_\lambda}(\lambda)B_{k_\phi}(\phi)B_{k_t}(t).
$$

经度采用周期三角 B 样条，纬度和时间采用端点插值的二次多项式 B 样条；与“空间球谐 + 时间 B 样条”的传统方案比较。

**结论**：在模拟的不均匀观测和数据缺口条件下，B 样条的误差主要局限在缺口邻域；球谐函数因全局支撑会产生远离缺口的伪振荡。文中例子中，B 样条方案的平均 RMS 明显低于球谐方案。

### 2.2 Nohutcu et al. (2010)：区域二维/三维 B 样条、鲁棒估计与病态处理

**二维模型**：在固定时间段内，将太阳固定经度与纬度归一化到 $[0,1]$，使用二次 B 样条张量积：

$$
\Delta\mathrm{VTEC}(x,y)
=\sum_{k_1}\sum_{k_2}d_{k_1k_2}
B_{k_1,2}(x)B_{k_2,2}(y).
$$

**三维模型**：在地固坐标系中增加时间维：

$$
\Delta\mathrm{VTEC}(x,y,z)
=\sum_{k_1,k_2,k_3}d_{k_1k_2k_3}
B_{k_1,2}(x)B_{k_2,2}(y)B_{k_3,2}(z).
$$

**估计策略**：使用 IRLS（bi-square 权函数）压低离群观测影响，并使用 LSQR 求解大型/病态系统。文中明确指出：无观测支撑区域对应的 B 样条系数不应任意估计；2D 更数值稳定，3D 更能描述窗口内时间变化。

### 2.3 Erdogan et al. (2021)：B 样条 + Kalman 滤波的实时全局扩展

**模型**：纬度二次多项式 B 样条与经度三角 B 样条的张量积；各 B 样条系数 $d_{k_1,k_2}(t)$ 随时间由自适应 Kalman 滤波递推，并与载波相位弧段偏差一同估计。

$$
\mathrm{VTEC}(\phi,\lambda,t)
=\sum_{k_1}\sum_{k_2}d_{k_1k_2}(t)
N_{k_1,2}(\phi)T_{k_2,3}(\lambda).
$$

**贡献**：局部 B 样条带来稀疏观测矩阵，适合实时大规模计算；数据空洞处使用由超快产品系数预报得到的补充信息稳定滤波。其 dSTEC 平均 RMS 报告为 1.05 TECU，略优于对比实时产品（1.12–1.22 TECU）。

---

## 3. 原始 IPP 宽带加权的一维纬向 B 样条

本方法**不先把 IPP 硬平均到 0.25° 格网**。每个 10 min 窗口内保留原始 IPP 的连续经纬度和 VTEC。设最终产品的经、纬输出间隔均为

$$
\Delta_\lambda=\Delta_\phi=0.25^\circ.
$$

第 $m$ 个输出列的中心经度为

$$
\lambda_m=\lambda_{min}+m\Delta_\lambda.
$$

对该列，建立一个只沿纬度拟合的一维模型；原始 IPP 的纳入范围可比输出格网宽。作为第一版方案，采用总宽 $W_\lambda=0.50^\circ$ 的经度条带，即半宽

$$
h_\lambda=\frac{W_\lambda}{2}=0.25^\circ,
$$

选择满足

$$
|\lambda_i-\lambda_m|\le h_\lambda
$$

的全部原始 IPP。以中心经度 $\lambda_m=120^\circ$ 为例，进入该条带拟合的 IPP 为

$$
119.75^\circ\le\lambda_i\le120.25^\circ.
$$

这些 IPP 的纬度保持原始值，不作纬向预平均。拟合得到的剖面只在该输出列的纬向格点上计算，并赋给以 $120^\circ$ 为中心、经度范围 $[119.875^\circ,120.125^\circ]$ 的 $0.25^\circ$ 输出格网单元。故 **0.25° 是输出采样间隔，0.50° 是拟合时借用观测的经向邻域宽度**。

相邻中心为 $120.25^\circ$ 的条带使用 $[120.00^\circ,120.50^\circ]$ 内的 IPP；两条带共同使用 $[120.00^\circ,120.25^\circ]$ 内的观测。该重叠使相邻列不会由互不相关的数据独立决定，从而减弱条带接缝。

#### 经向余弦权重与一维三次 B 样条

虽然进入条带的 IPP 来自不同经度，但模型自变量仍然**只有纬度**。经度只用于观测权重。对观测 $i$，定义经向余弦核权重：

$$
w_{im}^{(\lambda)}=
\begin{cases}
\dfrac12\left[1+\cos\left(\pi\dfrac{|\lambda_i-\lambda_m|}{h_\lambda}\right)\right],
&|\lambda_i-\lambda_m|\le h_\lambda,\\
0,&|\lambda_i-\lambda_m|>h_\lambda.
\end{cases}
$$

对 $\lambda_m=120^\circ$、$h_\lambda=0.25^\circ$，有

$$
w^{(\lambda)}(120^\circ)=1,\quad
w^{(\lambda)}(119.875^\circ)=w^{(\lambda)}(120.125^\circ)=0.5,\quad
w^{(\lambda)}(119.75^\circ)=w^{(\lambda)}(120.25^\circ)=0.
$$

因此，条带中心附近的 IPP 主导拟合，边缘 IPP 只提供平滑的辅助信息。将纬度标准化为

$$
x_i=\frac{\phi_i-\phi_{min}}{\phi_{max}-\phi_{min}}\in[0,1].
$$

固定 B 样条次数为三次（degree 3），第 $m$ 条带的纬向剖面写为

$$
\widehat V_m(\phi)=\sum_{j=1}^{P_m}c_{mj}B_{j,3}(x).
$$

其中 $P_m$ 是基函数数量，$c_{mj}$ 是待估系数。

对当前条带，以加权平滑最小二乘估计系数：

$$
\hat{\mathbf c}_m=
\arg\min_{\mathbf c_m}
\left\{
\sum_{i\in\mathcal I_m}
w_{im}^{(\lambda)}w_i^{(q)}
\left[V_i-\sum_jc_{mj}B_{j,3}(x_i)\right]^2
+\alpha\int_0^1\left[\widehat V_m''(x)\right]^2dx
\right\}.
$$

其中，$\mathcal I_m$ 为条带内 IPP 集合，$w_i^{(q)}$ 为观测质量权重（可综合仰角、VTEC 方差及距窗口中心的时间），$\alpha$ 为纬向曲率平滑参数。第一项保证贴合观测，第二项抑制无物理依据的锯齿和假梯度。

#### 最终格网化与相邻条带的平滑来源

当前原型中，各条带独立求解，平滑主要来自两点：

1. 相邻条带的观测域重叠，且对共享 IPP 使用连续的余弦权重；
2. 每条带内部的 B 样条曲率惩罚。

若试验后仍有明显接缝，再增加相邻剖面软约束；对一组公共纬度检查点 $\phi_r$，可在总目标函数中附加

$$
\beta\sum_m\sum_r
\left[\widehat V_m(\phi_r)-\widehat V_{m+1}(\phi_r)\right]^2,
$$

其中 $\beta$ 应从较小值开始；**该项是可选的后续升级**，不能过强，以免抹去真实经向差异。

#### 条带宽度的初始设置与比较及进一步优化建议

当前确定的起始值为 $W_\lambda=0.50^\circ$；它使相邻输出中心相差 $0.25^\circ$ 的两条带拥有 $0.25^\circ$ 的共同观测域。建议随后比较

$$
W_\lambda\in\{0.50^\circ,0.75^\circ,1.00^\circ\}.
$$

带宽扩大能提高每条带的有效样本量和列间平滑性，但会减弱经向细节。应以当前窗口的空间留出误差为主标准，选择预测性能接近最优时最窄的带宽；故不应把 0.25° 输出间隔表述为 0.25° 的真实经向有效分辨率。

建议 $w_i^{(q)}$ 至少综合仰角、VTEC 方差与时间距窗口中心；**例如令较低仰角权重更小，离 $t_0$ 越远权重越低**。

---

## 4. 动态选择内部节点数

### 4.1 选择对象与候选模型

这里不采用“动态阶次”这一名称，因为 B 样条次数固定为三次。动态选择的是内部节点数 $K$，即模型复杂度。$K$ 越大，剖面越灵活，也越可能把噪声拟合为虚假的局地梯度。第一版建议让一个 10 min 窗口内所有条带使用同一个 $K$，避免相邻条带复杂度突然变化；候选集可取

$$
\mathcal K=\{0,1,2,3,4,5,6,7\}.
$$

### 4.2 当前窗口的分组空间交叉验证

不能用训练 RMSE 选择 $K$，因为节点越多，训练误差几乎总会越小。应只使用当前 10 min 窗口内的观测，按接收机或站星弧段分组做 $G$ 折交叉验证。例如有 20 个接收机时可作 5 折：每次留出 4 个接收机的全部观测，用其余 16 个接收机拟合，再预测被留出的观测。

对候选 $K$，定义空间预测误差

$$
E_{CV}(K)=
\sqrt{\frac{1}{N_{test}}
\sum_{i\in test}
\left[V_i-\widehat V_{-g(i),K}(\lambda_i,\phi_i)\right]^2},
$$

其中 $g(i)$ 是观测 $i$ 所在的留出组，$\widehat V_{-g(i),K}$ 表示不使用该组数据拟合的模型。分组而非随机留点可避免同一站星弧段的相邻、高相关观测同时出现在训练集和测试集而导致误差虚低。

### 4.3 预测性能接近时选更简单、更平滑的节点数

令

$$
E_{min}=\min_{K\in\mathcal K}E_{CV}(K).
$$

记最优 CV-RMSE 的标准误为 $SE(E_{min})$，构造性能可接受的候选集合（one-standard-error rule）：

$$
\mathcal K_{good}=
\left\{K\in\mathcal K:
E_{CV}(K)\le E_{min}+SE(E_{min})\right\}.
$$

这一步的意思是：只要较少节点的预测性能与最优节点数没有显著差异，就不为极小的 RMSE 改善而增加复杂度。

然后使用全部当前窗口观测，对每个 $K\in\mathcal K_{good}$ 重新拟合，计算所有条带的曲率能量：

$$
J_{curv}(K)=
\sum_m\int_{\phi_{min}}^{\phi_{max}}
\left[\frac{d^2\widehat V_{m,K}(\phi)}{d\phi^2}\right]^2d\phi.
$$

最终选择

$$
K^*=\arg\min_{K\in\mathcal K_{good}}J_{curv}(K).
$$

若仍并列，则选较小的 $K$。例如，若 $K=5$ 的 CV-RMSE 最低为 0.80 TECU、标准误为 0.03 TECU，而 $K=4,5,6$ 的 CV-RMSE 分别为 0.81、0.80、0.81 TECU，则三者均进入 $\mathcal K_{good}$；若它们的曲率能量分别为 13、25、49，最终选择 $K^*=4$。

该选择的优先顺序为：

$$
\text{当前窗口空间预测能力}
\rightarrow
\text{避免不必要的模型复杂度}
\rightarrow
\text{最小化无物理依据的梯度振荡}.
$$

---



## 参考文献

- Schmidt, M., Dettmering, D., Mößmer, M., Wang, Y., & Zhang, J. (2011). Comparison of spherical harmonic and B spline models for the vertical total electron content. *Radio Science, 46*, RS0D11. https://doi.org/10.1029/2010RS004609
- Nohutcu, M., Karslioglu, M. O., Gucluer, B., Schmidt, M., Zeilhofer, C., Zhang, Z., & Ergintav, S. (2010). B-spline modeling of VTEC over Turkey using GPS observations. *Journal of Atmospheric and Solar-Terrestrial Physics, 72*, 617–624.
- Erdogan, E., Schmidt, M., Goss, A., Görres, B., & Seitz, F. (2021). Real-time monitoring of ionosphere VTEC using multi-GNSS carrier-phase observations and B-splines. *Space Weather, 19*, e2021SW002858. https://doi.org/10.1029/2021SW002858





