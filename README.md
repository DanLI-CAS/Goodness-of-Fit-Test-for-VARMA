### 弱误差情况下的向量自回归移动平均模型的拟合优度检验

本项目针对向量自回归移动平均模型（VARMA(p,q)）进行拟合优度检验。给定$d$维弱平稳时间序列数据$\{X_{t},t=1,\ldots,n\}$, 本项目检验
\begin{equation}\label{hypothesis2}
	%\begin{aligned}
	%&H_{0}:\left\{X_{t}\right\} \text{ 是VARMA(p,q)序列 vs.}\\
	%&H_{1}:\left\{X_{t}\right\}\text{ 不是VARMA(p,q)序列或者满足VARMA(P,Q)模型，且$P>p$或者$Q>q$} 
	%&H_{1}:\left\{X_{t}\right\}\text{ 不是VARMA(p,q)序列}
	%\end{aligned}
	H_{0}:\left\{X_{t}\right\} \text{ 是VARMA(p,q)序列\quad v.s}\quad H_{1}:\left\{X_{t}\right\}\text{ 不是VARMA(p,q)序列}
\end{equation}

针对宏观经济与金融数据中常见的弱误差（误差项不相关但不独立，比如存在条件异方差）情形，传统的拟合优度检验通常会产生严重的渐近分布偏移。本项目提出了一种新的拟合优度检验方法。 我们基于托普利茨矩阵构造统计量，并采用自助法来获得假设检验的临界值。本项目所提假设检验方法能够有效控制第一类错误概率，且具有良好的检验功效。

#### `zzr.varma.gof.boot1`与 `zzr.varma.gof.boot2` : 基于随机加权自助法的VARMA(p,q)拟合优度检验 

####  函数描述 
这两个函数用于对向量自回归移动平均模型VARMA(p,q)进行拟合优度检验, `zzr.varma.gof.boot1`使用标准指数分布$Exp(1)$作为随机权重, `zzr.varma.gof.boot2`使用Mammen分布作为随机权重。

#### 參數說明
* **`zt`**: 数值型矩阵, 输入的多元时间序列数据, 要求维度为 $d \times n$, 其中 $d$ 为数据维数, $n$ 为数据量.
* **`beta0`**: 数值型向量, VARMA 模型参数估计的迭代初始值向量, 为$d\times (p+q)$向量.
* **`p`**: 整数, 向量自回归部分的阶数.
* **`q`**: 整数, 向量移动平均部分的阶数.
* **`k_max`**: 整数, 检验统计量所需的最大滞后阶数.
* **`boot`**: 整数, 自助法的重抽样次数, 默认值为 1000.

#### 返回值
函数返回一个列表 (List), 包含以下四个元素：
* **'statistic'**: 计算得到的拟合优度检验统计量的值.
* **'p_value'**: 计算出的 $p$ 值.
* **'reject_5pct'**: 在 5% 的显著性水平下是否拒绝原假设, TRUE 表示拒绝, FALSE 表示不拒绝.
* **'reject_10pct'**:: 在 10% 的显著性水平下是否拒绝原假设.

#### 细节说明
* **'参数估计与残差计算'**：函数调用 zzr.varma.est 和底层的 zzrvarmaResiduals_cpp 对原始数据进行参数估计并提取残差.
* **'统计量计算'**: 函数调用 MahdiMcLeod 函数计算基于残差的多元拟合优度检验统计量.
* **'随机权重自助法'**：函数通过加权拟似然函数重估参数（zzr.varma.boot.est），并利用 C++ 拓展 (zzr_boot_cov_cpp) 计算自协方差矩阵.
* **'去除均值'**：本函数默认已进行数据的期望为0, 在使用时需首先除均值.

#### 依赖说明
运行此函数前需要确保已加载/编译对应的 C++ 函数(varmalikelihood.cpp), 并在环境中加载`zzr.varma.est`, `zzr.varma.boot.est`, `MahdiMcLeod`

#### 示例
```
data<-t(mvrnorm(500,rep(0,2),diag(2)) 

beta_init<-c(0,0,0,0,0,0,0,0)

test_results <- zzr.varma.gof.boot1(
  zt = my_data, 
  beta0 = beta_init, 
  p = 1, 
  q = 1, 
  k_max = 5, 
  boot = 1000
)

# 查看检验结果
print(test_results$p_value)
print(test_results$reject_5pct)
```


#### `MahdiMcLeod.wntest.boot1`与 `MahdiMcLeod.wntest.boot2` : 基于随机加权自助法的白噪声检验

####  函数描述 
这两个函数用于对数据进行白噪声检验, `MahdiMcLeod.wntest.boot1`使用标准指数分布$Exp(1)$作为随机权重, `MahdiMcLeod.wntest.boot2`使用Mammen分布作为随机权重。


#### 參數說明
* **`X`**: 数值型矩阵, 输入的多元时间序列数据, 要求维度为 $d \times n$, 其中 $d$ 为数据维数, $n$ 为数据量.
* **`k_max`**: 整数, 检验统计量所需的最大滞后阶数.
* **`boot`**: 整数, 自助法的重抽样次数, 默认值为 1000.

#### 返回值
函数返回一个列表 (List), 包含以下四个元素：
* **'statistic'**: 计算得到的拟合优度检验统计量的值.
* **'p_value'**: 计算出的 $p$ 值.
* **'reject_5pct'**: 在 5% 的显著性水平下是否拒绝原假设, TRUE 表示拒绝, FALSE 表示不拒绝.
* **'reject_10pct'**:: 在 10% 的显著性水平下是否拒绝原假设.

#### 示例
```
data<-t(mvrnorm(500,rep(0,2),diag(2)) 

test_results <- MahdiMcLeod.wntest.boot1(
  X = my_data, 
  k_max = 5, 
  boot = 1000
)

# 查看检验结果
print(test_results$p_value)
print(test_results$reject_5pct)
```
