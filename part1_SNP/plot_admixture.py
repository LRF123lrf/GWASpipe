import matplotlib
matplotlib.use('Agg')  # 使用 Agg 后端

import pandas as pd
import matplotlib.pyplot as plt

sort -k1,1n -k2,2n input.Q > sorted_input.Q

# 指定文件路径
Q_file = 'sorted_input.Q'

# 读取 Q 文件
data = pd.read_csv(Q_file, sep=' ', header=None)

# 每行代表一个个体，每列代表该个体分配到每个群体的概率
K = data.shape[1]  # 列数即为 K 的值

# 绘图
fig, ax = plt.subplots(figsize=(10, 6))

# 绘制堆叠条形图
data.plot(kind='bar', stacked=True, ax=ax, width=1, edgecolor='none')

# 设置标签
ax.set_xlabel('Individuals')
ax.set_ylabel('Ancestry Proportion')
ax.set_title(f'Ancestry Proportions for K={K}')

# 去掉 x 轴刻度
plt.xticks([])

# 保存图像到文件
plt.savefig('admixture_plot.png')

