print("hello world")
# 整数
x = 10

# 浮点数
y = 3.14

# 字符串
name = "Alice"

# 布尔值
is_student = True


# if语句
if x > 0:
    print("x是正数")

# for循环
for i in range(5):
    print(i)

# while循环
while x > 0:
    print(x)
    x -= 1
    
    
def greet(name):
    return "Hello, " + name + "!"

print(greet("Alice"))

# 列表
fruits = ["apple", "banana", "cherry"]

# 字典
person = {"name": "Alice", "age": 25}

import pandas as pd

# # 创建DataFrame
# data = {'Name': ["John", "Anna", "Peter", "Linda"],
#         'Location': ["New York", "Paris", "Berlin", "London"],
#         'Age': [24, 13, 53, 33]}
# 
# df = pd.DataFrame(data)

# 查看DataFrame
print(df)
library(reticulate)
# 安装AnnData和Scanpy
!pip install anndata scanpy

# 使用Scanpy进行单细胞数据分析
import scanpy as sc

# 加载你的单细胞数据集
adata = sc.read("your_data_file.h5ad")

# 基本的数据预处理
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
