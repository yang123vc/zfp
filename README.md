# Introduction 
基于中控指纹识别器的指纹识别课程项目（课程指导项目代码+个人实现）

# Content
  1. 创建项目和导入采集器SDK OK
  1. 接收和显示指纹图像 OK
  1. 指纹图像中值滤波 OK
  1. 指纹图像直方图均衡化 OK
  1. 指纹脊线方向计算 OK
  1. 指纹脊线频率计算 OK
  1. 指纹掩码计算 OK
  1. 指纹图像Gabor滤波增强 OK
  1. 指纹图像二值化 OK
  1. 指纹图像细化 OK
  1. 指纹特征提取 OK
  1. 指纹特征过滤 OK
  1. 指纹特征入库 OK
  1. 指纹特征匹配


# 改进
  1. 修改指纹边缘特征过滤算法
      - 修改边缘特征点过滤算法，并优化利于并行计算，原算法会误丢弃正常特征点
  2. 采用SQLite数据库存储指纹特征和登记信息
      - SQLite是一个轻量型的关系型数据库
      - 本系统的数据库直接存储BLOB二进制数据

