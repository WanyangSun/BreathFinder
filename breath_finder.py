# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 7:52:29 2018

@author: Thomas Sun
"""
import os
import sys
import pandas as pd
# 定义工作路径
workdir = r'D:\Desktop\Breath Coding\BreathFinder_v1.0\bin'
os.chdir(workdir)
sys.path.append(workdir)

# 导入自编程序包
import load_xml
import extract_signal as es
import batch_extract as be
# 设置初始变量
id_p, id_n = 0, 0
# 给出质谱数据名称
filename = 'data1.mzXML'
listname = 'Targeted ion list.csv'
###filename = input('Please input a ms data file (xxx.mzXML): ')


# In[1] 提交要处理的质谱数据文件名称
# 在be.gen_workspace()的括号中填写：'xxx.mzXML'（注：一定要有单引号）
in_path, out_path, list_path = be.gen_workspace(filename, listname)

# In[2] 读取mzXML质谱数据
ms, ms2 = load_xml.load_mzxml_file(in_path)
# 拆分ms正负离子数据，获取ms极性信息
ms_pos, ms_neg, mspola = es.sep_polarity(ms)
dic_msdata = {'+':ms_pos, '-':ms_neg}

# In[3] 测试某一离子获取呼气信号和背景信号的准确度
# target_mass输入格式：0 = 总离子流，xxx.xxxx = 目标离子

# 方法一：同时提取两个离子模式的数据：
dic_id = be.check_all_sig(dic_msdata, mspola, signal_weight=0.8)


# 方法二：单独提取一个离子模式的数据（推荐，可以分别设置参数）
#id_p, none = be.check_all_sig(ms1, '+', baseline_poly=1, signal_weight=0.8)
#
#none, id_n = be.check_all_sig(ms1, '-', baseline_poly=1, signal_weight=0.8, t_start=1.2, t_end=3.5)

# 详细参数
#id_p, id_n = be.check_all_sig(ms1, mspola, smooth_num=7, smooth_poly=1,
#                              baseline_poly=1, index_window=3, spec_num=8,
#                              signal_weight=0.9, count=0, t_start=0, t_end=0)
# 增加时间限定功能：
# 时间区间可以通过传入t_start和t_end两个参数限定。默认t_end=0，对应最大时间。


# In[4] 提取一个样品中强度高于一定阈值的全部m/z信号（“时间-强度矩阵”和“呼气信号&背景结果矩阵”），并生成数据提取结果.csv文件

###intensity_threshould = float(input('Please give a intensity threshould for signal extraction: '))

# 方法一：不传入id_p和id_n参数，提取时间-离子强度矩阵
#df_p, df_n = be.extract_all_mz(out_path, ms1, mspola, intensity_thre=1e7)
## 提取限定时间范围内的时间-离子强度矩阵
#df_p, df_n = be.extract_all_mz(out_path, ms1, mspola, intensity_thre=1e7,
#                               t_start=0.2, t_end=2.5)
#
# 方法二：传入id_p和id_n参数，提取时间-离子强度矩阵后，提取呼气信号&背景结果矩阵
###df_p, df_n = be.extract_all_mz(out_path, ms1, mspola, id_p, id_n,
###                               intensity_thre=intensity_threshould)

df_p, df_n, mzlist_unique = be.extract_all_mz(out_path, dic_msdata, '+', dic_id, intensity_thre=5e4, min_num=5)

# 不生成时间-强度矩阵&信号提取矩阵
df_p, df_n, mzlist_unique = be.extract_all_mz(out_path, dic_msdata, mspola, id_p, id_n, intensity_thre=1e5, ext_matrix=False, ext_signal=False)


###input('Please enter any key to exit.')

## 详细参数
#df_p, df_n = be.extract_all_mz(out_path, ms1, mspola, id_p, id_n,
#                               intensity_thre=1e5, mz_Da=0.001, mz_ppm=2, min_num=5, t_start=0, t_end=0, 
#                               ext_matrix=True, ext_signal=True)
# 1 增加时间限定功能：
# 时间区间可以通过传入t_start和t_end两个参数限定。默认t_end=0，对应最大时间。
# 2 增加某个离子最少连续在几张谱图中出现的限定
# min_num参数，默认=5，可自行设置
# 3 增加matrix和signal矩阵是否生成功能：
# ext_matrix表示是否生成提取离子文件，ext_signal表示是否生成提取信号文件。生成signal矩阵，必须先生成matrix


# In[5] 列表中m/z靶向批量提取






# In[5] 目标离子提取
#
#input_path_target, output_path_target = be.gen_workspace('Input_list_1_Targeted ion list.csv')
#
#
#
## In[]
#
#
#import matplotlib.pyplot as plt
#import seaborn as sns
#
#
#int_pos = df_intense_p.iloc[:, 1:-1]
#index_name = [round(i,4) for i in df_intense_p['m/z']]
#int_pos.index = index_name
#column_name = [round(i,2) for i in int_pos.columns]
#int_pos.columns = column_name
#
#sns.set()
#
#g= sns.clustermap(int_pos, fmt="d", cmap='YlGnBu', z_score=0, col_cluster=False, figsize=(12, 12))
#ax = g.ax_heatmap
#ax.set_title('Hierarchical Clustering', fontsize=16)
#ax.tick_params(axis='y',labelsize=12)
#ax.set_xticklabels([])
#ax.set_xlabel('time',fontsize=16)
#ax.set_ylabel('m/z',fontsize=16)
#plt.show()
#


