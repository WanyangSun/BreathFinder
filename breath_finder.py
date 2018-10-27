# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 7:52:29 2018

@author: Thomas Sun
"""
import os
import sys
# 定义工作路径
workdir = r'D:\Desktop\Breath Coding\BreathFinder_v1.0\bin'
os.chdir(workdir)
sys.path.append(workdir)

# 导入自编程序包
import load_xml
import extract_signal as es
import batch_extract as be

# 给出质谱数据名称
filename = 'data1.mzXML'
listname = 'Targeted ion list.csv'
#filename = input('Please input a ms data file (xxx.mzXML): ')


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

dic_id = be.check_all_sig(dic_msdata, mspola, signal_weight=0.8)


# In[4] 提取一个样品中强度高于一定阈值的全部m/z信号（“时间-强度矩阵”和“呼气信号&
# 背景结果矩阵”），并生成数据提取结果.csv文件

#intensity_threshould = float(input('Please give a intensity threshould for 
#                                   signal extraction: '))

df_dic, mzlist_dic = be.extract_all_mz(out_path, dic_msdata, mspola, dic_id, 
                                       intensity_thre=1e6, min_scan_num=5)


# In[5] 列表中m/z靶向批量提取

'''
使用m/z靶向提取功能，需要：
1. 在最上方提交文件名称（listname = 'Targeted ion list.csv'）。导入的target ion 
list，文件格式为.csv，至少需要包含Compound，m/z，Polarity (POS/NEG)三列。
2. 将dic_scr传入be.extract_all_mz函数（如下所示）
注：使用该功能时，intensity_thre和min_scan_num参数无效。详细阐述请运行“?(be.extract_all_mz)”查看。
'''

# 准备m/z靶向筛查数据
dic_scr = be.read_target_list(list_path)
# 执行m/z靶向筛查
df_dic, mzlist_dic = be.extract_all_mz(out_path, dic_msdata, mspola, dic_id, 
                                      dic_scr)

###input('Please enter any key to exit.')
