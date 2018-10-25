
import os
import time
import pandas as pd
import numpy as np
from interval import Interval
from progressbar import *
import extract_signal as es
import stat_signal as ss

dic_pola = {'+':'positive', '-':'negative'}

# 给出数据导入和导出的文件（夹）路径
def gen_workspace(filename, listname):
    data_input_path = r'..\data_input' + '\\' + filename
    list_input_path = r'..\data_input' + '\\' + listname
    data_output_path = r'..\data_output'
    try:
        os.makedirs(r'..\data_output')
        print('The output folder has been established.')
    except Exception as e:
        print('The output folder is already exist.')
    return data_input_path, data_output_path, list_input_path


## 批量提取离子流

## 获取唯一mz
def get_unique_mz(mz_list, mz_unique, mz_Da=0.001, mz_ppm=2):
    while len(mz_list) != 0:
        mz = mz_list[0]
        mz_unique.append(mz)
        if mz*(1e-6*mz_ppm) <= mz_Da:
            mz_interval = mz_list[(mz_list>mz-mz_Da) & (mz_list<mz+mz_Da)]
        else:
            mz_interval = mz_list[(mz_list>mz*(1-1e-6*mz_ppm)) & (mz_list<mz*(1+1e-6*mz_ppm))]
        mz_list = mz_list.drop(mz_interval)
    return mz_unique

# 获取ms数据中强度大于设定阈值的mzlist
def get_unique_mzlist(msdata, scan_index, intensity_threshold=1e3, mz_Da=0.001, mz_ppm=2):
    start_time = time.time()
    mzlist_all = []
    # 获取全部图谱中的mz
    for i in scan_index[0]:
        mzlist = msdata.data[i].peaks
        mzlist_all += [mzlist[i][0] for i in range(1,len(mzlist)) if
                       mzlist[i-1][1] < mzlist[i][1] > mzlist[i+1][1]
                       and mzlist[i][1] > intensity_threshold]
    mz_count = round(pd.Series(mzlist_all), 5).value_counts()
    mzlist = mz_count.index
    # print(len(mzlist_all))
    # 排序后，获取mz_Da和mz_ppm窗口以外的mz
    mz_unique = []
    mzlist_unique = get_unique_mz(mzlist, mz_unique, mz_Da, mz_ppm)
    end_time = time.time()
    print('{} unique m/z values were extracted. Running time:{:.4f} s.'.format(
            len(mzlist_unique), (end_time - start_time)))
    return mzlist_unique, mzlist_all

# 交互式判断分区是否合理
def check_sig(msdata, polarity, smooth_num=7, smooth_poly=1, baseline_poly=1,
                         index_window=3, spec_num=8, signal_weight=0.9, count=0, t_start=0, t_end=0):
    bg_id = 0
    sig_id = 0
    mz_for_sig = 0
    if count == 0:
        bg_id, sig_id = es.get_key_index(msdata, polarity, mz_for_sig, smooth_num, smooth_poly,
                                         baseline_poly, index_window, spec_num, signal_weight,
                                         t_start, t_end)
    count += 1
    check = input('Are all signals extracted? (Y/N):')
    if check in ['Y', 'y']:
        print('Detected signals and backgrounds are stored.')
    else:
        mz_for_sig = input('Please submit a m/z for signals searching (xxx.xxxx):')
        # 输入的m/z是0或者未输入，都认为是总离子流
        if mz_for_sig == '0' or mz_for_sig == '':
            mz_for_sig = 0
        else:
            mz_for_sig = float(mz_for_sig)
        bg_id, sig_id = es.get_key_index(msdata, polarity, mz_for_sig, smooth_num, smooth_poly,
                                             baseline_poly, index_window, spec_num, signal_weight,
                                             t_start, t_end)
        check_sig(msdata, polarity, count=count)
    return bg_id, sig_id


# 交互式判断分区是否合理
def check_all_sig(dic_msdata, pola, smooth_num=7, smooth_poly=1, baseline_poly=1,
                      index_window=3, spec_num=8, signal_weight=0.9,
                      t_start=0, t_end=0):
    id_p, id_n = 0, 0
    dic_id = {'+':id_p, '-':id_n}
    for i in pola:
        print('Please check the extraction result.')
        count = 0
        bg_id, sig_id = check_sig(dic_msdata[i], i, smooth_num, smooth_poly, baseline_poly,
                                      index_window, spec_num, signal_weight, count, t_start, t_end)
        dic_id[i] = [sig_id, bg_id]
    return dic_id


# 生成所有质荷比的时间×强度矩阵
def gen_intense_matrix(mz_intense_list, mzlist_uni, time_list, filepath_intense):
    df_mz_intense = pd.DataFrame(mz_intense_list, columns=time_list)
    df_mz_intense['Average'] = df_mz_intense.apply(lambda x: x.mean(), axis=1)
    df_mz_intense.insert(0, 'm/z', mzlist_uni)
    df_mz_intense.to_csv(filepath_intense, index=False)
    return df_mz_intense


# 检查某一个离子出现的次数是否符合要求
def check_continuous(eic, min_num=5, check=False):
    index = 0
    while index <= len(eic)-min_num:
        if sum(eic[index:index+min_num] != 0.01) == min_num:
            check = True
            break
        else:
            index += 1
    return check


# 生成时间强度矩阵，并保存到output文件夹
def get_intense_matrix(msdata, pola, time_list, mzlist_uni, filepath_intense, min_num=5, start_in_min=0, end_in_min=0):
    total = len(mzlist_uni)
    mzlist_uni_temp = mzlist_uni[:]
    mz_intense_list = []
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    print('\nStep 1 Establising %s intensity table.\n' % (dic_pola[pola]))
    for index, mz in enumerate(mzlist_uni):
        pbar.update(10 * index + 1)
        eic_temp = msdata.eic(mz, t_start=start_in_min, t_end=end_in_min)
        check = check_continuous(eic_temp, min_num)
        if check == True:
            mz_intense_list.append(eic_temp)
        else:
            mzlist_uni_temp.remove(mz)
#            print(len(mzlist_uni))
    pbar.finish()
    df_intense = gen_intense_matrix(mz_intense_list, mzlist_uni_temp, time_list, filepath_intense)
    print('\nIntensity table has been established.\n')
    return mzlist_uni_temp, mz_intense_list, df_intense


# 生成信号统计矩阵，并保存到output文件夹
def get_signals_matrix(msdata, pola, time_list, id_range, mzlist_uni, mz_intense_list, filepath):
    total = len(mzlist_uni)
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    extract_data = []
    print('\nStep 2 Extracting %s signals and backgrounds matrix.\n' % (dic_pola[pola]))
    for index, mz in enumerate(mzlist_uni):
        pbar.update(10 * index + 1)
        extract_data.append(ss.Stat_signal(id_range[0], id_range[1], msdata.eic(mz)))
    pbar.finish()
        # 生成呼气信号平均值&标准偏差矩阵
    ss.batch_statistic(mzlist_uni, time_list, mz_intense_list, extract_data, filepath)
    print('\nSignals and backgrounds matrix has been established.\n')


# 根据生成的信号区间，提取一个样品中的全部mz的对应信号
def extract_all_mz(output_path, dic_msdata, pola, dic_id,
                   intensity_thre=1e5, mz_Da=0.001, mz_ppm=2, min_num=5, t_start=0, t_end=0, 
                   ext_matrix=True, ext_signal=True):
    '''
    [summary]
        Extract time-intensity matrix and signal-background matrix
        in positive and negative ionization mode, respectively.

    [description]
        Step 1 Extract time-intensity matrix.
        Step 2 Extract signal-background matrix. (if id_p and id_n != 0)

    Arguments:
        output_path {[str]} -- [Output file path and name.]
        msdata {[class(Ioncurrent)]} -- [description]
        pola {[set]} -- [{'+', '-'}]

    Keyword Arguments:
        id_p {list} -- Signal and background ids in positive ionization mode.
                         (default: {0})
        id_n {list} -- Signal and background ids in negative ionization
                          mode] (default: {0})
        intensity_thre {num} -- Intensity threshold of a m/z (default: {1e5})
        mz_Da {num} -- m/z window in Da. This parameter is work with mz_ppm.
                       The bigger one will be excecute. (default: {0.001})
        mz_ppm {num} -- m/z window in ppm. This parameter is work with mz_Da.
                       The bigger one will be excecute.(default: {2})
        t_start {num} -- Time (min) to start extracting. (default: {0})
        t_end {num} -- Time (min) to start extracting. (default: {0})

    Returns:
        [DataFrame] -- [Time-intensity matrices.]
    '''
    msdata = list(dic_msdata.values())[0]
    filename = msdata.data[0].filename.split('.')[0]
    df_intense_p, df_intense_n = [], []
    dic_df_in = {'+':df_intense_p, '-':df_intense_p}
    for i in pola:
        msdata = dic_msdata.get(i)
        filepath = output_path + '\\' + filename + '_' + i + '_signals' + '.csv'
        filepath_intense = output_path + '\\' + filename + '_' + i + '_intensity_matrix' + '.csv'
        if t_end == 0:
            t_end = msdata.time[-1]
        # 将保留时间转换为scan id
        scan_index = np.where((msdata.time>=t_start) & (msdata.time<=t_end))
        time_list = msdata.time[scan_index]
        mzlist_uni, mzlist_all = get_unique_mzlist(msdata, scan_index, intensity_thre, mz_Da, mz_ppm)
        if ext_matrix:
            # 提取所有m/z的离子强度信息
            mzlist_uni, mz_intense_list, df_intense = get_intense_matrix(msdata, i, time_list,
                                                               mzlist_uni, filepath_intense,
                                                               min_num, t_start, t_end)
        dic_df_in[i], id_extract = df_intense, dic_id[i]
        if ext_signal:
            if id_extract != 0:
                # 判断id_p有无
                if id_extract[0] == []:
                    print('No signal has been detected.')
                else:
                    # 生成呼气信号平均值&标准偏差矩阵
                    get_signals_matrix(msdata, i, time_list, id_extract, mzlist_uni, mz_intense_list, filepath)
    # 输出离子时间×强度矩阵
    return df_intense_p, df_intense_n, mzlist_uni





