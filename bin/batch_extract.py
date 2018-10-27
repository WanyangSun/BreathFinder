
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
def get_unique_mz(mz_list, mz_Da=0.001, mz_ppm=2):
    mz_unique = []
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
    for i in scan_index:
        mzlist = msdata.data[i].peaks
        mzlist_all += [mzlist[i][0] for i in range(1,len(mzlist)) if
                       mzlist[i-1][1] < mzlist[i][1] > mzlist[i+1][1]
                       and mzlist[i][1] > intensity_threshold]
    mz_count = round(pd.Series(mzlist_all), 5).value_counts()
    mzlist = mz_count.index
    # print(len(mzlist_all))
    # 排序后，获取mz_Da和mz_ppm窗口以外的mz
    mzlist_unique = get_unique_mz(mzlist, mz_Da, mz_ppm)
    end_time = time.time()
    print('{} unique m/z values were extracted. Running time:{:.4f} s.'.format(
            len(mzlist_unique), (end_time - start_time)))
    return mzlist_unique


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
def check_all_sig(dic_msdata, pola, smooth_window=7, smooth_poly=1,
                  baseline_poly=1, index_window=3, spec_num=8,
                  signal_weight=0.8, t_start=0, t_end=0):
    '''
    Description
        Obtain signal and background scan ids.

    Parameters
    :param dic_msdata: dic
        Mass spectra data set.
    :param pola: dic
        Polarity ('+'/'-'/{'+','-'}).
    :param smooth_window: int
        The length of the filter window (i.e. the number of coefficients).
        smooth_window must be a positive odd integer. If mode is ‘interp’,
        window_length must be less than or equal to the size of x. Default is 7.
    :param smooth_poly: int
        The order of the polynomial used to fit the samples. smooth_poly must
        be less than window_length. Default is 1.
    :param baseline_poly: int
        Degree of the polynomial that will estimate the data baseline. A low
        degree may fail to detect all the baseline present, while a high degree
        may make the data too oscillatory, especially at the edges. Default is
        1.
    :param index_window: int
        The adjust window of scan id obtained by smoothed data. Default is 3.
    :param spec_num: int
        Minimum scan span over which the signal/background platform must be
        observed. Default is 8.
    :param signal_weight:
        The weight to combine the baselines to the bottom and top of the ion
        current. Default is 0.8.
    :param t_start: float
        Start time (min). Default is 0.
    :param t_end: float
        End time (min). Default is 0.

    Returns:
        Signal and background id intervals.
        :type: dic
    '''
    id_p, id_n = 0, 0
    dic_id = {'+':id_p, '-':id_n}
    for i in pola:
        print('Please check the extraction result.')
        count = 0
        bg_id, sig_id = check_sig(dic_msdata[i], i, smooth_window, smooth_poly, baseline_poly,
                                      index_window, spec_num, signal_weight, count, t_start, t_end)
        dic_id[i] = [sig_id, bg_id]
    return dic_id


# 生成所有质荷比的时间×强度矩阵
def gen_intense_matrix(mz_intense_list, mzlist_uni, time_list, output_path, df_screen=None):
    df_mz_intense = pd.DataFrame(mz_intense_list, columns=time_list)
    df_mz_intense['Average'] = df_mz_intense.apply(lambda x: x.mean(), axis=1)
    df_mz_intense.insert(0, 'm/z (exp)', mzlist_uni)
    if df_screen is not None:
        df_mz_intense = pd.merge(df_screen.reset_index(drop=True), df_mz_intense.iloc[:,1:], left_index=True, right_index=True)
    df_mz_intense.to_csv(output_path, index=False)
    return df_mz_intense


# 检查某一个离子出现的次数是否符合要求
def check_continuous(eic, min_scan_num=5, check=False):
    index = 0
    while index <= len(eic)-min_scan_num:
        if sum(eic[index:index+min_scan_num] != 0.01) == min_scan_num:
            check = True
            break
        else:
            index += 1
    return check


# 生成时间强度矩阵，并保存到output文件夹
def get_intense_matrix(msdata, pola, time_list, mzlist_uni, output_path, mz_Da,
                       min_scan_num=5, t_start=0, t_end=0, df_screen=None):
    filename = msdata.data[0].filename.split('.')[0]
    file_path = output_path + '\\' + filename + '_' + pola + '_intensity_matrix.csv'
    if t_end == 0:
        t_end = msdata.time[-1]
    total = len(mzlist_uni)
    mzlist_uni_temp = mzlist_uni[:]
    mz_intense_list = []
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    print('\nEstablising %s intensity table.\n' % (dic_pola[pola]))
    for index, mz in enumerate(mzlist_uni):
        pbar.update(10 * index + 1)
        eic_temp = msdata.eic(mz, mz_Da, False, t_start, t_end)
        check = check_continuous(eic_temp, min_scan_num)
        if check == True or df_screen is not None:
            mz_intense_list.append(eic_temp)
        else:
            mzlist_uni_temp.remove(mz)
    pbar.finish()
    df_intense = gen_intense_matrix(mz_intense_list, mzlist_uni_temp, time_list, file_path, df_screen)
    print('\nIntensity table has been established.\n')
    return mzlist_uni_temp, mz_intense_list, df_intense


# 生成信号统计矩阵，并保存到output文件夹
def get_signals_matrix(msdata, pola, scan_index, mzlist_uni, output_path, mz_Da, df_screen=None):
    filename = msdata.data[0].filename.split('.')[0]
    file_path = output_path + '\\' + filename + '_' + pola + '_signals' + '.csv'
    total = len(mzlist_uni)
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    extract_data = []
    print('\nExtracting %s signals and backgrounds matrix.\n' % (dic_pola[pola]))
    for index, mz in enumerate(mzlist_uni):
        pbar.update(10 * index + 1)
        extract_data.append(ss.Stat_signal(scan_index[0], scan_index[-1], msdata.eic(mz, mz_Da)))
    pbar.finish()
    # 生成呼气信号平均值&标准偏差矩阵
    ss.batch_statistic(mzlist_uni, extract_data, file_path, df_screen)
    print('\nSignals and backgrounds matrix has been established.\n')


# 导入target ion list（目标化合物离子列表）
def read_target_list(list_path):
    df_scr = pd.read_csv(list_path)
    dic_scr = {'+': df_scr[df_scr.Polarity == 'POS'], '-': df_scr[df_scr.Polarity == 'NEG']}
    return dic_scr



# 根据生成的信号区间，提取一个样品中的全部mz的对应信号
def extract_all_mz(output_path, dic_msdata, pola, dic_id, target_mz_df={},
                   intensity_thre=1e5, mz_Da=0.001, mz_ppm=2, min_scan_num=5,
                   t_start=0, t_end=0, ext_matrix=True, ext_signal=True):
    '''
    Description
        Extract time-intensity matrix and signal-background matrix.
        Step 1 Extract time-intensity matrix.
        Step 2 Extract signal-background matrix.

    Parameters
    :param output_path: str
        Output file path and name.
    :param dic_msdata: dic
        Mass spectra data set.
    :param pola: dic
        Polarity ('+'/'-'/{'+','-'}).
    :param dic_id: dic
        Signal and background id intervals.
    :param target_mz_df: dic
        Targeted ion list. Default is {}.
    :param intensity_thre: int
        Minimum intensity of the highest data point in the ion current. Default
        is 1e5.
    :param mz_Da: float
        Maximum m/z difference of data points in an extract ion current. This
        parameter is work with mz_ppm. The bigger one will be adopted. Default
        is 0.001.
    :param mz_ppm: float
        Maximum m/z difference of data points in an extract ion current. This
        parameter is work with mz_Da. The bigger one will be adopted. Default
        is 2.
    :param min_scan_num: int
        Minimum scan num over which the same ion must be observed in order to
        be recognized as a ion current. Default is 5.
    :param t_start: float
        Start time (min). Default is 0.
    :param t_end: float
        End time (min). Default is 0.
    :param ext_matrix: bool
        Generate time-intensity matrix or not. Default is True.
    :param ext_signal: bool
        Generate signal-background matrix. or not. Default is True.

    Returns
        Time-intensity matrix and signal-background matrix.
        :type: dic
        Extracted m/z list.
        :type: dic
    '''
    dic_df_in = {}
    dic_mzlist = {} if target_mz_df == {} else {'+':target_mz_df['+']['m/z'], '-':target_mz_df['-']['m/z']}
    for i in pola:
        msdata = dic_msdata.get(i)
        target_mz = target_mz_df.get(i, None)
        if t_end == 0:
            t_end = msdata.time[-1]
        # 将保留时间转换为scan id
        scan_index = np.where((msdata.time>=t_start) & (msdata.time<=t_end))
        time_list = msdata.time[scan_index]
        # 判断是否提交了target mz list
        if dic_mzlist.get(i, None) is None:
            mz_unique = get_unique_mzlist(msdata, scan_index[0], intensity_thre, mz_Da, mz_ppm)
        else:
            mz_unique = dic_mzlist[i].reset_index(drop=True)
        if ext_matrix:
            # 提取所有m/z的离子强度信息
            dic_mzlist[i], intense_list, df_intense = get_intense_matrix(msdata, i, time_list, mz_unique, output_path,
                                                                         mz_Da, min_scan_num, t_start, t_end, target_mz)
        dic_df_in[i], id_extract = df_intense, dic_id[i]
        if ext_signal == True and id_extract != 0:
                if id_extract[0] == []:     # 判断id有无
                    print('No signal has been detected.')
                else:   # 生成呼气信号平均值&标准偏差矩阵
                    get_signals_matrix(msdata, i, id_extract, dic_mzlist[i], output_path, mz_Da, target_mz)
    # 输出离子时间×强度矩阵
    return dic_df_in, dic_mzlist
