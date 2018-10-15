
import os
import time
import pandas as pd
from interval import Interval
from progressbar import *
import extract_signal as es
import stat_signal as ss


# 给出数据导入和导出的文件（夹）路径
def gen_workspace(filename):
    data_input_path = r'..\data_input' + '\\' + filename
    data_output_path = r'..\data_output'
    try:
        os.makedirs(r'..\data_output')
        print('The output folder has been established.')
    except Exception as e:
        print('The output folder is already exist.')
    return data_input_path, data_output_path


## 批量提取离子流
# 获取ms数据中强度大于设定阈值的mzlist
def get_unique_mzlist(msdata, intensity_threshold=5e4, mz_Da = 0.002, mz_ppm=2, start_id=0, end_id=0):
    start_time = time.time()
    mzlist_init = []
    mzlist_all = []
    # 获取全部图谱中的mz
    for i in range(start_id, end_id):
        mzlist = msdata[i].peaks
        mzlist_init = [mzlist[i][0] for i in range(1,len(mzlist)) if
                       mzlist[i-1][1] < mzlist[i][1] > mzlist[i+1][1]
                       and mzlist[i][1] > intensity_threshold]
        mzlist_all += mzlist_init
    mzlist_all.sort()
    # print(len(mzlist_all))
    # 排序后，获取mz_Da和mz_ppm窗口以外的mz
    mzlist_unique = [mzlist_all[i] for i in range(1, len(mzlist_all)) if
                     mzlist_all[i] not in Interval(mzlist_all[i-1] * (1 - mz_ppm * 1e-6),
                     mzlist_all[i-1] * (1 + mz_ppm * 1e-6)) and
                     mzlist_all[i] not in Interval(mzlist_all[i-1] - mz_Da,
                     mzlist_all[i-1] + mz_Da)]
    end_time = time.time()
    print('%s unique m/z values were extracted. Running time: %.2f s'
          % (len(mzlist_unique), (end_time - start_time)))
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
def check_all_sig(msdata, pola, smooth_num=7, smooth_poly=1, baseline_poly=1,
                      index_window=3, spec_num=8, signal_weight=0.9,
                      t_start=0, t_end=0):
    bg_id_p = 0
    sig_id_p = 0
    bg_id_n = 0
    sig_id_n = 0
    for i in pola:
        if i == '+':
            count = 0
            bg_id_p, sig_id_p = check_sig(msdata, i, smooth_num, smooth_poly, baseline_poly,
                                          index_window, spec_num, signal_weight, count, t_start, t_end)
        if i == '-':
            count = 0
            bg_id_n, sig_id_n = check_sig(msdata, i, smooth_num, smooth_poly, baseline_poly,
                                          index_window, spec_num, signal_weight, count, t_start, t_end)
    id_p = [sig_id_p, bg_id_p]
    id_n = [sig_id_n, bg_id_n]
    return id_p, id_n


# 生成所有质荷比的时间×强度矩阵
def gen_intense_matrix(mz_intense_list, mzlist_uni, time_list, filepath_intense):
    df_mz_intense = pd.DataFrame(mz_intense_list, columns=time_list)
    df_mz_intense['Average'] = df_mz_intense.apply(lambda x: x.mean(), axis=1)
    df_mz_intense.insert(0, 'm/z', mzlist_uni)
    df_mz_intense.to_csv(filepath_intense, index=False)
    return df_mz_intense


# 生成时间强度矩阵，并保存到output文件夹
def get_intense_matrix(msdata, pola, time_list, mzlist_uni, filepath_intense, start_in_min=0, end_in_min=0):
    total = len(mzlist_uni)
    mz_intense_list = []
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    mz_intense_list = []
    if pola == '+':
        print('\nStep 1 Establising positive intensity table.\n')
        for index, mz in enumerate(mzlist_uni):
            pbar.update(10 * index + 1)
            mz_intense_list.append(msdata.eic_p(mz, t_start=start_in_min, t_end=end_in_min))
    if pola == '-':
        print('\nStep 1 Establising negative intensity table.\n')
        for index, mz in enumerate(mzlist_uni):
            pbar.update(10 * index + 1)
            mz_intense_list.append(msdata.eic_n(mz, t_start=start_in_min, t_end=end_in_min))
    pbar.finish()
    df_intense = gen_intense_matrix(mz_intense_list, mzlist_uni, time_list, filepath_intense)
    print('\nIntensity table has been established.\n')
    return mz_intense_list, df_intense


# 生成信号统计矩阵，并保存到output文件夹
def get_signals_matrix(msdata, pola, time_list, id_range, mzlist_uni, mz_intense_list, filepath):
    total = len(mzlist_uni)
    widgets = ['Progress: ',Percentage(), ' ', Bar('#'),' ', Timer(),' ', ETA(),' ']
    pbar = ProgressBar(widgets=widgets, maxval=10*total).start()
    extract_data = []
    if pola == '+':
        print('\nStep 2 Extracting positive signals and backgrounds matrix.\n')
        for index, mz in enumerate(mzlist_uni):
            pbar.update(10 * index + 1)
            extract_data.append(ss.Stat_signal(id_range[0], id_range[1], msdata.eic_p(mz)))
    if pola == '-':
        print('\nStep 2 Extracting negative signals and backgrounds matrix.\n')
        for index, mz in enumerate(mzlist_uni):
            pbar.update(10 * index + 1)
            extract_data.append(ss.Stat_signal(id_range[0], id_range[1], msdata.eic_n(mz)))
    pbar.finish()
        # 生成呼气信号平均值&标准偏差矩阵
    ss.batch_statistic(mzlist_uni, time_list, mz_intense_list, extract_data, filepath)
    print('\nSignals and backgrounds matrix has been established.\n')


# 根据生成的信号区间，提取一个样品中的全部mz的对应信号
def extract_all_mz(output_path, msdata, pola, id_p=0, id_n=0,
                   intensity_thre=1e5, mz_Da=0.001, mz_ppm=2, t_start=0, t_end=0):
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
    filename = msdata.pos[0].filename.split('.')[0]
    df_intense_p = []
    df_intense_n = []
    for i in pola:
        filepath = output_path + '\\' + filename + '_' + i + '_signals' + '.csv'
        filepath_intense = output_path + '\\' + filename + '_' + i + '_intensity_matrix' + '.csv'
        if i == '+':
            if t_end == 0:
                t_end = msdata.time_p[-1]
            # 将保留时间转换为scan id
            start_id, end_id = es.time_to_scannum(msdata, i, t_start, t_end)
            time_list = [i for i in msdata.time_p if t_start <= i <= t_end]
            mzlist_uni = get_unique_mzlist(msdata.pos, intensity_thre, mz_Da, mz_ppm, start_id, end_id)
            # 提取所有m/z的离子强度信息
            mz_intense_list, df_intense_p = get_intense_matrix(msdata, '+', time_list, mzlist_uni,
                                                               filepath_intense, t_start, t_end)
            if id_p != 0:
                # 判断id_p有无
                if id_p[0] == []:
                    print('No signal has been detected.')
                else:
                    # 生成呼气信号平均值&标准偏差矩阵
                    get_signals_matrix(msdata, i, time_list, id_p, mzlist_uni, mz_intense_list, filepath)
        if i == '-':
            if t_end == 0:
                t_end = msdata.time_n[-1]
            # 将保留时间转换为scan id
            start_id, end_id = es.time_to_scannum(msdata, i, t_start, t_end)
            time_list = [i for i in msdata.time_n if t_start <= i <= t_end]
            mzlist_uni = get_unique_mzlist(msdata.neg, intensity_thre, mz_Da, mz_ppm, start_id, end_id)
            # 提取所有m/z的离子强度信息
            mz_intense_list, df_intense_n = get_intense_matrix(msdata, '-', time_list, mzlist_uni,
                                                               filepath_intense, t_start, t_end)
            if id_n != 0:
                # 判断id_p有无
                if id_n[0] == []:
                    print('No signal has been detected.')
                else:
                    # 生成呼气信号平均值&标准偏差矩阵
                    get_signals_matrix(msdata, i, time_list, id_p, mzlist_uni, mz_intense_list, filepath)
    # 输出离子时间×强度矩阵
    return df_intense_p, df_intense_n





