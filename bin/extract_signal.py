import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import peakutils


# 离子流相关操作
class Ioncurrent:

    def __init__(self, mspos, msneg, tic_p, tic_n, time_p, time_n):
        self.pos = mspos
        self.neg = msneg
        self.tic_p = tic_p
        self.tic_n = tic_n
        self.time_p = time_p
        self.time_n = time_n
        self.scan_p = len(mspos)
        self.scan_n = len(msneg)
#        self.dic_time_p = pd.Series(time_p)
#        self.dic_time_n = pd.Series(time_n)

    # 提取离子流
    def eic_p(self, target_mass=0, tole=0.001, start=0, end=5, plot=False, t_start=0, t_end=0):
        scan_start, scan_end = time_to_scannum(self, '+', t_start, t_end)
        if target_mass is 0:
            eic = np.array([intense for i, intense in enumerate(self.tic_p) if i in range(scan_start, scan_end)])
            if plot == True:
                plot_EIC(self.pos, self.time_p, target_mass, 1000,
                         start=self.time_p[scan_start], end=self.time_p[scan_end-1], title='TIC (+)')
        else:
            eic_with_zero = np.array(extract_ion_current(self.pos, target_mass, tole, scan_start, scan_end))
            # 将0强度点改成非0数字，防止除零错误出现
            eic = np.where(eic_with_zero == 0, 0.01, eic_with_zero)
            if plot == True:
                plot_EIC(self.pos, self.time_p, target_mass, tole,
                         start=self.time_p[scan_start], end=self.time_p[scan_end-1])
        return eic

    def eic_n(self, target_mass=0, tole=0.001, start=0, end=5, plot=False, t_start=0, t_end=0):
        scan_start, scan_end = time_to_scannum(self, '-', t_start, t_end)
        if target_mass is 0:
            eic = np.array([intense for i, intense in enumerate(self.tic_n) if i in range(scan_start, scan_end)])
            if plot == True:
                plot_EIC(self.neg, self.time_n, target_mass, 1000,
                         start=self.time_p[scan_start], end=self.time_p[scan_end-1], title='TIC (-)')
        else:
            eic_with_zero = np.array(extract_ion_current(self.neg, target_mass, tole, scan_start, scan_end))
            # 将0强度点改成非0数字，防止除零错误出现
            eic = np.where(eic_with_zero == 0, 0.01, eic_with_zero)
            if plot == True:
                plot_EIC(self.neg, self.time_n, target_mass, tole,
                         start=self.time_p[scan_start], end=self.time_p[scan_end-1])
        return eic


# 给定某一个时间区间，提取区间内数据
def time_to_scannum(msdata, pola, start_in_min, end_in_min):
    if pola == '+':
        if end_in_min == 0:
            end_in_min = msdata.time_p[-1]
        id_list = [i for i, t in enumerate(msdata.time_p) if
                   start_in_min <= t <= end_in_min]
    if pola == '-':
        if end_in_min == 0:
            end_in_min = msdata.time_n[-1]
        id_list = [i for i, t in enumerate(msdata.time_n) if
                   start_in_min <= t <= end_in_min]
    start_scan = id_list[0]
    end_scan = id_list[-1] + 1
    return start_scan, end_scan


# 拆分ms1正负离子数据集，获取polarity（极性）信息
def sep_polarity(ms):
    ms_pos = []
    time_pos= []
    ms_neg = []
    time_neg = []
    tic_pos = []
    tic_neg = []
    pola = set([i.polarity for i in ms[0:2]])
    if '+' in pola:
        ms_pos = [i for i in ms if i.polarity == '+']
        time_pos = np.array([round((ms_pos[i].rt_in_second / 60), 5) for i in range(0, len(ms_pos))])
        tic_pos = np.array([float(ms_pos[i].tic) for i in range(0, len(ms_pos))])
    if '-' in pola:
        ms_neg = [i for i in ms if i.polarity == '-']
        time_neg = np.array([round((ms_neg[i].rt_in_second / 60), 5) for i in range(0, len(ms_neg))])
        tic_neg = np.array([float(ms_neg[i].tic) for i in range(0, len(ms_neg))])
    output = Ioncurrent(ms_pos, ms_neg, tic_pos, tic_neg, time_pos, time_neg)
    print('Polarity: %s' % [i for i in pola])
    return output, pola


# 提取离子流
def extract_ion_current(mzfile, target_mass, tolerance=0.001, start_scan=0, end_scan=0):
    target_intense = []
    if end_scan == 0:
        end_scan = len(mzfile)
    for i in range(start_scan, end_scan):
        spec_ms1 = mzfile[i].peaks
        # 获取目标离子窗口中强度最高的离子
        try:
            intense_max = int(max([i[1] for i in spec_ms1 if (target_mass - tolerance)
                          < i[0] < (target_mass + tolerance)]))
        except:
            intense_max = 0
        target_intense.append(intense_max)
    return np.array(target_intense)


# 绘制离子流
def plot_EIC(ms1data, time, target_mass, tolerance=0.001, start=0, end=5, title='TIC'):
    intense = extract_ion_current(ms1data, target_mass, tolerance)
    plt.figure(figsize=(10.5, 5))
    plt.plot(time, intense)
    if target_mass is 0:
        plt.title(title)
    else:
        plt.title('m/z %s ' % (target_mass))
    plt.xlabel("time (min)")
    plt.ylabel("Intensity")
    plt.xlim(start, end)

# 绘制baseline
def plot_baseline(time, Y, baseline_para):
    # time_series = Y
    time_series = np.asarray(Y)
    baseline = peakutils.baseline(time_series, baseline_para)
    plt.plot(time, baseline, time, Y)


# 获得与分割线交叉的点
def cal_cross_index(intense, spec_num):
    index = [0]
    start = False
    start_index_temp = 0
    for i in range(intense.size-1):
        if i > start_index_temp + spec_num:
            if (intense[i] > 0 and intense[i+1] < 0):
                if start is False:
                    index.append(i+1)
                else:
                    index.append(i)
                start = not start
                start_index_temp = i
            if (intense[i] < 0 and intense[i+1] > 0):
                if start is False:
                    index.append(i+1)
                else:
                    index.append(i)
                start = not start
                start_index_temp = i
    index.append(intense.size-1)
    return index

# 获得与分割线交叉的点
def cal_cross_index_ori(intense):
    index = [0]
    start = False
    for i in range(intense.size-1):
        if i > 10:
            if (intense[i] > 0 and intense[i+1] < 0):
                if start is False:
                    index.append(i+1)
                else:
                    index.append(i)
                start = not start
            if (intense[i] < 0 and intense[i+1] > 0):
                if start is False:
                    index.append(i+1)
                else:
                    index.append(i)
                start = not start
    index.append(intense.size-1)
    index_dere = list(set(index))
    index_dere.sort(key=index.index)
    return index_dere


# 绘制基线、分割线&有效数据识别点
def plot_detect_result(intense, time, smooth, baseline_normal,
                       baseline_reversed, midline, title='TIC', pola=''):
    plt.figure(figsize=(10.5, 5))
    if title is 0:
        plt.title('TIC (%s)' % pola)
    else:
        plt.title('m/z %s (%s)' % (title, pola))
    p1, = plt.plot(time, intense, color='dimgrey', linestyle='--')
    p2, = plt.plot(time, smooth, color='firebrick', linestyle='-')
    p3, = plt.plot(time, baseline_normal, color='royalblue', linestyle='--')
    p4, = plt.plot(time, -baseline_reversed, color='royalblue', linestyle='--')
    p5, = plt.plot(time, midline, color='royalblue', linestyle='-')
    plt.legend([p1, p2, p3, p4, p5], ['intensity', 'intensity_smooth',
               'baseline_1', 'baseline_2', 'cutting line'], loc='upper right')


# 获取底部和顶部基线，取平均获得中间分割线
def get_baseline(intense, time, smooth_num=7, smooth_poly=1,
                 baseline_poly=1, signal_weight=0.9, title='TIC', pola=''):
    # Savitzky–Golay滤波平滑
    # 获得底部baseline
    smooth = signal.savgol_filter(intense, smooth_num,
                                  smooth_poly, mode='nearest')
    baseline_normal = peakutils.baseline(smooth, baseline_poly)
    # 获得顶部baseline
    intense_reverse = [-i for i in intense]
    smooth_reversed = signal.savgol_filter(intense_reverse, smooth_num, smooth_poly, mode='nearest')
    baseline_reversed = peakutils.baseline(smooth_reversed, baseline_poly)
    # 获得切割线，底部基线权重较大
    midline = (baseline_normal - signal_weight * baseline_reversed) / 2
    # 绘图
    plot_detect_result(intense, time, smooth, baseline_normal,
                       baseline_reversed, midline, title, pola)
    return smooth, midline


# 获得有效信号index
def get_key_index(msdata, pola, target_mass=0, smooth_num=7, smooth_poly=1, baseline_poly=1,
                  index_window=3, spec_num=8, signal_weight=0.9, start_in_min=0, end_in_min=0):
#    start_scan, end_scan = time_to_scannum(msdata, '+', start_time, end_time)
    if pola == '+':
        if end_in_min == 0:
            end_in_min = msdata.time_p[-1]
        intense = msdata.eic_p(target_mass, t_start=start_in_min, t_end=end_in_min)
        time = np.array([i for i in msdata.time_p if start_in_min <= i <= end_in_min])
    if pola == '-':
        if end_in_min == 0:
            end_in_min = msdata.time_n[-1]
        intense = msdata.eic_n(target_mass, t_start=start_in_min, t_end=end_in_min)
        time = np.array([i for i in msdata.time_n if start_in_min <= i <= end_in_min])
    smooth_sig, midline = get_baseline(intense, time, smooth_num, smooth_poly,
                                       baseline_poly, signal_weight, target_mass, pola)
    # 获取smooth切割点
    smooth_adjust = smooth_sig - midline
    smooth_index = cal_cross_index(smooth_adjust, spec_num)
    # 获取原始数据切割点
    intense_adjust = intense - midline
    ori_index = cal_cross_index_ori(intense_adjust)
    # 校正切割点id
    index = []
    for i in ori_index:
        for u in smooth_index:
            if (i - index_window) <= u <= (i + index_window):
                index.append(i)
    # index.append(intense.size-1)
    background_index = [index[i:i+2] for i in range(0, len(smooth_index), 2)]
    peak_index = [index[i:i+2] for i in range(1, len(smooth_index)-1, 2)]
    peak_time = [tuple(time[index[i:i+2]]) for i in range(1, len(smooth_index)-1, 2)]
    plt.scatter(time[index[1:-1]], intense[index[1:-1]], color='r')
    bg_id, bg_time = refine_signal(time, background_index)
    sig_id, sig_time = refine_signal(time, peak_index)
    plot_signal_level(intense, sig_id, sig_time)
    plt.show()
    print('A total of %s breath signals and %s background signals were detected.'
          % (len(peak_index), len(background_index)))
    print('Breath signals are between %s mins (%s scans).' % (peak_time, index[1:-1]))
    return bg_id, sig_id


# 获取中间80%有效信号
def refine_signal(time, signal_id, tail=0.1):
    refine_id = []
    time_temp = []
    # 获取每段信号中间80%的信号强度
    for i in signal_id:
        if len(i) == 2:
            rm_num = int(tail * (i[1] - i[0] + 1)) + 1
            start_id = i[0] + rm_num
            end_id = i[1] - rm_num
            refine_id.append([start_id, end_id])
            time_temp.append(time[start_id:end_id])
    return refine_id, time_temp


def plot_signal_level(intense, index, time):
    intense_temp = []
    for i, ind in enumerate(index):
        intense_temp = intense[ind[0]:ind[1]]
        ave_temp = np.average(intense_temp)
        x = time[i]
        y = ave_temp + x * 0
        plt.plot(x, y, color='firebrick', linewidth=5, linestyle='--')

