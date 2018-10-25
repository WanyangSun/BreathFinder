import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import peakutils


# 离子流相关操作
class Massdata:

    def __init__(self, ms, polarity):
        self.data = ms
        self.pola = polarity
        self.tic = np.array([float(ms[i].tic) for i in range(0, len(ms))])
        self.time = np.array([round((ms[i].rt_in_second / 60), 5) for i in range(0, len(ms))])
        self.scan = np.array([i for i in range(0, len(ms))])

    # 提取离子流
    def __extract_ion_current(self, target_mass, scan_index, tolerance=0.001):
        target_intense = []
        if target_mass == 0:
            target_intense = self.tic[scan_index]
        else:
            for i in scan_index[0]:
                spec_ms = self.data[i].peaks
                mzs = spec_ms[:, 0]
                intense = spec_ms[:, 1]
                # 获取目标离子窗口中强度最高的离子
                # 获取目标离子范围内的m/z
                index = np.argwhere((mzs > (target_mass - tolerance)) & (mzs < (target_mass + tolerance)))
                # 将0强度点改成非0数字，防止除零错误出现
                intense_max = int(max(intense[index])) if len(index) != 0 else 0.01
                target_intense.append(intense_max)
        return np.array(target_intense)

    # 绘制离子流
    def __plot_eic(self, eic, time, target_mass):
        plt.figure(figsize=(10.5, 5))
        plt.plot(time, eic)
        pola = '+' if self.pola == '+' else '-'
        if target_mass == 0:
            plt.title('TIC (%s)' % pola)
        else:
            plt.title('(%s) m/z %s ' % (self.pola, target_mass))
        plt.xlabel("time (min)")
        plt.ylabel("Intensity")
        plt.xlim(time[0], time[-1])
        plt.show()

    # 提取离子流
    def eic(self, target_mass=0, tole=0.001, plot=False, t_start=0, t_end=0):
        if t_end == 0:
            t_end = self.time[-1]
        scan_index = np.where((self.time>=t_start) & (self.time<=t_end))
        time_list = self.time[scan_index]
        eic = self.__extract_ion_current(target_mass, scan_index, tole)
        if plot == True:
            self.__plot_eic(eic, time_list, target_mass)
        return eic


# 拆分ms1正负离子数据集，获取polarity（极性）信息
def sep_polarity(ms):
    ms_pos = []
    ms_neg = []
    pola = set([i.polarity for i in ms[0:2]])
    if '+' in pola:
        ms_pos = [i for i in ms if i.polarity == '+']
    if '-' in pola:
        ms_neg = [i for i in ms if i.polarity == '-']
    output_pos = Massdata(ms_pos, '+')
    output_neg = Massdata(ms_neg, '-')
    print('Polarity of this data file: %s' % [i for i in pola])
    return output_pos, output_neg, pola


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
    if end_in_min == 0:
        end_in_min = msdata.time[-1]
    intense = msdata.eic(target_mass, t_start=start_in_min, t_end=end_in_min)
    time = np.array([i for i in msdata.time if start_in_min <= i <= end_in_min])
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

