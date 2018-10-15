import numpy as np
import pandas as pd


class Stat_signal:
    def __init__(self, signal_id, background_id, intense):
        self.sig_id = signal_id
        self.bg_id = background_id
        self.sig_data = [intense[i[0]:(i[1]+1)] for i in signal_id]
        self.bg_data = [intense[i[0]:(i[1]+1)] for i in background_id]
        self.sig_detail = [round(np.average(i), 2) for i in self.sig_data]
        self.bg_detail = [round(np.average(i), 2) for i in self.bg_data]
        self.sig = round(np.average(self.sig_detail), 2)
        self.bg = round(np.average(self.bg_detail), 2)
        self.sig_sd = round(np.std(self.sig_detail), 2)
        self.bg_sd = round(np.std(self.bg_detail), 2)


# 批量处理
def batch_statistic(mzlist, time, mz_intense_list, signal_data, output_path):
    mzlist_round = [round(i, 5) for i in mzlist]
    signals = []
    signals_detail = []
    signals_sd = []
    backgrounds = []
    backgrounds_detail = []
    backgrounds_sd = []
    for i in signal_data:
        signals.append(i.sig)
        signals_detail.append(i.sig_detail)
        signals_sd.append(i.sig_sd)
        backgrounds.append(i.bg)
        backgrounds_detail.append(i.bg_detail)
        backgrounds_sd.append(i.bg_sd)
#    signals = [i.sig for i in signal_data]
#    signals_detail = [i.sig_detail for i in signal_data]
#    signals_sd = [i.sig_sd for i in signal_data]
#    backgrounds = [i.bg for i in signal_data]
#    backgrounds_detail = [i.bg_detail for i in signal_data]
#    backgrounds_sd = [i.bg_sd for i in signal_data]
    column = [mzlist_round, signals, signals_sd, backgrounds, backgrounds_sd, signals_detail, backgrounds_detail]
    column_name = ['m/z', 'Signals', 'Signals_sd', 'Backgrounds', 'Backgrounds_sd', 'Signals_detail', 'Backgrounds_detail']
    dfresult = pd.DataFrame(dict(zip(column_name, column)), columns = column_name)
    dfresult.to_csv(output_path, index=False)



#    
#    # 提取呼气信号的时间区间
#    def sig_time(self):
#        time = []
#        for i in self.signal_id:
#            time_singal = [self.time[i[0]], self.time[i[1]]]
#            time.append(time_singal)
#        return time
#    # 提取一个呼气信号强度
#    def sig_intense(self):
#        sig_list = []
#        for i in self.signal_id:
#            intense_temp = self.intense[i[0]:(i[1]+1)]
#            sig_list.append(intense_temp)
#        return sig_list
#
#
#intense = ms1.tic_p
#signal_id = sig_p_id
#
#sig_list = [self.intense[i[0]:(i[1]+1)] for i in self.signal_id]
#
#sig_list = [intense[i[0]:(i[1]+1)] for i in signal_id]
#sig = [np.average(i) for i in sig_list]
#sig_ave = np.average(sig)
#
#    def bg_intense(self):
#        bg_list = []
#        for i in self.signal_id:
#            intense_temp = self.intense[i[0]:(i[1]+1)]
#            bg_list.append(intense_temp)
#        return bg_list
#    
#    def stat(self):
#        bg_list = []
#        for i in self.signal_id:
#            intense_temp = self.intense[i[0]:(i[1]+1)]
#            bg_list.append(intense_temp)
#        return bg_list
#
#    def sig(self):
#        sig_single = []
#        for i in sig_intense(self):
#            ave_temp = np.average(i)
#            sig_single.append(ave_temp)
#        return sig_single
#
#
#
#
#
#
#
#
#time = ms1.time_p
#time_list = []
#for i in sig_p_id:
#    time_temp = [time[i[0]], time[i[1]]]
#    time_list.append(time_temp)



