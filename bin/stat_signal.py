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
def batch_statistic(mzlist, signal_data, output_path, df_screen=None):
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
    column_name = ['m/z (exp)', 'Signals', 'Signals_sd', 'Backgrounds', 'Backgrounds_sd', 'Signals_detail', 'Backgrounds_detail']
    dfresult = pd.DataFrame(dict(zip(column_name, column)), columns = column_name)
    if df_screen is not None:
        dfresult = pd.merge(df_screen.reset_index(drop=True), dfresult.iloc[:,1:], left_index=True, right_index=True)
    dfresult.to_csv(output_path, index=False)


