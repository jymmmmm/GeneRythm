# -*- coding: utf-8 -*-
"""
@File   ： generate_f_info.py
@Time   ： 2023/10/17 17:06
@Author ： Jia Yiming
"""
import pandas as pd
from scipy.fftpack import fft
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

trajectory_info = pd.read_csv("/Users/jiayiming/Desktop/麦吉尔暑研/summer_intern/human_fetal_immune.csv") #dataset1

a = pd.Index(trajectory_info.loc["path"] ==5)

y = np.array(trajectory_info.iloc[:3000,a])

x = np.array(trajectory_info.loc["time",a])

x.sort()
x_smooth = np.linspace(start=x.min(), stop=x.max(), num=100)
y_smooth = make_interp_spline(x, y.transpose(), k =3)(x_smooth)

yf_smooth = fft(y_smooth.transpose())

yf_info = np.array([abs(yf_smooth[i,:50]) for i in range(yf_smooth.shape[0])])

for i in range(len(yf_info)):
    if yf_info[i].max() == 0:
        continue
    yf_info[i] = yf_info[i]/yf_info[i].max()

y_info = y_smooth.transpose()
for i in range(len(y_info)):
    if y_info[i].max() == 0:
        continue
    y_info[i] = y_info[i]/y_info[i].max()
res = np.array([np.append(y_info[i],yf_info[i]) for i in range(yf_smooth.shape[0])])
np.save('ML/res3.npy', np.array(res[:,:]))
# np.save('ML/res_t.npy', np.array(res[non_zero_index,50:]))
np.save('ML/res_t3.npy', np.array(y_info[:,:]))
np.save('ML/non_zero_index3.npy', np.array(non_zero_index))

