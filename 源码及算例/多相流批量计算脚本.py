# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 11:02:05 2021

@author: 舒晋
"""

from beggsbrill import pressuremethod
from properties import fluid
from iterator2 import iterator 
import pandas as pd
import numpy as np
# =============================================================================
# 数据导入模块
# =============================================================================
#总数据
data0=pd.read_excel('多相流计算案例.xlsx')#,encoding='gb18030')
def pwf(data):
    number=data.shape[0]
    #油压
    pt=data.iloc[:,3]*10**6
    #日液
    ql=data.iloc[:,4]
    #日油
    qo=data.iloc[:,5]
    #含水
    fw=data.iloc[:,6]/100
    #生产气油比
    Rp=data.iloc[:,7]
    print(pt,ql,qo,fw,Rp)
    pwf=[]
    #基础所需数据
    for i in range(number):
        H=1000
        # t_bottom=136
        # print(p0,H,Rp,fw,ql,t_bottom)
        try:
            t0=25
            GT=3
            angle=0
            pwf.append(iterator(pt[i],H,int(Rp[i]),fw[i],ql[i],t0,GT,angle)[-1]/10**6)
        except:
            print('第{}行数据故障'.format(i))
            pwf.append(np.nan)
    return pwf

if __name__=='__main__':
    data0_pwf=pwf(data0)
