# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 21:32:37 2021

@author: 舒晋
"""

import numpy as np
from properties import temperature
# =============================================================================
# 迭代器函数
# =============================================================================

def iterator(p0,H,Rp,fw,qo,t0,GT,angle):
    """
    Parameters
    ----------
    p0 : float
        初始压力，Pa.
    H : float
        油井深度，m.
    Rp : float
        生产气油比.
    fw : float
        含水率，小数.
    t0 : float
        恒温层温度，℃.
    GT : float
        地温梯度
    angle : float
        井斜角

    Returns
    -------
    p : float
        DESCRIPTION.

    """
    #确定起始点压力及计算段深度和分段节点数
    DZ=100
    n=int(H/DZ)+1
    
    #初设计算段的压力降△P1，并计算下端压力p21=p1+△p1
    p=[0]*n#压力分布列表
    p2=[10**6]*n#迭代辅助压力分布列表
    pa = [0] * (n - 1)#平均压力分布列表
    p[0]=p0
    p2[0]=p[0]
    
    #井筒温度分布(近似看作线性变化，且将地层温度看作井筒最下端温度)
    # t=[i for i in np.linspace(25,t_bottom,n)]
    # ta=[0]*(n-1)
    t=[]
    for i in range(n):
        L=DZ*i
        
        i=temperature(fw,qo,GT,t0,H,L)
        t.append(i)
    ta=[0]*(n-1)
    
    #开始迭代
    for i in range(1,n):
        while np.abs(p2[i]-p[i])>0.1:
            p[i]=p2[i]
            pa[i-1]=(p[i-1]+p[i])/2
            ta[i-1]=(t[i-1]+t[i])/2
            dp=np.abs(cal(pa[i-1],ta[i-1],fw,Rp,qo,angle)[1])*DZ
            p2[i]=p2[i-1]+dp
        
        # print('第{}段:'.format(i-1))
        # print('p2({})='.format(i-1),p2[i-1],'p({})='.format(i-1),p[i-1])
        # print('压降:{}MPa'.format(round(dp/10**6,3)))
        # print('平均压力:{}MPa'.format(round(pa[i-1]/10**6,3)))
        # print('平均温度:{}℃'.format(round(ta[i-1],3)))
        # print('流型:{}'.format(cal(pa[i-1],ta[i-1],fw,Rp,qo)[0]))   
        # print('----------------------------') 
        
    return p

# =============================================================================
# beggs-brill方法计算函数
# =============================================================================

def cal(p,t,fw,Rp,qo,angle):
    """
    Parameters
    ----------
    p : float
        压力，Pa.
    t : float
        温度，℃.
    fw : float
        含水率，小数.
    Rp : float
        生产气油比.

    Returns
    -------
    res : list
        beggsbrill方法所得参数(其中1号位参数dp是压力梯度，单位：Pa/m).

    """
    fu=fluid(p,t,Rp,fw)
    pm = pressuremethod()
    
    p0=0.101#标况下大气压
    t0=293#标况下温度
    pm.Ql = qo*fu.Bo/86400#液体流量
    pm.Pa = p#平均压力
    pm.rL = fu.rol#液体密度
    pm.IDT = 0.062#油管半径
    pm.Qg = p0*t*fu.Z*qo*(Rp-fu.Rs)/86400/p/t0#气体流量
    pm.angle = angle#管段倾斜角度
    pm.rg = 40.05#气体密度
    pm.Sm = fu.sigmal#1.83*10**-2#表面张力
    pm.ug = fu.ug#1#气体粘度
    pm.uL = fu.ul#9#液体粘度
    # print('液体流量=',pm.Ql,'平均压力=',pm.Pa,'液体密度=',pm.rL,'表面张力',pm.Sm)
    # print('气体流量=',pm.Qg,'气体密度=',pm.rg,'气体粘度=',pm.ug,'液体粘度',pm.uL)
    res = pm.beggsbrill()
    return res
    
    
        
if __name__=='__main__':
    from beggsbrill import pressuremethod
    from properties import fluid 
    #基础所需数据
    qo=30.8
    fw=9.1/100
    Rp=425
    p0=45*10**6
    H=3500
    angle=0
    
    t0=25
    GT=3
    
    # p=iterator(p0,H,Rp,fw,qo,t_bottom)
    p=iterator(p0,H,Rp,fw,qo,t0,GT,angle)
    
else:
    from beggsbrill import pressuremethod
    from properties import fluid   