# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 16:23:47 2021

@author: 舒晋
"""
import math
import numpy as np

# =============================================================================
# 流体性质参数以及流动参数
# =============================================================================
class fluid:
    def __init__(self,P,t,Rp,fw):
        #初始值
        gamao=0.85
        gamag=0.55
        # fw=0.8
        # Rp=0.5
        self.row=1
        self.ug=1.22*10**-5
        
        #计算值
        self.yAPI=self.yAPI(gamao)
        self.Rs=self.Rs(self.yAPI, t, gamag, gamao, P, Rp)
        self.Bo=self.Bo(self.Rs, gamag, gamao, t)
        self.roo=self.roo(gamao, gamag, self.Rs, self.Bo)
        self.rol=self.rol(fw, self.roo, self.row)
        self.uo=self.uo(self.yAPI, t, self.Rs)
        self.uw=self.uw(t)
        self.ul=self.ul(fw, self.uo, self.uw)
        self.sigmaog=self.sigmaog(self.yAPI,P,t)
        self.sigmawg=self.sigmawg(P,t)
        self.sigmal=self.sigmal(fw,self.sigmaog,self.sigmawg)
        self.Z=self.Z(gamag,P,t)
        self.rog=self.rog(gamag,self.Z,P,t)
        
    def roo(self,gamao,gamag,Rs,Bo):#原油密度
        """
        Parameters
        ----------
        gamao : float
            地面条件下的原油相对密度.
        gamag : float
            地面条件下的气体相对密度.
        Rs : float
            溶解气油比.
        Bo : float
            当前温压条件下的原油体积系数.
    
        Returns
        -------
        roo : float
            原油密度,kg/m^3.
    
        """
        roo=1000*(gamao+1.206*10**-3*Rs*gamag)/Bo
        return roo
    def yAPI(self,gamao):#API度（相对密度）
        """
        Parameters
        ----------
        gamao : float
            地面条件下的原油相对密度.
    
        Returns
        -------
        yAPI : float
            原有的API度.
    
        """
        yAPI=141.5/gamao-131.5
        return yAPI
    def Bo(self,Rs,gamag,gamao,t):
        """
        Parameters
        ----------
        Rs : float
            溶解气油比.
        gamag : float
            地面条件下气体相对密度.
        gamao :float
            地面条件下原油相对密度.
        t :float
            温度，℃.

        Returns
        -------
        Bo : float
            原油体积系数.

        """
        F=5.615*Rs*math.sqrt(gamag/gamao)+2.25*t+40#这里的T在教程里是小写
        Bo=0.972+0.000147*F**1.175
        return Bo
    def Rs(self,yAPI,t,gamag,gamao,P,Rp):#溶解气油比
        """
        Parameters
        ----------
        yAPI : float
            API度（原油相对密度）.
        t : TYPE
            温度，℃.
        gamag : float
            地面条件下气体相对密度.
        gamao : float
            地面条件下原油相对密度.
        P : float
            压力，Pa.

        Returns
        -------
        Rs : float
            溶解气油比.

        """
        if yAPI<15:
            # print('API<15')
            A=0.0125*yAPI-0.00091*(1.8*t+32)
            Rs=0.17812*gamag*(8.0558*P*10**A)**1.2048
        if yAPI>=15:
            if yAPI>=38.3:
                # print('API>=38.3')
                mo=10**(0.6631*math.log((1346/(yAPI-2.1)),math.e))
            if yAPI<38.3:
                # print('API<38.3')
                mo=(61.933-yAPI)/0.0943
            xg=8.0558*P*gamag/(10**5*(t+273.15))
            if xg>=3.448:
                yng=0.3531*np.log(xg/0.5967)
            if 0.7<xg<3.448:
                yng=0.2401*np.log(xg/0.27)
            if xg<=0.7:
                # print('xg',xg)
                yng=0.1236*np.log(xg/0.1223)
            Rs=23650*gamao/mo*yng/(1-yng)
        if Rs>Rp:
            Rs=Rp
        return Rs
    def rol(self,fw,roo,row):#液体密度
        """
        Parameters
        ----------
        fw : float
            含水率，小数.
        roo : float
            油密度，kg/^3.
        row : float
            水密度，kg/^3.

        Returns
        -------
        rol : float
            混合液体密度,kg/m^3.

        """
        rol=roo*(1-fw)+row*fw
        return rol
    def uo(self,yAPI,t,Rs):#原油粘度
        """
        Parameters
        ----------
        yAPI : float
            API度，相对密度.
        t : float
            温度，℃.
        Rs : float
            溶解气油比.

        Returns
        -------
        uo : float
            粘度，Pa.s.

        """
        z=3.0324-0.02023*yAPI
        y=10**z
        x=y*(32+1.8*t)**-1.163

        uod=(np.power(10,x)-1)/1000
        B=5.44*(5.615*Rs+150)**-0.338
        A=10.715*(5.615*Rs+100)**-0.515
        uo=A*(1000*uod)**B/1000
        return uo
    def uw(self,t):#水粘度
        """
        Parameters
        ----------
        t : float
            温度，℃.

        Returns
        -------
        uw : float
            水的粘度，pa.s.

        """
        uw=(math.e**(1.003-1.479*10**-2*(32+1.8*t)+1.982*10**-5*(32+1.8*t)**2))/1000
        return uw
    def ul(self,fw,uo,uw):#流体粘度
        """
        Parameters
        ----------
        fw : float
            含水率，小数.
        uo : float
            油的粘度,Pa.s.
        uw : float
            水的粘度,Pa.s.

        Returns
        -------
        ul : float
            液体粘度,Pa.s.

        """
        ul=uo*(1-fw)+uw*fw
        return ul
    def sigmaog(self,yAPI,P,t):#油气表面张力
        """
        Parameters
        ----------
        P : float
            压力，MPa.
        t : float
            温度，℃.

        Returns
        -------
        sigmaog : float
            油气表面张力，N/m.

        """
        sigmaog=(42.4-0.047*(1.8*t+32)-0.267*yAPI)*math.e**(-0.015*10**-7*P)/1000
        return sigmaog
    def sigmawg(self,P,t):#水气表面张力
        """
        Parameters
        ----------
        P : float
            压力，Pa.
        t : float
            温度，℃.

        Returns
        -------
        sigmawg : float
            水气表面张力，N/m.

        """
        
        sigma2333=76*math.e**(-3.62575*10**-7*P)/1000
        sigma13778=(52.5-8.7018*10**-7*P)/1000
        sigmawg=((248-1.8*t)/206*(sigma2333-sigma13778)+sigma13778)
        return sigmawg
    def sigmal(self,fw,sigmaog,sigmawg):#油水混合物和天然气的表面张力
        """
        Parameters
        ----------
        fw : float
            含水率，小数.
        sigmaog : float
            油气表面张力.
        sigmawg : float
            水气表面张力.

        Returns
        -------
        sigmal : float
            油水混合物和天然气的表面张力.

        """
        sigmal=sigmaog*(1-fw)+sigmawg*fw
        return sigmal
    def Z(self,gamag,P,t):#天然气压缩因子
        """
        Parameters
        ----------
        gamag : float
            地面条件下气体相对密度.
        P : float
            压力,Pa.
        t : float
            温度，℃.

        Returns
        -------
        Z : float
            天然气压缩因子.

        """
        Tc=9.22+176.67*gamag
        if gamag>=0.7:
            Pc=10**6*(4.88-0.39*gamag)
        if gamag<0.7:
            Pc=10**6*(4.78-0.25*gamag)
        Tr=(273.15+t)/Tc
        Pr=P/Pc
        Z=1
        PR=0.27*Pr/Z/Tr
        Z=1+(0.31506+(-1.0467/1.55)+(-0.5783/1.55**3))*PR+(0.5353-0.023/1.55)*PR**2+0.6815/1.55**3*PR**2
        return Z
    def rog(self,gamag,Z,P,t):#天然气密度
        """
        Parameters
        ----------
        gamag : float
            地面条件下气体相对密度.
        Z : float
            天然气压缩因子.
        P : float
            压力，Pa.
        t : float
            温度，℃.

        Returns
        -------
        rog : float
            天然气密度，Kg/m^3.

        """
        rog=3.4844*10**-3*gamag*P/Z/(t+273.15)
        return rog

# def temperature(roo,row,fw,QL,GT,t0,H,L):
def temperature(fw,QL,GT,t0,H,L):
    """
    Parameters
    ----------
    roo : float
        原油密度，kg/m^3.
    row : float
        水的密度,kg/m^3.
    fw : float
        含水率，小数.
    QL : float
        油井产液量，t/d.
    GT : float
        地温梯度,/m.
    t0 : float
        恒温层温度，℃.
    H : float
        油层中深，m.
    L : float
        井筒中任一点深度.

    Returns
    -------
    t : float
        DESCRIPTION.

    """
    QL=QL*1.1357
    row=1000
    roo=800
    FW=row*fw/(row*fw+roo+(1-fw))
    G=QL*1000/24
    KP=1/(1.1573+5.4246*np.power(math.e,-G/1000))
    BATA=2*np.pi*KP/G/(1+FW)
    tr=t0+GT*L/100
    t=t0+(tr-t0)/(BATA*H)*(BATA*L+1-np.power(math.e,-BATA*(H-L)))
    return t

if __name__=='__main__':
    P,t,Rp,fw=5*10**6,1,400,0.5
    
    fluid=fluid(P,t,Rp,fw)#实例化测试
    roo=fluid.roo
    row=fluid.row
    yAPI=fluid.yAPI
    Bo=fluid.Bo
    Rs=fluid.Rs
    rol=fluid.rol
    uo=fluid.uo
    uw=fluid.uw
    ul=fluid.ul
    sigmaog=fluid.sigmaog
    sigmawg=fluid.sigmawg
    sigmal=fluid.sigmal
    Z=fluid.Z
    rog=fluid.rog
    
    QL=50
    GT=3
    t0=50
    H=2000
    L=10
    tem=temperature(fw,QL,GT,t0,H,L)