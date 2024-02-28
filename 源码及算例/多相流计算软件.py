# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:20:18 2021

@author: 舒晋
"""

from tkinter import *
from tkinter.ttk import *
from beggsbrill import pressuremethod
from properties import fluid
from iterator2 import iterator
from properties import temperature

class mainform():
    def __init__(self):
        self.root=Tk()
        self.root.title('706多相流计算器(基于Beggs-Brill方法)')
        self.root.geometry('500x470')
        
        self.para()
        self.menu()
        self.layout()
        self.root.mainloop()
    def para(self):
        self.input_name=['油压','井深','生产气油比',
                         '含水率','油井产液量','恒温层温度',
                         '地温梯度','井斜角']
        self.input_value=[DoubleVar(),DoubleVar(),DoubleVar(),
                          DoubleVar(),DoubleVar(),DoubleVar(),
                          DoubleVar(),DoubleVar()]
        self.input_unit=['MPa','m','m^3/m^3',
                         '小数','m^3/d','℃',
                         '℃/100m','°']

        value_list=[45,3500,425,0.091,30.8,25,3,0]
        for i in range(8):
            self.input_value[i].set(value_list[i])
        
        self.output_name=['井底流压']
        self.output_value=[DoubleVar()]
        self.output_unit=['MPa']
    
        
    def menu(self):
        self.menubar=Menu(self.root)
        
        self.menubar.add_command(label='归零',command=self.zero)
        self.menubar.add_command(label='退出',command=self.root.destroy)
        
        self.root.config(menu=self.menubar)
    def layout(self):
        #总框架布局
        self.lf=[[],[]]
        self.lf[0]=LabelFrame(self.root,text='输入区')
        self.lf[1]=LabelFrame(self.root,text='输出区')
        self.btn=[[]]
        self.btn[0]=Button(self.root,text='计算',command=self.cal)
        
        for i in range(2):
            self.lf[i].pack()
        self.btn[0].pack()
        
        
        #输入框布局
        self.input_name_lb=[[],[],[],[],[],[],[],[],[]]
        self.input_value_et=[[],[],[],[],[],[],[],[],[]]
        self.input_unit_lb=[[],[],[],[],[],[],[],[],[]]
        for i in range(8):
            self.input_name_lb[i]=Label(self.lf[0],text=self.input_name[i])
            self.input_value_et[i]=Entry(self.lf[0],textvariable=self.input_value[i])
            self.input_unit_lb[i]=Label(self.lf[0],text=self.input_unit[i])
            
            self.input_name_lb[i].grid(row=i,column=0,padx=10,pady=10,sticky=W)
            self.input_value_et[i].grid(row=i,column=1,padx=10,pady=10,sticky=W)
            self.input_unit_lb[i].grid(row=i,column=2,padx=10,pady=10,sticky=W)
        #输出框布局
        self.output_name_lb=[[],[]]
        self.output_value_et=[[],[]]
        self.output_unit_lb=[[],[]]
        for i in range(1):
            self.output_name_lb[i]=Label(self.lf[1],text=self.output_name[i])
            self.output_value_et[i]=Entry(self.lf[1],textvariable=self.output_value[i])
            self.output_unit_lb[i]=Label(self.lf[1],text=self.output_unit[i])
            
            self.output_name_lb[i].grid(row=i,column=0,padx=10,pady=10,sticky=W)
            self.output_value_et[i].grid(row=i,column=1,padx=10,pady=10,sticky=W)
            self.output_unit_lb[i].grid(row=i,column=2,padx=10,pady=10,sticky=W)
    
    def zero(self):
        value_list=[45,3500,425,0.091,30.8,25,3,0]
        for i in range(8):
            self.input_value[i].set(value_list[i])
        self.output_value[0].set(0)
    def cal(self):
        pt=self.input_value[0].get()
        H=self.input_value[1].get()
        Rp=self.input_value[2].get()
        fw=self.input_value[3].get()
        ql=self.input_value[4].get()
        t0=self.input_value[5].get()
        GT=self.input_value[6].get()
        angle=self.input_value[7].get()
        pwf=iterator(pt*10**6,H,int(Rp),fw,ql,t0,GT,angle)[-1]/10**6
        self.output_value[0].set(round(pwf,2))

if __name__=='__main__':
    m=mainform()