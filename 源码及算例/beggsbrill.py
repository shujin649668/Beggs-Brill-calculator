import math


class pressurevalias:
    Qg = 0
    Ql = 0
    Pa = 0
    Sm = 0
    rL = 0
    rg = 0
    ug = 0
    uL = 0
    IDT = 0
    angle = 0
    FlowType = ''


class pressuremethod(pressurevalias):
    def beggsbrill(self):
        Qg = self.Qg
        Ql = self.Ql
        IDT = self.IDT
        rg = self.rg
        rL = self.rL
        angle = self.angle
        Sm = self.Sm
        Pa = self.Pa
        ug = self.ug
        uL = self.uL
        FlowType = self.FlowType

        Qm = Qg + Ql
        A = 3.1415926 * IDT ** 2 / 4
        vsg = Qg / A
        vsl = Ql / A
        vm = vsl + vsg
        Gm = Qg * rg + Ql * rL
        El = Ql / Qm
        Flag = 0
        Hl = Hl1 = Hl2 = 0
        # 将度转化为弧度
        angle = (90 - angle) * 3.1415926 / 180
        # 计算弗鲁德数NFr
        NFr = vm * vm / (9.8 * IDT)
        # 计算液相速度准数
        Nvl = vsl * (rL / (9.8 * Sm)) ** 0.25
        # 计算分区线方程
        l1 = 316 * El ** 0.302
        l2 = 0.0009252 * El ** (-2.4684)
        l3 = 0.1 * El ** (-1.4516)
        l4 = 0.5 * El ** (-6.738)

        # 判断当前流型
        if (El < 0.01 and NFr < l1) or (El >= 0.01 and NFr < l2):
            Flag = 1  # 分离流动
            FlowType = "Bubble"
            if (El >= 0.01 and NFr > l2) and NFr <= l3:
                Flag = 2  # 过渡流动
                FlowType = "Slug"

        if (((0.01 <= El < 0.4) and NFr > l3) and NFr < l1) or ((El >= 0.4 and NFr > l3) and NFr <= l4):
            Flag = 3  # 间歇流动
            FlowType = "Churn"

        if (El < 0.4 and NFr >= l1) or (El >= 0.4 and NFr > l4):
            Flag = 4  # 分散流动
            FlowType = "Annulus"

        # 计算真实持液率
        if Flag == 1 or Flag == 2:  # 计算上坡流动的系数 （流体从低处向高处流动）
            C0 = (1 - El) * math.log(0.011 * math.pow(El, -3.768) * math.pow(Nvl, 3.539) * math.pow(NFr, -1.614))
            Hl = (1 + C0 * (math.sin(1.8 * angle) - 0.3333 * math.pow(math.sin(1.8 * angle), 3))) * \
                 9.8 * math.pow(El, 0.4846) / math.pow(NFr, 0.0868)
            if Flag == 2:
                Hl1 = Hl

        if Flag == 3 or Flag == 2:
            C0 = (1 - El) * math.log(2.96 * math.pow(El, 0.305) * math.pow(Nvl, -0.4473) * math.pow(NFr, 0.0978))
            # 计算上坡流动的系数 （流体从低处向高处流动）
            Hl = (1 + C0 * (math.sin(1.8 * angle) - 0.3333 * math.pow(math.sin(1.8 * angle), 3))) * 9.8 * math.pow(El,
                                                                                                                   0.5351) \
                 / math.pow(NFr, 0.0173)
            if Flag == 2:
                Hl2 = Hl

        if Flag == 4:
            Hl = 1.065 * math.pow(El, 0.5929) / math.pow(NFr, 0.0609)  # 计算分散流水平持液率

        if Flag == 2:  # 过渡流持液率计算
            Hl = (l3 - NFr) / (l3 - l2) * Hl1 + (1 - (l3 - NFr) / (l3 - l2)) * Hl2

        if Hl > 1 or Hl < El:
            if Hl < El:
                Hl = El
            else:
                Hl = 1
        # 计算比持液率
        y = El / (Hl * Hl)
        # 计算系数
        if 1 < y < 1.2:
            S = math.log(2.2 * y - 1.2)
        else:
            a = math.log(y)
            S = a / (-0.0523 + 3.182 * a - 0.8725 * a * a + 0.01853 * math.pow(a, 4))

        # 计算两相流动的雷诺数
        rm = rg * (1 - El) + rL * El  # 混合流体密度
        Nre = rm * vm * IDT / (ug * (1 - El) + uL * El)  # 两相流雷诺数
        # 计算无滑脱气液两相流阻力系数
        lamb = (0.0056 + 0.5 / math.pow(Nre, 0.32)) * math.exp(S)
        dPh = rm * 9.8 * math.sin(angle)
        dPf = lamb * vm * Gm / (2 * IDT * A)
        dP = (dPh + dPf) / (1 - rm * vm * vsg / Pa)
        dPa = dP * (rm * vm * vsg) / Pa
        vl = vsl / Hl
        vg = vsg / (1 - Hl)
        return FlowType, dP, dPa, vl, vg
