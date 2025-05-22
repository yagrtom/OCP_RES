# -*- coding: cp1251 -*-

from copy import deepcopy
from math import exp
from re import A
import matplotlib.pyplot as plt
import numpy as np
import math
import kptau1
import kptau2
def All_modeling_CR_RC(R1, R2, C1, C2, K):
    tau1 = R1*C1
    tau2 = R2*C2
    # tau1 = round(R1*C1, 4)
    # tau2 = round(R2*C2, 4)
    print("\u03C41 = ",tau1, "\u03C41 = ", tau2)
    if (tau1 > tau2):
        T_max = 5*tau1
        T_d = tau2/5
    else:
        T_max = 5*tau2
        T_d = tau1/5
        
    t = np.arange(-2*T_d, T_max, T_d)

    delta_func = np.where(t == 0, 1, 0)  
    unit_step = np.where(t >= 0, 1, 0)

    PH = kptau1.modeling_invariant_Impulse_response(unit_step, t, T_d, tau1, tau2, K)
    PH2 = kptau1.modeling_invariant_transient_response(unit_step, t, T_d, tau1, tau2, K)
    PH3 = kptau1.modeling_consistent_z_transformation(unit_step, t, T_d, tau1, tau2, K)
    PH4 = kptau1.modeling_bilinear_z_transformation(unit_step, t, T_d, tau1, tau2, K)

    plt.plot(t, PH, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH2, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH3, markersize=4, label='����Z')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH4, markersize=4, label='�����Z')  # 'o-' �������� �����, ����������� �������  
    plt.xlabel('t')  # ������� ��� x
    plt.ylabel('')  # ������� ��� y
    plt.title('')  # ��������� �������
    plt.legend()  # ���������� �������
    plt.grid(True)  
    plt.show()  # ����������� �������

def All_modeling_RL_LR(R1, R2, L1, L2, K):
    tau1 = L1/R1
    tau2 = L2/R2
    # tau1 = round(R1*L1, 4)
    # tau2 = round(R2*L2, 4)
    print("\u03C41 = ",tau1, "\u03C42 = ", tau2)
    if (tau1 > tau2):
        T_max = 5*tau1
        T_d = tau2/5
    else:
        T_max = 5*tau2
        T_d = tau1/5
    t = np.arange(-2*T_d, T_max, T_d)

    delta_func = np.where(t == 0, 1, 0)  
    unit_step = np.where(t >= 0, 1, 0)

    PH = kptau1.modeling_invariant_Impulse_response(unit_step, t, T_d, tau1, tau2, K)
    PH2 = kptau1.modeling_invariant_transient_response(unit_step, t, T_d, tau1, tau2, K)
    PH3 = kptau1.modeling_consistent_z_transformation(unit_step, t, T_d, tau1, tau2, K)
    PH4 = kptau1.modeling_bilinear_z_transformation(unit_step, t, T_d, tau1, tau2, K)

    plt.plot(t, PH, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH2, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH3, markersize=4, label='����Z')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH4, markersize=4, label='�����Z')  # 'o-' �������� �����, ����������� �������  
    plt.xlabel('t')  # ������� ��� x
    plt.ylabel('')  # ������� ��� y
    plt.title('')  # ��������� �������
    plt.legend()  # ���������� �������
    plt.grid(True)  
    plt.show()  # ����������� �������
    
def All_modeling_RC_CR(R1, R2, C1, C2, K):
    tau1 = R1*C1
    tau2 = R2*C2
    # tau1 = round(R1*C1, 4)
    # tau2 = round(R2*C2, 4)
    print("\u03C41 = ",tau1, "\u03C41 = ", tau2)
    if (tau1 > tau2):
        T_max = 5*tau1
        T_d = tau2/5
    else:
        T_max = 5*tau2
        T_d = tau1/5
        
    t = np.arange(-2*T_d, T_max, T_d)

    delta_func = np.where(t == 0, 1, 0)  
    unit_step = np.where(t >= 0, 1, 0)

    PH = kptau2.modeling_invariant_Impulse_response(unit_step, t, T_d, tau1, tau2, K)
    PH2 = kptau2.modeling_invariant_transient_response(unit_step, t, T_d, tau1, tau2, K)
    PH3 = kptau2.modeling_consistent_z_transformation(unit_step, t, T_d, tau1, tau2, K)
    PH4 = kptau2.modeling_bilinear_z_transformation(unit_step, t, T_d, tau1, tau2, K)

    plt.plot(t, PH, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH2, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH3, markersize=4, label='����Z')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH4, markersize=4, label='�����Z')  # 'o-' �������� �����, ����������� �������  
    plt.xlabel('t')  # ������� ��� x
    plt.ylabel('')  # ������� ��� y
    plt.title('')  # ��������� �������
    plt.legend()  # ���������� �������
    plt.grid(True)  
    plt.show()  # ����������� �������
    
def All_modeling_LR_RL(R1, R2, L1, L2, K):
    tau1 = L1/R1
    tau2 = L2/R2
    # tau1 = round(R1*L1, 4)
    # tau2 = round(R2*L2, 4)
    print("\u03C41 = ",tau1, "\u03C42 = ", tau2)
    if (tau1 > tau2):
        T_max = 5*tau1
        T_d = tau2/5
    else:
        T_max = 5*tau2
        T_d = tau1/5
    t = np.arange(-2*T_d, T_max, T_d)

    delta_func = np.where(t == 0, 1, 0)  
    unit_step = np.where(t >= 0, 1, 0)

    PH = kptau2.modeling_invariant_Impulse_response(unit_step, t, T_d, tau1, tau2, K)
    PH2 = kptau2.modeling_invariant_transient_response(unit_step, t, T_d, tau1, tau2, K)
    PH3 = kptau2.modeling_consistent_z_transformation(unit_step, t, T_d, tau1, tau2, K)
    PH4 = kptau2.modeling_bilinear_z_transformation(unit_step, t, T_d, tau1, tau2, K)

    plt.plot(t, PH, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH2, markersize=4, label='���')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH3, markersize=4, label='����Z')  # 'o-' �������� �����, ����������� �������  
    plt.plot(t, PH4, markersize=4, label='�����Z')  # 'o-' �������� �����, ����������� �������  
    plt.xlabel('t')  # ������� ��� x
    plt.ylabel('')  # ������� ��� y
    plt.title('')  # ��������� �������
    plt.legend()  # ���������� �������
    plt.grid(True)  
    plt.show()  # ����������� �������
    
All_modeling_CR_RC(1e3, 2e3, 1e-6, 1e-7, 100) #��� �������
# All_modeling_RL_LR(1e3, 2e3, 1e-3, 0.5e-3, 10) # ������� ����
# All_modeling_LR_RL(1e3, 1e3, 1e-3, 0.5e-3, 1) # ������� ������
# All_modeling_CR_RC(2e3, 1e3, 1e-6, 0.5e-6, 10) # ������ ����
# All_modeling_RL_LR(1e3, 2e3, 1e-3, 0.5e-3, 10) # ��������� ������� �� ����