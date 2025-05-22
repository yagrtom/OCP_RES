# -*- coding: cp1251 -*-
# передаточная функция K*p*tau1/((1+p*tau1)(p*tau2+1))
from copy import deepcopy
from math import exp
from re import A
import numpy as np
import math
# Моделирование методом инвариантной импульсной характеристики
def modeling_invariant_Impulse_response(y_input, t, T_d, tau1, tau2, K):
    Am1 = K/(tau2-tau1)
    Am2 = K*tau1/(tau1*tau2-tau2*tau2)
    # print (Am1, Am2)
    a0 = T_d*(Am1 + Am2)
    a1 = T_d*(-Am1*exp(-T_d/tau2)-Am2*exp(-T_d/tau1))
    b1 = exp(-T_d/tau1)+exp(-T_d/tau2)
    b2 = -exp(-T_d/tau1)*exp(-T_d/tau2)
    print("_____________________________________________________")
    print("Метод инвариантной импульсной характеристики:")
    print('a0 = ', a0)
    print('a1 = ', a1)
    print('b1 = ', b1)
    print('b2 = ', b2)
    print("_____________________________________________________")
    y_output1 = np.zeros_like(t)
    for i in range(2, np.size(y_input), 1):
        y_output1[i] = y_input[i]*a0 + y_input[i-1]*a1 + y_output1[i-1]*b1 + y_output1[i-2]*b2
    return (y_output1)

# Моделирование методом инвариантной переходной характеристики
def modeling_invariant_transient_response(y_input, t, T_d, tau1, tau2, K):
    Am1 = K/(tau2-tau1)
    Am2 = K*tau1/(tau1*tau2-tau2*tau2)
    a0 = -Am1*tau1-Am2*tau2
    a1 = Am1*tau1+Am2*tau2+Am1*tau1*exp(-T_d/tau2)+Am2*tau2*exp(-T_d/tau1)
    a2 = -Am1*tau1*exp(-T_d/tau2) - Am2*tau2*exp(-T_d/tau1)
    b1 = exp(-T_d/tau1)+exp(-T_d/tau2)
    b2 = -exp(-T_d/tau1)*exp(-T_d/tau2)
    print("_____________________________________________________")
    print("Метод инвариантной переходной характеристики:")
    print('a0 = ', a0)
    print('a1 = ', a1)
    print('a2 = ', a2)
    print('b1 = ', b1)
    print('b2 = ', b2)
    print("_____________________________________________________")
    y_output1 = np.zeros_like(t)
    for i in range(2, np.size(y_input), 1):
        y_output1[i] = y_input[i]*a0 + y_input[i-1]*a1 + y_input[i-2]*a2 + y_output1[i-1]*b1 + y_output1[i-2]*b2
    return (y_output1)

# Моделирование методом согласованного z
def modeling_consistent_z_transformation(y_input, t, T_d, tau1, tau2, K):
    a0 = 1*K*T_d/tau2
    a1 = -1*K*T_d/tau2
    b1 = exp(-T_d/tau1)+exp(-T_d/tau2)
    b2 = -exp(-T_d/tau1)*exp(-T_d/tau2)
    print("_____________________________________________________")
    print("Метод согласованное z-преобразование:")
    print('a0 = ', a0)
    print('a1 = ', a1)
    print('b1 = ', b1)
    print('b2 = ', b2)
    print("_____________________________________________________")
    y_output1 = np.zeros_like(t)
    for i in range(2, np.size(y_input), 1):
        y_output1[i] = y_input[i]*a0 + y_input[i-1]*a1 + y_output1[i-1]*b1 + y_output1[i-2]*b2
    return (y_output1)

# Моделирование методом ,билинейного z
def modeling_bilinear_z_transformation(y_input, t, T_d, tau1, tau2, K):
    a0 = K*tau1*2*T_d/((T_d+2*tau1)*(T_d+2*tau2))
    a1 = 0
    a2 = -a0
    b1 = (-2*T_d*T_d+8*tau1*tau2)/((T_d+2*tau1)*(T_d+2*tau2))
    b2 = (-T_d+2*tau2)*(T_d-2*tau1)/((T_d+2*tau1)*(T_d+2*tau2))
    print("_____________________________________________________")
    print("Метод билинейного z-преобразования:")
    print('a0 = ', a0)
    print('a1 = ', a1)
    print('a2 = ', a2)
    print('b1 = ', b1)
    print('b2 = ', b2)
    print("_____________________________________________________")
    y_output1 = np.zeros_like(t)
    for i in range(2, np.size(y_input), 1):
        y_output1[i] = y_input[i]*a0 + y_input[i-1]*a1 + y_input[i-2]*a2 + y_output1[i-1]*b1 + y_output1[i-2]*b2
    return (y_output1)
