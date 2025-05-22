from copy import deepcopy
from math import exp
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal

#Расчет погрешностижж
def get_sigma(g, arr):
    M = np.size(g)
    diff = g - arr  # Разность между g и y_output1
    squared_diff = diff ** 2  # Квадрат разности
    mean_squared_diff = np.sum(squared_diff) / (M + 1)  # Среднее значение квадратов разностей
    sigma = np.sqrt(mean_squared_diff)  # Среднеквадратическое отклонение
    return sigma

# Функция для округления до N значащих цифр
def round_to_significant(x, sig = 2):
    if x == 0:
        return 0  # Если число равно нулю, возвращаем 0
    # Определяем порядок числа (степень 10)
    order = int(np.floor(np.log10(abs(x))))
    # Вычисляем множитель для округления
    factor = 10 ** (sig - 1 - order)
    factor = Decimal(factor)
    x = Decimal(x)
    # Округляем число
    return float(round(x * factor) / factor)

#ИХ инерциального RC
def IH_inerc_RC(t, tau):
    g = np.where(t >= 0, np.exp(-t/tau)*(1/tau), 0) # ИХ = 0 для t < 0
    g = g/np.max(g)
    return g

# Моделирование RC звено методом инвариантной импульсной характеристики
def RC_modeling_invariant_Impulse_response(y_input, T_d, tau):
    #RC звено методом инвариантной импульсной характеристики
    N = np.size(y_input)
    if N<1:
        pass
    y_output1 = np.zeros(N)
    sumator = np.zeros(N)
    b1 = round_to_significant(exp(-T_d/tau), 2)
    a0 = round_to_significant(T_d/tau, 2)
    for i in range(1, N, 1):
        y_output1[i] = y_input[i]*a0 + y_output1[i-1]*b1
    y_output1 = y_output1/np.max(y_output1)
    return (y_output1)


# Моделирование RC звено методом инвариантной переходной характеристики
def RC_modeling_invariant_Transitional_response(y_input, T_d, tau):
    N = np.size(y_input)
    y_output2 = np.zeros(N)
    sumator1 = np.zeros(N)
    a1 = round_to_significant(1-exp(-T_d/tau))
    b1 = round_to_significant(exp(-T_d/tau))
    # 
    for i in range(1, N, 1):
        y_output2[i] = y_input[i-1]*a1 + y_output2[i-1]*b1
    y_output2 = y_output2/np.max(y_output2)
    return y_output2

# Моделирование RC звена методом билинейного z-преобразования
def RC_modeling_bilinear_z_transformation(y_input, T_d, tau):
    N = np.size(y_input)
    y_output3 = np.zeros(N)
    sumator1 = np.zeros(N)
    a0 = round_to_significant(T_d/(T_d+2*tau))
    a1 = round_to_significant(T_d/(T_d+2*tau))
    b1 = round_to_significant((2*tau - T_d)/(2*tau+T_d))
    for i in range(1, N, 1):
        y_output3[i] = y_input[i]*a0 + y_input[i-1]*a1 + y_output3[i-1]*b1
    y_output3 = y_output3/np.max(y_output3)
    return y_output3

def graf_sigma_from_T_d(T_max, T_min, tau):
    #Построение графиков
    T_d_min = tau/1000
    T_d_max = tau/10
    
    T_d_delta = tau/10000
    
    T_d_values = [ 0.01, 0.009, 0.008, 0.0075, 0.005, 0.003, 0.001, 0.0009, 0.00075]

    sigma1 = np.zeros(np.size(T_d_values))
    sigma2 = np.zeros(np.size(T_d_values))    
    sigma3 = np.zeros(np.size(T_d_values))
    
    for j in range(np.size(T_d_values)):
        
        t = np.arange(T_min, T_max + T_d_values[j], T_d_values[j])  # +T_d чтобы включить T_max
        N = np.size(t)
        # Входной сигнал - дельта функция
        y_input = np.zeros(N) 
        y_input[0] = 1

        #сдвигаем массив входненого сигнала на 1
        y_shifted = np.roll(y_input, 1)
        y_shifted[0] = 0  # Заменяем первый элемент на 0
        y_input = y_shifted

        #сдвигаем массив времени на 1
        t_shifted = np.roll(t, 1)
        t_shifted[0] = -T_d_values[j]  # Заменяем первый элемент на -T_d
        t = t_shifted

        y_output1 = RC_modeling_invariant_Impulse_response(y_input, T_d_values[j], tau)
        y_output2 = RC_modeling_invariant_Transitional_response(y_input, T_d_values[j], tau)
        y_output3 = RC_modeling_bilinear_z_transformation(y_input, T_d_values[j], tau)
        
        g = IH_inerc_RC(t, tau)
        
        sigma1[j] = get_sigma(g, y_output1)
        sigma2[j] = get_sigma(g, y_output2)
        sigma3[j] = get_sigma(g, y_output3)
    plt.plot(T_d_values, sigma1, 'o', markersize=2, label='IIH')  # 'o-' означает точки, соединенные линиями#
    plt.plot(T_d_values, sigma3, 'o', markersize=2, label='BZ')  # 'o-' означает точки, соединенные линиями
    plt.xlabel('T')  # подпись оси x
    plt.ylabel('')  # подпись оси y
    plt.title('')  # заголовок графика
    plt.legend()  # добавление легенды
    plt.grid(True)  
    plt.show()  # отображение графика

# Параметры моделирования
R = 1e5
C = 1e-6
sig = 2 #округление до 2 значащих чисел
T_max = 0.1 # Конец времени моделирования в 'с' 0.1


T_min = 0 # Начало времени моделирования в 'с'
tau = R*C
tau = 0.1
T_d = 0.005  # Период дискретизации (шаг по времени)
# Создание массива времени
t = np.arange(T_min, T_max + 2*T_d, T_d)  # +T_d чтобы включить T_max
N = np.size(t)
# Входной сигнал - дельта функция
y_input = np.zeros(N) 
y_input[0] = 1

#сдвигаем массив входненого сигнала на 1
y_shifted = np.roll(y_input, 1)
y_shifted[0] = 0  # Заменяем первый элемент на 0
y_input = y_shifted

#сдвигаем массив времени на 1
t_shifted = np.roll(t, 1)
t_shifted[0] = -T_d  # Заменяем первый элемент на -T_d
t = t_shifted

#инерциальное RC звено методом инвариантной импульсной характеристики
y_output1 = RC_modeling_invariant_Impulse_response(y_input, T_d, tau)

#инерциальное RC звено методом инвариантной переходной характеристики    
y_output2 = RC_modeling_invariant_Transitional_response(y_input, T_d, tau)

#инерциальное RC звено методом билинейного z-преобразования
y_output3 = RC_modeling_bilinear_z_transformation(y_input, T_d, tau)

#ИХ инерциального RC
g = IH_inerc_RC(t, tau)
sigma1 = get_sigma(g, y_output1)
sigma2 = get_sigma(g, y_output2)
sigma3 = get_sigma(g, y_output3)
print(sigma1, sigma2, sigma3)

arr = [1, 2, 3, 4]
g = g[1:]
t = t[1:]
y_output1 = y_output1[1:]
y_output3 = y_output3[1:]
# построение графика
plt.plot(t, g, 'o-', markersize=4, label='g(t)', color = 'red')  # 'o-' означает точки, соединенные линиями    
plt.plot(t, y_output1, 'o-', markersize=1, label='IIH')  # 'o-' означает точки, соединенные линиями#
plt.plot(t, y_output3, 'o-', markersize=1, label='Bz')  # 'o-' означает точки, соединенные линиями
plt.xlabel('t')  # подпись оси x
plt.ylabel('')  # подпись оси y
plt.title('')  # заголовок графика
plt.legend()  # добавление легенды
plt.grid(True)  
plt.show()  # отображение графика

graf_sigma_from_T_d(T_max, T_min, tau)