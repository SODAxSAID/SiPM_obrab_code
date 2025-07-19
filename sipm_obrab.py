import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import scipy.signal
from scipy.optimize import curve_fit
import math
import os


def find_m_and_s_all_file(x, y):
    m = 0
    num_of_sred = 0
    maximum = max(y)
    for i in range(len(x)):
        m += (x[i]*y[i]/sum(y))
    for i in range(len(x)):
        if x[i] < m:
            num_of_sred = i
    for i in range(0, num_of_sred):
        if y[i] > maximum/4:
            num_start = i
            break
    #print(num_start)
    return num_of_sred, num_start


def start_of_graph(x, y, data_file):
    voltage = data_file.split("\\")[-1].split(".")[0].split("_")[-1]
    graph_info = {'x': [], 'y': []}
    start_num = 0
    for i in range(len(x)):
        if y[i] < y[i+1] < y[i+2] < y[i+3] < y[i+4]:
            start_num = i
            break
    num_seredina, num_start = find_m_and_s_all_file(x, y)
    stop_num = num_seredina
    for i in range(7, -30, -1):
        #print(x[num_seredina+i], x[num_seredina+i-1], x[num_seredina])
        if y[num_seredina+i] < y[num_seredina+i-1] < y[num_seredina+i-2] < y[num_seredina+i-3]:
            stop_num = num_seredina+i
            #print(stop_num)
            break
    for j in range(num_start, stop_num):
        graph_info['x'].append(x[j])
        graph_info['y'].append(y[j])
    return graph_info, num_start


def smooth_plot(x, y):
    y_new = []
    x_new = []
    for i in range(2, len(x)-2, 1):
        a = int((y[i-2]+2*y[i-1]+2*y[i]+2*y[i+1]+y[i-2])/8)
        y_new.append(a)
        x_new.append(x[i])
    smooth_grath = {'x': x_new, 'y': y_new}
    return smooth_grath


def find_m_s(peaks, graph):
    # grath is smooth_grath
    m_1 = []
    sigma = []
    range_for_x = []
    range_for_y = []
    right_difs = []
    left_difs = [int(round((peaks[1]-peaks[0])/3, 0))]
    for i in range(len(peaks)-1):
        right_difs.append(int(round((peaks[i+1]-peaks[i])/1.7, 0)))
    right_difs.append(right_difs[len(right_difs)-1])
    for i in range(1, len(peaks)):
        left_difs.append(int(round((peaks[i]-peaks[i-1])/3, 0)))
    for i in range(len(peaks)):
        y_ranges = []
        x_ranges = []
        summ_for_m_sigma = 0
        summ_for_m = 0
        left = min(left_difs[i], 4)
        right = min(right_difs[i], 3)
        if peaks[i] + right > len(graph['x']):
            right = len(graph['x']) - peaks[i]
        if peaks[i] - left < 0:
            left = len(graph['x']) + peaks[i]
        for k in range(peaks[i] - left, peaks[i] + right):
            summ_for_m += graph['x'][k]*graph['y'][k]
            summ_for_m_sigma += graph['x'][k]*graph['x'][k]*graph['y'][k]
            y_ranges.append(graph['y'][k])
            x_ranges.append(graph['x'][k])
        m_1.append(round(summ_for_m/sum(y_ranges), 3))
        sigma.append(round(math.sqrt((summ_for_m_sigma/sum(y_ranges)) -
                               (summ_for_m/sum(y_ranges))*(summ_for_m/sum(y_ranges))), 2))
        #print(peaks[i], round(summ_for_m/sum(y_ranges), 3))
        range_for_x.append(x_ranges)
        range_for_y.append(y_ranges)
    return m_1, sigma, range_for_x, range_for_y


def multi_fit(x, a, m, s):
    func = a*np.exp(-(x - m)**2/(2*s**2))
    return func


def codes_to_Q(x_range, data_file):
    a = data_file.split("\\")[-1].split(".")[0].split("_")[-3]
    if a == 'sipm':
        chan_info = data_file.split("\\")[-1].split(".")[0].split("_")[-4]
    else:
        chan_info = data_file.split("\\")[-1].split(".")[0].split("_")[-3]
    x_range_new = []
    for i in x_range:
        x_range_new.append(channels[chan_info][1]*i + channels[chan_info][0])
    #print(x_range_new)
    return x_range_new


def analyse_spectrum(data_file, sum_of_sigmas):
    x = []
    y = []
    print(data_file)
    fig, ax = plt.subplots(figsize=[12, 9])
    sipm_num = data_file.split("\\")[-1].split(".")[0].split("_")[-2]
    voltage = data_file.split("\\")[-1].split(".")[0].split("_")[-1]
    a = data_file.split("\\")[-1].split(".")[0].split("_")[-3]
    if a == 'sipm':
        chan_info = data_file.split("\\")[-1].split(".")[0].split("_")[-4]
    else:
        chan_info = data_file.split("\\")[-1].split(".")[0].split("_")[-3]
    with open(data_file) as f:
        for lines in f:
            buff = lines.split('\t')
            buff = [int(i) for i in buff]
            x.append(buff[0])
            y.append(buff[1])
    smooth_graph = smooth_plot(x, y)
    smooth_graph2, num_start = start_of_graph(smooth_graph['x'], smooth_graph['y'], data_file)
    plt.plot(codes_to_Q(x, data_file), y, '.-', c='black', lw=2.12, label='Spectrum of SiPM')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.plot(codes_to_Q(smooth_graph2['x'], data_file), smooth_graph2['y'], '.-', c='green')
    #plt.plot(codes_to_Q(smooth_graph2['x'], data_file), smooth_graph2['y'], '.-', c='red')
    #plt.show()
    peaks, _ = scipy.signal.find_peaks(smooth_graph2['y'], prominence=prominences[voltage], distance=distances[voltage])
    print(peaks)
    #print(num_start)
    k = 0
    peaks = np.delete(peaks, 0)
    while (smooth_graph2['x'][peaks[k]] < num_start):
        peaks = np.delete(peaks, k)
    if (smooth_graph2['y'][peaks[0]] < smooth_graph2['y'][peaks[0]-15]):
        peaks = np.delete(peaks, 0)
    x2 = [smooth_graph2['x'][i] for i in peaks]
    y2 = [smooth_graph2['y'][i] for i in peaks]
    plt.plot(codes_to_Q(x2, data_file), y2, '.', color='red', markersize=15)
    #plt.show()
    found_means[voltage] = []
    found_sigmas[voltage] = []
    means, sigmas, x_need, y_need = find_m_s(peaks, smooth_graph2)
    schet = 0
    for i in range(len(means)):
        x_range = np.arange(min(smooth_graph2['x']), max(smooth_graph2['y']), 0.01)
        params = [y2[i], means[i], sigmas[i]]
        popt, pcov = curve_fit(multi_fit, x_need[i], y_need[i], p0=params)
        print(popt)
        print(channels[chan_info][1]*popt[1] + channels[chan_info][0])
        if popt[2] > 27:
            schet -= 1
            continue
        sum_of_sigmas += (popt[2]/popt[1])**2
        y_fit = multi_fit(x_range, *popt)
        found_means[voltage].append(round((channels[chan_info][1]*popt[1] + channels[chan_info][0]), 3))
        found_sigmas[voltage].append(round((channels[chan_info][1] * popt[0] + channels[chan_info][0]), 3))
        plt.plot(codes_to_Q(x_range, data_file), y_fit, lw=2.2, label=f'GaussFit {str(i+schet+1)} peak')
        plt.plot(channels[chan_info][1]*popt[1] + channels[chan_info][0], popt[0], '.', c='red', markersize=15)
    plt.grid(True)
    plt.xlim(min(codes_to_Q(x, data_file)), 3)
    plt.xlabel('Заряд, пКл', fontsize=25)
    plt.ylabel('Количество событий', fontsize=25)
    plt.title(f'Спектр SiPM номер '+sipm_num+' при напряжении питания ' + str(voltages_for_graph[voltage])+' В',
              fontsize=30)
    plt.legend(fontsize=20)
    plt.savefig(fname=f'..sipm{sipm_num}_{voltage}.png', dpi=250) # Put here your path, where you want to see resulting pictures
    fig.clear()
    plt.close()
    return sum_of_sigmas, chan_info


def find_M():
    mean_of_differences = {}
    for i in found_means:
        dif_of_means = []
        for j in range(1, len(found_means[i])):
            dif_of_means.append(found_means[i][j] - found_means[i][j-1])
        mean_of_differences[i] = round((sum(dif_of_means)/(len(dif_of_means)*0.16)), 3)
    return mean_of_differences


if __name__ == "__main__":
    global found_means
    global channels
    global distances
    global found_sigmas
    global prominences
    global voltages_for_graph
    voltages_for_graph = {'16V': 25.8, '32V': 26.0, '48V': 26.3, '64V': 26.6, '80V': 26.8, '96V': 27.2, '112V': 27.3}     # code works only with this voltages, for other voltages u need to rewrite code (112V - 112 - codes ADC, not Volts)
    file_with_M = f'' # Put here your path, where you want to see resulting file with Gain
    V_for_file = 'N_sipm\t' 
    V = [26.3, 26.8, 27.3]
    count = 1
    for volts in V:
        if count != len(V):
            V_for_file += str(volts) + '\t'
            count += 1
        else:
            V_for_file += str(volts) + '\n'
    prominences = {'16V': 4,  '32V': 4, '48V': 1, '64V': 1, '80V': 3, '96V': 4, '112V': 2}
    distances = {'16V': 1, '32V': 3, '48V': 8, '64V': 10, '80V': 11, '96V': 11, '112V': 14}
    channels = {'0': [0.0215, 0.00113], '1': [0.32, 0.01], '2': [0.07173, 0.011875], '3': [0.6166, 0.00831],
                '4': [0.0477, 0.0109], '5': [0.5596, 0.0092], '6': [0.1, 0.1], '7': [0.5898, 0.009112],
                '8': [0.1, 0.1], '9': [0.1742, 0.01025], '10': [0.1811, 0.009893], '11': [0.24222, 0.00951],
                '12': [0.1433, 0.01011], '13': [0.2119, 0.00950], '14': [0.1815, 0.00977], '15': [0.2236, 0.00981],
                '16': [0.1507, 0.01034], '17': [0.1479, 0.01062], '18': [0.1819, 0.00950], '19': [0.1628, 0.00999],
                '20': [0.5481, 0.00952], '21': [0.1603, 0.01011], '22': [0.5737, 0.00915], '23': [0.1697, 0.01004],
                '24': [-0.1378, 0.01272], '25': [0.1556, 0.01031], '26': [0.1383, 0.01437], '27': [0.1450, 0.01035],
                '28': [0.5033, 0.00971], '29': [0.0886, 0.01079], '30': [0.5547, 0.00907], '31': [0.0948, 0.01116]}


    found_means = {}
    found_sigmas = {}
    sipm_path = ""    # Put here your path to sipm files 
    sipm_list = os.listdir(sipm_path)
    data_path = []
    Y_all = []
    names_all = []
    for sipm in sipm_list:
        data_path.append(sipm_path+sipm+'\\')
    for j in data_path:
        files = [j + i for i in os.listdir(j)]
        deltas = {'32V': 0, '48V': 0, '64V': 0, '80V': 0, '96V': 0, '112V': 0}
        for file in files:
            voltage = file.split("\\")[-1].split(".")[0].split("_")[-1]    # change this if your file has another name
            sum_of_sigmas = 0
            sum_of_sigmas, chan_info = analyse_spectrum(file, sum_of_sigmas)
            #print(sum_of_sigmas)
            sum_of_sigmas += ((0.1/voltages_for_graph[voltage])**2 + (10**(-4)/channels[chan_info][0])**2 +
                              (10**(-5)/channels[chan_info][1])**2)
            deltas[voltage] = round(math.sqrt(sum_of_sigmas), 6)
        #print(deltas)
        Aga = find_M()
        #print(Aga)
        sipm_num = file.split("\\")[-1].split(".")[0].split("_")[-2]    # change this if your file has another name
        names_all.append('sipm '+sipm_num)
        V_for_file += sipm_num+'\t'
        Coefs = []
        M_delta = []
        for i in ['48V', '80V', '112V']:
            Coefs.append(Aga[i])
            M_delta.append(Aga[i]*deltas[i])
        print('ошибка', M_delta)
        Y_all.append(Coefs)
        count = 1
        for i in range(len(Coefs)):
            if count != len(Coefs):
                V_for_file += str(Coefs[i]) + '\t'
                count += 1
            else:
                V_for_file += str(Coefs[i]) + '\n'
        #print(V_for_file)
        fig, ax = plt.subplots(figsize=[12, 9])
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.plot(V, Coefs, '.-', markersize=1, color='black')
        plt.errorbar(V, Coefs, yerr=M_delta, fmt='o', markersize=10, lw=2.5, color='red', ecolor='red', capsize=8)
        plt.grid(True)
        plt.ylim(0.45, 1.71)
        plt.xlabel('Напряжение питания, В', fontsize=24)
        plt.ylabel('Коэффициент усиления, 10\u2076', fontsize=24)
        plt.title(f'Зависимость коэффициента усиления от напряжения для SiPM №' + sipm_num, fontsize=22)
        plt.savefig(fname=f'M{sipm_num}.png', dpi=400)
        print(Coefs)
        plt.close()
    fig, ax = plt.subplots(figsize=[12, 9])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for coef in Y_all:
        plt.plot(V, coef, '.-', markersize=10)
    plt.grid(True)
    plt.legend(names_all, fontsize=10, ncol=6)
    plt.ylim(0.45, 1.71)
    plt.xlabel('Напряжение питания, В', fontsize=24)
    plt.ylabel('Коэффициент усиления, 10\u2076', fontsize=24)
    plt.title(f'Зависимость коэффициента усиления от напряжения для всех SiPM', fontsize=22)
    plt.savefig(fname=f'..\\results\\M_all.png', dpi=400) # Put here your path, where you want to see resulting pictures
    plt.close()

    with open(file_with_M, 'w') as f:
        f.write(V_for_file)
