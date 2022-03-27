import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import math


#Читаем файл с настройками
def file_settings_file():

    file_itp_name = ""; file_trajectory = ""; step = 0.2; working_mode = 1; atom_name = "si"; G = 0; core = 0; points = 0

    file_settings = open("setting_density.txt", 'r')

    for line in file_settings:
        
        if "#" in line or '\n' == line: continue

        if "file_itp_name" in line: file_itp_name = line.split()[2]

        if "file_trajectory_name" in line: file_trajectory = line.split()[2]

        if "step" in line:  step = float(line.split()[2])

        if "working_mode" in line: working_mode = int(line.split()[2])

        if "atom_name" in line: atom_name = line.split()[2]

        if "G" in line: G = int(line.split()[2])

        if "core" in line: core = int(line.split()[2])

        if "points" in line:  points = int(line.split()[2])

    file_settings.close()
    return file_itp_name, file_trajectory, step, working_mode, atom_name, G, core, points

#читаем файл itp и берём оттуда массу каждого атома
def file_itp_read(name_file):

    data_atom = [[], [], [], []]
    file_itp = open(name_file, 'r')

    Flag = False
    for line in file_itp: 

# Начинаем записывать массы как только натыкаемся на эту строку; после этой строки идут характеристики каждого атома 
        if "; nr type resnr resid atom cgnr q m" in line:
            Flag = True    
            continue

# То когды мы перестаём записывать
        if Flag and line == "\n":
            break

#Запись масс
        if Flag:
            data_atom[3].append(float(line.split()[-1]))

    file_itp.close()
    return data_atom

#Тут должна быть функция которая выделяет отдельные виды атомов
def atom_sort(data_atom_sistem, atom_name):  
#   Создаём массив в который будет записаны координаты атомов определённого сорта, их масса и id
    data_atom_sort = [[], [], [], [], []]

#   Определённый типы атомов отличаются по массе и исходя из названия мы получаем массу
#   Затем сортируем массив всех атомов по массе
    data_atom_mass = 0
    if "si" in atom_name: data_atom_mass = 28
    if "C2" in atom_name: data_atom_mass = 14
    if "C3" in atom_name: data_atom_mass = 15
 
    for i in range(len(data_atom_sistem[3])):
        if data_atom_sistem[3][i] == data_atom_mass:
            data_atom_sort[3].append(data_atom_sistem[3][i])
            data_atom_sort[4].append(i+1)

    return data_atom_sort

#выделение только атомов кремния по генерационным слоям
def atom_si(data_atom_sistem, G, core, points):
#   Создаём массив с координатами и массой кождого si и его номер; создаётся отдельный слот для каждого слоя точек ветвления  
    data_atom_si = [[[], [], [], [], []] for i in range(G+1)] 
#   Создаём массив в которой хранится число атомов кремния в отдельном генрационном слое и заполняем его
#   Для удобства сразу в него кладём ядро (т.е. добавляем единичку)
    number_atom_si_in_layers = [1]
    for i in range(G):
        number_atom_si_in_layers.append(core)
#       Числа атомов в генерационном слое = число атомов слое до * функциональность точек ветвления - 1
#       Тут так же учтено что ядро может иметь другую функциональность
        core *= points - 1

    number_layers = 0 # Номер слоя точек ветления. Отсчёт начниаем от ядра (ядро - нулевой слой)
    number_atom = 1 #   Колличество атомов в слое
    for i in range(len(data_atom_sistem[3])):
        if data_atom_sistem[3][i] == 28: 
            data_atom_si[number_layers][3].append(data_atom_sistem[3][i])
            data_atom_si[number_layers][4].append(i+1)

#           Это всё нужно чтобы соотвествующие слои точек ветления расположились в соотвествующие слои в массиве            
            if number_atom < number_atom_si_in_layers[number_layers]: #   Пока счёткик атомов в слое меньше колличества атом в слое мы делаем запись в слой
                number_atom += 1
            else: # как только счётчик принял значение атомов в слое мы перескакиваем на следующий слой, обновляя счётчик
#               Как только дошли до крайненго атома в крайнем слое оканчиваем запись в массив
                if number_atom == number_atom_si_in_layers[-1]: 
                    break
                else:
                    number_atom = 1 # Обновляем счётчик 
                    number_layers += 1 # переходим на следующий слой

    return data_atom_si





#тут будет происходить вызов функции 
def density_profile(): 
#   Читаем файл с настройками и заполняем соотвествующие параметры
    file_itp_name, file_trajectory_name, step, working_mode, atom_name, G, core, points = file_settings_file()

#   Получаем массив с массами всей системы
    data_atom_sistem = file_itp_read(file_itp_name)


    if working_mode == 1:
#       Получаем массив с массами всей системы
        pass

    if working_mode == 2:
#       Получаем массив с отсортированными атомами        
        data_atom_sort = atom_sort(data_atom_sistem, atom_name)


    if working_mode == 3:
#       Получаем массив с атомами кремния
        data_atom_si = atom_si(data_atom_sistem, G, core, points)


density_profile()