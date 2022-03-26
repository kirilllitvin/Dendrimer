from enum import Flag
from itertools import count
import matplotlib as plt
import numpy as np


#функция которая читает файл itp и выделяет нужную нам часть ветки
def file_itp_read(name_file):
#Создаём массив куда будем записывать номера нужных нам атомов
#Начинаем счёт именно с 2 ниже будет добавлен 1
    data_number = np.array([2])
#Открываем файл itp
    file_itp = open(name_file, 'r')

    Flag = False
    ai = 0; aj = 0; aj_1 = 0    #   aj_1 - это a_(j-1)
    for line in file_itp: 

#Как только натыкаемся на эту строку то начинаем анализировать файл. Это область bonds
        if "; ai aj funct c0 c1" in line:
            Flag = True    
            continue

#Как только область bonds заканчивается (попадается пустая строка) то мы перестаём читать файл 
        if Flag and line == "\n":
            break

        if Flag:
            ai = int(line.split()[0])
            aj = int(line.split()[1])

#Требется для выделения 1 дендрона и 1 "прямой линии" в нём, но тут остаётся 1 метильная группа и одна CH2 от следующего спейсера
#Допустим у нас есть точка ветвления с 3 функциональностями 
            if ai in data_number and not(aj in data_number) and aj_1 <= aj:
                data_number = np.append(data_number, aj)
                aj_1 = aj 

#Здесь мы убираем метильную группу и CH2 от спейсера. Удаляя из масивая лишение номера
            if data_number.size >= 3 and data_number[-1] > data_number[-2] + 1 and data_number[-2] > data_number[-3] + 1:
                data_number = np.delete(data_number, [-3, -2])

#Добавляем ядро  
    data_number = np.append(data_number, 1)
    data_number.sort()
           
    
    file_itp.close()
    return data_number

#Анализируем файл с траекторией дендримера. Тут мы выделяем движение нужной нам части дендрона
def file_gro_read(name_file_itp, name_file_gro, name_file_trr):

#В этот массив будут записаны номера нужных нам атомов. 
    data_number = np.array([])
#Для этого мы вызываем функцию file_itp_read
    data_number = file_itp_read(name_file_itp)

#Открываем .gro файл с иходной конформацией 
#Его открываем для того чтобы получить 2 и последнею строку 
    file_gro = open(name_file_gro, 'r')
    line_file_gro = file_gro.readlines()
#   one_line = line_file_gro[0]
    two_line = line_file_gro[1]
    last_line = line_file_gro[-1]

#2 строка содержит колличество атомов в системе, выделяя нужную нам часть дендрона изменилось колличество атомов по сравнению с целым дендроном 
#Мы тут создаём "новоую" 2 строку с нужным нам колличеством атомов 
    two_line_new = " " + str(data_number.size) + '\n'

#открываем файл .gro в котором записана траектория
    file_trr = open(name_file_trr, 'r')

#создаём новый файл в который будем записывать новую траекторию
    dendron_trr = open("dendron_trr.gro", 'w')

    Flag = False
#x0, y0, z0 - начальные координаты
    x0 = 0 ; y0 = 0; z0 = 0
#считываем файл построчно
    for line in file_trr:

#если встречаем строчку со временем то без изменений записываем её в файл
        if "t=" in line:
            dendron_trr.write(line)
            continue

#после строчки со временем нужно записть строчку с колличеством атомов в системе, которая была создана выше
#как только это строка была создана то я буду записывать траекторию
        if two_line in line:
            dendron_trr.write(two_line_new)
            Flag = True
            continue

#       как только наткунились на конец траектории в данный момент времени то записываем эту строку без изменений
        if last_line in line:
            dendron_trr.write(line)
            continue

#       Анализ 1 строки в траектории которая отвечает за кооординату ядра дендримеры 
#       Наша цель сделать ядро неподвижным и разместить его в начало координат. Точку (0, 0, 0)
        if Flag:
            line0 = line.split()

#           Запоминаем точку в которой находилось ядро
            x0 = float(line0[-3]) 
            y0 = float(line0[-2])
            z0 = float(line0[-1])

#           Создаём строку отвечаюзую за ядро и записываем её в файл
#           первые 3 элемента это параметры атома, последние 3 это его координата
            line0 = '{0:>7}{1:>8}{2:>5} {3:7.3f} {4:7.3f} {5:7.3f}'.format(line0[0], line0[1], line0[2], 0, 0, 0) + '\n'
            dendron_trr.write(line0)
            Flag = False
            continue

#       Запись траекторий остальных атомов
        if int(line.split()[2]) in data_number:
            n_line = line.split()

#       Если мы перенесли ядро в другую точку то чтобы всё было ок то мы должны перенести и остальные атомы системы
#       То есть совершить перенос системы 
#       Это достигается как x - x0; y - y0; z - z0 для каждого атома
            x = float(n_line[-3]) - x0
            y = float(n_line[-2]) - y0
            z = float(n_line[-1]) - z0

#Создаем строку отвечающую за каждый атом, анологично с тем как было сделано для ядра
            N_line = '{0:>8}{1:>7}{2:>5} {3:7.3f} {4:7.3f} {5:7.3f}'.format(n_line[0], n_line[1], n_line[2], x, y, z) + '\n'
            dendron_trr.write(N_line)

    file_trr.close()
    dendron_trr.close()
    file_gro.close()

    
#    for i in range(2, len(line_file_gro) - 2):
#        line = line_file_gro[i]
#        if int(line.split()[2]) in data_number:
#            dendron_gro.write(line)

file_gro_read("c43g6.itp", "c43g6_300.gro", "c43g6_trr.gro")


