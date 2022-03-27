from copy import copy
from operator import le
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import math


#функция которая читает файл itp и берёт оттуда массу каждого атома и записывает их в отдельный массив
def file_itp_read(name_file):
    data_massa = []
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
            data_massa.append(float(line.split()[-1]))

    file_itp.close()
    return data_massa

#функция которая читает файл itp и берёт оттуда номер каждого выбранного нами сорта атома и записывает их в отдельный массив (так же добавляет массу)
def file_itp_read_atom(name_file, atom_name):
    data_massa_id = [[], []]
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

#Запись масс и номеров выбранного атома
        if Flag and atom_name in line:
            data_massa_id[0].append(float(line.split()[-1])) # Добавляем массы 
            data_massa_id[1].append(int(line.split()[0])) # Добавляем номера атомов
    file_itp.close()
    return data_massa_id

#функция которая читает файл itp и берёт оттуда номер каждого выбранного нами сорта атома и записывает их в отдельный массив (так же добавляет массу)
def file_itp_read_si(name_file, data_atom, G, Nсore, Npoints):

#   Создаём массив в которой хранится число атомов кремния в отдельном генрационном слое и заполняем его
#   Для удобства сразу в него кладём ядро (т.е. добавляем единичку)
    number_atom_si_in_layers = [1]
    for i in range(G):
        number_atom_si_in_layers.append(Nсore)
#       Числа атомов в генерационном слое = число атомов слое до * функциональность точек ветвления - 1
#       Тут так же учтено что ядро может иметь другую функциональность
        Nсore *= Npoints - 1

    file_itp = open(name_file, 'r')

    number_layers = 0 # Номер слоя точек ветления. Отсчёт начниаем от ядра (ядро - нулевой слой)
    number_atom = 1 #   Колличество атомов в слое
    Flag = False
    for line in file_itp: 

# Начинаем записывать массы как только натыкаемся на эту строку; после этой строки идут характеристики каждого атома 
        if "; nr type resnr resid atom cgnr q m" in line:
            Flag = True    
            continue

#       То когды мы перестаём записывать
        if Flag and line == "\n":
            break

#       Запись масс и номеров выбранного атома в слой
        if Flag and "si" in line:
            data_atom[number_layers][-2].append(float(line.split()[-1])) # Добавляем массы 
            data_atom[number_layers][-1].append(int(line.split()[0])) # Добавляем номера атомов

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

    file_itp.close()
    return data_atom

#функция которая высчитывает отдельную координату центра масс
#data_k - координата частиц по определённой оси
def coordinate_center_mass(data_k0, data_mass0):
    data_k = data_k0.copy(); data_mass = data_mass0.copy()

#считываем сумму mi * ri для каждой отдельной координаты и всю массу системы 
    mi_xi = 0
    M = 0
    for i in range(len(data_k)):
        mi_xi += data_k[i] * data_mass[i]
        M += data_mass[i]
    
    return mi_xi / M 


#функция которая находит координаты центра масс
def center_mass(data_atom0):
    data_atom = data_atom0.copy()
    x0 = coordinate_center_mass(data_atom[0], data_atom[3])
    y0 = coordinate_center_mass(data_atom[1], data_atom[3])
    z0 = coordinate_center_mass(data_atom[2], data_atom[3])
    return x0, y0, z0
    
#функция которая возвращает r - расстояние от частицы до центра масс
#data_atom - атомы одного сорта если просисходит их разделение; в таком случае data_atom_all - атомы всей системы
def distance_to_center_of_mass(data_atom_all0, data_atom_sort0 = []):
    data_atom_all = data_atom_all0.copy(); data_atom_sort = data_atom_sort0.copy()

#   Если не пустой массив где записаны атомы всей системы то работаем с ним иначе с атомами определённого сорта
    if data_atom_sort == []:
        data_atom_cm = data_atom_all
    else:
        data_atom_cm = data_atom_sort
    
    x0, y0, z0 = center_mass(data_atom_all)
#   print(x0, y0, z0)

#Массив с растояниями
    distance = []

#Тут считается расстояние от центра масс до каждой частицы
    for i in range(0, len(data_atom_cm[0])):
        ri = ((data_atom_cm[0][i] - x0) ** 2 + (data_atom_cm[1][i] - y0) ** 2 + (data_atom_cm[2][i] - z0) ** 2 ) ** 0.5
        distance.append(ri)

    return distance
 
#функция которя считает объём шарового слоя, на вход падётся два радиуса "начала" и "конца" шарового слоя
def spherical_layer_volume(r1, r2):
    return  4 * math.pi / 3 * (r2 ** 3 - r1 ** 3)

#Создаём массив массивов из нулей в который будем добавлять плотности каждого слоя в соотвествующий массив
def data_density(file_name_g96, step):
    #   Открываем файл с траекторией
    g96 = open(file_name_g96, 'r')
    a = 0 # колличество элементов в массиве плотности
    Flag = False
    for line in g96:
        if "BOX" in line: #     После этой строки на следующей будет инфа про длину сторон массива      
            Flag = True
            continue
        if Flag:
            box = [float(line.split()[i]) for i in range(0, 3)] #   находим длину сторон нашего массива
            #Тут мы находим нибольшую сторону бокса так как размер дендримера не может быть больше его больше и делим его на на шаг по радиусу чтобы получить сколько эллементов в массиве плотности должно быть
            #P.S. длина бокса измеряется в нанометрах, step в ангстремах 
            a = int(max(box) * 10 / step) #  Наибольшая сторона бокса. Делим на два так как как-будто наш дендример не выходит за пределы бокса и его плотности там не может быть (максимальная оценка). 
            break
    g96.close()
    return [[] for i in range(0, a)]


#Добавляем массив элементов с плотностью в соотвустующие ячейки массива с общей плотностью
def sum_density(density_time0, density):
    density_time = density_time0.copy()
    for i in range(0, len(density_time)):
        density[i].append(density_time[i])
    return density

#ищем первое ненулевое вхождение в массив 
def first_none_zero(data):
    min_index = 0
    for i in range(len(data)):
        if data[i] != 0:
            if i == 0:
                min_index = 0
                break
            else:
                min_index = i 
                break
    return min_index

def density_mean(density, step):
    density_final = [] # тут будет записана средняя плотность по времени каждого слоя 

    for i in range(len(density)):
        densityj = density[i].copy() # тут записаны плотности в различные моменты времени сферы номера j
        if densityj == []:
            break

        sumj = 0 # сумма плотностей в за каждый момент времени
        for element in densityj: 
            sumj += element
        
        if densityj.count(0) == len(densityj):
            density_final.append(0)
        else:
            density_final.append(sumj / (len(densityj) - densityj.count(0))) # считаем среднее

#    j_min = first_none_zero(density_final)

    spher_radius = [step * (i + 0.5) for i in range(len(density_final))]

    return spher_radius, density_final 


#функция которая строит график, на вход подаётся массив иксов и игреков
#k - сколько графиков на одном листе мы хотим 
def make_graph(datax, datay, k = 1):

    fig = plt.figure()
    gs = GridSpec(ncols=1, nrows=1, figure=fig)

    ax = fig.add_subplot(gs[0,0])
    ax.set_xlabel(r'R, $\AA$', fontsize=20)
    ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)

    if k == 1:
#        max_x = ((max(datax) // 5) + 1) * 5
#        ax.set_xlim([0, max_x])
#        ax.set_xlim([0, 30])

#        max_y = ((max(datay) // 0.5) + 1) * 0.5
#        ax.set_ylim([0, max_y])
#        ax.set_ylim([0, 2.5])
        ax.plot(datax, datay, color='k')
    else:
        for i in range(k):
            ax.scatter(datax[i], datay[i])

    ax.grid()
    plt.show()

#функция которая записывает массив файл

#Функция которая считает профиль плотности в какой-то момент времни
#data_atom - атомы одного сорта если просисходит их разделение; в таком случае data_atom_all - атомы всей системы
#normalization - делаем ли мы нормировку ? Если True, то да 
def density_profile_one_time(data_atom_all0, step, data_atom_sort0 = [], normalization = False):
    data_atom_all = data_atom_all0.copy(); data_atom_sort = data_atom_sort0.copy()

#   Тут мы запишем расстояние от каждой частицы до центра масс
    disctance_to_cm = distance_to_center_of_mass(data_atom_all, data_atom_sort)

    if data_atom_sort == []:
        data_atom_sort = data_atom_all 

#   Делим всё прострнаство на слои, и присваиваем им номер j; и определяем какой атом какому слою принадлежит (если что слои нумеруются c 0, а их колличество на 1 больше чем макс. значение j в массиве)
#  Чтобы определить какой атом какому слою принадлежит мы делим целочисленно его расстояние от центра масс на наш выбранный шаг step
    dataj = [] # тут будут храниться то какой атом какому слою принадлежит
    for i in range(len(disctance_to_cm)):
        j = disctance_to_cm[i] * 10  // step # disctance_to_cm[i] имеет размерность нанометра, а step - ангстрем
        dataj.append(j)

        if False:
            x = data_atom_sort[0][i]
            y = data_atom_sort[1][i]
            z = data_atom_sort[2][i] 
            m = data_atom_sort[3][i]       
            print(x, y, z, m, data_atom_sort[4][i], j, disctance_to_cm[i])
    
    #Исходная формула для профиля плотности выглядит как:
    #плотность(R) = (сумма по i (mi)) /  Vi(R), где mi - частица попавшая в шаровой слой Vi(R); R - расстояние от центра масс

    #Созданим массив заполненный нулями в который будет помещена плотность(R) 
    density = [0.] * (int(max(dataj)) + 1)

    # умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
    # P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3
    if normalization:
#       Заполняем соотвуствующие ячейки колличеством атомов в слое
        for i in range(0, len(data_atom_sort[0])):
            j = int(dataj[i]) # то какому слою принадлежит атом
            density[j] += 1 # 1 столбик

    else:
#       Заполняем соотвуствующие ячейки атомной массой
        for i in range(len(data_atom_sort[3])):
            j = int(dataj[i]) # то какому слою принадлежит атом
            density[j] += data_atom_sort[3][i] # добавляем массу
        

#       умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
#       P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3
        for i in range(0, len(density)):
            valume = spherical_layer_volume(step * i , step * (i + 1)) # считаем объём шарового слоя
            density[i] *= 1.67 / valume # считаем плотность в сферическом слое

    
    spher_radius = [step * (i + 0.5) for i in range(len(density))]
#    make_graph(spher_radius , density)
    return density


#Функция в которой просиходит считываение профилей плотности в разные моменты времени для одного дендримера, которые потом записываются в массив
#Для всех атомов системы
def density_profile__for_all_the_time(density, data_atom, file_name_g96, step):

#   Открываем файл
    g96 = open(file_name_g96, 'r')

    Flag = False
    for line in g96:

        if "POSITIONRED" in line: 
            Flag = True
            continue
        if "END" in line:
            continue

        if Flag and "BOX" in line:
            density_time = density_profile_one_time(data_atom, step)
            sum_density(density_time, density)            
            
            [data_atom[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми
            Flag = False
            continue

        if Flag:
            data_xyz = line.split() 
            [data_atom[i].append(float(data_xyz[i])) for i in range(3)] # добвляем координаты x, y, z атома в определённый момент времени в соотвествующие ячейки
   
    return density_mean(density, step)


#Функция в которой просиходит считываение профилей плотности в разные моменты времени для всех атомов одного сорта дендримера, которые потом записываются в массив
#Атомы которые мы выбрали
#data_atom_sort - атомы одного сорта если просисходит их разделение; в таком случае data_atom_all - атомы всей системы
#normalization - делаем или нет нормировку
def density_profile__for_all_the_time_atom(density0, data_all0, data_atom_sort0, file_name_g96, step, normalization = False):
    density = density0.copy(); data_all = data_all0.copy(); data_atom_sort = data_atom_sort0.copy()

#   Открываем файл
    g96 = open(file_name_g96, 'r')
    
    count_line = 0
    Flag = False
    for line in g96:

        if "POSITIONRED" in line: 
            Flag = True
            continue
        
        if "END" in line:
            continue

        if Flag and "BOX" in line:
            density_time = density_profile_one_time(data_all, step, data_atom_sort, normalization)        
            sum_density(density_time, density)

            [data_all[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми
            [data_atom_sort[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми
            Flag = False
            count_line = 0
            continue

#       Заполняем массив с координатами для атомов всей системы
        if Flag:
            data_xyz = line.split() 
            [data_all[i].append(float(data_xyz[i])) for i in range(3)] # добвляем координаты x, y, z атома в определённый момент времени в соотвествующие ячейки
            count_line += 1
                                      
#       Заполняем массив с координатами для выбранных нами атомов
        if Flag and count_line in data_atom_sort[-1]:
            data_xyz = line.split() 
            [data_atom_sort[i].append(float(data_xyz[i])) for i in range(3)] # добвляем координаты x, y, z атома в определённый момент времени в соотвествующие ячейки


    g96.close()
    return density_mean(density, step)

#функция которая считает радиальный профиль плотности во все моменты времени:
#берём все атомы
def density_profile(file_name_itp, file_name_g96, step):
#   Создаём массив с координатами и массой каждый частцы
    data_atom = [[], [], [], []]

#   Добавляем в массив массу каждого атома (она измеряется в атомных массах)
    data_atom[3] = file_itp_read(file_name_itp)

#   Сюда будет записан профиль плотности
    density = data_density(file_name_g96, step) 
    
#   Получаем массив с плотностями 
    datax, datay = density_profile__for_all_the_time(density ,data_atom, file_name_g96, step)

#   Вызываем функцию которая строит график
    make_graph(datax, datay)

#функция которая считает радиальный профиль плотности во все моменты времени, но только для атомов одного сорта 
def density_profile_atom(file_name_itp, file_name_g96, step, atom_name):
#   Создаём массив с координатами и массой кождого выбранного нами атома и его номер
    data_atom_sort = [[], [], [], [], []]

#   Создаём массив в котором будет масса и координаты всех атомов (это нужно чтобы посчитать координаты центра масс системы)
    data_atom = [[], [], [], []]

#   Добавляем в массив массу каждого атома (она измеряется в атомных массах)
    data_atom_sort[3:5] = file_itp_read_atom(file_name_itp, atom_name)

#   Заполняем общий массив с массами каждого атома
    data_atom[3] = file_itp_read(file_name_itp)

#   Сюда будет записан профиль плотности
    density = data_density(file_name_g96, step) 
    
#   Получаем файл с плотностями 
    #   Получаем массив с плотностями 
    datax, datay = density_profile__for_all_the_time_atom(density ,data_atom, data_atom_sort, file_name_g96, step)
    

#   Вызываем функцию которая строит график
    make_graph(datax, datay)
    


#Фунция которая делает нормировку (чтобы площадь под графиком была равно колличеству атомов в слое который мы рассматривали)   
# суть нормировки в том что мы каждый столбик в нашей гистограмме умножаем на коэффициент x 
# x = k / (суммирает длину всех столбцов * шаг нашей гистограммы)
def normalization(datay, step, k):
    sum_h = 0
    for element in datay: # ищем сумму наших столбцов
        sum_h += element
    x = k / (step * sum_h) # ищем коэффициент x
    for i in range(len(datay)): # умножаем каждый столбик на сотвествующий коэффициент
        datay[i] *= x
    return datay
        
#функция которая считает радиальный профиль плотности во все моменты времени, для атомов кремения отдельно для каждой точки ветвления
#G - номер генирации дендримера; Ncore - функциональность ядра; Npoints - функциональность точек ветвления
def si(file_name_itp, file_name_g96, step, G, Ncore, Npoints):
#   Создаём массив с координатами и массой кождого si и его номер; создаётся отдельный слот для каждого слоя точек ветвления
    data_atom_si = [ [[], [], [], [], []] for i in range(G+1)]
     
    
#   Добавляем в массив массу каждого атома (она измеряется в атомных массах)
    file_itp_read_si(file_name_itp, data_atom_si, G, Ncore, Npoints)

#   Созадём файл в котором будет храниться информация о всех атомах в системе (нужно для подсчёта центра масс системы)
    data_atom = [[], [], [], []]

#   Заполняем его массами 
    data_atom[3] = file_itp_read(file_name_itp)

#   Сюда будет записан профиль плотности
    density = data_density(file_name_g96, step) 

    density = data_density(file_name_g96, step).copy()

#    datax, datay = density_profile__for_all_the_time_atom(density, data_atom, data_atom_si[1], file_name_g96, step)
    

#   Получаем файл с плотностями
    density_final = []
    density_final_r = []
    for i in range(G+1):
        density = data_density(file_name_g96, step).copy()

        datax, datay = density_profile__for_all_the_time_atom(density, data_atom, data_atom_si[i], file_name_g96, step, True)
        normalization(datay, step, len(data_atom_si[i]))

        density_final.append(datay)
        density_final_r.append(datax)


    make_graph(density_final_r, density_final, G+1)





#density_profile("c43g6.itp", "trg6_300.g96", 0.2)
#density_profile_atom("c43g5.itp", "tr_g5_300_1500.g96", 0.2, "si")
#density_profile_si("c43g6.itp", "trg6_300.g96", 0.2, 6, 4, 3)

#data_atom_sistem = file_itp_read("c43g6.itp")
#deta_atom_sistem_si = file_itp_read_atom("c43g6.itp", "si")

file_itp_read_si("c43g6.itp", )


