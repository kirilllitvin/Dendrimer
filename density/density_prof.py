import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import math


#Читаем файл с настройками
def file_settings_file():

    file_itp_name = ""; file_trajectory = ""; step = 0.2; working_mode = 1; atom_name = "si"; system_center_core = False; G = 0; core = 0; points = 0; file_name = 'density.txt'

    file_settings = open("setting_density.txt", 'r')

    for line in file_settings:
        
        if "#" in line or '\n' == line: continue

        if "file_itp_name" in line: file_itp_name = line.split()[2]

        if "file_trajectory_name" in line: file_trajectory = line.split()[2]

        if "step" in line:  step = float(line.split()[2])

        if "working_mode" in line: working_mode = int(line.split()[2])

        if "atom_name" in line: atom_name = line.split()[2]

        if "system_center" in line and "2" == line.split()[2]: system_center_core = True

        if "G" in line: G = int(line.split()[2])

        if "core" in line: core = int(line.split()[2])

        if "points" in line:  points = int(line.split()[2])

        if "file_name" in line: file_name = line.split()[2]

    file_settings.close()
    return file_itp_name, file_trajectory, step, working_mode, atom_name, system_center_core, G, core, points, file_name

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
def atom_sort(data_atom_sistem, atom_name, system_center_core):  
#   Создаём массив в который будет записаны координаты атомов определённого сорта, их масса и id
    data_atom_sort = [[], [], [], [], []]

#   Определённый типы атомов отличаются по массе и исходя из названия мы получаем массу
#   Затем сортируем массив всех атомов по массе
    data_atom_mass = 0
    if "si" in atom_name: data_atom_mass = 28
    if "C2" in atom_name: data_atom_mass = 14
    if "C3" in atom_name: data_atom_mass = 15
    

#   Если мы мы считаем сферические слои от ядра то для построения графиков мы его не берём  
    k = 0
    if system_center_core:
        k = 1
        
    for i in range(k, len(data_atom_sistem[3])):
        if data_atom_sistem[3][i] == data_atom_mass:
            data_atom_sort[3].append(data_atom_sistem[3][i])
            data_atom_sort[4].append(i+1)

    return data_atom_sort

# функция которая вызывает генерационные слои

#Запись в файл
def writing_to_file(datax, datay, file_name, mode_writing = 'w', G = 0):
    
    file = open(file_name, mode_writing)

#   Нужно когда мы будем производить запись генерационних слоя    
    if mode_writing == 'a':
        line0 = '{0}'.format(G) + '\n'
        file.write(line0)

    for i in range(len(datay)):
        line = '{0:7.3f}  {1:7.3f}'.format(datax[i], datay[i]) + '\n'
        file.write(line)

    file.close()

#Создаём массив с плотностями
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
    return [0. for i in range(0, a)]

#функция которая считает центр масс
def coordinate_center_mass(data_atom_sistem):
#   Сюда мы запишем координаты
    data_coordinate = []
    for k in range(3):
#   считываем сумму mi * ri для каждой отдельной координаты и всю массу системы 
        mi_xi = 0
        M = 0
        for i in range(len(data_atom_sistem[3])):
            mi_xi += data_atom_sistem[k][i] * data_atom_sistem[3][i]
            M += data_atom_sistem[3][i]
        data_coordinate.append(mi_xi / M)
    return data_coordinate[0], data_coordinate[1], data_coordinate[2]

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

def spherical_layer_volume(r1, r2):
    return  4 * math.pi / 3 * (r2 ** 3 - r1 ** 3)
   

#Функция которая считает профиль плотности в какой-то момент времни
#data_atom - атомы одного сорта если просисходит их разделение; в таком случае data_atom_all - атомы всей системы
#normalization - делаем ли мы нормировку ? Если True, то да 
def density_profile_one_time(density, step, data_atom_sistem, data_atom_sort = [], data_core = [], normalization = False):

#   Нужно если мы строим профиль плотности не всей системы    
    if  data_atom_sort == []:
        data_atom_sort = data_atom_sistem

    if data_core == []:        
#   Получаем расстояение до центра масс:
        x0, y0, z0 = coordinate_center_mass(data_atom_sistem)   
    else:
        x0, y0, z0 = data_core[0], data_core[1], data_core[2] 

#   Делим всё прострнаство на слои, и присваиваем им номер j; и определяем какой атом какому слою принадлежит (если что слои нумеруются c 0, а их колличество на 1 больше чем макс. значение j в массиве)
#   Чтобы определить какой атом какому слою принадлежит мы делим целочисленно его расстояние от центра масс на наш выбранный шаг step
    dataj = [] # тут будут храниться то какой атом какому слою принадлежит
    for i in range(len(data_atom_sort[3])):
        rcm = ((data_atom_sort[0][i] - x0) ** 2 + (data_atom_sort[1][i] - y0) ** 2 + (data_atom_sort[2][i] - z0) ** 2) ** 0.5 #  Считаем расстояние до центра масс
        j = rcm * 10  // step # disctance_to_cm[i] имеет размерность нанометра, а step - ангстрем
        dataj.append(j)
    
#   Если необходимо делать нормировку
    if normalization:
#       Заполняем соотвуствующие ячейки колличеством атомов в слое
        for i in range(len(data_atom_sort[0])):
            j = int(dataj[i]) # то какому слою принадлежит атом
            density[j] += 1 # 1 столбик

#   Считаем профиль плотности
    else:
#   Исходная формула для профиля плотности выглядит как:
#   плотность(R) = (сумма по i (mi)) /  Vi(R), где mi - частица попавшая в шаровой слой Vi(R); R - расстояние от центра масс
#   умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
#   P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3

#       Заполняем соотвуствующие ячейки атомной массой
        for i in range(len(data_atom_sort[3])):
            j = int(dataj[i]) # то какому слою принадлежит ато
            density[j] += data_atom_sort[3][i] # добавляем массу
    
def density_profile__for_all_the_time_atom(file_name_g96, step, data_atom_sistem, data_atom_sort = [], system_center_core = False  ,normalization = False):

    #сюда будет записана плотность 
    density = data_density(file_name_g96, step)

    Flag_sort_atom = False
#   Смотрим работаем ли мы с определёнными атомами системы или нет
    if data_atom_sort != []:
        Flag_sort_atom = True

#   Открываем файл
    g96 = open(file_name_g96, 'r')
    data_core = []
    count_line = 0 # Счётчик строчек 
    Flag = False
    time_step = 0 # Смотрим сколько шагов по времени мы сделалаи; это нужно для усреднения

    for line in g96:

        if "POSITIONRED" in line: 
            time_step += 1
            print(time_step)
            Flag = True
            continue
        
        if "END" in line:
            continue

        if Flag and "BOX" in line:
            density_profile_one_time(density, step, data_atom_sistem, data_atom_sort, data_core, normalization)       
            [data_atom_sistem[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми

            if Flag_sort_atom:
                [data_atom_sort[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми
            
            if system_center_core:
                data_core.clear()

            Flag = False
            count_line = 0
            continue

#       Заполняем массив с координатами для атомов всей системы
        if Flag:
            data_xyz = line.split() 
            [data_atom_sistem[i].append(float(data_xyz[i])) for i in range(3)] # добвляем координаты x, y, z атома в определённый момент времени в соотвествующие ячейки
            count_line += 1

        if system_center_core and Flag and count_line == 1:
            cm_xyz = line.split()
            [data_core.append(float(cm_xyz[i])) for i in range(3)] 

#       Заполняем массив с координатами для выбранных нами атомов если они есть
        if Flag_sort_atom and Flag and count_line in data_atom_sort[-1]:
            data_xyz = line.split() 
            [data_atom_sort[i].append(float(data_xyz[i])) for i in range(3)] # добвляем координаты x, y, z атома в определённый момент времени в соотвествующие ячейки


#   Если мы не делаем нормировку, то делаем усреднение по времени 
#   Производим усреднение по времени и умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
#   P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3
    if normalization: # Если есть нормировка то просто считаем среднее
        for i in range(len(density)):
            density[i] /= time_step   
    else:
        for i in range(len(density)):
            volume = spherical_layer_volume(step * i, step * (i+1))
            density[i] *= 1.66 / (volume * time_step)
#   Ищем последний ненулевой элемент (всё дедалось с расчётом что он будет и тут не учтено что последний не может быть ненулевой)
    last_none_zero = 0
    for i in range(1, len(density) + 1):
        if density[-i] != 0:
            last_none_zero = -i + 1
            break   

    datar= [step * (i + 0.5) for i in range(len(density[:last_none_zero]))]
    return datar, density[:last_none_zero]


#тут будет происходить вызов функции 
def density_profile(): 
#   Читаем файл с настройками и заполняем соотвествующие параметры
    file_itp_name, file_trajectory_name, step, working_mode, atom_name, system_center_core, G, core, points, file_name = file_settings_file()

#   Получаем массив с массами всей системы
    data_atom_sistem = file_itp_read(file_itp_name)

    

#   Необходимо для построения графиков
    fig = plt.figure()
    gs = GridSpec(ncols=1, nrows=1, figure=fig)

    ax = fig.add_subplot(gs[0,0])
    ax.set_xlabel(r'R, $\AA$', fontsize=20)

#    if True:
    if working_mode == 1:
        datar1, density1 = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem)
        writing_to_file(datar1, density1, file_name)

        ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)
        ax.plot(datar1, density1)

    if working_mode == 2:
#       Получаем массив в котором только атомы кремния         
        data_atom_sort = atom_sort(data_atom_sistem, atom_name, system_center_core)
        datar2, density2 = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem, data_atom_sort, system_center_core, False)
        writing_to_file(datar2, density2, file_name)

        ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)
        ax.plot(datar2, density2)
        

    if working_mode == 3:
#       Получаем массив с атомами кремния
        data_atom_si = atom_si(data_atom_sistem, G, core, points)
        file = open(file_name, 'w') # очищаем файл перед записью
        file.close()
        for i in range(1, G+1):
            datar3, density3 = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem, data_atom_si[i], system_center_core, True)
            writing_to_file(datar3, density3, file_name,'a', i)
            ax.plot(datar3, density3)

            
    if working_mode == 4:
        pass

    ax.grid()
    plt.show()

density_profile()