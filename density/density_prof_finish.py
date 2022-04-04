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

#   характеристики каждого атома
#   x y z координата масса id и в какой номер массива по плотности мы его записываем
#   если последний параметр равен -1 то это значит не куда по итогу мы его не записываем
    data_atom = [[], [], [], [], []]
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
            data_atom[4].append(-1)

    file_itp.close()
    return data_atom

#Функция которая выделяет отдельные виды атомов
def atom_sort(data_atom_sistem, atom_name, system_center_core):  

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
            data_atom_sistem[4][i] = 0    
            
#выделение только атомов кремния по генерационным слоям (не берём при этом ядро)
def atom_si(data_atom_sistem, G, core, points):

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
    for i in range(len(data_atom_sistem[4])):
        if data_atom_sistem[3][i] == 28: 
            data_atom_sistem[4][i] = number_layers - 1  #   считаем слои с 0, а не с 1          

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
    

# функция которая разделяет атомы по генерационным слоям
def generation_layer(data_atom_sistem, G, core, points):

#   Создаём массив в которой хранится число атомов кремния в отдельном генрационном слое и заполняем его
#   При этом ядро не учитывается
    number_atom_si_in_layers = []
    for i in range(G):
        number_atom_si_in_layers.append(core)
#       Числа атомов в генерационном слое = число атомов слое до * функциональность точек ветвления - 1
#       Тут так же учтено что ядро может иметь другую функциональность
        core *= points - 1
    
#   Отдельно добовляем ядро в первые генрационный слой(отсчёт начинается с нуля)
    data_atom_sistem[4][0] = 0

    number_layers = 0 # Номер слоя точек ветления. Отсчёт начниаем от ядра (ядро - нулевой слой)
    number_atom = 1 #   Колличество атомов si в слое
    Flag = False # После атомов si идут сразу по номеру атомы CH3 которые нужно сразу записать в следующий ген слой
    for i in range(1, len(data_atom_sistem[4])):

#       Если попадется атом кремния то мы его  записываем в следующий генерационный слой 
        if data_atom_sistem[3][i] == 28: 
            data_atom_sistem[4][i] = number_layers + 1
            
#           Записываем CH3(которые идут сразу после кремния)
            data_atom_sistem[4][i+1] = number_layers + 1
            Flag = True
    
#           Это всё нужно чтобы соотвествующие слои точек ветления расположились в соотвествующем генерационном слое             
            if number_atom < number_atom_si_in_layers[number_layers]: #   Пока счёткик атомов в слое меньше колличества атом si в слое мы делаем запись в слой
                number_atom += 1
            else: # как только счётчик принял значение атомов в слое мы перескакиваем на следующий слой, обновляя счётчик
                    number_atom = 1 # Обновляем счётчик si
                    number_layers += 1 # переходим на следующий слой
                    print(number_layers)                  
            continue

#       Пропуск из-за записи CH3
        if Flag:
            Flag = False
            continue

#       Записывам остальные атомы
        data_atom_sistem[4][i] = number_layers
            
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
def data_density(file_name_g96, step, number_data_density):
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
    return [[0. for i in range(0, a)] for i in range(number_data_density)] 

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

def spherical_layer_volume(r1, r2):
    return  4 * math.pi / 3 * (r2 ** 3 - r1 ** 3)
   
#Функция которая считает профиль плотности в какой-то момент времни
#data_atom - атомы одного сорта если просисходит их разделение; в таком случае data_atom_all - атомы всей системы
#normalization - делаем ли мы нормировку ? Если True, то да 
def density_profile_one_time(density, step, data_atom_sistem, data_core = [], normalization = False):

    if data_core == []:        
#       Получаем расстояение до центра масс:
        x0, y0, z0 = coordinate_center_mass(data_atom_sistem)   
    else:
        x0, y0, z0 = data_core[0], data_core[1], data_core[2] 

#   Исходная формула для профиля плотности выглядит как:
#   плотность(R) = (сумма по i (mi)) /  Vi(R), где mi - частица попавшая в шаровой слой Vi(R); R - расстояние от центра масс
#   умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
#   P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3

#   Делим всё прострнаство на слои, и присваиваем им номер j; и определяем какой атом какому слою принадлежит (если что слои нумеруются c 0, а их колличество на 1 больше чем макс. значение j в массиве)
#   Чтобы определить какой атом какому слою принадлежит мы делим целочисленно его расстояние от центра масс на наш выбранный шаг step
      
#   Заполняем соотвуствующие ячейки колличеством атомов в слое
    for i in range(len(data_atom_sistem[0])):
        k = data_atom_sistem[4][i] # то какому слою принадлежит атом
        if k != -1:
            rcm = ((data_atom_sistem[0][i] - x0) ** 2 + (data_atom_sistem[1][i] - y0) ** 2 + (data_atom_sistem[2][i] - z0) ** 2) ** 0.5 #  Считаем расстояние до центра масс
            j = int(rcm * 10  // step)  # то какому шаровому слою принадлежит атом
            if normalization: #Если необходимо делать нормировку
                density[k][j] += 1 # 1 столбик
            else:
                density[k][j] += data_atom_sistem[3][i] # добавляем массу

   
def density_profile__for_all_the_time_atom(file_name_g96, step, data_atom_sistem, number_data_density, system_center_core = False, normalization = False):

    #сюда будет записана плотность 
    density = data_density(file_name_g96, step, number_data_density)
  
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
            density_profile_one_time(density, step, data_atom_sistem, data_core, normalization)       
            [data_atom_sistem[i].clear() for i in range(3)] # делаем массивы содержащие координаты атомов пустыми
                        
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


#   Если мы не делаем нормировку, то делаем усреднение по времени 
#   Производим усреднение по времени и умножаем каждую ячейку в которой хранится атомная масса на соотвествующий массовый коэффициент и делим на объём сверического слоя
#   P.S: масса 1 нуклона 1,67 * 10 ^ (-27) кг; размерность step - ангстрем; плотность (R) - имеет размерность г/см ^ 3
    data_density_finish = []
    datar_finish = []
    for j in range(len(density)):

        for i in range(len(density[j])):
            if normalization:  # Если есть нормировка то просто считаем среднее
                density[j][i] /= time_step
            else:
                volume = spherical_layer_volume(step * i, step * (i+1))
                density[j][i] *= 1.66 / (volume * time_step)          

#       Ищем последний ненулевой элемент (всё дедалось с расчётом что он будет и тут не учтено что последний не может быть ненулевой)
        last_none_zero = 0
        for i in range(1, len(density[j]) + 1):
            if density[j][-i] != 0:
                last_none_zero = -i + 1
                break
        datar = [step * (i + 0.5) for i in range(len(density[j][:last_none_zero]))]

        data_density_finish.append(density[j][:last_none_zero])
        datar_finish.append(datar)   

    return datar_finish, data_density_finish


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
    
#   Куда будет записана плотность и расстояние до шаровых слоёв    
    datar = []
    density = []
    writing_mode = 'w' # как мы открываем файл

    if working_mode == 1:
        for i in range(len(data_atom_sistem[3])):
            data_atom_sistem[4][i] = 0

        number_data_density = 1
        datar, density = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem, number_data_density)
        ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)

    if working_mode == 2:
        number_data_density = 1
#       Выделяем нужные нам атомы        
        atom_sort(data_atom_sistem, atom_name, system_center_core)
        datar, density = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem, number_data_density, system_center_core, False)
        ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)

        
    if working_mode == 3 or working_mode == 4: 
#       Делаем ли нормировку 
        norm = True
        if working_mode == 3:
            number_data_density = G
#           Выделяем атомы кремния            
            atom_si(data_atom_sistem, G, core, points)
            norm = True

        else:
            number_data_density = G + 1
#           Получаем массив с атомами разделённые по генерационным слоя            
            generation_layer(data_atom_sistem, G, core, points)
            norm = False
            ax.set_ylabel(r'$\rho$, г/$см^3$', fontsize=20)

        datar, density = density_profile__for_all_the_time_atom(file_trajectory_name, step, data_atom_sistem, number_data_density, system_center_core, norm)

        file = open(file_name, 'w') # очищаем файл перед записью
        file.close()
        writing_mode = 'a'

    [writing_to_file(datar[i], density[i], file_name, writing_mode, i+1) for i in range(len(density))]
    [ax.plot(datar[i], density[i]) for i in range(len(density))]
    ax.grid()
    plt.show()

density_profile()

