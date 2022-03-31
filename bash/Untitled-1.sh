#/bin/bash

function_grompp(){
    gmx grompp -f mdrun.mdp -p c43g${1}.top -c c43g${1}_600.gro
}

function_mdrun(){
    gmx mdrun -v -cpt 5 -c c43g${1}_${2}_${3}.gro -nt 6
}

function_gyrate_begin(){
    echo 0 | gmx gyrate -o g${1}_${2}.xvg -p -b ${3}
}

# G; температура; длительность моделирования; начало в gyrate; 
function_start(){
    function_grompp $1 
    function_mdrun $1 $2 $3 
    function_gyrate_begin $1 $2 $4 
}

# 1 - G; 2 - начальная температура; 3 - конечная температура; 4 - шаг по температуре
function_MD(){
    cd G${1}
    Temperature=$2

    while [ $Temperature -le $3 ]
    do
        cd $Temperature
        function_start $1 $Temperature 2000 1000

        cd 100
        function_start $1 $Temperature 20 0

        let Temperature+=$4 
        cd ../..
        echo $Temperature
    done 
    cd ..
}

function_MD 5 350 350 50 
#function_MD 6 250 550 50