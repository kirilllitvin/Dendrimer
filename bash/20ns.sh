#/bin/bash


function_grompp(){
    gmx grompp -f mdrun.mdp -p c43g${1}.top -c c43g${1}_b.gro
}

function_mdrun(){
    gmx mdrun -v -cpt 5 -c c43g${1}_${2}.gro -nt 6
}

function_gyrate_begin(){
    echo 0 | gmx gyrate -o g${1}_${2}.xvg -p -b ${3}
}

# G; температура; длительность моделирования; начало в gyrate; 
function_start(){
    function_grompp $1 
    function_mdrun $1 $2 $3 
#    function_gyrate_begin $1 $2 $4 
}

# 1 - G; 2 - начальная температура; 3 - конечная температура; 4 - шаг по температуре
function_MD(){
    cd ../..
    cd G6
    Temperature=$2
    ls
    while [ $Temperature -le $3 ]
    do  
        cd ${Temperature}k_51ns
        function_start $1 $Temperature 51ns 
        ls 
        cd ..
        
        echo $Temperature
        let Temperature+=$4 
        
    done 
    cd ..
}

function_MD 6 300 600 300
#function_MD 6 250 550 50