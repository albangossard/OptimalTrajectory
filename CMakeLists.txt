#Déclaration du projet
project("OptimalTrajectory")


#Déclaration de l'exécutable
add_executable(
        OptimalTrajectory
        grid.h
        grid.cxx
        set.h
        set.cxx
        pathFinding.cxx
)