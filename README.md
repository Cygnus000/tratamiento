# Tratamiento

Codigos fortran 90 que implementan el metodo Runge Kutta de orden 2 para resolver Sistemas de Ecuaciones Diferenciales Ordinarias (EDOs) de diferentes modelos de tratamiento de tumores. Los modelos resueltos son:

+ Modelo de Efecto Proporcional (logkill.f90)
+ Modelo de Efecto Proporcional con ecuación de Hill (logkill-Hill.f90)
+ Modelo de Efecto Proporcional con Resistencia a los Fármacos (logkillResistencia.f90)
+ Modelo de Norton-Simon (norton.f90)
+ Modelo Antiogenesis ---
+ Modelo de Stepanova (stepanova.f90)

Modelos de Tratamiento de tumores con Runge Kutta de orden 4:

+ Modelo de Efecto Proporcional (logkill4.f90)
+ Modelo de Efecto Proporcional con ecuación de Hill (logkill-Hill4.f90)
+ Modelo de Efecto Proporcional con Resistencia a los Fármacos (logkillResistencia4.f90)
+ Modelo de Norton-Simon (norton4.f90)
+ Modelo Antiangiogenesis ang4.f90
+ Modelo de Stepanova (stepanova4.f90)

Modelos de Tratamiento de tumores con Dormand-Prince:
+ Modelo de Efecto Proporcional (logkillDP54.f90)
+ Modelo de Efecto Proporcional con ecuación de Hill (logkill-HillDP54.f90)
+ Modelo de Efecto Proporcional con Resistencia a los Fármacos (logkillResistenciaDP54.f90)
+ Modelo de Norton-Simon (nortonDP54.f90)
+ Modelo Antiogenesis ---
+ Modelo de Stepanova (stepanovaDP54.f90)

Codigos gnuplot para generar graficos en formato eps y latex automaticamente a partir de los archivos de datos (*.dat) generados en fortran:

+ Modelo de Efecto Proporcional (logkill.gplot)
+ Modelo de Efecto Proporcional con ecuación de Hill (logkill-H.gplot)
+ Modelo de Efecto Proporcional con Resistencia a los Fármacos (resistencia.gplot)
+ Modelo de Norton-Simon (norton.gplot)
+ Modelo Antiangiogenesis (angiogenesis.gplot)
+ Modelo de Stepanova (nova.gplot)

![Modelo log-kill](https://github.com/Cygnus000/tratamiento/blob/main/logkill.png)
![Modelo log-kill con ecuacion de Hill](https://github.com/Cygnus000/tratamiento/blob/main/logkill-H.png)
![Modelo log-kill con resistencia al farmaco](https://github.com/Cygnus000/tratamiento/blob/main/resistencia.png)
![Modelo norton-simon](https://github.com/Cygnus000/tratamiento/blob/main/norton.png)
![Modelo angiogenesis](https://github.com/Cygnus000/tratamiento/blob/main/angiogenesis.png)
![Modelo stepanova](https://github.com/Cygnus000/tratamiento/blob/main/nova.png)



