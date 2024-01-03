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
+ Modelo Antiangiogenesis ---
+ Modelo de Stepanova (stepanova4.f90)

Modelos de Tratamiento de tumores con Dormand-Prince:
+ Modelo de Norton-Simon (nortonDP45.f90)

Codigos gnuplot para generar graficos en formato png automaticamente a partir de los archivos de datos (*.dat) generados en fortran, estos son los mismos para RK2 y RK4:

+ Modelo de Efecto Proporcional (logkill.p)
+ Modelo de Efecto Proporcional con ecuación de Hill (logkill-H.p)
+ Modelo de Efecto Proporcional con Resistencia a los Fármacos (resistencia.p)
+ Modelo de Norton-Simon (norton.p)
+ Modelo Antiangiogenesis ---
+ Modelo de Stepanova (nova.p)

![Modelo log-kill](https://github.com/Cygnus000/tratamiento/blob/main/logkill.png)
![Modelo log-kill con ecuacion de Hill](https://github.com/Cygnus000/tratamiento/blob/main/logkill-H.png)
![Modelo log-kill con resistencia al farmaco](https://github.com/Cygnus000/tratamiento/blob/main/resistencia.png)
![Modelo norton-simon](https://github.com/Cygnus000/tratamiento/blob/main/norton.png)
![Modelo stepanova](https://github.com/Cygnus000/tratamiento/blob/main/nova.png)



