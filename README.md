# MCOC2020-P1

# Integración de ecuaciones diferenciales

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega1/balistica.png)

# Primeras predicciones con la EDM básica del satélite

*Con el sátelite a 700 km de la superficie de la Tierra, se fue probando para distintas velocidades tangeciales en la dirección y, y se llego a una velocidad de 7000 m/s, en la cual el satélite se estabiliza sobre la atmósfera para varias órbitas completas:

*Con una velocidad de 5000 m/s se tuvo:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/v_5000.png)


*Con una velocidad de 6000 m/s se tuvo:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/v_6000.png)


*Con una velocidad de 6500 m/s se tuvo:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/v_6500.png)

*Finalmente, con una velocidad de 7000 m/s se tuvo:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/v_7000.png)

*En cuánto a la posición x(t) y(t) z(t), se observó:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/distancia_satelite.png)

*Y viendo su distancia hacia la Tierra y la atmósfera:

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega2/distancia_tierra.png)

# Estudio de convergencia Método de Euler

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega4/oscilador_armonico.png)

# Mejoras al modelo y estudio de convergencia

## Posición predicha versus la posición real

* Con un modelo básico usando la funcion ODEINT, se compara con la posición real del satélite  y su respectiva deriva, donde la deriva final es alrededor de los 877 km, es decir muy lejos todavía de la posición real.

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/1.posicion_real_pred.png)

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/1.deriva_ode_real.png)


## Comparación ODEINT y EULERINT

* Se presenta una comparación del modelo ocupando odeint versus otro modelo que ocupa la función EULERINT con 1 subdivisión.

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/2.deriva_ode_euler.png)

* El tiempo que tardo cada función en producir resultados fue de 0.205 segundos para ODEINT y 0.732 segundos para EULERINT.
* La deriva máxima del modelo con ODEINT es la anteriormente mencionada de 877 km y del modelo con EULERINT de 19367 km.

## Predicción con EULERINT

### 1 Subdivisión
![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/3.deriva_euler_real_1.png)

* Tiempo de ejecución: 1.01 s

### 200 Subdivisiones
![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/3.deriva_euler_real_200.png)

* Tiempo de ejecución: 177.54 s

### 500 Subdivisiones
![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/3.deriva_euler_real_500.png)

* Tiempo de ejecución: 297.88 s

## Mejoras al modelo

* Finalmente se le agregan correcciones al modelo, y se presenta la prediccion_edm_avanzada y se replican los mismos gráficos anteriores, donde se puede apreciar que la deriva del modelo disminuye considerablemente agregando las correcciones J2 y J3 del modelo geopotencial gravitatorio.

* La deriva del modelo usando ODEINT con las correcciones fue alrededor de 80 km

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/4.posicion_real_predJ2J3.png)

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/4.deriva_ode_realJ2J3.png)

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/4.deriva_ode_eulerJ2J3.png)

![Alt Text](https://github.com/raimolid/MCOC2020-P1/blob/master/Entrega5/4.deriva_euler_real_1J2J3.png)

* Para la presentación final del código se presenta la carpeta Entrega_final con el código completo a implementar

