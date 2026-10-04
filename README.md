# Actividad 1 — Problemas de contorno unidimensionales: métodos de disparo

## Descripción

En esta actividad se resuelve un **problema de contorno unidimensional no lineal de segundo orden con condiciones de tipo Robin** mediante el **método de disparo**.

El parámetro de disparo se define como el valor de la solución en el extremo izquierdo y se determina utilizando dos métodos iterativos:

- Método de la **secante**.
- Método de **Newton**.

En ambos casos, el problema de valor inicial asociado se integra mediante el método de **Runge-Kutta de orden 4 (RK4)**.

## Parámetros utilizados

- Número de subintervalos: `N = 10`
- Tolerancia: `10^-6`
- Integrador: Runge-Kutta de orden 4

## Resultados

Ambos métodos convergen al mismo valor del parámetro de disparo:

`t* ≈ 0.2710333105`

La **secante** alcanza la solución en **7 iteraciones**, mientras que **Newton** necesita **5 iteraciones**.

Los órdenes de convergencia estimados son próximos a los valores teóricos. Sin embargo, el menor número de iteraciones del método de Newton no implica necesariamente un menor coste computacional, ya que requiere integrar un sistema de dimensión 4, frente al sistema de dimensión 2 utilizado con la secante.

## Análisis de convergencia

Los resultados muestran una elevada sensibilidad de ambos métodos a los valores iniciales. Esto se debe principalmente a la forma de la función de disparo y a la existencia de regiones en las que la integración numérica desborda.

Además, una discretización demasiado gruesa puede producir **raíces espurias**, que se identifican al refinar la malla.

## Comparación con `bvp5c`

Como referencia se utiliza el comando `bvp5c` de MATLAB.

Las soluciones obtenidas mediante el método de disparo presentan diferencias del orden de `10^-5` respecto a la solución proporcionada por `bvp5c`.

Esta diferencia se atribuye principalmente al **error de discretización del método RK4** y no a la tolerancia empleada en los métodos iterativos.

## Conclusión

Los métodos de la secante y Newton permiten resolver satisfactoriamente el problema de contorno mediante disparo y conducen al mismo parámetro de disparo. No obstante, los resultados ponen de manifiesto la importancia de la elección de los valores iniciales y del refinamiento de la discretización para garantizar una solución numérica fiable.
