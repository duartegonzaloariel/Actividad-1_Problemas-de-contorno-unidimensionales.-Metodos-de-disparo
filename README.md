# Actividad 1 — Problemas de contorno unidimensionales: métodos de disparo

## Descripción

En esta actividad se resuelve un problema de contorno unidimensional no lineal de segundo orden con condiciones de tipo Robin mediante el método de disparo.

El parámetro de disparo se define como el valor de la solución en el extremo izquierdo y se determina utilizando dos métodos iterativos:

- Método de la secante.
- Método de Newton.

En ambos casos, el problema de valor inicial asociado se integra mediante el método de Runge-Kutta de orden 4 (RK4).

## Parámetros utilizados

- Número de subintervalos: `N = 10`
- Tolerancia: `10^-6`
- Criterio de parada: `|F(t_k)| < tol`, es decir, la condición de contorno en el extremo derecho se cumple con un residuo menor que la tolerancia
- Número máximo de iteraciones: `50`
- Integrador: Runge-Kutta de orden 4

## Resultados

Ambos métodos convergen a la misma solución dentro de la tolerancia, con `t ≈ 0.271033`:

| Método  | Valores iniciales    | Iteraciones | Último `t`     |
|---------|----------------------|-------------|----------------|
| Secante | `t0 = 0`, `t1 = 0.2` | 6           | `0.2710333001` |
| Newton  | `t0 = 0.2`           | 4           | `0.2710333108` |

Los órdenes de convergencia estimados son compatibles con los valores teóricos. Sin embargo, el menor número de iteraciones del método de Newton no implica un menor coste computacional, ya que requiere integrar un sistema de dimensión 4, frente al sistema de dimensión 2 utilizado con la secante.

## Análisis de convergencia

Los resultados muestran una elevada sensibilidad de ambos métodos a los valores iniciales. Esto se debe principalmente a la forma de la función de disparo y a la existencia de regiones en las que la integración numérica desborda.

Además, una discretización demasiado gruesa puede producir raíces espurias, que se identifican al refinar la malla.

## Comparación con `bvp5c`

Como referencia se utiliza el comando `bvp5c` de MATLAB.

Las soluciones obtenidas mediante el método de disparo presentan diferencias del orden de `10^-5` respecto a la solución proporcionada por `bvp5c`.

Esta diferencia se atribuye principalmente al error de discretización del método RK4 y no a la tolerancia empleada en los métodos iterativos.

## Archivos

| Archivo                  | Contenido                                                             |
|--------------------------|-----------------------------------------------------------------------|
| `rk4_sistema.m`          | Runge-Kutta de orden 4 para sistemas de primer orden                  |
| `f_problema.m`           | Sistema de primer orden asociado a la ecuación diferencial            |
| `disparo_secante.m`      | Método de disparo con el esquema de la secante                        |
| `Act1_DisparoSecante.m`  | Aproximación con secante y análisis de convergencia según `(t0, t1)`  |
| `disparo_newton.m`       | Aproximación con el método de disparo con Newton                      |
| `convergencia_newton.m`  | Análisis de convergencia de Newton según `t0`                         |
| `comparacion.m`          | Comparación entre secante, Newton y `bvp5c`                           |

## Conclusión

Los métodos de la secante y Newton permiten resolver satisfactoriamente el problema de contorno mediante disparo y conducen a la misma solución dentro de la tolerancia. No obstante, los resultados ponen de manifiesto la importancia de la elección de los valores iniciales y del refinamiento de la discretización para garantizar una solución numérica fiable.
