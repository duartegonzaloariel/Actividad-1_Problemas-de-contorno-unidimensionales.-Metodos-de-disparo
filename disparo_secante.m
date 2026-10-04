function [x, Y, t, iter, historial, conv] = disparo_secante(f, a, b, CI, Fcont, t0, t1, N, tol, maxiter)
% Metodo de disparo con el esquema de la secante
% Entradas:
%   f       - sistema de primer orden Y' = f(x, Y)
%   [a, b]  - intervalo
%   CI      - condiciones iniciales en funcion del parametro: Y(a) = CI(t)
%   Fcont   - condicion de contorno en b: F(t) = Fcont(y(b), y'(b))
%   t0, t1  - valores iniciales de la recurrencia
%   N       - numero de subintervalos de RK4
%   tol     - tolerancia en |F(t_k)|
%   maxiter - numero maximo de iteraciones
% Salidas:
%   x, Y      - nodos y solucion aproximada para el ultimo t
%   t         - ultimo valor del parametro
%   iter      - numero de iteraciones realizadas
%   historial - filas [k, t_k, F(t_k), |t_k - t_{k-1}|]
%   conv      - 1 si converge, 0 en caso contrario

[~, Y] = rk4_sistema(f, a, b, CI(t0), N);
F0 = Fcont(Y(end,1), Y(end,2));
[x, Y] = rk4_sistema(f, a, b, CI(t1), N);
F1 = Fcont(Y(end,1), Y(end,2));

historial = [0, t0, F0, NaN;
             1, t1, F1, abs(t1 - t0)];
conv = 0;
iter = 0;
t = t1;

while iter < maxiter
    iter = iter + 1;
    % Recurrencia de la secante
    t2 = t1 - F1*(t1 - t0)/(F1 - F0);
    % Nuevo disparo
    [x, Y] = rk4_sistema(f, a, b, CI(t2), N);
    F2 = Fcont(Y(end,1), Y(end,2));
    historial = [historial; iter+1, t2, F2, abs(t2 - t1)];
    t = t2;
    if ~isfinite(t2) || ~isfinite(F2)
        break;              % el disparo desborda: no converge
    end
    % Criterio de parada: la condicion de contorno en b se cumple
    if abs(F2) < tol
        conv = 1;
        break;
    end
    t0 = t1;  F0 = F1;
    t1 = t2;  F1 = F2;
end
end