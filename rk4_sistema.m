function [x, Y] = rk4_sistema(f, a, b, Y0, N)
% Metodo de Runge-Kutta de orden 4 para sistemas Y' = f(x, Y)
% Entradas: f (funcion), [a, b] (intervalo), Y0 (cond. inicial), N (subintervalos)
% Salidas:  x (nodos, (N+1)x1), Y (solucion en los nodos, (N+1)x m)
h = (b - a) / N;
x = (a:h:b)';
Y = zeros(N+1, length(Y0));
Y(1,:) = Y0(:)';
for j = 1:N
    yj = Y(j,:)';
    k1 = f(x(j),       yj);
    k2 = f(x(j) + h/2, yj + h/2*k1);
    k3 = f(x(j) + h/2, yj + h/2*k2);
    k4 = f(x(j) + h,   yj + h*k3);
    Y(j+1,:) = (yj + h/6*(k1 + 2*k2 + 2*k3 + k4))';
end
end