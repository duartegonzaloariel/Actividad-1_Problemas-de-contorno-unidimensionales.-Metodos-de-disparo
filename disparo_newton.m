% Laboratorio 1 - Metodo de disparo con Newton
% y'' + (y')^2 + y = -1/2 cos(x) y',  x en [0, pi]
% 2y(0) + y'(0) = 1/4,   -4y(pi) + 2y'(pi) = 2
clear; clc; close all;

a = 0;  b = pi;
N = 10;             % subintervalos
tol = 1e-6;         % tolerancia
maxiter = 50;
t0 = 0.2;           % valor inicial del parametro

t = t0;
T = [];  Fv = [];  dFv = [];
disparos = {};
convergio = false;

for k = 1:maxiter
    % Sistema de dimension 4: Y = [y, y', z, z']
    [x, Y] = rk4_sistema(@sistema, a, b, N, [t; 1/4 - 2*t; 1; -2]);
    F  = -4*Y(end,1) + 2*Y(end,2) - 2;
    dF = -4*Y(end,3) + 2*Y(end,4);

    T(end+1) = t;  Fv(end+1) = F;  dFv(end+1) = dF;
    disparos{end+1} = Y(:,1);

    if ~isfinite(F) || ~isfinite(dF) || dF == 0
        break
    end

    tn = t - F/dF;                       % paso de Newton

    if abs(tn - t) < tol
        convergio = true;
        t = tn;
        [x, Y] = rk4_sistema(@sistema, a, b, N, [t; 1/4 - 2*t; 1; -2]);
        T(end+1)   = t;
        Fv(end+1)  = -4*Y(end,1) + 2*Y(end,2) - 2;
        dFv(end+1) = -4*Y(end,3) + 2*Y(end,4);
        disparos{end+1} = Y(:,1);
        break
    end
    t = tn;
end

n = length(T);

% ---- TABLA DE ITERACIONES ----
fprintf('\nDisparo con Newton: t0 = %.4f, N = %d, tol = %.0e\n\n', t0, N, tol);
fprintf('%3s %14s %15s %13s %15s %7s\n', 'k', 't_k', 'F(t_k)', 'dF(t_k)', '|t_k-t_{k-1}|', 'p');
for k = 1:n
    if k == 1
        fprintf('%3d %14.10f %15.6e %13.6f %15s %7s\n', 0, T(1), Fv(1), dFv(1), '---', '---');
    else
        d1 = abs(T(k) - T(k-1));
        if k >= 4
            d2 = abs(T(k-1) - T(k-2));
            d3 = abs(T(k-2) - T(k-3));
            p = log(d1/d2) / log(d2/d3);
            fprintf('%3d %14.10f %15.6e %13.6f %15.6e %7.3f\n', k-1, T(k), Fv(k), dFv(k), d1, p);
        else
            fprintf('%3d %14.10f %15.6e %13.6f %15.6e %7s\n', k-1, T(k), Fv(k), dFv(k), d1, '---');
        end
    end
end

if convergio
    fprintf('\nConverge en %d iteraciones. Ultimo t = %.10f\n', n-1, T(end));
else
    fprintf('\nNO converge. Ultimo t = %.10f\n', T(end));
end

% ---- COCIENTES |e_{k+1}| / |e_k|^2 (convergencia cuadratica) ----
fprintf('\nCocientes |t_{k+1}-t*| / |t_k-t*|^2 (t* = ultimo t):\n');
for k = 2:n-2
    fprintf('  k = %d: %.4f\n', k-1, abs(T(k+1) - T(end)) / abs(T(k) - T(end))^2);
end

% ---- SOLUCION EN LOS NODOS ----
fprintf('\n%8s %12s %12s %12s %12s\n', 'x_j', 'y', 'y''', 'z', 'z''');
for j = 1:N+1
    fprintf('%8.4f %12.6f %12.6f %12.6f %12.6f\n', x(j), Y(j,1), Y(j,2), Y(j,3), Y(j,4));
end

% ---- CONDICIONES DE CONTORNO ----
fprintf('\n2y(0) + y''(0)     = %.6f\n', 2*Y(1,1) + Y(1,2));
fprintf('-4y(pi) + 2y''(pi) = %.6f\n', -4*Y(end,1) + 2*Y(end,2));
fprintf('F''(t*)            = %.6f\n', dFv(end));

% ---- GRAFICO ----
figure;
hold on;
leyenda = cell(1, n);
for k = 1:n-1
    plot(x, disparos{k}, '--', 'LineWidth', 1);
    leyenda{k} = sprintf('t_%d = %.4f', k-1, T(k));
end
plot(x, Y(:,1), 'k-o', 'LineWidth', 1.8, 'MarkerFaceColor', 'k', 'MarkerSize', 5);
leyenda{n} = sprintf('Solucion (t = %.6f)', T(end));
xlabel('x');
ylabel('y(x)');
title(sprintf('Disparo con Newton: disparos sucesivos (N = %d)', N));
legend(leyenda, 'Location', 'southwest');
grid on;
hold off;
print('-dpng', '-r300', 'disparo_newton.png');

% ================= FUNCIONES LOCALES =================

function dY = sistema(x, Y)
    % Y = [y; y'; z; z']
    dY = zeros(4,1);
    dY(1) = Y(2);
    dY(2) = -Y(2)^2 - Y(1) - 0.5*cos(x)*Y(2);
    dY(3) = Y(4);
    dY(4) = -Y(3) - (2*Y(2) + 0.5*cos(x))*Y(4);
end

function [x, Y] = rk4_sistema(f, a, b, N, Y0)
    % Runge-Kutta de orden 4 para sistemas, N subintervalos
    h = (b - a)/N;
    x = linspace(a, b, N+1)';
    Y = zeros(N+1, length(Y0));
    Y(1,:) = Y0';
    for j = 1:N
        yj = Y(j,:)';
        k1 = f(x(j), yj);
        k2 = f(x(j) + h/2, yj + h/2*k1);
        k3 = f(x(j) + h/2, yj + h/2*k2);
        k4 = f(x(j) + h, yj + h*k3);
        Y(j+1,:) = (yj + h/6*(k1 + 2*k2 + 2*k3 + k4))';
    end
end