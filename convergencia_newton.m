% Laboratorio 1 - Disparo con Newton: convergencia segun t0
clear; clc; close all;

a = 0;  b = pi;
N = 10;
tol = 1e-6;
maxiter = 50;
tstar = 0.2710333108;      % raiz obtenida con disparo_newton.m

% ---- 1) TABLA: valores iniciales seleccionados ----
t0s = [0.25 0.30 0.20 0.35 0.10 0.40 0.06 0.05 0.00 -1.00 -1.50 -2.00 -2.70 -3.00];
fprintf('%8s %16s %6s %18s   %s\n', 't0', 't1', 'Iter', 'Ultimo t', 'Resultado');
for i = 1:length(t0s)
    [t, it, conv, t1] = newton_disparo(t0s(i), a, b, N, tol, maxiter);
    fprintf('%8.2f %16.6f %6d %18.10f   %s\n', t0s(i), t1, it, t, clasifica(conv, t, tstar));
end

% ---- 2) BARRIDO: 81 valores de t0 en [-3.5, 0.5] ----
G = linspace(-3.5, 0.5, 81);
[res, its] = barrido(G, a, b, N, tol, maxiter, tstar);
n = numel(G);
fprintf('\nBarrido de %d valores de t0 en [-3.5, 0.5]:\n', n);
fprintf('  Converge a t*:      %3d (%.1f %%)\n', sum(res == 1), 100*sum(res == 1)/n);
fprintf('  Raiz espuria:       %3d (%.1f %%)\n', sum(res == 2), 100*sum(res == 2)/n);
fprintf('  No converge:        %3d (%.1f %%)\n', sum(res == 0), 100*sum(res == 0)/n);
fprintf('  Iteraciones (t*):   min %d, max %d, mediana %g\n', ...
        min(its(res == 1)), max(its(res == 1)), median(its(res == 1)));
fprintf('  t0 que convergen a t*:      %s\n', mat2str(G(res == 1), 4));
fprintf('  t0 que van a raiz espuria:  %s\n', mat2str(G(res == 2), 4));

% ---- 3) BUSQUEDA FINA del intervalo que converge a t* ----
G2 = 0:1e-4:0.45;
res2 = barrido(G2, a, b, N, tol, maxiter, tstar);
i1 = find(res2 == 1, 1, 'first');
i2 = find(res2 == 1, 1, 'last');
fprintf('\nBusqueda fina (paso 1e-4) en [0, 0.45]:\n');
fprintf('  Converge a t* para t0 en [%.4f, %.4f]\n', G2(i1), G2(i2));
fprintf('  Intervalo sin huecos: %d\n', all(res2(i1:i2) == 1));

% ---- 4) DATOS para explicar el caso t0 = 0 ----
[~, Y] = rk4_sistema(@sistema, a, b, N, [0; 1/4; 1; -2]);
F0  = -4*Y(end,1) + 2*Y(end,2) - 2;
dF0 = -4*Y(end,3) + 2*Y(end,4);
fprintf('\nt0 = 0:  F(0) = %.6f,  dF(0) = %.6f,  t1 = %.6f\n', F0, dF0, -F0/dF0);

% Mayor t para el que RK4 devuelve un valor finito de F
G4 = 0.40:1e-4:0.50;
fin = false(size(G4));
for i = 1:numel(G4)
    [~, Y] = rk4_sistema(@sistema, a, b, N, [G4(i); 1/4 - 2*G4(i); 1; -2]);
    fin(i) = all(isfinite(Y(end,:)));
end
fprintf('Mayor t con F finita (paso 1e-4): %.4f\n', G4(find(fin, 1, 'last')));

% Derivada en la zona plana
for tt = [0 0.20 0.40 -1.00 -1.47 -1.50]
    [~, Y] = rk4_sistema(@sistema, a, b, N, [tt; 1/4 - 2*tt; 1; -2]);
    fprintf('dF(%6.2f) = %12.6f\n', tt, -4*Y(end,3) + 2*Y(end,4));
end

% ---- 5) GRAFICO: iteraciones segun t0 ----
G3 = linspace(-3.5, 0.5, 401);
[res3, its3] = barrido(G3, a, b, N, tol, maxiter, tstar);

figure;
hold on;
plot(G3(res3 == 0), its3(res3 == 0), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 10);
plot(G3(res3 == 1), its3(res3 == 1), 'o', 'Color', [0 0.45 0.74], ...
     'MarkerFaceColor', [0 0.45 0.74], 'MarkerSize', 5);
plot(G3(res3 == 2), its3(res3 == 2), 's', 'Color', [0.85 0.33 0.10], ...
     'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5);
xlabel('t_0');
ylabel('Iteraciones');
title(sprintf('Disparo con Newton: resultado segun t_0 (N = %d)', N));
legend('No converge', 'Converge a t*', 'Converge a raiz espuria', 'Location', 'north');
grid on;
hold off;
print('-dpng', '-r300', 'convergencia_newton.png');

% ================= FUNCIONES LOCALES =================

function [t, it, conv, t1] = newton_disparo(t0, a, b, N, tol, maxiter)
    % Devuelve el ultimo t, el numero de iteraciones, si converge y el primer iterado
    t = t0;  it = 0;  conv = false;  t1 = NaN;
    for k = 1:maxiter
        [~, Y] = rk4_sistema(@sistema, a, b, N, [t; 1/4 - 2*t; 1; -2]);
        F  = -4*Y(end,1) + 2*Y(end,2) - 2;
        dF = -4*Y(end,3) + 2*Y(end,4);
        if ~isfinite(F) || ~isfinite(dF) || dF == 0
            return
        end
        % Criterio de parada: la condicion de contorno en b se cumple
        if abs(F) < tol
            conv = true;
            return
        end
        t = t - F/dF;                    % paso de Newton
        it = k;
        if k == 1
            t1 = t;
        end
    end
end

function [res, its] = barrido(G, a, b, N, tol, maxiter, tstar)
    % res: 1 = converge a t*, 2 = converge a raiz espuria, 0 = no converge
    res = zeros(size(G));
    its = zeros(size(G));
    for i = 1:numel(G)
        [t, it, conv] = newton_disparo(G(i), a, b, N, tol, maxiter);
        its(i) = it;
        if conv && abs(t - tstar) < 1e-4
            res(i) = 1;
        elseif conv
            res(i) = 2;
        end
    end
end

function s = clasifica(conv, t, tstar)
    if conv && abs(t - tstar) < 1e-4
        s = 'Converge a t*';
    elseif conv
        s = 'Converge a raiz espuria';
    else
        s = 'No converge';
    end
end

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