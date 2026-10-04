% Laboratorio 1 - Comparacion: disparo con secante, disparo con Newton y bvp5c
clear; clc; close all;

a = 0;  b = pi;
N = 10;
tol = 1e-6;
maxiter = 50;
x = linspace(a, b, N+1)';

% ---- Disparo con secante y con Newton ----
[tS, itS, YS] = disparo_secante(0, 0.2, a, b, N, tol, maxiter);
[tN, itN, YN] = disparo_newton(0.2, a, b, N, tol, maxiter);
yS = YS(:,1);
yN = YN(:,1);

% ---- bvp5c ----
solinit = bvpinit(linspace(a, b, N+1), [0 0]);
solD = bvp5c(@sistema2, @contorno, solinit);                 % tolerancias por defecto
opts = bvpset('RelTol', 1e-8, 'AbsTol', 1e-10);
solF = bvp5c(@sistema2, @contorno, solinit, opts);           % tolerancias exigentes
YD = deval(solD, x')';
YF = deval(solF, x')';
yD = YD(:,1);
yF = YF(:,1);

% ---- RESUMEN ----
fprintf('Secante: %d iteraciones, t = %.10f\n', itS, tS);
fprintf('Newton:  %d iteraciones, t = %.10f\n', itN, tN);
fprintf('|t_secante - t_newton| = %.3e\n', abs(tS - tN));
fprintf('bvp5c por defecto:  y(0) = %.10f, %d puntos de malla\n', yD(1), numel(solD.x));
fprintf('bvp5c exigente:     y(0) = %.10f, %d puntos de malla\n', yF(1), numel(solF.x));

% ---- TABLA EN LOS NODOS ----
fprintf('\n%8s %12s %12s %12s %11s %11s %11s\n', ...
        'x_j', 'Secante', 'Newton', 'bvp5c', '|S-N|', '|S-B|', '|N-B|');
for j = 1:N+1
    fprintf('%8.4f %12.6f %12.6f %12.6f %11.2e %11.2e %11.2e\n', ...
            x(j), yS(j), yN(j), yF(j), abs(yS(j)-yN(j)), abs(yS(j)-yF(j)), abs(yN(j)-yF(j)));
end

fprintf('\nDiferencias maximas en los nodos:\n');
fprintf('  Secante - Newton:                 %.3e\n', max(abs(yS - yN)));
fprintf('  Secante - bvp5c (exigente):       %.3e\n', max(abs(yS - yF)));
fprintf('  Newton  - bvp5c (exigente):       %.3e\n', max(abs(yN - yF)));
fprintf('  Secante - bvp5c (por defecto):    %.3e\n', max(abs(yS - yD)));
fprintf('  bvp5c por defecto - exigente:     %.3e\n', max(abs(yD - yF)));

% ---- ERROR DEL DISPARO AL REFINAR LA MALLA ----
fprintf('\n%6s %16s %14s %10s\n', 'N', 't (Newton)', '|t - y_B(0)|', 'cociente');
errAnt = NaN;
for NN = [10 20 40]
    tNN = disparo_newton(0.2, a, b, NN, tol, maxiter);
    err = abs(tNN - yF(1));
    fprintf('%6d %16.10f %14.3e %10.2f\n', NN, tNN, err, errAnt/err);
    errAnt = err;
end

% ---- GRAFICO ----
xf = linspace(a, b, 200);
Yf = deval(solF, xf);

figure('Position', [100 100 1000 400]);

subplot(1,2,1);
hold on;
plot(xf, Yf(1,:), 'k-', 'LineWidth', 1.5);
plot(x, yS, 'o', 'Color', [0 0.45 0.74], 'MarkerSize', 8, 'LineWidth', 1.5);
plot(x, yN, 'x', 'Color', [0.85 0.33 0.10], 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('x');
ylabel('y(x)');
title('Soluciones aproximadas');
legend('bvp5c', 'Disparo con secante', 'Disparo con Newton', 'Location', 'north');
grid on;
hold off;

subplot(1,2,2);
semilogy(x, max(abs(yS - yF), eps), 'o-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
hold on;
semilogy(x, max(abs(yN - yF), eps), 'x--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
semilogy(x, max(abs(yS - yN), eps), 's-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
xlabel('x');
ylabel('Diferencia absoluta');
title('Diferencias en los nodos');
legend('|Secante - bvp5c|', '|Newton - bvp5c|', '|Secante - Newton|', 'Location', 'east');
grid on;
hold off;

print('-dpng', '-r300', 'comparacion.png');

% ================= FUNCIONES LOCALES =================

function [t, it, Y] = disparo_secante(t0, t1, a, b, N, tol, maxiter)
    [~, Y] = rk4_sistema(@sistema2, a, b, N, [t0; 1/4 - 2*t0]);
    F0 = -4*Y(end,1) + 2*Y(end,2) - 2;
    [~, Y] = rk4_sistema(@sistema2, a, b, N, [t1; 1/4 - 2*t1]);
    F1 = -4*Y(end,1) + 2*Y(end,2) - 2;
    t = t1;  it = 0;
    for k = 1:maxiter
        t = t1 - F1*(t1 - t0)/(F1 - F0);
        [~, Y] = rk4_sistema(@sistema2, a, b, N, [t; 1/4 - 2*t]);
        F = -4*Y(end,1) + 2*Y(end,2) - 2;
        it = k;
        if abs(t - t1) < tol
            break
        end
        t0 = t1;  F0 = F1;
        t1 = t;   F1 = F;
    end
end

function [t, it, Y] = disparo_newton(t0, a, b, N, tol, maxiter)
    t = t0;  it = 0;
    for k = 1:maxiter
        [~, Y] = rk4_sistema(@sistema4, a, b, N, [t; 1/4 - 2*t; 1; -2]);
        F  = -4*Y(end,1) + 2*Y(end,2) - 2;
        dF = -4*Y(end,3) + 2*Y(end,4);
        tn = t - F/dF;
        it = k;
        if abs(tn - t) < tol
            t = tn;
            [~, Y] = rk4_sistema(@sistema4, a, b, N, [t; 1/4 - 2*t; 1; -2]);
            break
        end
        t = tn;
    end
end

function dY = sistema2(x, Y)
    % Y = [y; y']
    dY = zeros(2,1);
    dY(1) = Y(2);
    dY(2) = -Y(2)^2 - Y(1) - 0.5*cos(x)*Y(2);
end

function dY = sistema4(x, Y)
    % Y = [y; y'; z; z']
    dY = zeros(4,1);
    dY(1) = Y(2);
    dY(2) = -Y(2)^2 - Y(1) - 0.5*cos(x)*Y(2);
    dY(3) = Y(4);
    dY(4) = -Y(3) - (2*Y(2) + 0.5*cos(x))*Y(4);
end

function res = contorno(ya, yb)
    % 2y(0) + y'(0) = 1/4,   -4y(pi) + 2y'(pi) = 2
    res = [ 2*ya(1) + ya(2) - 1/4;
           -4*yb(1) + 2*yb(2) - 2 ];
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