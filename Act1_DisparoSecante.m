% Actividad 1 - Codigo 2
% Metodo de disparo con el esquema de la secante
% y'' + (y')^2 + y = -1/2 cos(x) y',  x en [0, pi]
% 2y(0) + y'(0) = 1/4,   -4y(pi) + 2y'(pi) = 2

clear; clc; close all;

% ---- Guardar la salida de la consola en un archivo de texto ----
if exist('Act1_Secante_resultados.txt', 'file'), delete('Act1_Secante_resultados.txt'); end
diary('Act1_Secante_resultados.txt');

% ---- Datos del problema ----
a = 0;  b = pi;
N = 10;                 % subintervalos de RK4
tol = 1e-6;
maxiter = 50;

CI    = @(t) [t; 1/4 - 2*t];                 % y(0) = t, y'(0) = 1/4 - 2t
Fcont = @(yb, dyb) -4*yb + 2*dyb - 2;        % F(t) = -4y(pi) + 2y'(pi) - 2

% ---- Disparo con secante (valores iniciales de referencia) ----
t0 = 0;  t1 = 0.2;
[x, Y, t, iter, hist, conv] = disparo_secante(@f_problema, a, b, CI, Fcont, ...
                                              t0, t1, N, tol, maxiter);

fprintf('Disparo con secante: t0 = %.2f, t1 = %.2f\n', t0, t1);
% Estimacion del orden de convergencia (ACOC)
e = hist(:,4);
p = nan(size(e));
for i = 4:length(e)
    p(i) = log(e(i)/e(i-1)) / log(e(i-1)/e(i-2));
end
fprintf('%-4s %-16s %-16s %-14s %-8s\n', 'k', 't_k', 'F(t_k)', '|t_k - t_k-1|', 'p');
for i = 1:size(hist, 1)
    fprintf('%-4d %-16.10f %-16.6e %-14.6e %-8.3f\n', hist(i,1), hist(i,2), hist(i,3), hist(i,4), p(i));
end
fprintf('Iteraciones: %d   Ultimo t: %.10f\n\n', iter, t);

fprintf('%-6s %-12s %-12s\n', 'x_j', 'y(x_j)', "y'(x_j)");
for j = 1:N+1
    fprintf('%-6.4f %-12.6f %-12.6f\n', x(j), Y(j,1), Y(j,2));
end
fprintf('Comprobacion CC en 0:  2y(0) + y''(0)      = %.10f\n', 2*Y(1,1) + Y(1,2));
fprintf('Comprobacion CC en pi: -4y(pi) + 2y''(pi) = %.10f\n\n', -4*Y(end,1) + 2*Y(end,2));

% Raices de F con N = 10: t* y las dos raices espurias obtenidas con
% los valores iniciales (-2.5, -2.7) en secante y t0 = 0.05 en Newton
raices = [t, -2.7633317317, 0.4224562914];

% ---- GRAFICO 1: disparos sucesivos y solucion ----
figure;
hold on;
colores = lines(size(hist,1));
for i = 1:size(hist,1)-1
    [xs, Ys] = rk4_sistema(@f_problema, a, b, CI(hist(i,2)), N);
    plot(xs, Ys(:,1), '--', 'Color', colores(i,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('t_{%d} = %.4f', hist(i,1), hist(i,2)));
end
plot(x, Y(:,1), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k', ...
     'DisplayName', sprintf('Solucion (t = %.6f)', t));
xlabel('x');
ylabel('y(x)');
title('Disparo con secante: disparos sucesivos (N = 10)');
legend('Location', 'eastoutside');
grid on;
hold off;
saveas(gcf, 'Secante_solucion.png');

% ---- GRAFICO 2: funcion F(t) discreta ----
t_vals = linspace(-3, 0.45, 400);
F_vals = nan(size(t_vals));
for i = 1:length(t_vals)
    [~, Ys] = rk4_sistema(@f_problema, a, b, CI(t_vals(i)), N);
    Fi = Fcont(Ys(end,1), Ys(end,2));
    if isfinite(Fi) && abs(Fi) < 50
        F_vals(i) = Fi;
    end
end
figure;
hold on;
plot(t_vals, F_vals, 'b-', 'LineWidth', 1.8, 'DisplayName', 'F(t)');
yline(0, 'k--', 'HandleVisibility', 'off');
plot(raices(1), 0, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.47 0.67 0.19], ...
     'MarkerEdgeColor', 'k', 'DisplayName', sprintf('t* = %.6f', raices(1)));
plot(raices(2:3), [0 0], 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.85 0.33 0.10], ...
     'MarkerEdgeColor', 'k', 'DisplayName', 'Raices espurias');
xlabel('t');
ylabel('F(t)');
title('Funcion de disparo F(t) discretizada con RK4 (N = 10)');
legend('Location', 'north');
ylim([-6 8]);
grid on;
hold off;
saveas(gcf, 'Funcion_F.png');

% ---- COMPROBACION DE RAICES ESPURIAS: refinamiento de la malla ----
fprintf('%-6s %-14s %-14s %-14s\n', 'N', 'F(t*)', 'F(0.4224563)', 'F(-2.7633317)');
for Nr = [10 20 40 80]
    Fr = zeros(1, 3);
    for i = 1:3
        [~, Ys] = rk4_sistema(@f_problema, a, b, CI(raices(i)), Nr);
        Fr(i) = Fcont(Ys(end,1), Ys(end,2));
    end
    fprintf('%-6d %-14.4e %-14.4e %-14.4e\n', Nr, Fr(1), Fr(3), Fr(2));
end
fprintf('\n');

% ---- ANALISIS DE CONVERGENCIA: pares (t0, t1) ----
t_ref  = t;                  % raiz principal
pares = [ 0.2  0.3;  0.25 0.3;  0.3  0.4;  0    0.2;
          0.1  0.4;  0    0.1; -0.5  0;   -1    0.2;
         -1    0;   -2   -1;   -2.5 -2.7; -3   -2.8];
fprintf('%-8s %-8s %-6s %-16s %-10s\n', 't0', 't1', 'iter', 'ultimo t', 'resultado');
for i = 1:size(pares,1)
    [~, ~, tf, it, ~, cv] = disparo_secante(@f_problema, a, b, CI, Fcont, ...
                                            pares(i,1), pares(i,2), N, tol, maxiter);
    if cv && abs(tf - t_ref) < 1e-4
        res = 'raiz t*';
    elseif cv
        res = 'raiz espuria';
    else
        res = 'no converge';
    end
    fprintf('%-8.2f %-8.2f %-6d %-16.10f %-10s\n', pares(i,1), pares(i,2), it, tf, res);
end

% ---- GRAFICO 3: mapa de convergencia en el plano (t0, t1) ----
v = linspace(-3.5, 0.5, 81);
mapa = zeros(length(v));
iters = nan(length(v));
for i = 1:length(v)
    for j = 1:length(v)
        if i == j, continue; end
        [~, ~, tf, it, ~, cv] = disparo_secante(@f_problema, a, b, CI, Fcont, ...
                                                v(j), v(i), N, tol, maxiter);
        if cv && abs(tf - t_ref) < 1e-4
            mapa(i,j) = 1;  iters(i,j) = it;
        elseif cv
            mapa(i,j) = 2;  iters(i,j) = it;
        end
    end
end
figure;
imagesc(v, v, mapa);
set(gca, 'YDir', 'normal');
colormap([0.85 0.85 0.85; 0.20 0.45 0.85; 0.90 0.45 0.15]);
caxis([-0.5 2.5]);
cb = colorbar('Ticks', [0 1 2], 'TickLabels', {'No converge', 'Raiz t*', 'Raiz espuria'});
xlabel('t_0');
ylabel('t_1');
title('Convergencia del disparo con secante segun (t_0, t_1)');
axis square;
saveas(gcf, 'Secante_mapa.png');

fprintf('\nPorcentaje de pares que convergen a t*: %.1f %%\n', 100*mean(mapa(:)==1));
fprintf('Iteraciones (convergencia a t*): min = %d, max = %d, mediana = %g\n', ...
        min(iters(mapa==1)), max(iters(mapa==1)), median(iters(mapa==1)));

% ---- Cerrar el archivo de resultados ----
diary off;