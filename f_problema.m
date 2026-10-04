function dY = f_problema(x, Y)
% Sistema de primer orden asociado a la EDO
%   y'' = -(y')^2 - y - (1/2) cos(x) y'
% con Y = [y1; y2] = [y; y']
dY = [ Y(2);
      -Y(2)^2 - Y(1) - 0.5*cos(x)*Y(2) ];
end