function r = sep2(xt,yt) % Функция, вычисляющая расстояние между частицами.
% m-m distance
dx = xt(2) - xt(1);
dy = yt(2) - yt(1);
dxy = dx.*dx + dy.*dy;
r = sqrt(dxy);