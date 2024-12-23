function f=forc2(r)
% Функция, вычисляющая силу взаимодействия между частицами по потенциалу Леннарда-Джонса.
%  m-m force & Len-Jones potential
%  used in mo22de.m // tested in tstforc2.m
ri = 1./r;
ri3 = ri.*ri.*ri;
ri6 = ri3.*ri3;
g = ri.*ri6.*(2*ri6-1);
f(1,:) = g./r;
% 1/r canceled in Fx Fy by dx - dy
% f(2) - L-J potential
f(2,:) = 4*ri6.*(ri6-1);