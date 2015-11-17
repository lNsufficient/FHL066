function x = bar_def(ec,ed)
%x = bar_def(ec,ed)
% x den deformerade st�ngen representerad som vektor
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ed f�rskjutningar [a1 a2 ... a6]
a = ec(:) + ed(:);
x = a(4:6)-a(1:3);