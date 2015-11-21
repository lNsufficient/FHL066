function x = bar_def(ec,ed)
%x = bar_def(ec,ed)
% x den deformerade stången representerad som vektor
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ed förskjutningar [a1 a2 ... a6]
a = ec + [ed(1:3), ed(4:6)];
x = a(:,2)-a(:,1);
