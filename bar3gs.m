function [es,ee] = bar3gs(ec,ep,ed)
%[es,ee] = bar3gs(ec,ep,ed)
% es normalkraft
% ee Greentöjning
% ---
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
% ep materialparametrar [E A0]
% ed förskjutningar [a1 a2 ... a6]
E=ep(1); A0 = ep(2);
L0 = norm(bar(ec));

%X = ec(:) + ed(:);
L = norm(bar_def(ec,ed));

ee = (L^2-L0^2)/(2*L0^2);

es=E*A0*ee;