function x0 = bar(ec)
% Skriver stången som en vektor givet
% ec odeformerade koordinater [x1 x2; y1 y2 ; z1 z2]
x0 = ec(:,2)-ec(:,1);