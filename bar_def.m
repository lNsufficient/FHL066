function x = bar_def(ec,ed)

a = ec(:) + ed(:);
x = a(4:6)-a(1:3);