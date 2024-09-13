syms x y X Y;

u = 2 + sin((2*pi*x)/X)*sin((4*pi*y)/Y);
h1 = diff(u, x);
h2 = diff(u, y);
