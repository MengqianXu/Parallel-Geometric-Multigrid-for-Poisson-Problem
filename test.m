syms x y a b c d;

u = (x^a)*((1 - x)^b)*(y^c)*((1 - y)^d);
h1 = diff(u, x, 2);
h2 = diff(u, y, 2);
f = -(h1 + h2);