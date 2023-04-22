function [f,g] = RastriginFcn(x)
n = length(x);
f = 0;
g = zeros(n,1);
for i = 1:n
    f = f + x(i)^2 - 10*cos(2*pi*x(i));
end
f = 10*n + f;
if nargout > 1
   for i = 1:n
    g(i,1) = g(i,1) + 2*x(i) + 20*pi*sin(2*pi*x(i));
   end
end