function f = RosenbrockFcn(x)
d = length(x);
sum = 0;
for i = 1:(d-1)
    xi = x(i);
    xnext = x(i+1);
    new = 100*(xnext - xi^2)^2 + (xi - 1)^2;
    sum = sum + new;
end
    f = sum;
end
