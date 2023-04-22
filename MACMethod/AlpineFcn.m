function f = AlpineFcn(xx)
d = length(xx);
sum = 0;
for ii = 1:d
    xi = xx(ii);
    sum = sum + abs(xi*sin(xi)+0.1*xi);
end
f = sum;
end