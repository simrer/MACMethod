function cost = ZakharovFcn(x)
d = length(x);
sum1 = 0;
sum2 = 0;
for i = 1:d;
    sum1 = sum1 + x(i)^2;
    sum2 = sum2 + 0.5*i*x(i);
end
cost = sum1 + sum2^2 + sum2^4;
end