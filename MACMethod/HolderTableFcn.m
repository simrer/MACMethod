function f = HolderTableFcn(x)
f = -abs(sin(x(1))*cos(x(2))*exp(abs(1-(sqrt(x(1)^2+x(2)^2)/pi))));
end