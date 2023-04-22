function y = Layeb12(x)
% multimodal, non-separable, scalable
% very hard to optimize and most algorithm fails
% optimal value = -(e+1)(n-1)and it occurs at xi = 2
% Search space [-5 5]
n = length(x);
indx = zeros(1,n-1);
for j = 1:n-1
    indx(j) = cos(pi/2*x(j)-pi/4*x(j+1)-pi/2)*exp(cos(2*pi*x(j)*x(j+1)))+1;
end
y = -sum(indx);