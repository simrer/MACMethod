function y = Layeb11(x)
% multimodal, non-separable, scalable
% very hard to optimize 
% optimal value is -n+1 and it occurs at an alternation of -1 and 0
% Search space [-10 10]
n = length(x);
indx = zeros(1,n-1);
for j = 1:n-1
    indx(j) = (cos(x(j)*x(j+1)+pi))/((100*abs(x(j)^2-x(j+1)-1))^2+1);
end
y = sum(indx);