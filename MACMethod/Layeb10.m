function y = Layeb10(x)
% multimodal, non-separable, scalable
% very hard to optimize 
% optimal value is 0 and it occurs at xi = 0.5
% Search space [-100 100]
n = length(x);
indx = zeros(1,n-1);
for j = 1:n-1
    indx(j) = (log(x(j)^2+x(j+1)^2+0.5))^2 + abs(100*sin(x(j)-x(j+1)));
end
y = sum(indx);