function y = Layeb03(x)
% multimodal, non-separable, scalable
% very hard to optimizers in higher dimentions
% optimal value is -n+1(not reproduceble it seems wrong) and it occurs at kpi
%search space (-10,10)
n = length(x);
indx = zeros(1,n-1);
for i = 1:n-1
    indx(i) = abs(sin(x(i))*(exp(abs(100-((sqrt(x(i)^2+x(i+1)^2))/pi))))...
        + sin(x(i+1))+1)^-(0.1);
end
y = -sum(indx);