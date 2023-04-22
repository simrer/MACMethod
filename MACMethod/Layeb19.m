function f = Layeb19(x)
% Noiselog function
% unimodal, separable, scalable
% very hard to optimize and it contains a random component
% optimal value 0 and it occurs at xi = 1
% search space [-5 5]
n = length(x);
indx = zeros(1,n);
rng('default')
for i = 1:n
    indx(i) = 100*rand^i*(log((x(i)-1)^2+1))^2;
end
f = sum(indx);