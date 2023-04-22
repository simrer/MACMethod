function y = Layeb01(x)
% unimodal, separable and scalable
% optimal value is 0 and it occurs at x = 1
% domain [-100, 100]
n = length(x);
exp1 = zeros(1,n);
for i = 1:n
    exp1(i) = sqrt(abs(exp(x(i)-1)^2-1));
end
y = sum(100*exp1);

    