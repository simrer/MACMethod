function f = Layeb15(x)
% multimodal, non-separable, scalable
% very hard to optimize. most methods fail
% optimal value 0 and it occurs at an alternation of 1 and -1
% search space [-100 100]
n = length(x);
indx = zeros(1,n-1);
for i = 1:n-1
    %indx(i) = 10*sqrt(tanh(2*abs(x(i))-x(i+1)^2-1))+ abs(exp(x(i)*x(i+1)+1)-1);
    indx(i) = (tanh(2*abs(x(i))-x(i+1)^2-1))^(1/10)+ abs(exp(x(i)*x(i+1)+1)-1);
end
f = sum(indx);