function y = Layeb02(x)
% uniimodal, separable, scalable
% has high value aroung the point [0,0,0,...]
% optimal value is 0 and it occurs at x = 1;
% search space (-10,10)
n = length(x);
indx = zeros(1,n);
for i = 1:n
    indx(i) = abs(exp((100*(x(i)-1)^2)/(exp(x(i)+1)))-1);
end
y = sum(indx);