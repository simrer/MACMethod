function f = Layeb17(x)
% Wings function
% multimodal, non-separable, scalable
% optimal value 0 and it occurs at an alternation of -1 and 0
% search space [-100 100]
n = length(x);
indx = zeros(1,n-1);
for i = 1:n-1
    indx(i) = 10*abs(log((x(i)+x(i+1)+2)^2))-1/((1000*(x(i)^2-x(i+1)-1))^2+1)+1;
end
f = sum(indx);