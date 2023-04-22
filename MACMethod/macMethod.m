function [xmin,fmin,y,out] = macMethod(tfun,nvars,LB,UB,U_0,u_0,alpha_0,N,delta,gamma_0)

    % New stochastic optimization algorithm proposed by Dr. Attila Nagy
    % test function is the function need to optimize (default test function
    % s for this test are given in this script function

    % INPUT ARGUMENTS
    % "nvars" are the number of decision variables
    % "LB" and "UB" are the lower and upper bounds of the decision variables
    % respectively.
    % "U_0" positive definite matrix of size nvars (in this case I take the 
    % identity matrix.
    % "u_0" any initial vector of size nvars within the domain of the
    % objective function.
    % "alpha_0" initial positive integer whose value updated at every step
    % (it is taken 1 in this case).
    % "N" any positive integer whose value increases at every step (it is 2
    % in this case)
    % "delta" is the tolerance in the objective function value.
    % "gamma_0" initial positive real number whose value increases at every
    % step (it is 0.001 in this case).

    % OUTPUT ARGUMENTS
    % "xmin" is the optimizer point
    % "fmin" is the optimum value of the objective function.
    % "y" is the fmin value at every iteration
    % "out" is a field that contains information like the number of function
    % evaluation, time,...
    
    tStart = tic;
    if ~ nargin
       tfun = @demo1;   % Demo
       nvars = 3;
       u_0 = 5*ones(nvars,1);
       U_0 = eye(nvars);
       LB = -10*ones(nvars,1);
       UB = -LB;  
    end

% Default parameters

    if nargin > 10    
       error('function requires at most 10 input arguments')
    end

    if nargin < 10
       gamma_0 = 0.001;
    end

    if nargin < 9
       delta = 1e-6;
    end

    if nargin < 8
       N = 4;
    end

    if nargin < 7
       alpha_0 = 1;
    end

    if nargin && nargin < 6 
       u_0 = [];
    end

    if nargin && nargin < 5
        U_0 = [];
    end
    
    if nargin && nargin < 4
        UB = [];
    end

    if nargin && nargin < 3
        LB = [];
    end

    if nargin && nargin < 2
        disp('MAC method requires at least two input arguments')
    end
    
    if ischar(tfun)
       tfun = str2func(tfun);
    end

    if isempty(LB)
        LB = -100*ones(nvars,1);
    end

    if isempty(UB)
        UB = 100*ones(nvars,1);
    end

    if isempty(U_0)
        U_0 = eye(nvars);
    end

    if isempty(u_0)
        u_0 = 5*ones(nvars,1);
    end

    % Generating random numbers from basic distribution with mean 0 and identity standard deviation (lognormal
    % distribution in this case)
        rng(123)
        mu = zeros(nvars,1);
        sigma = eye(nvars);
        zeta = mvnrnd(mu,sigma,alpha_0)'; % multivariate normal with mean mu and standard deviation sigma.
        zeta = (0.5-1./(1+exp(-zeta)))/min(std(zeta)); % multivariate logit-normal (by transforming the multivariate normal to logistic function)
        numeratorCovSum = zeros(nvars,nvars);   
        WeightedP = zeros(nvars,alpha_0); 
        
        % Random points pi for i = 1 to alpha_n

        p = u_0 + U_0*zeta;

    % Calculating the weights for every points pi, emperical expected value
    % and emperical covariance matrix.

        if all (p >= LB) && all(p <= UB) % if pi is in K (domain of the objective function)  
           weight1 = weightFcnOfGamma(gamma_0, tfun(transpose(p))); % weight/penality
        else
            weight1 = 0;
        end

        if weight1 ~= 0 && (weight1 > 1 || weight1 < 1e-9)
            weight1 = rand;
        end
     
       sum1 = sum(weight1);

      if sum1 == 0
          u_n = u_0/2;
          U_n = U_0;
      else
        u_n = weight1*p/2*sum1; % Emperical expected value
        weightedCovMatrix = (weight1*(p - u_n)*transpose(p - u_n))/2*sum1;
        numeratorCovSum = numeratorCovSum + weightedCovMatrix;  
        U_n = sqrt(numeratorCovSum); % Emperical covariance matrix
     end

    % Calculating objective function value at several random points pi and
    % storing the minimum value

        p_all = [p u_n u_0];
        [xmin,fmin] = getMin(p_all,tfun);
        store_fmin = fmin;
%         u_n = xmin;
        store_xmin = xmin;
        xout = store_xmin; 
        fvec = store_fmin;
        iter = 0;
        n = 0;
        tolX = 1e-4;
        countf = size(p_all,2); % Initial fevals
        MaxIter = 100;
        maxEvals = 10000;

fprintf('Result of iteration number: %d is %5.4f\n',n,store_fmin);

% Main Loop

while n < MaxIter && (norm(u_n - u_0) >= delta)% || norm(U_n-U_0) >= delta)

    if countf > maxEvals
        disp('The method stopped because it reaches maximum number of function evaluation');
        break 
     end

    if n == MaxIter
       disp('The method stoped because it reaches the maximum number of iteration')
       break
    end

% Updating initial values
        gamma_0 = 4*gamma_0; 
        u_0 = u_n;
        U_0= U_n;
        n = n + 1;
        N = N + 1;
        s_p = size(p,2) + 1;
        alpha_nN = n*N+s_p-1;
%         alpha_nN = n*N+size(p,2);
        rng(123)
        % Generating random numbers zeta_i from basic distribution for i = alpha_n-1 to alpha_nN 
        len = alpha_nN-s_p+1;
        zeta(:,s_p:alpha_nN) = mvnrnd(mu,sigma,len)';
        zeta(:,s_p:alpha_nN) = (0.5-1./(1+exp(-zeta(:,s_p:alpha_nN))))/min(std(zeta(s_p:alpha_nN)));
       
    for ii = s_p:alpha_nN
        p(:,ii) = u_0 + U_0*zeta(:,ii);%
    end

    for kk = s_p:size(p,2)
        if all (p(:,kk) >= LB) && all(p(:,kk) <= UB) % if pi is in K  
           weight1(kk) = weightFcnOfGamma(gamma_0, tfun(transpose(p(:,kk))));
           WeightedP(:,kk) = weight1(kk)*p(:,kk);%
        else
            weight1(kk) = 0;
        end

        if weight1(kk) ~= 0 && (weight1(kk) > 1 || weight1(kk) < 1e-9)
            weight1(kk) = rand;
            WeightedP(:,kk) = weight1(kk)*p(:,kk);%
        end     
    end

       sum1 = sum1 + sum(weight1(s_p:end)); % Sum of weights
       p_cut = p(:,s_p:end);
       u_n = sum(WeightedP,2)/sum1;
       p_all = [p_cut u_n];
       [xmin,fmin] = getMin(p_all,tfun);
       store_fmin = min([fmin,store_fmin]);
       
       if store_fmin == fmin && all(xmin == u_n)
          store_xmin = u_n;
     
       else 
          store_xmin = xmin;
          u_n = u_0/2;
       end

    if sum1 ~= 0
       for ii = s_p: size(p,2) %alpha_nN
          weightedCovMatrix = (weight1(ii)*(p(:,ii) - u_n)*transpose(p(:,ii) - u_n))/sum1;
          numeratorCovSum = numeratorCovSum + weightedCovMatrix;
       end 
    end
        U_n = sqrt(numeratorCovSum); 

        xout = [xout store_xmin];
        fprintf('Result of iteration number: %d is %5.4f\n',n,store_fmin);
        s = size(p(:,s_p:end),2);
        countf = countf + s;
        iter = [iter n];
        fvec = [fvec store_fmin];

        if norm(xout(:,end-1)-xout(:,end)) < tolX
        disp('The method stopped because it reaches optimality tolerance');
        break 
       end
end

    [fmin,idx] = min(fvec);
    xmin = xout(:,idx);
    out.iter = iter(:);
    out.fvec = fvec(:);
    out.fcount = countf;
    out.xvec = xout;
    fvec = sort(fvec,'descend');
    y = [iter;fvec]';
    fid = fopen('logitFromLogisticMVNDResult.txt','w');
    fprintf(fid,'it     f\n');
    fprintf(fid,'%d  %8.4f\n',y');
    fclose(fid);
    plot(y(:,1),y(:,2),'bo-','LineWidth',2)
    xlabel('Iteration number')
    ylabel('function value')
    title('Minimization of an objective function by new method')
    tEnd = toc(tStart);
    out.tEnd = tEnd;
    out.store_xmin = store_xmin;
    out.u_0 = u_0;

 
function f = demo1(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
f  = (x1-4)^2 + 5*(x2-4)^4 +3*(x3-4)^6;

function g = weightFcnOfGamma(gamma, ObjFcn)
         g = exp(-gamma*ObjFcn);

function [x_min,f_min] = getMin(p,ObjFcn)
        n = size(p,2);
        f_min = inf;
        x_min = [];
        for ii = 1:n
          a = ObjFcn(p(:,ii)');
          if a < f_min
             x_min = p(:,ii);
            f_min = a;
         end
       end
    

