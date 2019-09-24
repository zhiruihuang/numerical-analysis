function p = Bisection(f, a, b, TOL, N_0)
%% Algorithm 2.1 in page 49 of Numerical Analysis (10E)
% To find a solution to f(x)=0 given the continuous function f 
%   on the interval [a, b], where f(a) and f(b) have opposite signs: 
% INPUT:    function f; 
%           endpoints a, b; 
%           tolerance TOL; 
%           maximum number of iterations N_0. 
% OUTPUT:   approximate solution p or a message of failure. 

% Example: 
% Bisection(@(x) x^3+4*x^2-10, 1, 2, 1e-4, 15);

% GMT+8 2019/9/24 17:08 By Rex HUANG
% Github: github.com/zhiruihuang

%% Step 1
i = 1;
FA = f(a);

%% Step 2
while true
    if i > N_0
        fprintf('The method failed after N0 iterations, N0 = %d. \n', N_0);
        break;
    end
    
    % Step 3
    p = (a+b)/2; % p = a + (b-a)/2;
    FP = f(p);
    
    % Step 4. Stopping condition.
    % We also can use other stopping conditions, such as
        % abs(p-a)/abs(a) < TOL
        % abs(f(p)) < TOL
    if FP==0 || (b-a)/2<TOL 
        fprintf(['Approximate solution p = %12.15f \n', ..., 
            ' with relative erorr %12.15f \n after %d iterations. \n'], ...,
            p, (b-a)/2, i);
        break;
    end
    
    % Step 5
    i = i+1;
    
    % Step 6
    if FA*FP > 0
        a = p;
        FA = FP;
    else
        b = p;
    end
end
