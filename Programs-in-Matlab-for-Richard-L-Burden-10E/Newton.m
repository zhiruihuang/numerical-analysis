function p = Newton(f, p_0, TOL, N_0)
%% Algorithm 2.3 in page 67 of Numerical Analysis (10E)
% To find a solution to f(x)=0 given an initial approximation p_0: 
% INPUT:    function f; 
%           initial approximation p_0; 
%           tolerance TOL; 
%           maximum number of iterations N_0. 
% OUTPUT:   approximate solution p or a message of failure. 

% Example: 
% Newton(@(x) cos(x)-x, pi/4, 1e-10, 10);

% GMT+8 2019/9/24 17:11 By Rex HUANG
% Github: github.com/zhiruihuang

%% Compute the 1st derivative function of function f
syms x;
df = matlabFunction(diff(f(x)));

%% Step 1
i = 1;

%% Step 2
while true
    if i > N_0
        fprintf('The method failed after N0 iterations, N0 = %d. \n', N_0);
        break;
    end
    
    % Step 3
    if df(p_0) == 0
        fprintf("Failure, zero derivative. \n");
        break;
    end
    p = p_0 - f(p_0)/df(p_0);
    
    % Step 4. Stopping condition.
    % We also can use other stopping conditions, such as
        % abs(p-p0)/abs(p0) < TOL
        % abs(f(p)) < TOL
    if abs(p - p_0) < TOL
        fprintf(['Approximate solution p = %12.15f \n', ..., 
            ' with relative erorr %12.15f \n after %d iterations. \n'], ...,
            p, abs(p-p_0), i);
        break;
    end
    
    % Step 5
    i = i+1;
    
    % Step 6
    p_0 = p;
end

