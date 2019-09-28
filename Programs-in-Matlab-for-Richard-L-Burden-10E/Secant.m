function p = Secant(f, p_0, p_1, TOL, N_0)
%% Algorithm 2.4 in page 71 of Numerical Analysis (10E)
% To find a solution to f(x)=0 given an initial approximation p_0: 
% INPUT:    function f; 
%           initial approximation p_0, p_1; 
%           tolerance TOL; 
%           maximum number of iterations N_0. 
% OUTPUT:   approximate solution p or a message of failure. 

% Example: 
% Secant(@(x) cos(x)-x, 1/2, pi/4, 1e-7, 10);

% Matlab R2017b
% GMT+8 2019/9/28 18:24 By Rex HUANG
% Github: github.com/zhiruihuang

%% Step 1
i = 2;
q_0 = f(p_0);
q_1 = f(p_1);

%% Step 2
while true
    if i > N_0
        fprintf('The method failed after N0 iterations, N0 = %d. \n', N_0);
        break;
    end
    
    % Step 3
    p = p_1 - q_1*(p_1 - p_0)/(q_1 - q_0);
    
    % Step 4. Stopping condition.
    % We also can use other stopping conditions, such as
        % abs(p-p_1)/abs(p_1) < TOL
        % abs(f(p)) < TOL
    if abs(p - p_1) < TOL
        fprintf(['Approximate solution p = %12.15f \n', ..., 
            ' with relative erorr %12.15f \n after %d iterations. \n'], ...,
            p, abs(p-p_1), i);
        break;
    end
    
    % Step 5
    i = i+1;
    
    % Step 6
    p_0 = p_1;
    q_0 = q_1;
    p_1 = p;
    q_1 = f(p);
    
end

