function p = FixedPoint(g, p_0, TOL, N_0)
%% Algorithm 2.2 in page 59 of Numerical Analysis (10E)
% To find a solution to p=g(x) given an initial approximation p_0:
% INPUT:    function f; 
%           initial approximation p_0;
%           tolerance TOL; 
%           maximum number of iterations N_0. 
% OUTPUT:   approximate solution p or a message of failure. 

% Example: 
% FixedPoint(@(x) cos(x), pi/4, 1e-4, 20);

% Matlab R2017b
% GMT+8 2019/9/28 18:24 By Rex HUANG
% Github: github.com/zhiruihuang

%% Step 1
i = 1;

%% Step 2
while true
    if i > N_0
        fprintf('The method failed after N0 iterations, N0 = %d. \n', N_0);
        break;
    end
    
    % Step 3
	p = g(p_0);

    % Step 4. Stopping condition.
    % We also can use other stopping conditions, such as
        % abs(p-p_0)/abs(p_0) < TOL
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
