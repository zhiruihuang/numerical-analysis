function p = NewtonDD(x, f_x)
%% Algorithm 3.2 in page 126 of Numerical Analysis (10E)
% To obtain the divided-difference coefficients 
% of the interpolatory polynomial P 
% on the (n+1) distinct numbers x_0, x_1, \ldots, x_n for the function f:

% INPUT:    column vector x, including numbers x_0, x_1, \ldots, x_n; 
%           column vector f_x, including values f(x_0), f(x_1), \ldots, 
%               f(x_n) as F_{0,0}, F_{1,0}, \ldots, F_{n,0}. 
% OUTPUT:   the numbers F_{0,0}, F_{1,1}, \ldots, F_{n,n} where 
%               P_n(x) = F_{0,0} 
%               + \sum_{i=1}^n F_{i,i} \prod_{j=0}^{i-1} (x-x_j) .
%               (F_{i,i} is f[x_0, x_1, \ldots, x_i]. )

% Note that some subscript indices are different from the corresponding
% pseudocode, since the subscript index of matrices/arrays starts at 1 
% in Matlab. 

% Example: 
% x = [1; 1.3; 1.6; 1.9; 2.2];
% f_x = [0.7651997; 0.6200860; 0.4554022; 0.2818186; 0.1103623];
% NewtonDD(x, f_x);

% Matlab R2017b
% GMT+8 2019/9/30 22:46 By Rex HUANG
% Github: github.com/zhiruihuang

%% Step 1
n = length(x)-1;
T = zeros(n+1);
T(:, 1) = f_x;
for i=2:(n+1)
    for j=2:i
        T(i, j) = (T(i, j-1)-T(i-1, j-1)) / (x(i)-x(i-j+1));
    end
end
T = [(0:n)', x, T];
F = T(2:end, 4:end);
%% Step 2
p = [T(1, 3); diag(F)];
