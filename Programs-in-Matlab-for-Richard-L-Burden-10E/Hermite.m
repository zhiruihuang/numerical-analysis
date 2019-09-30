function p = Hermite(x, f_x, f_prime_x)
%% Algorithm 3.3 in page 139 of Numerical Analysis (10E)
% To obtain the coefficients of the Hermite interpolatory polynomial H(x)  
% on the (n+1) distinct numbers x_0, x_1, \ldots, x_n for the function f:

% INPUT:    column vector x, including numbers x_0, x_1, \ldots, x_n; 
%           column vector f_x, including values f(x_0), f(x_1), 
%               \ldots, f(x_n); 
%           column vector f_prime_x, including values f'(x_0), f'(x_1), 
%               \ldots, f'(x_n).
% OUTPUT:   the numbers Q_{0,0}, Q_{1,1}, \ldots, Q_{2n+1,2n+1} where 
%               H(x) = Q_{0,0} + Q_{1,1}(x-x_0) + Q_{2,2}(x-x_0)^2 
%                   + Q_{3,3}(x-x_0)^2 (x-x_1) 
%                   + Q_{4,4}(x-x_0)^2 (x-x_1)^2 + \ldots 
%                   + Q_{2n+1,2n+1}(x-x_0)^2 (x-x_1)^2 
%                       \cdots (x-x_{n-1})^2 (x-x_n) .

% Note that some subscript indices are different from the corresponding
% pseudocode, since the subscript index of matrices/arrays starts at 1 
% in Matlab. 

% Example: 
% x = [1.3; 1.6; 1.9];
% f_x = [0.6200860; 0.4554022; 0.2818186];
% f_prime_x = [-0.5220232; -0.5698959; -0.5811571];
% Hermite(x, f_x, f_prime_x);

% Matlab R2017b
% GMT+8 2019/9/30 22:47 By Rex HUANG
% Github: github.com/zhiruihuang


%% Step 1
n = length(x)-1;
T = zeros(2*n+2);
z = kron(x, [1; 1]);
T(:, 1) = kron(f_x, [1; 1]);
T(:, 2) = kron(f_prime_x, [0; 1]);
% i>1
for i=2:(n+1)
    T(2*i-1, 2) = (T(2*i-1, 1)-T(2*i-2, 1)) / (z(2*i-1)-z(2*i-2));
end

% z = zeros(1, 2*n+2)';
% for i=1:(n+1)
%     z(2*i-1) = x(i);
%     z(2*i) = x(i);
%     T(2*i-1, 1) = f_x(i);
%     T(2*i, 1) = f_x(i);
%     T(2*i, 2) = f_prime_x(i);
%     if i>1
%         T(2*i-1, 2) = (T(2*i-1, 1)-T(2*i-2, 1)) / (z(2*i-1)-z(2*i-2));
%     end
% end

%% Step 4
for i=3:(2*n+2)
    for j=3:i
        T(i, j) = (T(i,j-1)-T(i-1,j-1)) / (z(i)-z(i-j+1));
    end
end
T = [(0:(2*n+1))', z, T];
Q = T(2:end, 4:end);
%% Step 5
p = [T(1, 3); diag(Q)];