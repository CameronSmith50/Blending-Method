function [ A ] = A_mat( P_lambda, K_mat, r )
%A_mat Calculates the A matrix for the p function
%   Takes in P_lambda, the kernel matrix K_mat and the mesh r, and outputs 
%   the matrix A.

% Find N and M, the number of meshpoints in [0,1] and [0,a]
h1 = r(2)-r(1);
h2 = r(end) - r(end-1);
N = 1/h1;
M = length(r)-1;

% Matrix of pre-multipliers
v = zeros(M+1,1);
nums = 1:(M+1);
v(nums<N+1 & mod(nums,2) == 0) = 4*h1/3*(1-P_lambda);
v(nums<N+1 & mod(nums,2) == 1) = 2*h1/3*(1-P_lambda);
v(nums>N+1 & nums<M+1 & mod(nums,2) == 0) = 4*h2/3;
v(nums>N+1 & nums<M+1 & mod(nums,2) == 1) = 2*h2/3;
v(1) = h1/3*(1-P_lambda);
v(N+1) = 1/3*(h1-h1*P_lambda+h2);
v(M+1) = h2/3;
X = repmat(v,1,M+1)';
A = K_mat.*X;

end

