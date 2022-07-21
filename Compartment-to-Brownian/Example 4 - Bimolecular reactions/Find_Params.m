function [ P_lambda ] = Find_Params( gam, kappa, K_mat, r )
%Find_Params Finds the two unknown parameters Plambda and sigma.
%   Uses the model parameters that are known in order to find the values of
%   the probability of reacting in a time-step for the bimolecular
%   reaction (Plambda) and the dissociation radius (dissoc).

% STEP 1: Find h1, N and M
h1 = r(2)-r(1);
N = 1/h1;
M = length(r) - 1;

% STEP 2: Generate the matrix A for the g function as a function of Plambda
A = @(Plambda) A_mat(Plambda,K_mat,r);

% STEP 5: Generate the RHS vector for the g function and turn it into a
% column vector
Jg = 0.5*(2-(erf((r(end)-r)./(gam*sqrt(2))) + erf((r(end)+r)./(gam*sqrt(2))))+2*gam^2/r(end).*K_fun(r,r(end),gam))';

% STEP 6: Solve the numerical systems for both p and g
g = @(Plambda) (eye(M+1)-A(Plambda))\Jg;

% STEP 8: Calculate the vector b such that g.b = kappa
v = zeros(M+1,1);
nums = 1:(M+1);
v(mod(nums,2) == 1) = 2;  % Even and odd switched due to matlab indexing starting at 1 and not 0
v(mod(nums,2) == 0) = 4;
v([1,N+1,M+1]) = 1;
D = diag(v);
b = @(Plambda) 2*pi*h1*Plambda/3*D*[r(r<=1).^2,r(r>1)*0]';

% STEP 9: Solve g.b = kappa using fsolve
% options = optimoptions('fsolve','Display','iter');
% P_lambda = fsolve(@(Plambda) g(Plambda,dissoc_fun(Plambda))'*b(Plambda)-kappa,1,options);
options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',1e-16,'OptimalityTolerance',1e-16);
P_lambda = lsqnonlin(@(Plambda) g(Plambda)'*b(Plambda)-kappa,0.5,0,1,options);

end

