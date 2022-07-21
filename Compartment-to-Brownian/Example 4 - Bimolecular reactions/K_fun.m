function [ K ] = K_fun( x,y,gam )
%K_fun Evaluates the kernel from Lipkova et al. (2011)

K = y./(x.*gam*sqrt(2*pi)).*(exp(-(x-y).^2./(2*gam^2)) - exp(-(x+y).^2./(2*gam^2)));

% Set NaN's to 0
K(isnan(K)) = 0;

% If x = 0, we use de L'Hopital's rule in order to find that K(0,y) =
% sqrt(2/pi)*y^2*exp(-y^2/(2*gam^2))/gam^3
if abs(x) < 1e-10
    K = sqrt(2/pi)*y^2/gam^3*exp(-(y^2)/(2*gam^2));
end


end

