function [ Dist ] = Pairwise_Dist( P )
%Pairwise_Dist Efficiently calculates the pairwise distances between
%particles
%   Particle positions are given in the variable P (an N*3 matrix where N
%   is the number of particles and 3 is the dimension of the space), and
%   outputs the pairwise distances in a matrix Dist. The code also makes
%   the matrix upper triangular.

% Find the number of particles and dimension
[N,~] = size(P);

% Find the positions of each particle in their coordinate directions
x_pos = ones(N,1)*P(:,1)';
y_pos = ones(N,1)*P(:,2)';
z_pos = ones(N,1)*P(:,3)';

% Calculate the pairwise distances in each direction
x_dist = abs(x_pos - x_pos.');
y_dist = abs(y_pos - y_pos.');
z_dist = abs(z_pos - z_pos.');

% Find the distance matrix
Dist = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2);

% Remove the lower triangular part
Dist = Dist - tril(Dist,-1);

end

