function [ B ] = Bin_Parts( x, P )
%Bin_Parts Bins particles according to their x coordinate.
%   Takes in the boundaries of the binning mesh (x) and the particle list
%   P, where P is an n x 3 matrix, with each column containing the x, y and
%   z coordinate, and each row an individual particles location.

% Calculate some constants
[n,~] = size(P);  % Number of particles
m = length(x)-1;  % Number of bins
B = zeros(m,1);  % Pre-define the output vector

if ~isempty(P)
    
    % Error message if any particle lies outside the x range
    if min(P(:,1)) < x(1) || max(P(:,1)) > x(end)
        error('At least one of the particles is outside of the binning width')
    end
    
    % Calculate the bins
    for ii = 1:n
        
        Loc_Log1 = x(2:end) > P(ii,1);  % Find the points in the mesh which are bigger than the location (e.g [ | | | |*|*|*|*|*])
        Loc_Log2 = x(2:end) > P(ii,1) + x(2) - x(1);  % Find all points as above, but shifted by a width (e.g [ | | | | |*|*|*|*])
        B = B + (double(Loc_Log1-Loc_Log2))';  % Add on the difference between the two                   (e.g [ | | | |*| | | | ])
        
    end
end
    
end

