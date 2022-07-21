function [ Q ] = Remove_Duplicates( P, N )
%Remove_Duplicates Removes duplicated particles that have reacted at random
%   Takes in a matrix P, a Rx2 matrix, containing pairs of indices deemed
%   to react, and outputs Q, the same matrix, but with rows removed so that
%   each index appears at most once in Q. If rows are to be removed, they
%   are removed at random. The value of N is the total number of particles
%   in the system.
%
% Cameron Smith
% Date created: 27/12/17
% Last modified: 27/12/17

% Extract the value of R
[R,~] = size(P);

% Firstly, check if there are any duplicates
uniqueP = unique(P(:));
if length(P(:)) == length(uniqueP)
    Q = P;

% Otherwise, we remove those repeated at random
else

    % Create a cell array to contain index locations and a vector for the
    % number of times the vector appears
    C = cell(1,N);
    v = zeros(1,N);
    
    % Loop through the matrix P to store indices and number of times each
    % particle appears
    for ii = 1:R
        
        C{P(ii,1)} = [C{P(ii,1)},ii];
        C{P(ii,2)} = [C{P(ii,2)},ii];
        v(P(ii,:)) = v(P(ii,:)) + 1;
        
    end
    
    % Run through a loop to remove the duplicates
    To_Remove = [];  % Vector containing row indices to be deleted
    for jj = 1:N
        
        % Firstly check if there is a duplicated index here
        if v(jj) > 1
            
            % Extract the vector of indices
            x = C{jj};
            x_len = length(x);
            
            % Choose an index to keep
            keep = randi(x_len);
            
            % Update both v and C to reflect the reaction that has
            % occurred. First find all of the reaction indices that can no
            % longer occur
            Remove = [];  % These are the row indices
            y = P(x(keep),:);  % These are the two particles that have reacted
            z1 = C{y(1)};  % Look up all reaction pathways connected to the first particle
            Remove = [Remove, z1(z1 ~= x(keep))]; %#ok
            z2 = C{y(2)};  % Look up all reaction pathways connected to the second particle
            Remove = [Remove, z2(z2 ~= x(keep))]; %#ok
            Remove = sort(Remove);  % Sort the vector
            
            % Now loop through these indices and remove them from the count
            % vector v and the cell array C
            for kk = 1:length(Remove)
                w = P(Remove(kk),:);  % Extract the particles on this pathway
                v(w) = v(w) - 1;  % Remove one from the count
                vec1 = C{w(1)};  % Extract reaction indices related to the first particle
                vec1(vec1 == Remove(kk)) = [];  % Remove the pathway
                C{w(1)} = vec1;
                vec2 = C{w(2)};  % Extract reaction indices related to the second particle
                vec2(vec2 == Remove(kk)) = [];  % Remove the pathway
                C{w(2)} = vec2;
            end
            
            % Add the indices to the remove list
            To_Remove = [To_Remove,Remove];  %#ok
        end
        
    end
    
    % Remove the duplicates
    To_Remove = unique(To_Remove);
    
    % Generate a vector which extracts the correct indices
    keeping = ones(1,R);
    keeping(To_Remove) = 0;
    keeping = (1:R).*keeping;
    keeping(keeping == 0) = [];
    
    % Finally, remove those that have been removed
    Q = P(keeping,:);
    
end
            

end

