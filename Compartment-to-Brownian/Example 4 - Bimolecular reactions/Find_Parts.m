function [ Ords ] = Find_Parts( L )
%Find_Parts Turns a logical matrix into a t column matrix containing the
%non-zero coordinates
%   The function takes in a logical matrix L and outputs the co-ordinates
%   of the non-zero entries (i.e which row and column).
%
% Cameron Smith
% Date created: 30/10/17
% Last modified: 30/10/17

% Find the size of the matrix L
[r,~] = size(L);

% Find the total number of elements and pre-define the output matrix
num_els = sum(sum(L));
Ords = zeros(num_els,2);

% Create a dummy index for the recording
rec = 1;

% Loop through the rows
for ii = 1:r
        
    % Check to see if there are any ones in the row
    rows = sum(L(ii,:));
    
    % If this is non-zero, then go through
    if rows > 0
        
        % Create dummy variables for looping through the columns, once all
        % have been found, breaks from loop
        sums = 0;
        col = 1;
        while sums < rows
            
            % If the column is one then place into the output matrix,
            % otherwise don't. Update all dummy variables appropriately
            if L(ii,col) == 1
                Ords(rec,:) = [ii,col];
                col = col+1;
                sums = sums+1;
                rec = rec+1;
            else
                col = col+1;
            end
        end
    end
end
            
            


end

