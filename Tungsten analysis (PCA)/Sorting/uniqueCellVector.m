function [C,ia,ic] = uniqueCellVector(A)
    % Written by Benjamin Kraus
    % MathWorks Technical Support Department
    % June 17, 2014
    
    n = numel(A);
    
    ic = (1:n)';
    for cur = 2:n
        % Compare current item with all previous items.
        for prev = 1:(cur-1)
            % If we find a match, then we are done, update the list of
            % unique items to show that current item = previous item.
            if(isequal(A(cur),A(prev)))
                ic(cur) = prev;
                break
            end
        end
    end
    
    % Now that we have numerical vectors, use 'unique' to make the unique
    % integers continuous (instead of skipping over numbers). This will
    % also generate 'ia' and regenerate 'ic' automatically.
    [~,ia,ic] = unique(ic);
    
    % Use 'ia' to resample the original matrix and return the unique
    % elements of the original matrix.
    C = A(ia);
end
