function W = RowsPartition( X, partitions)
% Segments the X matrix into blocks
%
% INPUT :
% X = Tall data matrix
% partitions= Number of segments
%
% OUTPUT :
% W = Array of cells containing blocks of partitioned data
%
    rows = size(X,1);
    stride = floor( rows/partitions);
    remaining = mod( rows, partitions);
    
    count = 1;
    W = cell( 1, partitions );
    
    for i = 1:partitions
        step = count + stride -1;
        if ( i == partitions )
            step = step + remaining;
        end
        W{i} = X(count:step,:);
        count = count + stride;
    end
end