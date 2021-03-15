function W = ColumnsPartition(X, partitions)
    cols = size(X,2);
    stride = floor( cols/partitions);
    remaining = mod( cols, partitions);
    
    count = 1;
    W = cell( 1, partitions );
    
    for i = 1:partitions
        step = count + stride -1;
        if ( i == partitions )
            step = step + remaining;
        end
        W{i} = X(:,count:step);
        count = count + stride;
    end
end