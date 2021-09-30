function m = checkpts(m)
    [rows, cols] = size(m);
    if rows > cols
        m = m';
        [cols, rows] = deal(rows, cols); % swap
    end
    
    if rows == 2
        m = [m;
             ones(1,cols)];
    elseif rows == 3
        m = m./m(3,:);
    else
        error('input points must be 2xN or Nx2 or 3xN or Nx3');
    end
    
end