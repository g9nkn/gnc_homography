function err = homoDist_reproj(x, m1, m2)
    n = size(m1,2);

    k = 9;
    x1 = [x( k:k+n-1)'
          x(k+n:k+2*n-1)'];
    x2 = [x(k+2*n:k+3*n-1)'
          x(k+3*n:k+4*n-1)'];
    
%     x1 = [x(  9:9+n-1)'
%           x(9+n:9+2*n-1)'];
%     x2 = [x(9+2*n:9+3*n-1)'
%           x(9+3*n:9+4*n-1)'];
    
    err1 = ( x1 - m1(1:2,:) )';
    err2 = ( x2 - m2(1:2,:) )';
    err = [err1, err2];
end
