% Reference http://www.leet.it/home/giusti/teaching/matlab_sessions/stitching/stitch.html

function result_image = stitch(image, H, target_image, alpha, invert)
    if nargin < 4, alpha = 0.5; end
    if nargin < 5, invert = false; end
    
    if invert
        H = inv(H);
        [target_image, image] = deal(image, target_image);
    end
    H = H/H(3,3);

    if isinteger(image),        image        = double(image); end
    if isinteger(target_image), target_image = double(target_image); end
    
    [h, w, ~] = size(image);
    [h2, w2, ~] = size(target_image);
    max_size = max([2*w, 2*h]);
    toolarge = 1;
    itr = 0;
    while toolarge
        
    m = [1 w 1 w
         1 1 h h
         1 1 1 1];
        x  = H(1,:)*m./(H(3,:)*m);
        y  = H(2,:)*m./(H(3,:)*m);
        xdata = [min(x), max(x)];
        ydata = [min(y), max(y)];

        xdata_out=[min(1,xdata(1)) max(w2, xdata(2))];
        ydata_out=[min(1,ydata(1)) max(h2, ydata(2))];

        xsize = round(sum(abs(xdata_out)));
        ysize = round(sum(abs(ydata_out)));
        if xsize <= max_size && ysize <= max_size
            toolarge = 0;
        else
            H = diag([0.5, 0.5, 1]) * H * inv(diag([0.5, 0.5, 1]));
            w = w/2;
            h = h/2;
            w2 = w2/2;
            h2 = h2/2;
            itr = itr + 1;
        end
    end
    image = imresize(image, 0.5^itr);
    target_image = imresize(target_image, 0.5^itr);
    Rout  = imref2d([ysize, xsize], xdata_out, ydata_out);
    
    result1 = imwarp(       image, projective2d(H'), 'OutputView', Rout);
    result2 = imwarp(target_image, affine2d(eye(3)), 'OutputView', Rout);
    
    result_image = result1 + result2;
    overlap      = (result1 > 0) & (result2 > 0);
    result_avg   = alpha*result1 + (1-alpha)*result2;
    
    result_image(overlap) = result_avg(overlap);
    result_image = uint8(result_image);
    
end