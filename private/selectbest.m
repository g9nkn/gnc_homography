function [H, inliers, beststats, best_idx] = selectbest(stats)
    
    thresh = [stats.thresh];
    ratio  = [stats.ratio];
    
    diff_r = [NaN, diff(ratio)];
    diff_r( ~(thresh > min(thresh) & thresh < 30) ) = NaN;
    
    if all(isnan(diff_r))
        beststats = stats(end);
        best_idx  = length(stats);
    else
        [~, minidx] = min( abs(diff_r) );
        best_idx  = minidx;
        beststats = stats(best_idx);
        
    end
    

    H = beststats.model;
    inliers = beststats.weight > 0.5;
    
    
end
    