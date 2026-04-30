function result = regulaRootHelper(gi, gi_1)
    %  a helper function for computing regulaRoot 
    %  gi and gi_1 are two equal_size 2d matrix
    difference = gi - gi_1;
    abs_dif = abs(difference);
    result = difference./abs_dif;
    result(~isfinite(result))=0;
end