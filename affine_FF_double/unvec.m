%convert vector index to matrix index, assuming the matrix dim to be nr*nc
function [x,y] = unvec(k, ncol)
    y = mod(k-1, ncol)+1;
    x  = (k-y)/ncol + 1;
end