%convert vector index to matrix index, assuming the matrix dim to be nr*nc
function [x,y] = unvec(k, ncol)
    x = mod(k-1, ncol)+1;
    y  = (k-x)/ncol + 1;
end