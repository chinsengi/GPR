

% rng(1)
% randi([1,50], 3,3)

tmp = 10;
x4 = [1,2];
[y, dy] = dlfeval(@fun, dlarray([1 2 3 4 5 6]), x4);
dy

function [y,dy1] = fun(x1, x4)
x4 = 2;
x2 = 3;
x3 = anotherfun(x1);
tmp = x3([1 3 5])
tmp = tmp.*x3;
y = sum(tmp,'all') + x2 + mean(x4);
% tmp  = x3;
% y = y+tmp;

% dy = dlgradient(y,x2); % Error: x2 is untraced
dy1 = dlgradient(y,x1); % No error even though y has an untraced portion
end

function [ret] = anotherfun(x)
    ret = reshape(x, 2, 3);
end
function [x,y] = unvec(k, ncol)
    y = mod(k-1, ncol)+1;
    x  = (k-y)/ncol + 1;
end