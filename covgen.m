function [cov, cov2] = covgen(sigma, l, cov, cov2) %#codegen
    for i = 0:2499
        for j = 0:2499
            x1 = mod(i,50);
            y1 = floor((i-x1)/50);
            x2 = mod(j,50);
            y2 = floor((j-x2)/50);
            dist = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
            cov(i+1,j+1) = sigma*sigma*exp(-l*dist/2);
            cov2(i+1,j+1) = -sigma*sigma*dist*exp(-l*dist/2)/2;
        end
    end
end