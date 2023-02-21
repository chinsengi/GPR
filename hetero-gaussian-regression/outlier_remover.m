function ret = outlier_remover(data)
    ret = data;
    dim = size(data);
    for i = 1: numel(data)
        [x,y] = ind2sub(size(data), i);
        
        surronding = [];
        for j = x-1:x+1
            for k = y-1:y+1
                if j>0 && k>0 &&j <= dim(1) && k <= dim(2) && ~(j== x&& k==y)
                    surronding(end+1) = data(j,k);
                end
            end
        end
            
        medi = median(surronding);
        variance = var(surronding);
        if abs(ret(x,y)- medi)>variance*2
            ret(x,y) = medi;
        end
    end
end