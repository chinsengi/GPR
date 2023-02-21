function Abar = svd_complete(A, orig)
    Abar = zeros(size(A));
    while norm(A-Abar)>1e-6
        [U,S,V] = svd(A,'econ');
        s = diag(S);
%         tau = getThreshold(s);
%         tau = 5;
%         s(tau+1:end) = 0;
        s(s<1) = 0;
        S = diag(s);
        Abar = A;
        A = U*S*V';
        A(orig~=0) = orig(orig~=0);
%         fprintf('error: %1.2e\n',norm(A-Abar));
    end
end

function tau = getThreshold(s)
    best = 0;
    for i = 2:length(s)-1
        if abs(s(i+1)-2*s(i)+s(i-1))<1
            tau = i;
            break;
        end
    end
end