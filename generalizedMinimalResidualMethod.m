function x = generalizedMinimalResidualMethod(A, b)
    if condest(A)>1e6
        tol = 1e-11;
        maxit = 200;
        [P,R,C] = equilibrate(A);
        B = R*P*A*C;
        d = R*P*b;
        [x,fl,rr,it,rv] = gmres(B,d,[],tol,maxit);
        x = C*x;
        if fl ~= 0
%             warning(['the algorithm does not converge', num2str(fl)]);
            x = A\b;
        end
    else 
        x = A\b;
    end
end