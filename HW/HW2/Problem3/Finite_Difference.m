function [x] = Finite_Difference(n)
    DMain = ones(1,n)*2;
    DUp = ones(1,n)*-1;
    DLow = ones(1,n)*-1;
    D = [DUp; DMain; DLow];
    A = spdiags(D',-1:1,n,n);
    p = 0:1/(n+1):1;
    p = p(2:n+1);
    f = -(2*pi)^2 * sin(2*pi*p);
    b = -(1/(n+1))^2 *f;
    x = inv(A)*b';
end

