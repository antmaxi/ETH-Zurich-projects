function [w_x, w_y, E] = tps_model(Xunwarped,Y,lambda)
    n = size(Y,1);
    P = [ones(n,1), Xunwarped];
    K = U(Xunwarped);
    Kl = K + lambda.*eye(n);
    A = [Kl, P;
         P', zeros(3,3)];
    bx = [Y(:,1); zeros(3, 1)];
    by = [Y(:,2); zeros(3, 1)];
    w_x = (A\bx);
    w_y = (A\by);
    E = w_x(1:n)'*K*w_x(1:n) + w_y(1:n)'*K*w_y(1:n);
    disp(E);
end



function K = U(X)
t = squareform(pdist(X));
K = zeros(size(t,1), size(t,1));
K(t ~= 0) = 2*t(t ~= 0).^2.*log(t(t ~= 0));
end