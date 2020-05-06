function [LL,FF,MM,A]=Linear(L_e,M_e,f_e,n,dx,x)
    %% use Gauss-Legendre quadrature, for Q = 2, w = 1, so we ignore w in implement
    Je = dx/2;
    xi = [1/sqrt(3),-1/sqrt(3)];
    phi = [1+xi;1-xi]/2;        %%construct phi used in mass matrix
    dphi = [-1,1]/2;            %%construct dphi used in Laplacian matrix
    
    %%Equation in page 4-17 in the lecture note
    for i = 1 : n               %%construct mass matrix using numerical integration with Gaussian Quadrature
        M_e(:,i) = {[phi(1,1)*phi(1,1)+phi(1,2)*phi(1,2),phi(1,1)*phi(2,1)+phi(1,2)*phi(2,2);...
                     phi(2,1)*phi(1,1)+phi(2,2)*phi(1,2),phi(2,1)*phi(2,1)+phi(2,2)*phi(2,2)]*Je};%
    end
    
    for i = 1 : n               %%construct Laplacian matrix using numerical integration with Gaussian Quadrature
        L_e(:,i) = {[dphi(1)*dphi(1),dphi(1)*dphi(2);...
                     dphi(2)*dphi(1),dphi(2)*dphi(2)]*2/Je};%
    end
    
    for i = 1 : n               %%construct forcing matrix using numerical integration with Gaussian Quadrature
        f  = (4*pi*pi+1)*cos(2*pi*(dx/2*(ones(1,2)+xi)+x(i)*ones(1,2)));
        f_e(:,i) = {[phi(1,1)*f(1)+phi(1,2)*f(2);...
                     phi(2,1)*f(1)+phi(2,2)*f(2)]*Je};
    end
    
    % There is an alternative way(more efficient way is introduced in the report to do these two step
    %%matrix A
    A = sparse(2*n,n+1);
    A(1,1) = 1;
    A(end,end) = 1;
    A(2:end-1 , 2:end-1) = kron(eye(n-1),[1;1]);
    
    %%extend matrix L and f and M
    LL = zeros(2*n,2*n);
    FF = zeros(2*n,1);
    MM = zeros(2*n,2*n);
    for i = 1 : n
        LL(2*i-1:2*i,2*i-1:2*i) = cell2mat(L_e(i));
        FF(2*i-1:2*i,1) = cell2mat(f_e(i));
        MM(2*i-1:2*i,2*i-1:2*i) = cell2mat(M_e(i));
    end
end