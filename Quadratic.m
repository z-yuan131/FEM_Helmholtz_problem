function [LL,FF,MM,A]=Quadratic(L_e,M_e,f_e,n,dx,x)
    %% use Gauss-Lobatto quadrature, for Q = 3, w = [1/3,4/3,1/3].
    Je = dx/2;
    xi = [-1,0,1];
    w  = [1/3,4/3,1/3];
    
%%alternatively Gauss-Legendre quadrature can be used here as an option
%     xi = [-sqrt(3/5),0,sqrt(3/5)];        
%     w  = [5/9,8/9,5/9];

    phi = [xi.*(xi-1)/2;(1-xi).*(1+xi);xi.*(1+xi)/2];                   %%construct phi used in mass matrix
    dphi = [xi-ones(1,length(xi))/2;-2*xi;xi+ones(1,length(xi))/2];     %%construct dphi used in Laplacian matrix
    
    %%Equation in page 4-17 in the lecture note
    for i = 1 : n            %%construct mass matrix using numerical integration with Gaussian Quadrature
        M_e(:,i) = {[phi(1,1)*phi(1,1)*w(1)+phi(1,2)*phi(1,2)*w(2)+phi(1,3)*phi(1,3)*w(3),...
                     phi(1,1)*phi(2,1)*w(1)+phi(1,2)*phi(2,2)*w(2)+phi(1,3)*phi(2,3)*w(3),...
                     phi(1,1)*phi(3,1)*w(1)+phi(1,2)*phi(3,2)*w(2)+phi(1,3)*phi(3,3)*w(3);...
                     phi(2,1)*phi(1,1)*w(1)+phi(2,2)*phi(1,2)*w(2)+phi(2,3)*phi(1,3)*w(3),...
                     phi(2,1)*phi(2,1)*w(1)+phi(2,2)*phi(2,2)*w(2)+phi(2,3)*phi(2,3)*w(3),...
                     phi(2,1)*phi(3,1)*w(1)+phi(2,2)*phi(3,2)*w(2)+phi(2,3)*phi(3,3)*w(3);...
                     phi(3,1)*phi(1,1)*w(1)+phi(3,2)*phi(1,2)*w(2)+phi(3,3)*phi(1,3)*w(3),...
                     phi(3,1)*phi(2,1)*w(1)+phi(3,2)*phi(2,2)*w(2)+phi(3,3)*phi(2,3)*w(3),...
                     phi(3,1)*phi(3,1)*w(1)+phi(3,2)*phi(3,2)*w(2)+phi(3,3)*phi(3,3)*w(3)]*Je};%
    end
    
    for i = 1 : n           %%construct Laplacian matrix using numerical integration with Gaussian Quadrature
        L_e(:,i) = {[dphi(1,1)*dphi(1,1)*w(1)+dphi(1,2)*dphi(1,2)*w(2)+dphi(1,3)*dphi(1,3)*w(3),...
                     dphi(1,1)*dphi(2,1)*w(1)+dphi(1,2)*dphi(2,2)*w(2)+dphi(1,3)*dphi(2,3)*w(3),...
                     dphi(1,1)*dphi(3,1)*w(1)+dphi(1,2)*dphi(3,2)*w(2)+dphi(1,3)*dphi(3,3)*w(3);...
                     dphi(2,1)*dphi(1,1)*w(1)+dphi(2,2)*dphi(1,2)*w(2)+dphi(2,3)*dphi(1,3)*w(3),...
                     dphi(2,1)*dphi(2,1)*w(1)+dphi(2,2)*dphi(2,2)*w(2)+dphi(2,3)*dphi(2,3)*w(3),...
                     dphi(2,1)*dphi(3,1)*w(1)+dphi(2,2)*dphi(3,2)*w(2)+dphi(2,3)*dphi(3,3)*w(3);...
                     dphi(3,1)*dphi(1,1)*w(1)+dphi(3,2)*dphi(1,2)*w(2)+dphi(3,3)*dphi(1,3)*w(3),...
                     dphi(3,1)*dphi(2,1)*w(1)+dphi(3,2)*dphi(2,2)*w(2)+dphi(3,3)*dphi(2,3)*w(3),...
                     dphi(3,1)*dphi(3,1)*w(1)+dphi(3,2)*dphi(3,2)*w(2)+dphi(3,3)*dphi(3,3)*w(3)]/Je};%
    end
    
    for i = 1 : n           %%construct forcing matrix using numerical integration with Gaussian Quadrature
        f  = (4*pi*pi+1)*cos(2*pi*(dx/2*(ones(1,3)+xi)+x(i)*ones(1,3)));
        f_e(:,i) = {[phi(1,1)*f(1)*w(1)+phi(1,2)*f(2)*w(2)+phi(1,3)*f(3)*w(3);...
                     phi(2,1)*f(1)*w(1)+phi(2,2)*f(2)*w(2)+phi(2,3)*f(3)*w(3);...
                     phi(3,1)*f(1)*w(1)+phi(3,2)*f(2)*w(2)+phi(3,3)*f(3)*w(3)]*Je};
    end
   
    % There is an alternative way(more efficient way is introduced in the report to do these two step
    %%matrix A
    A_e = [1,0;1,0;0,1];
    A = sparse(kron(eye(n),A_e));
    A(1,:)=[];
    A(3*n,2*n+1) = 1;
    
    %%extend matrix L and f and M
    LL = zeros(3*n,3*n);
    FF = zeros(3*n,1);
    MM = zeros(3*n,3*n);
    for i = 1 : n
        LL(3*i-2:3*i,3*i-2:3*i) = cell2mat(L_e(i));
        FF(3*i-2:3*i,1) = cell2mat(f_e(i));
        MM(3*i-2:3*i,3*i-2:3*i) = cell2mat(M_e(i));
    end
end