%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Matlab Code for Helmholtz Problem %
%           Zhenyang Yuan                 %
%           04/04/2020                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code can calculate the whole assigment at a time by applying both
%linear and quadratic expansion. Function documents 'Linear' and
%'Quadratic' are used to achieve the goal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Initialization
nel = [5,10,20,50,100];        %%element number appeared int he assignment
L = 1;                         %%domain length
er = zeros(length(nel),2);     %%preallocate for error analysis

for method = 1 : 2
    if method == 1             %%Linear expansion 
        Q = 2;
    elseif method == 2         %%Quadratic expansion
        Q = 3;
    end
        
    for count = 1:length(nel)
        n = nel(count);         %%number of elements
        x = linspace(0,L,n+1);  %%construct element nodes
        dx = L/n;               %%element length

        %% Construct Laplacian, Mass and RHS matrices
        L_e = cell(1,n);
        M_e = cell(1,n);
        f_e = cell(1,n);

        if (Q == 2)             %%using linear expansion in this case

            [LL,FF,MM,A]=Linear(L_e,M_e,f_e,n,dx,x);

        elseif (Q == 3)         %%using quadratic expansion in this case

            [LL,FF,MM,A]=Quadratic(L_e,M_e,f_e,n,dx,x);

        else            
            disp('error: Only Q = 2 or Q = 3 is accepted in this algrithm')
        end   


        %% Build global matrices %%alternative way(more efficient) is discussed in the report
        L_g = A'*LL*A;
        F_g = A'*FF;
        M_g = A'*MM*A;

        %% Apply boundary conditions(lower boundary is direchlet BC where u = 1; Upper boundary with numann with du = 0)
        F_g = F_g(2:end);                             %%on direchelt BC, the value is fixed, so elimated row and column respect to this node
        F_g = F_g - (M_g(2:end,1) + L_g(2:end,1))*1;  %%apply direchlet BC
        L_g = L_g(2:end,2:end);                         
        M_g = M_g(2:end,2:end);  

        %% Solver
        u = (M_g + L_g) \ F_g;                        %%invert to have homogenous solution, matlab back slash is fast enough to do this job                          
        if Q == 3
            n = 2*n;
            x = linspace(0,L,n+1);
        end

        u(2:n+1) = u;
        u(1) = 1;                                 %%reconstruct final solution by applying homogenous solution and boundary
        
        %% Plot
        uex = cos(2*pi*x');                           %%exact solution
        figure(method)
        plot(x,u,'Linewidth',2);                      %%plot result of simulation, figure will be rewrite when element number is changing
        hold on
        scatter(x,uex);
        hold off
        legend('Simulation','Exact solution','Location','southwest')
        xlabel('x')
        ylabel('u')
        if Q == 2
            title(['Linear finite element expansion, N_{el} = ',num2str(n)])
        else
            title(['Quadratic finite element expansion N_{el} = ',num2str(n/2)])
        end
        drawnow

        %% Error analysis
        du = uex - u;                                           
        L2 = norm(du)/sqrt(n+1);                      %%matlab norm2 is calculated as sum(abs(X).^2)^(1/2)
        er(count,method) = L2;
    end
end

figure(3)                                             %%plot error in both methids
loglog(1./nel,er(:,1),'b-',1./nel,er(:,2),'r--','Linewidth',2)
xlabel('Size of element:{1/N_{el}}')
ylabel('L_2 norm of error')
grid on
legend('Linear expension','Quadratic expension','Location','northwest')