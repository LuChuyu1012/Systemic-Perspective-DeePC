

function [g_opt, y_pred,u_opt] = solve_DeePC(U_p, Y_p, U_f, Y_f, u_ini, y_ini, r, Q, R,lambda_g)



    % Sizes
    n_uf = size(U_f, 1); % Number of future inputs
    n_yf = size(Y_f, 1); % Number of future outputs
    n_g = size(U_p, 2); % Number of selector elements
    ny_p = size(Y_p, 1); % Past-output stacked length




 cvx_begin quiet

        % Declare the optimization variables
        variable u(n_uf,1)
        variable y(n_yf,1)
        variable g(n_g)

        % Cost function initialization
        J = 0;
        
        for ii=1:length(r)
            
            e = y(ii) - r(ii);
            J = J + e'*Q*e + u(ii)'*R*u(ii);

        end
            J = J + lambda_g * sum_square_abs(g); 



        minimize(J)
        subject to
            % Model constraints
            U_p*g==u_ini;
            U_f*g==u;
            Y_p*g == y_ini;
            Y_f*g==y;
            u >= -3;
            u <= 3;
    cvx_end

    % Extract the optimal control action
    g_opt = g;
    u_opt = u;
    y_pred = y;    
    fprintf('g = [');
    fprintf(' %.6e', g_opt);
    fprintf(' ]\n');

    fprintf('||g||_2 = %.6e\n', norm(g_opt, 2));



end
