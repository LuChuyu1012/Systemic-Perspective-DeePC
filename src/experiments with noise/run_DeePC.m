

function [u_seq, y_seq,x_seq,y_pred_seq] = run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R,x0, T,sys,lambda_g,noise,randseed) 
rng(randseed)
    % Runs the DeePC algorithm in a receding horizon manner.
    % Inputs:
    %   U_p, Y_p, U_f, Y_f - Hankel matrices
    %   desired_trajectory - Full reference trajectory
    %   T_ini - Initial condition length
    %   N - Prediction horizon
    %   Q, R - Performance weight matrices
    %   lambda_g, lambda_y - Regularization parameters
    %   T           - Total simulation steps
    % Output:
    %   u_seq - Optimized control sequence over time
    
    % Initialize input and output sequences
    u_seq = zeros(T, 1);  % Stores the control inputs applied
    nx=size(x0,1);

    y_seq = zeros(T, 1);  % Stores the measured system outputs
    x_seq = zeros(nx, T+1);  % Stores the system states
    y_pred_seq = zeros(T, N);    
    y_ini_seq = zeros(T, T_ini);    
    u_ini_seq = zeros(T, T_ini);    

    
    % 1. Initialize past input and output data

     u_ini=4*rand(1, T_ini)' -2;

    y_ini = zeros(T_ini, 1);  

    % Compute initial y_ini from system
    x_tmp = x0;                 
    for i = 1:T_ini
        [y_ini(i), x_tmp] = system_dynamics_global(u_ini(i), x_tmp,sys,i,noise);
    end
    x_seq(:,1) = x_tmp; 



    % 2. Control loop (rolling horizon optimization)
    for t = 1:T
        if t > 1
            u_ini = [u_ini(2:end); u_seq(t-1)];
            y_ini = [y_ini(2:end); y_seq(t-1)];
        end

        %  Define reference trajectory for next N steps
        r = desired_trajectory(t:t+N-1);

        % Solve the DeePC optimization problem
        [g_opt, y_pred,u_opt] = solve_DeePC(U_p, Y_p, U_f, Y_f, u_ini, y_ini, r, Q, R,lambda_g);
        y_pred_seq(t,:) = y_pred';
        y_ini_seq(t,:) = y_ini';    
        u_ini_seq(t,:) = u_ini';


        % 5. Compute the optimal input sequence
        if ~isempty(g_opt)
            u_seq(t) = u_opt(1); % Apply only the first input u_0*
        else
            keyboard;
            warning('DeePC optimization failed at step %d, keeping previous control input.', t);
            if t > 1
                u_seq(t) = u_seq(t-1); % Keep previous control input if optimization fails
            end
        end
        

        % 6. Update measurements
        [ y_seq(t), x_next ] = system_dynamics_global(u_seq(t), x_seq(:,t),sys,t+T_ini,noise);
       



        x_seq(:,t+1) = x_next; 

        
        
        fprintf('Step %d | u = %.4f | y = %.4f\n', t, u_seq(t), y_seq(t));  

    end
    y_seq = y_seq';

    % Return the computed control input sequence
    disp('DeePC optimization completed.');
end






