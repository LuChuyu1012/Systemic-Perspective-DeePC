function  [y_seq,x_seq] = system_dynamics_no_noise(u_seq,x0,sys_id)
    % Discrete-time system simulation (handles input as a sequence)
    % Inputs:
    %   x0 - Initial state (n×1)
    %   u_seq - Control input sequence (1×T)
    % Outputs:
    %   x_seq - State sequence (n×T)
    %   y_seq - Output sequence (p×T)
    % System matrices
      switch sys_id
        case 1
            % system 1
            A = [0.7326, -0.0861;
                 0.1722,  0.9909];
            B = [0.0609;
                 0.0064];
            C = [0, 1.4142];

        case 2
            % system 2
            A = [0.85, 0.30;
                 -0.20, 0.75];
            B = [1.0;
                 0.2];
            C = [0.0, 1.0];

        case 3
            % system 3
            A = [0.9000, -0.4000;
                 0.6000,  0.7000];
            B = [0.0000;
                 1.0000];
            C = [1.0000, 0.0000];
            
        case 4
            % system 4
            A = [  2.20632892, -2.29641076,  1.49610313, -0.5184;
                    1           0           0           0;
                    0           1           0           0;
                    0           0           1           0 ];
            B = [1; 0; 0; 0];
            C = [ 0, 1, -1.17557050, 1 ];
            

            
       case 5
            % system 5
            A = [  1.01192665, -0.31215078,  0.33767358, -0.50765625;
                     1           0           0           0;
                     0           1           0           0;
                     0           0           1           0 ];
            B = [1; 0; 0; 0];
            C = [ 0, 1, -0.21900825, 0.49 ];
            

        otherwise
            error('invalid sys_id');
     end

    % Get the length of the input sequence
    T = length(u_seq);
    
    % Initialize state and output sequences
    n = size(A,1); % State dimension
    p = size(C,1); % Output dimension
    x_seq = zeros(n, T);
    y_seq = zeros(p, T);

    % Initial state
    x = x0;

    % Iterate through the time steps
    for t = 1:T

        % Compute the next state
        x_next = A * x + B * u_seq(t);

        % Compute the output
        y = C * x;

        % Store results
        x_seq(:, t) = x;
        y_seq(:, t) = y;

        % Update state
        x = x_next;
    end
end
