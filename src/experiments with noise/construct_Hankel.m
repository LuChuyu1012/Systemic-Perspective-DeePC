function [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_data, y_data, T_ini, N)
    % Construct past and future Hankel matrices with L = T_ini + N.
    % Inputs:
    %   u_data - Input sequence (1×T)
    %   y_data - Output sequence (1×T)
    %   T_ini - Initial state estimation window
    %   N - Prediction horizon
    % Outputs:
    %   U_p, U_f - Past and future Hankel blocks for input
    %   Y_p, Y_f - Past and future Hankel blocks for output

    % Compute L
    L = T_ini + N;

    % Get total sequence length
    T = length(u_data);

    % Check if L is valid
    if L > T
        error('L must be less than or equal to the length of data sequences');
    end

    % Construct full Hankel matrices
    H_u = hankel(u_data(1:L), u_data(L:T));
    H_y = hankel(y_data(1:L), y_data(L:T));

    % Split into past (T_ini) and future (N) blocks
    U_p = H_u(1:T_ini, :);
    U_f = H_u(T_ini+1:end, :);
    Y_p = H_y(1:T_ini, :);
    Y_f = H_y(T_ini+1:end, :);
%      Combine all into a matrix for rank check
%     M = [U_p; U_f];
    M = H_u;

    rM = rank(M);
    if rM < size(M,1)
        fprintf("Rank deficiency detected: rank = %d < %d. Data does NOT span full behavioral space.", rM, size(M,1));
    else
        fprintf("[OK] Input matrix full rank: %d / %d\n", rM, size(M,1));
    end

end


