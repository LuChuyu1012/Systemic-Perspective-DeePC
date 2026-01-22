



clc; clear; close all;

T_D   = 216;         
Ts    = 0.01;
Fs    = 1/Ts;
f_min = 0;           % Hz
f_max = Fs/2;        % 50 Hz

f_target = 2;       
T_ini = 5;           
N     = 15;          
L     = T_ini + N;

T  = 2000;         
TT = 100;           
nT = 0:(T-1);        % time index

Q = 100;             
R = 1;              

for sysnum = 1:6

    if ismember(sysnum, [1, 2, 3])
        dim = 2;
    else
        dim = 4;
    end
    x0 = zeros(dim,1);

    for si = 0:3

        fname_big = sprintf('sf%02d_r%02d_case%02d.mat', sysnum, f_target, si);
        big = struct(); 

        for i = 1:20
            rng(i); 

            switch si
                case 0
                    u_seq = randn(1, T);

                case 1
                    n_mults = 20;
                    u_seq = zeros(size(nT));
                    for k = 0:n_mults
                        phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder 
                        u_seq = u_seq + sin(((0 + 1*k) * 2*pi / Fs) * nT + phase_k);
                    end
                    u_seq = u_seq / (n_mults + 1);

                case 2
                    n_mults = 18;
                    u_seq = zeros(size(nT));
                    for k = 0:n_mults
                        phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder 
                        u_seq = u_seq + sin(((0 + 1*k) * 2*pi / Fs) * nT + phase_k);
                    end
                    u_seq = u_seq / (n_mults + 1);

                case 3
                    n_mults = 18;
                    u_seq = zeros(size(nT));
                    for k = 0:n_mults
                        phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder 
                        u_seq = u_seq + sin(((20 + 1*k) * 2*pi / Fs) * nT + phase_k);
                    end
                    u_seq = u_seq / (n_mults + 1);


              
              

                otherwise
                    error('unknown si = %d', si);
            end

            [y_seq, ~] = system_dynamics_no_noise(u_seq, x0, sysnum);

            u_seq_cut = u_seq(end - T_D + 1 : end);
            y_seq_cut = y_seq(end - T_D + 1 : end);
            
            desired_trajectory = sin( f_target * 2*pi / Fs * (0:(T-1)) );


            [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N);

            PE_rank          = rank([U_p; U_f]);
            condition_number = cond([U_p; U_f; Y_p; Y_f]);
            disp(rank([U_p; U_f;Y_p; Y_f]));

            % ---------- DeePC ----------
            tic;
            [u_seq_opt, y_seq_opt, x_seq, y_pred_seq, coverage, g_seq, y_ini_seq, u_ini_seq, y_seq_all] = ...
                run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, ...
                           x0, TT, sysnum,i);
            elapsed_time = toc;
            fprintf('[sys=%d si=%d] Run %d | execution time: %.6f s\n', sysnum, si, i, elapsed_time);

            desired_traj_TT = desired_trajectory(1:TT);
            error = desired_traj_TT - y_seq_opt;
            rms_error = sqrt(mean(error.^2));
            fprintf('[sys=%d si=%d] Run %d | Final RMS error: %.4f\n', sysnum, si, i, rms_error);

            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            big.(varname) = struct( ...
                'seed',          i, ...
                'PE_rank',       PE_rank, ...
                'cond_num',      condition_number, ...
                'rms_error',     rms_error, ...
                'coverage',      coverage(:), ...
                'error',         error(:), ...
                'y_seq_opt',     y_seq_opt(:), ...
                'u_seq_opt',     u_seq_opt(:), ...
                'elapsed_time',  elapsed_time ...
            );

        end % for i = 1:20

        save(fname_big, '-struct', 'big');
        fprintf('Saved: %s\n', fname_big);

    end % for si = 0:3

end % for sysnum = 1:6
