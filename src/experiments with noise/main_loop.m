clc; clear; close all;
T_D   = 200;         
Ts    = 0.01;
Fs    = 1/Ts;
f_min = 0;           % Hz
f_max = Fs/2;        % 50 Hz

f_target = 2;       
T_ini = 8;           
N     = 25;          
L     = T_ini + N;



lambda_g=0.01;     
sigma_e  = 0.005;
T  = 20000;         
TT = 100;           
nT = 0:(T-1);        % time index

Q = 100;             
R = 1; 


rng(114514);
noise=sigma_e*randn(1, T);


for gidx = 0:3                        
    lambda_g_curr = lambda_g * 10^(gidx-1);

    for sysnum = 1:6

        if ismember(sysnum, [1, 2, 3])
            dim = 2;
        else
            dim = 4;
        end 
        x0 = zeros(dim,1);





        for si = 0:3
            sigma_tag = sprintf('se%s', erase(sprintf('%.0e', sigma_e), '+'));  % 0.01 -> 'se1e-02'
fname_big = sprintf('sb%02d_r%02d_case%02d_g%d_%s.mat', ...
                    sysnum, f_target, si, gidx, sigma_tag);            big = struct(); 

            for i = 1:20
                rng(i+1); 

                switch si
                    case 0
                        u_seq = randn(1, T);

                    case 1
                        n_mults = 60;
                        u_seq = zeros(size(nT));
                        for k = 0:n_mults
                            % phase_k = rand * 2 * pi;
                            phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder

                            u_seq = u_seq + 5*sin(((0 + 0.5*k) * 2*pi / Fs) * nT + phase_k);
                        end
                        u_seq = u_seq / (n_mults + 1);

                    case 2
                        n_mults = 55;
                        u_seq = zeros(size(nT));
                        for k = 0:n_mults
                            % phase_k = rand * 2 * pi;
                            phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder

                            u_seq = u_seq + 5*sin(((0 + 0.5*k) * 2*pi / Fs) * nT + phase_k);
                        end
                        u_seq = u_seq / (n_mults + 1);

                    case 3
                        n_mults = 55;
                        u_seq = zeros(size(nT));
                        for k = 0:n_mults
                            % phase_k = rand * 2 * pi;
                            phase_k = -pi * k*(k-1) / n_mults; % ← Schroeder

                            u_seq = u_seq + 5*sin(((20 + 0.5*k) * 2*pi / Fs) * nT + phase_k);
                        end
                        u_seq = u_seq / (n_mults + 1);

                    otherwise
                        error('unknown si = %d', si);
                end

                [y_seq, ~] = system_dynamics(u_seq,x0, sysnum,noise);
                [y_seq_clean, ~] = system_dynamics_no_noise(u_seq,x0, sysnum);

                u_seq_cut = u_seq(end - T_D + 1 : end);
                y_seq_cut = y_seq(end - T_D + 1 : end);
                y_seq_clean = y_seq_clean(end - T_D + 1 : end);




                SNRdB_yn = snr(y_seq_clean, noise(1:T_D)); 
                SNRdB_un = snr(u_seq_cut, noise(1:T_D)); 

                desired_trajectory = sin( f_target * 2*pi / Fs * (0:(T-1)) );

                [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N);

                PE_rank          = rank([U_p; U_f]);
                condition_number = cond([U_p; U_f; Y_p; Y_f]);



                % ---------- DeePC ----------
                tic;
                [u_seq_opt, y_seq_opt, x_seq, y_pred_seq, coverage, g_seq, y_ini_seq, u_ini_seq, y_seq_all] = ...
                    run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, ...
                                x0, TT, sysnum,lambda_g_curr,noise,i);   
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
                    'SNRYN',     SNRdB_yn, ...
                    'SNRUN',     SNRdB_un, ...
                    'u_seq_cut',     u_seq_cut, ...
                    'y_seq_cut',     y_seq_cut, ...
                    'y_seq_clean',     y_seq_clean, ...
                    'noise',  noise ...
                );

            end % for i = 1:20

            save(fname_big, '-struct', 'big');
            fprintf('Saved: %s\n', fname_big);

        end % for si = 0:3

    end % for sysnum = 1:6
end % for gidx
