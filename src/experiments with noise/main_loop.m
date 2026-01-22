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

lambda_g = 0.01;
sigma_e  = 0.005;
T  = 10*T_D;
TT = 100;

Q = 100;
R = 1;

rng(114514);
noise = sigma_e*randn(1, T);

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

            sigma_tag = sprintf('se%s', erase(sprintf('%.0e', sigma_e), '+'));
            fname_big = sprintf('s%02d_r%02d_case%02d_g%d_%s.mat', ...
                                sysnum, f_target, si, gidx, sigma_tag);
            big = struct();

            for i = 1:20
                rng(i);

                % ======= (only used by multisine) =======
                k  = 0:(T-1);       % time index
                Nper = T_D;
                f0   = Fs / Nper;   % f0 = fs/N

                switch si
                    case 0
                        u_seq = randn(1, T);

                    case 1
                        Jmax = 60;
                        u_seq = zeros(size(k));

                        for j = 0:Jmax
                            phi_j = -pi * j*(j-1) / Jmax;
                            A_j   = 5;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                        u_seq = u_seq / Jmax;

                    case 2
                        % 0:0.5:27.5 Hz -> Jmax = 27.5 / f0 = 55
                        Jmax = 55;
                        u_seq = zeros(size(k));

                        for j = 0:Jmax
                            phi_j = -pi * j*(j-1) / Jmax;
                            A_j   = 5;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                        u_seq = u_seq / Jmax;

                    case 3
                        % 20:0.5:47.5 Hz
                        % j_min = 20/f0 = 40, j_max = 47.5/f0 = 95
                        j_min = 40;
                        j_max = 95;
                        Jnum  = j_max - j_min + 1;

                        u_seq = zeros(size(k));
                        idx = 0;
                        for j = j_min:j_max
                            idx = idx + 1;
                            phi_j = -pi * idx*(idx-1) / Jnum;
                            A_j   = 5;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                        u_seq = u_seq / Jnum;

                    otherwise
                        error('unknown si = %d', si);
                end

                [y_seq, ~] = system_dynamics(u_seq, x0, sysnum, noise);

                u_seq_cut = u_seq(end - T_D + 1 : end);
                y_seq_cut = y_seq(end - T_D + 1 : end);

                desired_trajectory = sin( f_target * 2*pi / Fs * (0:(T-1)) );

                [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N);

                PE_rank          = rank([U_p; U_f]);
                condition_number = cond([U_p; U_f; Y_p; Y_f]);

                % ---------- DeePC ----------
                tic;
                [u_seq_opt, y_seq_opt, x_seq, y_pred_seq] = ...
                    run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, ...
                              x0, TT, sysnum, lambda_g_curr, noise, i);
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
                    'rms_error',     rms_error, ...
                    'error',         error(:), ...
                    'y_seq_opt',     y_seq_opt(:), ...
                    'u_seq_opt',     u_seq_opt(:), ...
                    'noise',         noise ...
                );

            end % for i = 1:20

            save(fname_big, '-struct', 'big');
            fprintf('Saved: %s\n', fname_big);

        end % for si = 0:3

    end % for sysnum = 1:6
end % for gidx
