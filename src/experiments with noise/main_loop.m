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
sigma_e  = 0.01;          % noise scale is FIXED (do not change across cases)
T  = 10*T_D;
TT = 100;

Q = 100;
R = 1;

% ======= Enforce same output SNR across cases (EXACT formula) by scaling INPUT, noise fixed =======
SNRdB_target = 30;        % desired SNR in dB (based on clean output vs noise)
rng(10);
noise = sigma_e * randn(1, T);   % FIXED noise sequence + FIXED scale

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


            snr_tag  = sprintf('snr%02ddB', round(SNRdB_target));
            fname_big = sprintf('s%02d_r%02d_case%02d_g%d_%s.mat', ...
                                sysnum, f_target, si, gidx, snr_tag);

            big = struct();

            for i = 1:20
                rng(1);   % write-dead: inputs identical across runs

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
                        for j = 1:Jmax
                            phi_j = -pi * j*(j-1) / Jmax;
                            A_j   = 1;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                    case 2
                        Jmax = 55;
                        u_seq = zeros(size(k));
                        for j = 1:Jmax
                            phi_j = -pi * j*(j-1) / Jmax;
                            A_j   = 1;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                    case 3
                        j_min = 40;
                        j_max = 95;
                        Jnum  = j_max - j_min + 1;

                        u_seq = zeros(size(k));
                        idx = 0;
                        for j = j_min:j_max
                            idx = idx + 1;
                            phi_j = -pi * idx*(idx-1) / Jnum;
                            A_j   = 1;
                            u_seq = u_seq + A_j * sin( 2*pi*j*(f0/Fs)*k + phi_j );
                        end

                    otherwise
                        error('unknown si = %d', si);
                end

                % ======= match SNR by scaling INPUT (noise fixed) =======
                idx_ss = (T - T_D + 1) : T;

                % 1) noiseless output for current (unscaled) input
                [y0_no_noise, ~] = system_dynamics_no_noise(u_seq, x0, sysnum);

                s0  = y0_no_noise(idx_ss);          % clean output segment (signal)
                Es0 = sum(abs(s0).^2);              % sum |s_n|^2

                v   = noise(idx_ss);                % noise segment (fixed)
                Ev  = sum(abs(v).^2);               % sum |v_n|^2 (fixed)

                % 2) target signal energy from exact SNR formula:
                % SNRdB = 10*log10(Es/Ev)  =>  Es_target = Ev * 10^(SNRdB/10)
                Es_target = Ev * (10^(SNRdB_target/10));

                % 3) scale entire input so that Es matches Es_target (LTI scaling)
                scale_u = sqrt(Es_target / Es0);
                u_seq = u_seq * scale_u;

                % 4) re-simulate with fixed noise
                [y_seq_no_noise, ~] = system_dynamics_no_noise(u_seq, x0, sysnum);
                [y_seq, ~]         = system_dynamics(u_seq, x0, sysnum, noise);

                % ======= cut steady-state parts =======
                u_seq_cut = u_seq(idx_ss);
                y_seq_cut = y_seq(idx_ss);
                y_seq_cut_no_noise = y_seq_no_noise(idx_ss);

                % ======= report SNR using MATLAB snr() on clean output vs noise =======
                snr_dB_report = snr(y_seq_cut_no_noise, noise(idx_ss));
                fprintf('[sys=%d si=%d] Run %d | SNR(clean output vs noise) = %.3f dB\n', ...
                        sysnum, si, i, snr_dB_report);

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
