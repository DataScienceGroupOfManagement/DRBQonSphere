clc; clear; close all; fclose all; format long; format compact;
addpath(genpath('algs_tools'));

% add log infos
start_run_time = datestr(now,'mmmm_dd_yyyy_HH:MM:SS')
run_file_name = mfilename

%% param settings
func_type = 'SD15';
k_type = 2;
points_type = 'random';
NumExps = 50; % the repeated times
nfolds = 5;
noise_level = 1e-3;

%% the hyper-param list
if strcmp(func_type,'func2')
    q_value = 0.64;
else
    q_value = 0.74;
end

if k_type == 2
    switch func_type
        case 'func2'
            q_value = 0.61;
        case 'circles'
            q_value = 0.79;
        case 'SD15'
            q_value = 6/7;
    end
end

% the range of lambdas
mm = 100;
qq = 1.0/q_value;
lambda_seq = Lambda_q(qq, qq, mm);
lambda_seq(lambda_seq<1e-16) = [];

%% choose sim function
switch func_type
    case 'func2'
        fun = @(x) sim_func2(x);
        true_value =  6.6961822200736179523;
        disp('using func2')

    case 'circles'
        fun = @(x) sim_func_circles(x);
        exact_integral_value = -2*pi*(cos(pi/16) - 1) + -2*pi*(cos((3*pi)/16) - (2^(1/2) + 2)^(1/2)/2) + -2*pi*(cos((5*pi)/16) - 2^(1/2)/2);
        true_value =  exact_integral_value;
        disp('using circles')

    case 'SD15'
        fun = @(x) spherical_cap_SD15_func(x);
        true_value =    -2*pi*(cos(pi/32) - 1) * size(SD(15),1);
        disp('using SD15')
end

%% choose kernel function
switch k_type
    case 2
        KerPara.KernelType = 2;  % 使用 (1-r)^4*(4*r+1),    if  0<r<=1
        phi0_value = pi/7;
        disp('using (1-r)^4*(4*r+1), if  0<r<=1')
end

num_samples_list = [ 16000: -1000: 2000 ];
num_local_machines_list =  [3:3:30];
opt_lambda_index = 70;
mse_list = zeros(length(num_samples_list),  length(num_local_machines_list));
time_cost_list = zeros(length(num_samples_list),   length(num_local_machines_list));

for num_idx = 1:1:length(num_samples_list)
    n_samples = num_samples_list(num_idx);

    for m_idx = 1:1: length(num_local_machines_list)
        NumMj = num_local_machines_list(m_idx);
        tmp_avg_mc_list =  zeros( NumExps , 1);
        tmp_avg_drbq_list =  zeros( NumExps , 1);
        tmp_times_list = zeros( NumExps , 1);

        for exp_idx = 1:1:NumExps
            rng('shuffle')

            dim = 3;
            xtr = unifsphere(n_samples, dim);
            Ntr = size(xtr,1);
            yptr = fun(xtr);
            ytr = yptr + randn(Ntr,1) * noise_level ;

            % baseline monte-carlo
            equal_weights = ones(Ntr,1) * (pi*4)/Ntr;
            mc_integrate_value= ytr'*equal_weights;
            mc_error = abs(mc_integrate_value - true_value)^2;

            % DRBQ
            opt_lambda_value = lambda_seq(opt_lambda_index);
            m = NumMj;
            [opt_quadrature_value, mean_time] = DKRR_bayesian_quadrature_with_lambda(xtr,ytr, KerPara, opt_lambda_value, m,  phi0_value);
            opt_drbq_error =  abs(opt_quadrature_value - true_value)^2;
            flag = 0;
            if  opt_drbq_error <= mc_error
                flag = 1;
            end
            disp([  'n_samples = ' num2str(n_samples) ', NumMj = ' num2str( NumMj) ,' , repeated_times = ' num2str(exp_idx ), ', mc_error = ' num2str(mc_error) , ' , drbq_error = ' num2str(opt_drbq_error) ,' , flag = ' ,num2str(flag) ])
            tmp_avg_mc_list (exp_idx) = mc_error;
            tmp_avg_drbq_list (exp_idx) = opt_drbq_error;
            tmp_times_list (exp_idx) = mean_time;
        end
        avg_mc = mean(tmp_avg_mc_list) ;
        avg_drbq = mean(tmp_avg_drbq_list);
        avg_time = mean(tmp_times_list);
        if avg_drbq  <= avg_mc
            mse_list ( num_idx,  m_idx) = 1;
        else
            mse_list ( num_idx,  m_idx) = 0;
        end
        time_cost_list( num_idx,  m_idx) = avg_time;
    end
end

heatmap_data_matrix = mse_list;
heatmap_times_matrix = time_cost_list;
heatmap_samples = num_samples_list;
heatmap_machines = num_samples_list;

fig_plot_acc = figure,
heatmap(num_local_machines_list, num_samples_list, heatmap_data_matrix )
colormap('jet');
xlabel('m')
ylabel('|D|')
title('1 if SquareError(DRBQ) \leq SquareError(MC)')
set(gca, 'FontName', 'Times New Roman','FontSize',12);
save_figs_name = [ 'sim_figures_new_dkrr/figs/' 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)  '_noise_level' num2str(noise_level)  '_min_idx' num2str(opt_lambda_index)  '_dkrr_translation_map_accuracy'];
saveas(fig_plot_acc, [save_figs_name '.fig']);
save_pngs_name = [ 'sim_figures_new_dkrr/pngs/' 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps) '_noise_level' num2str(noise_level)  '_min_idx' num2str(opt_lambda_index)  '_dkrr_translation_map_accuracy'];
saveas(fig_plot_acc, [save_pngs_name '.png']);

time_threshod = 0.3;
heatmap_times_matrix(heatmap_times_matrix>time_threshod) = time_threshod;
fig_plot_time = figure,
heatmap(num_local_machines_list, num_samples_list, heatmap_times_matrix)
colormap('jet');
xlabel('m')
ylabel('|D|')
title([' Computation Time  \leq ' num2str(time_threshod), ' second'])
save_figs_name = [ 'sim_figures_new_dkrr/figs/' 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)  '_noise_level' num2str(noise_level) '_min_idx' num2str(opt_lambda_index)   '_dkrr_translation_map_time'];
saveas(fig_plot_time, [save_figs_name '.fig']);
save_pngs_name = [ 'sim_figures_new_dkrr/pngs/' 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps) '_noise_level' num2str(noise_level) '_min_idx' num2str(opt_lambda_index)   '_dkrr_translation_map_time'];
saveas(fig_plot_time, [save_pngs_name '.png']);
