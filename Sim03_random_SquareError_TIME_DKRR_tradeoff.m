clc; clear; close all; fclose all; format long; format compact;
addpath(genpath('algs_tools'));

% add log infos
start_run_time = datestr(now,'mmmm_dd_yyyy_HH:MM:SS')
run_file_name = mfilename

%% param settings
func_type = 'func2'; %  func2, circles, SD15
k_type = 2;   
points_type = 'random';   
NumExps = 50;  % the repeated times
nfolds = 5;
n_samples = 101^2;
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

num_local_machines_list = [6:3:60];

mse_list = zeros(length(num_local_machines_list),NumExps);
time_cost_list = zeros(length(num_local_machines_list),NumExps);

for num_idx = 1:1:length(num_local_machines_list)
    NumMj = num_local_machines_list(num_idx)

    for exp_idx = 1:1:NumExps
        rng(num_idx + exp_idx )

        dim = 3;
        xtr = unifsphere(n_samples, dim);
        Ntr = size(xtr,1);
        yptr = fun(xtr);
        ytr = yptr + randn(Ntr,1) * noise_level ;

        lambda_idx = 70; % such as: 30 for func2, 70 for SD15, 45 for circles
        opt_lambda_value=lambda_seq(lambda_idx);
        m = NumMj;
        [opt_quadrature_value, mean_time] = DKRR_bayesian_quadrature_with_lambda(xtr,ytr, KerPara, opt_lambda_value, m,  phi0_value);
        opt_param_error =  abs(opt_quadrature_value - true_value).^2

        mse_list( num_idx ,exp_idx ) = opt_param_error;
        time_cost_list( num_idx ,exp_idx ) = mean_time;
    end
end

mse_results = mean(mse_list,2);
time_results = mean(time_cost_list,2);

% save images
fig_plot = figure
hold on
mse_color = [0.6350 0.0780 0.1840];
time_color =[0 0.4470 0.7410]	;
marker_size = 8;
ax = gca;
ax.XAxis.Color = 'k';

yyaxis left
loglog(num_local_machines_list, mse_results,'o', ...
    'color', mse_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',mse_color,  ...
    'MarkerFaceColor',mse_color, ...
    'MarkerSize',marker_size);
ylabel('Square Error')
ax.YColor = mse_color;

marker_size = 8;
yyaxis right
loglog(num_local_machines_list, time_results,'o', ...
    'color', time_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',time_color,  ...
    'MarkerFaceColor',time_color, ...
    'MarkerSize',marker_size);
ylabel('Second')
ax.YColor = time_color;

legend( 'Integral Error' ,   ...
    'Computation Cost' , ...
    'Location' , 'north','TextColor','k','FontSize',14)

set(gca, 'FontName', 'Times New Roman','FontSize',12);
xlabel('m')
xlim([6, 60])
save_figs_name = [ 'sim_figures_new_dkrr/figs/' 'YspaceNoise_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps) '_noise_level' num2str(noise_level)  '_n_samples'  num2str(n_samples)  '_Errror_Time_tradeoff'];
saveas(fig_plot, [save_figs_name '.fig']);
save_pngs_name = [ 'sim_figures_new_dkrr/pngs/' 'YspaceNoise_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)  '_noise_level' num2str(noise_level)  '_n_samples'  num2str(n_samples)  '_Errror_Time_tradeoff'];
saveas(fig_plot, [save_pngs_name '.png']);
