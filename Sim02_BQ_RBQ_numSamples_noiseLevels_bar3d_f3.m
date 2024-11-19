clc; clear; close all; fclose all; format long; format compact;
addpath(genpath('algs_tools'));

% add log infos
start_run_time = datestr(now,'mmmm_dd_yyyy_HH:MM:SS')
run_file_name = mfilename
 
% param settings
func_type = 'SD15';  
k_type =  2 ;  
NumExps = 50; % the repeated times

% the save path
save_path = './sim02_save_results_new/';
if ~ exist(save_path)
    mkdir( save_path )
end

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

% choose sim function
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

% choose kernel function
switch k_type
    case 2
        KerPara.KernelType = 2;  % 使用 (1-r)^4*(4*r+1),    if  0<r<=1
        phi0_value = pi/7;   
        disp('using (1-r)^4*(4*r+1), if  0<r<=1')
end
 
noise_level_list = [ 1e-6,  1e-5,  1e-4,  1e-3, 5e-3,  1e-2, 5e-2, 1e-1,  0.3,  0.5 ];
num_noises  = length(noise_level_list);

points_type = 'random';
num_samples_list = 500: 500: 3000;

inv_bar3_matrix = zeros( length(num_samples_list),   num_noises);
krr_bar3_matrix = zeros( length(num_samples_list),   num_noises);

for n_idx = 1:1:length(num_samples_list)
    n_samples = num_samples_list(n_idx);

    for l_idx= 1:1:length(noise_level_list)
        noise_level = noise_level_list(l_idx);
        tmp_inv = zeros( NumExps, 1);
        tmp_krr =  zeros( NumExps, 1);

        for exp_idx =1:1:NumExps
            rng(n_idx+l_idx+ exp_idx)

            dim = 3;
            xtr = unifsphere(n_samples, dim);
            yptr = fun(xtr);
            Ntr = size(xtr,1);
            ytr = yptr + randn(Ntr, 1) * noise_level;

            % BQ
            Ktr = KernelComputation(xtr, xtr, KerPara);
            Ktr = (Ktr + Ktr')/2;
            Vphi0 = ones(Ntr,1)*phi0_value;
            inv_weights  = Ktr\Vphi0;
            inv_integral_value = ytr'*inv_weights;
            inv_integral_error_value = abs(inv_integral_value - true_value).^2;

            % RBQ
            KPhi = KernelComputation(xtr, xtr, KerPara); % the kernel matrix
            KPhi = (KPhi + KPhi')/2;
            Vphi0 = ones(Ntr,1)*phi0_value;
            KRR_weights_seq = bayesian_quadrature_weights(KPhi, Vphi0, lambda_seq); % the weights
            KRR_integrate_value_seq = bayesian_quadrature_rules(ytr, KRR_weights_seq); % compute integration value
            KRR_error_seq = abs(KRR_integrate_value_seq - true_value).^2;
            KRR_baye_quad_error_value = min(KRR_error_seq);

            disp( [ 'points_type = ' points_type,  ' , n_samples = ' , num2str(Ntr), ...
                ', noise_level =' num2str(noise_level) , ' , exp_idx = ' num2str(exp_idx), ...
                ' , inv_result = ' num2str(inv_integral_error_value) , ...
                ' , krr_results = ' num2str(KRR_baye_quad_error_value)  ])

            tmp_inv(exp_idx) = inv_integral_error_value;
            tmp_krr(exp_idx) = KRR_baye_quad_error_value;
        end
        %% collect results
        inv_bar3_matrix(n_idx,  l_idx) =  mean(tmp_inv);
        krr_bar3_matrix (n_idx,  l_idx) = mean(tmp_krr );
    end
end

max_value = max(  max(max( inv_bar3_matrix) ) ,  max(max( krr_bar3_matrix) ) );

if ~ exist([save_path, 'figs/' ])
    mkdir( [save_path, 'figs/' ]  )
end

if ~ exist([save_path, 'pngs/' ])
    mkdir( [save_path, 'pngs/' ]  )
end

% figure 1
fig_plot_inv = figure,
hold on
b = bar3( inv_bar3_matrix)

set(gca,'XTick',   1:1:length(noise_level_list) )
set(gca, 'XTickLabel',  noise_level_list,'FontSize',12 )
set(gca, 'YTick',  1:1:length(num_samples_list)  )
set(gca, 'YTickLabel',  num_samples_list,'FontSize',12 )
set(gca, 'ZScale','log' )

zlim([0, max_value])
colorbar
caxis([0, max_value]);  
for k =1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor ='interp';
end

ylabel('|D|','FontSize',14 )
xlabel('\delta of N( 0, \delta^{2})','FontSize',14 )
zlabel('Square Error','FontSize',14 )
title('BQ','FontSize',14 )

save_figs_name = [save_path, 'figs/' , 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)     '_INV_bar3dV2_results'];
saveas(fig_plot_inv, [save_figs_name '.fig']);
save_pngs_name = [ save_path, 'pngs/' , 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)   '_NumExps'  num2str(NumExps)    '_INV_bar3dV2_results'];
saveas(fig_plot_inv, [save_pngs_name '.png']);

save_matrix_name = [ save_path, 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)   '_NumExps'  num2str(NumExps) '_INV_bar3dV2_results'];
writematrix(inv_bar3_matrix,[save_matrix_name,'.xls'])

% figure 2
fig_plot_krr = figure,
hold on
b = bar3( krr_bar3_matrix)
zlim([0, max_value])

set(gca,'XTick',   1:1:length(noise_level_list) )
set(gca, 'XTickLabel',  noise_level_list,'FontSize',12 )
set(gca, 'YTick',  1:1:length(num_samples_list)  )
set(gca, 'YTickLabel',  num_samples_list,'FontSize',12 )
set(gca, 'ZScale','log' )

colorbar
caxis([0, max_value]);  
for k =1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor ='interp';
end
 
ylabel('|D|','FontSize',14 )
xlabel('\delta of N( 0, \delta^{2})','FontSize',14 )
zlabel('Square Error','FontSize',14 )
title('RBQ','FontSize',14 )

save_figs_name = [save_path, 'figs/' , 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)     '_KRR_bar3dV2_results'];
saveas(fig_plot_krr, [save_figs_name '.fig']);
save_pngs_name = [ save_path, 'pngs/' , 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)   '_NumExps'  num2str(NumExps)    '_KRR_bar3dV2_results'];
saveas(fig_plot_krr, [save_pngs_name '.png']);

save_matrix_name = [ save_path, 'YnoiseGaussian_' points_type '_' func_type '_ktype' num2str(k_type)   '_NumExps'  num2str(NumExps)    '_KRR_bar3dV2_results'];
writematrix(krr_bar3_matrix,[save_matrix_name,'.xls'])

 