clc; clear; close all; fclose all; format long; format compact;
addpath(genpath('algs_tools'));

% add log infos
start_run_time = datestr(now,'mmmm_dd_yyyy_HH:MM:SS')
run_file_name = mfilename

%% param settings
func_type = 'func2'; %  func2,  circles,  SD15
k_type =  2 ;   
points_type = 'random';   

degree  = 5; % the poly-degree of spherical harmonic
NumExps = 10;

save_path = './sim01_save_results_new/' ;
if ~ exist(save_path)
    mkdir( save_path )
end

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

num_samples_list = 200:200:5000;
lsq_integral_error_value_list = zeros(length(num_samples_list),NumExps);
mc_integral_error_value_list = zeros(length(num_samples_list),NumExps);
bq_integral_error_value_list = zeros(length(num_samples_list),NumExps);
rbq_integral_error_value_list  = zeros(length(num_samples_list),NumExps);

for num_idx= 1:1:length(num_samples_list)
    n_samples = num_samples_list(num_idx)

    for exp_idx =1:1:NumExps

        rng(num_idx+exp_idx)
        dim = 3;
        xtr = unifsphere(n_samples, dim);
        Ntr = size(xtr,1);
        yptr = fun(xtr);
        ytr = yptr; % noise free

        %% QMC
        equal_weights = ones(Ntr,1) * (pi*4)/Ntr;
        integrate_value= ytr'*equal_weights;
        mc_integral_error_value = abs(integrate_value - true_value)^2; 

        %% BQ
        Ktr = KernelComputation(xtr, xtr, KerPara);
        Ktr = (Ktr + Ktr')/2;
        Vphi0 = ones(Ntr,1)*phi0_value;
        inv_weights  = Ktr\Vphi0;
        inv_integral_value = ytr'*inv_weights;
        bq_error_value = abs(inv_integral_value - true_value)^2;

        %% RBQ
        KM = KernelComputation(xtr, xtr, KerPara); % the kernel matrix
        KM = (KM + KM')/2;
        Vphi0 = ones(Ntr,1)*phi0_value;
        weights_seq = bayesian_quadrature_weights(KM, Vphi0, lambda_seq); % the weights
        integrate_value_seq = bayesian_quadrature_rules(ytr, weights_seq); % compute integration value
        error_seq = abs(integrate_value_seq - true_value).^2;
        rbq_quad_error_value = min(error_seq);

        %% starting LSQ
        PN = (degree + 1)^2; % the number of basis functions
        Y_matrix = zeros(PN, Ntr);
        for n_idx = 1:1:Ntr
            P = zeros(PN, 1);
            c1 = xtr(n_idx,1);
            c2 = xtr(n_idx,2);
            c3 = xtr(n_idx,3);
            for l =0:1: degree
                for m=-l:1:l
                    p_idx = l * (l + 1) + m;
                    sh_value = sphHarm(l,m,c1,c2,c3);
                    P (int64(p_idx) +1 ) = sh_value;
                end
            end
            Y_matrix(:, n_idx) =  P ;
        end
        v = ones(Ntr,1) * (1/Ntr) ;
        ry =  zeros(PN, 1);
        Y_matrix(1,:) = 1;
        ry(1) = 1;
        b =  (Y_matrix * diag(v)* Y_matrix' ) \ ry;
        % b =  pinv(Y_matrix * diag(v)* Y_matrix') *ry;
        w_lsq = ( Y_matrix' * b) .* v;
        w_lsq = w_lsq *   4 * pi ; 
        integral_value = ytr' * w_lsq ;
        lsq_integral_error = abs( integral_value - true_value)^2; 

        %% collect results
        mc_integral_error_value_list (num_idx,exp_idx) = mc_integral_error_value;
        lsq_integral_error_value_list(num_idx,exp_idx) = lsq_integral_error;
        bq_integral_error_value_list(num_idx,exp_idx) =  bq_error_value; 
        rbq_integral_error_value_list(num_idx,exp_idx )= rbq_quad_error_value;

    end
end

% save images
fig_plot = figure
hold on
mc_color = [0.4660 0.6740 0.1880];
lsq_color =  [0  0 0.7410];
bq_color = [0 0 0]; 
rbq_color = [0.6350 0.0780 0.1840];
marker_size = 8; 
loglog(num_samples_list, mean(mc_integral_error_value_list,2),'o', ...
    'color',mc_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',mc_color, ...
    'MarkerFaceColor',mc_color,...
    'MarkerSize',10);hold on
loglog(num_samples_list, mean(lsq_integral_error_value_list,2),'^', ...
    'color', lsq_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',lsq_color,  ...
    'MarkerFaceColor',lsq_color, ...
    'MarkerSize',marker_size);hold on
loglog(num_samples_list, mean(bq_integral_error_value_list,2),'square', ...
    'color', bq_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',bq_color,  ...
    'MarkerFaceColor',bq_color, ...
    'MarkerSize', 12 );hold on
loglog(num_samples_list, mean(rbq_integral_error_value_list,2),'diamond', ...
    'color', rbq_color, 'Linewidth', 1.5, 'linestyle', ':',...
    'MarkerEdgeColor',rbq_color,  ...
    'MarkerFaceColor',rbq_color, ...
    'MarkerSize',marker_size);hold on
legend('MC', 'LSQ' ,  'BQ',  'RBQ', 'FontSize', 14)
xlabel('|D|','FontSize', 14)
ylabel('Square Error','FontSize', 14)
xlim([0, 5000 ])

set(gca,'FontName', 'Times New Roman', 'FontSize',12,  'YGrid','on',...
     'YMinorTick','on','YScale','log');

if ~ exist([save_path, 'figs/' ])
    mkdir( [save_path, 'figs/' ]  )
end
if ~ exist([save_path, 'pngs/' ])
    mkdir( [save_path, 'pngs/' ]  )
end
save_figs_name = [ save_path, 'figs/' 'noisefree_' points_type '_' func_type '_ktype' num2str(k_type)   '_NumExps'  num2str(NumExps)   '_QMC_LSQ_BQ_RBQ_results'];
saveas(fig_plot, [save_figs_name '.fig']);
save_pngs_name = [save_path,  'pngs/' ,'noisefree_' points_type '_' func_type '_ktype' num2str(k_type)  '_NumExps'  num2str(NumExps)   '_QMC_LSQ_BQ_RBQ_results'];
saveas(fig_plot, [save_pngs_name '.png']);


