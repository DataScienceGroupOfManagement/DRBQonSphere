clear, clc,close all; 

currentFilePath = mfilename('fullpath');
currentFolder = fileparts(currentFilePath);
save_path = [currentFolder,'/'];

noise_level_list = [ 1e-6,  1e-5,  1e-4,  1e-3, 5e-3,  1e-2, 5e-2, 1e-1,  0.3,  0.5 ];
num_samples_list = 500: 500: 3000;

func_type = 'SD15';

tick_size = 16;
label_size = 18;
title_size = 22;

inv_bar3_matrix = readmatrix('YnoiseGaussian_random_SD15_ktype2_NumExps50_INV_bar3dV2_results.xls');
krr_bar3_matrix = readmatrix('YnoiseGaussian_random_SD15_ktype2_NumExps50_KRR_bar3dV2_results.xls');
advantage_KRR = inv_bar3_matrix - krr_bar3_matrix;

max_value = max(  max(max( inv_bar3_matrix) ) ,  max(max( krr_bar3_matrix) ) );

% figure 1
fig_plot_inv = figure,
hold on
b = bar3( inv_bar3_matrix)
% zlim([0, max_value]);
colorbar
caxis([0, max_value]);  
for k =1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor ='interp';
end
set(gca,'XTick', 1:1:length(noise_level_list) )
set(gca,'XTickLabel',  noise_level_list,'FontSize',tick_size )
set(gca,'YTick',  1:1:length(num_samples_list)  )
set(gca,'YTickLabel',  num_samples_list,'FontSize',tick_size )
set(gcf,'Position', [100, 100, 800, 600]);
set(gca,'FontName', 'Times New Roman', 'FontSize',tick_size,  'ZGrid','on',...
     'ZMinorTick','on');
ylabel('|D|','FontSize',label_size ,'FontWeight', 'bold' )
xlabel('\delta of N( 0, \delta^{2})','FontSize',label_size,'FontWeight', 'bold' )
zlabel('Square Error','FontSize',label_size ,'FontWeight', 'bold')
title('BQ','FontSize',title_size ,'FontWeight', 'bold')
azimuth = -135;
elevation = 20;
view(azimuth, elevation);
grid on; 

% figure 2
fig_plot_krr = figure,
hold on
b = bar3( krr_bar3_matrix)
% zlim([0, max_value])
colorbar
caxis([0, max_value]);  
for k =1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor ='interp';
end
set(gca,'XTick', 1:1:length(noise_level_list) )
set(gca,'XTickLabel',  noise_level_list,'FontSize',tick_size )
set(gca,'YTick',  1:1:length(num_samples_list)  )
set(gca,'YTickLabel',  num_samples_list,'FontSize',tick_size )
set(gcf,'Position', [100, 100, 800, 600]);
set(gca,'FontName', 'Times New Roman', 'FontSize',tick_size,  'ZGrid','on',...
     'ZMinorTick','on');
ylabel('|D|','FontSize',label_size,'FontWeight', 'bold'  )
xlabel('\delta of N( 0, \delta^{2})','FontSize',label_size,'FontWeight', 'bold' )
zlabel('Square Error','FontSize',label_size ,'FontWeight', 'bold')
title('RBQ','FontSize',title_size,'FontWeight', 'bold' )
azimuth = -135;
elevation = 20;
view(azimuth, elevation);
grid on; 

% figure 3
fig_plot_advantage = figure,
hold on
b = bar3( advantage_KRR)
% zlim([0, max_value])
colorbar
caxis([0, max_value]);  
for k =1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor ='interp';
end
set(gca,'XTick', 1:1:length(noise_level_list) )
set(gca,'XTickLabel',  noise_level_list,'FontSize',tick_size )
set(gca,'YTick',  1:1:length(num_samples_list)  )
set(gca,'YTickLabel',  num_samples_list,'FontSize',tick_size )
set(gcf,'Position', [100, 100, 800, 600]);
set(gca,'FontName', 'Times New Roman', 'FontSize',tick_size,  'ZGrid','on',...
     'ZMinorTick','on');
ylabel('|D|','FontSize',label_size,'FontWeight', 'bold' )
xlabel('\delta of N( 0, \delta^{2})','FontSize',label_size,'FontWeight', 'bold' )
zlabel('Square Error','FontSize',label_size,'FontWeight', 'bold' )
title('Superiority of RBQ = BQ-RBQ','FontSize',title_size,'FontWeight', 'bold' )
azimuth = -135;
elevation = 20;
view(azimuth, elevation);
grid on; 

save_pngs_name = [save_path,func_type,'_BQ_results_V2'];
saveas(fig_plot_inv, [save_pngs_name '.png']);

save_pngs_name = [save_path,func_type,'_RBQ_results_V2'];
saveas(fig_plot_krr, [save_pngs_name '.png']);

save_pngs_name = [save_path,func_type,'_advantage_results_V2'];
saveas(fig_plot_advantage, [save_pngs_name '.png']);
