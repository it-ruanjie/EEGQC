clear;clc;
dataset = 200;
addpath 'E:\BCI\fieldtrip-master\connectivity'
% for num = 1:dataset
%     load(['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\samp.mat']);
%     [clean_pro(num), ~] = my_qcspectra(samp,3,250,13);
%     for noise_t = 1:41
%         EEGOUT = pop_loadset(['noise' num2str(noise_t) '.set'],['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num)]);
%         [noise_pro(num,noise_t), ~] = my_qcspectra(EEGOUT.data,3,250,13);
%     end
%     for debrain_t = 1:5
%         EEGOUT = pop_loadset(['debrain' num2str(debrain_t) '.set'],['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num)]);
%         [debrain_pro(num,debrain_t), ~] = my_qcspectra(EEGOUT.data,3,250,13);
%     end
%     
% end
% %----PaLos
% ft_data_result(clean_pro,noise_pro,debrain_pro,'changes of PaLos in different denoising level','PaLos','E:\BCI\my_workplace\icoh\ft_test0','PaLos');

method = {'coh','dtf','granger','pdc','plv','psi'}; 
method = {'granger','pdc','plv','psi'};
for method_i = 1:length(method)
    for num = 1:dataset
        disp(['--------->dataset:' num2str(num) '<------------']);
        clean_datapath = ['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\samp.set'];
        mfreq = ft_return_freq1(clean_datapath);
        clean_matrix = ft_multi_method(method{method_i},mfreq);
        [clean.efficiency_global(num),clean.efficiency_local(num,:),clean.transitity(num),clean.clustering_coef(num,:),clean.degrees(num,:),clean.mean_distance(num),~,clean.edg(num)] = bct_measures(clean_matrix);
        clean.global_clustering_coef(num) = clustCoeff(clean_matrix);
    
        for noise_t = 1:41
            
            noise_datapath = ['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\noise' num2str(noise_t) '.set'];
            mfreq = ft_return_freq1(noise_datapath);
            noise_matrix = ft_multi_method(method{method_i},mfreq);
            noise_DeVI(num,noise_t) = my_DeVI(clean_matrix,noise_matrix);
    
            [noise.efficiency_global(num,noise_t),noise.efficiency_local(num,noise_t,:),noise.transitity(num,noise_t),noise.clustering_coef(num,noise_t,:),noise.degrees(num,noise_t,:),noise.mean_distance(num,noise_t),~,noise.edg(num,noise_t)] = bct_measures(noise_matrix);
            noise.global_clustering_coef(num,noise_t) = clustCoeff(noise_matrix);
        end
    
        for debrain_t = 1:5
    
            debrain_datapath = ['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\debrain' num2str(debrain_t) '.set'];
            mfreq = ft_return_freq1(debrain_datapath);
            debrain_matrix = ft_multi_method(method{method_i},mfreq);
            debrain_DeVI(num,debrain_t) = my_DeVI(clean_matrix,debrain_matrix);
    
            [debrain.efficiency_global(num,debrain_t),debrain.efficiency_local(num,debrain_t,:),debrain.transitity(num,debrain_t),debrain.clustering_coef(num,debrain_t,:),debrain.degrees(num,debrain_t,:),debrain.mean_distance(num,debrain_t),~,debrain.edg(num,debrain_t)] = bct_measures(debrain_matrix);
            debrain.global_clustering_coef(num,debrain_t) = clustCoeff(debrain_matrix);
        end
    end
    mkdir(['E:\BCI\my_workplace\icoh\ft_test0\' method{method_i}]);
    result_path = ['E:\BCI\my_workplace\icoh\ft_test0\' method{method_i}];
    %----DeVI
    ft_data_result(zeros(200,1),noise_DeVI,debrain_DeVI,'changes of DeVI in different denoising level','DeVI',result_path,'DeVI');
    %----brain network
    ft_data_result(clean.edg,noise.edg,debrain.edg,'changes of edge in different denoising level','edge',result_path);
    ft_data_result(clean.efficiency_global,noise.efficiency_global,debrain.efficiency_global,'changes of global efficiency in different denoising level','global efficiency',result_path);
    ft_data_result(clean.transitity,noise.transitity,debrain.transitity,'changes of transitivity in different denoising level','transitivity',result_path);
    ft_data_result(clean.mean_distance,noise.mean_distance,debrain.mean_distance,'changes of average path length in different denoising level','average path length',result_path);
    ft_data_result(clean.global_clustering_coef,noise.global_clustering_coef,debrain.global_clustering_coef,'changes of global clustering coefficient in different denoising level','global clustering coefficient',result_path);
    ft_data_result(mean(clean.efficiency_local,2),mean(noise.efficiency_local,3),mean(debrain.efficiency_local,3),'changes of local efficiency in different denoising level','local efficiency',result_path);
    ft_data_result(mean(clean.clustering_coef,2),mean(noise.clustering_coef,3),mean(debrain.clustering_coef,3),'changes of mean clustering coefficient in different denoising level','mean clustering coefficient',result_path);
    %%-------------------------------------------------%%
end