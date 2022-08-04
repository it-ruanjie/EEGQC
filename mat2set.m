clear;clc;
dataset = 200;
for num = 1:dataset
    disp(['--------->dataset:' num2str(num) '<------------']);
    datapath = ['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\'];
    load([datapath,'samp.mat']);
    EEG = pop_importdata('data',samp,'chanlocs','E:\BCI\my_workplace\icoh\ft_test0\loc15.xyz','srate',250,'setname','clean_eeg');
    pop_saveset(EEG,'samp',datapath);
end 