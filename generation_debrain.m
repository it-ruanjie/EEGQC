clear;clc;
dataset = 200;
mu = zeros(15,1);
sigma = eye(15,15);
for num = 1:dataset
    clear samp;
    disp(['--------->dataset:' num2str(num) '<------------']);
    datapath = ['E:\BCI\my_workplace\icoh\ft_test0\data200\data' num2str(num) '\'];
    load([datapath,'samp.mat']);

    R = mvnrnd(mu,sigma,2500)';
    for t = 1:5
        debrain = samp - t*R/norm(R);
        EEG = pop_importdata('data',debrain,'chanlocs','E:\BCI\my_workplace\icoh\ft_test0\loc15.xyz','srate',250,'setname','debrain');
        pop_saveset(EEG,['debrain' num2str(t)],datapath);
    end
end