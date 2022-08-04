function [mfreq] = ft_return_freq1(datapath)

sub_cfg = [];
sub_cfg.dataset = datapath;
data = ft_preprocessing(sub_cfg);

sub_cfg = [];
sub_cfg.method = 'mtmfft';
sub_cfg.taper = 'dpss';
sub_cfg.output = 'fourier';
sub_cfg.tapsmofrq = 2;
mfreq = ft_freqanalysis(sub_cfg, data);
end
