function [mfreq] = ft_return_freq(datapath)

sub_cfg = [];
sub_cfg.dataset = datapath;
data = ft_preprocessing(sub_cfg);

sub_cfg = [];
sub_cfg.order = 5;
sub_cfg.method = 'bsmart';
mdata = ft_mvaranalysis(sub_cfg,data);

sub_cfg = [];
sub_cfg.method = 'mvar';
mfreq = ft_freqanalysis(sub_cfg,mdata);
end

