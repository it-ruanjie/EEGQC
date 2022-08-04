function [pro, rkr, f]= my_qcspectra(EEG,nw,fs,fm)  %varargin为元胞数组
% Quality check of cross spectra  交叉光谱的质量检查
% Input
%       EEG: data (nc-nt-ns / nc-nt)
%       nw: time-halfbandwidth product
%       fs: sampling rate
%       fmax: maximum frequency stored in varargin{1}
%       fmin: minimum frequency  stored in varargin{2}
%       loc --- chanlocs varargin{3}
%       svpath --- varargin{4}
% Output
%         pro ---- proportion
%         rkr ---- [rank of the EEG waveforms, the power spectra and
%    crossspectr at 10Hz, mean correlation interchannels powers topomaps]
%    EEG波形的等级，10Hz时的功率谱和交叉谱，平均相关通道间功率地形图
%
% The sum of power spectra equals to the sum of explained variance/engergy
% (eigenvalues) 功率谱之和等于解释方差/能量之和（特征值）


% Shiang Hu, Jan. 9, 2020

% check if cross spectra failes 
if size(EEG,1)<2
    pro=0; 
    return; 
end

% fmax = fm;   % filtering
% fmin = 0.5;   % filtering
fmax = fm;   % filtering
fmin = 8;   % filtering
loc = 'E:\BrainStrom\brainstorm_db\Protocol02\data\Subject01\@rawsub-CBM00001_task-protmap_eeg\Standard-10-20-Cap19.locs';
%[svfd,nm]=fileparts(varargin{4});
[S, f, nss] = xspt(EEG,nw,fs,fmax,fmin);
fbd = [fmin fmax];

% idxing and referencing
nf = size(S,3);
n = nss*ones(1,nf);
pmax = 10;
lmax = 50;

% CPC
[lmd,Q] = CPCstepwise1(S,n,pmax,lmax);
Q = Q(:,1:pmax);
lmd = abs(lmd(1:pmax,:)');  % freq by CPs
psd = abs(tdiag(S));        % get the multichannel psd
ssd = sum(psd(:));          % sum of powers (explained variance)
% ssl = sum(lmd(:));        % sum of eigenvalues

% palos index
profd = lmd(:,1)/ssd;               % frequency dependent palosi of 1st cp
procd = sum(lmd)/ssd;               % component dependent palosi
% proch = Q./repmat(sum(Q.^2),[nc,1]);% channal weights in the CPC
pro = sum(profd);                   % total palosi

% check if plot results
if nargin==7
    figure
    % eigenvalues
    subplot(231), semilogy(real(lmd)),
    xlim(fbd); xlabel('Freq'); title('Log engergy');
    
    subplot(234), plot(real(lmd)),
    xlim(fbd); xlabel('Freq'); title('Natural engergy');
    
    % proportion of components
    subplot(232), pareto(procd);
    title('Proportion of cumulated > 95% variance');
    
    subplot(235),
    plot(1:pmax,cumsum(procd),'.',1:pmax,procd,'^');
    legend({'CumPro','VarPro'}); xlim([1 pmax]); xlabel('Ordered CPC');
    title('Proportion of each CPC');
    
    % components
    subplot(233), imagesc(real(Q)), colorbar;
    xlabel('Ordered CPC'), ylabel('Channels'); title('CPC patterns');
    
    subplot(236),
    topoplot(real(Q(:,1)),loc,'style','map'); title('CPC1: map');
    fg=gcf; fg.Position = [24 317 1105 585];
    saveas(gcf,fullfile(svfd,['Par_',nm]),'svg');
    
    figure,
    for i=1:8
        subplot(4,2,i)
        topoplot(real(Q(:,i)),loc,'style','map'); title(['CPC:',num2str(i)]);
    end
    fg=gcf; fg.Position = [1135 206 380 771];
    saveas(gcf,fullfile(svfd,['Map_',nm]),'svg');
end
close all;

% ouput other qc measures
rou = triu(corr(log10(psd)),1);
mr = sum(rou(:))./(nf*(nf-1)/2);
fid10 = (abs(f-10)==min(abs(f-10)));
rkr = [rank(EEG), rank(psd), rank(S(:, :, fid10)), mr];
end

