function [matrix] = ft_multi_method(method,mfreq)

   switch method
       case 'coh'
            cfg = [];
            cfg.method = 'coh';
            cfg.complex = 'imag';
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.cohspctrm(:,:,stat.freq == 10);
            matrix = abs(matrix);
            matrix(find(matrix<0.05)) = 0;
            matrix(find(matrix~=0)) = 1;
       case 'dtf'
            cfg = [];
            cfg.method = 'dtf';
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.dtfspctrm(:,:,stat.freq == 10);
            matrix = abs(matrix);
            matrix=matrix-diag(diag(matrix));%对角线置0
            matrix = (matrix + matrix.')/2;
            matrix(find(matrix<0.2)) = 0;
            matrix(find(matrix~=0)) = 1;
       case 'granger'
            cfg = [];
            cfg.method = 'granger';
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.grangerspctrm(:,:,stat.freq == 10);
            matrix = matrix/norm(matrix);%归一化
            matrix = abs(matrix);
            matrix = (matrix + matrix.')/2;
            matrix(find(matrix<0.06)) = 0;
            matrix(find(matrix~=0)) = 1;
       case 'pdc'
            cfg = [];
            cfg.method = 'pdc';
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.pdcspctrm(:,:,stat.freq == 10);
            matrix = abs(matrix);
            matrix=matrix-diag(diag(matrix));%对角线置0
            matrix = (matrix + matrix.')/2;
            matrix(find(matrix<0.2)) = 0;
            matrix(find(matrix~=0)) = 1;
       case 'plv'
            cfg = [];
            cfg.method = 'plv';
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.plvspctrm(:,:,stat.freq == 10);
            matrix = abs(matrix);
            matrix=matrix-diag(diag(matrix));%对角线置0
            matrix(find(matrix<0.4)) = 0;
            matrix(find(matrix~=0)) = 1;
            
       case 'psi'
            cfg = [];
            cfg.method = 'psi';
            cfg.bandwidth = 2;
            stat = ft_connectivityanalysis(cfg, mfreq);
            matrix = stat.psispctrm(:,:,stat.freq == 10);
            matrix = abs(matrix);
            matrix=matrix-diag(diag(matrix));%对角线置0
            matrix(find(matrix<0.005)) = 0;
            matrix(find(matrix~=0)) = 1;
   end
end

