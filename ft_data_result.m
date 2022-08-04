function [] = ft_data_result(clean_data,noise_data,debrain_data,result_title,y_text,result_path,varargin)

   i=1;
   for snr=-10:5:30
       com(:,i) = noise_data(:,snr+11);
       i = i+1;
   end
   com(:,10) = clean_data;
   com(:,11:15) = debrain_data;
   figure(1);
   boxplot(com,'notch','on');
   set(gca,'xticklabels','')
   ylabel(y_text);
   h = findobj(gca,'Tag','Box');
   title(result_title);
   set(gca, 'FontName', 'Times New Roman','FontSize',14);
   for j=1:length(h) 
        if j<=5 
            h(j).Color = 'b';
        elseif j==6 
            h(j).Color = 'g';
        else 
            h(j).Color = 'm';
        end
   end
   hold on;
   ch = colorbar('horiz');
   set(ch,'ticklabels','');
%    set(get(ch,'title'),'string','Denoising   Level','position',[600 -20],'FontSize',14,'FontName','Times New Roman','FontWeight','bold');
   set(get(ch,'title'),'string','Denoising   Level','position',[160 -20],'FontSize',14,'FontName','Times New Roman','FontWeight','bold');
   saveas(gcf,[result_path '\' y_text '.png']);
   hold off;
   
   if nargin>6  %--------palos  DeVi
        clear com;
        com = noise_data;
        com(:,42) = clean_data;
        com(:,43:47) = debrain_data;
        figure(2);
        imagesc([-10 36],[1 200],com);colorbar;axis xy;
        set(gca, 'FontName', 'Times New Roman','FontSize',14);
        set(gca,'xticklabels','')
        ylabel('simulation');
        title([y_text ' : 200 simulations of different SNR']);
        hold on;
        ch = colorbar('horiz');
        set(ch,'ticklabels','');
        set(get(ch,'title'),'string','Denoising   Level','position',[160 -20],'FontSize',14,'FontName','Times New Roman','FontWeight','bold');
        saveas(gcf,[result_path '\' y_text '0.png']);
   end
end

