function [efficiency_global,efficiency_local,transitity,clustering_coef,degrees,mean_distance,global_diffusion_efficiency,edg] = bct_measures(icoh)
    
    efficiency_global = efficiency_bin(icoh);
    efficiency_local = efficiency_bin(icoh,1)';%1X15
    transitity = transitivity_bu(icoh);
    clustering_coef = clustering_coef_bu(icoh)';%1X15
    degrees = degrees_und(icoh);%1X15
    mean_distance = distance_bin(icoh);
    mean_distance = mean(mean(mean_distance));
    %[global_diffusion_efficiency,~] = diffusion_efficiency(icoh);
    global_diffusion_efficiency = 0;
    edg = length(find(icoh~=0))/2;%∂‘≥∆’Û
    
end

