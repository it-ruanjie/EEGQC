function [DeVI] = my_DeVI(matrice1, matrice2)
   
   [sim,~] = my_siminet(matrice1, matrice2);
   DeVI = 1/sim-1;
end

