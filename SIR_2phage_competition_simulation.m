clear all;
close all;


A1_value = [40];
A2_value = [100];
for i = 1:length(A1_value)   
    
 for j = 1:length(A2_value)
   SIR_2phage_competition_function(strcat(sprintf('Smesh_L=20_init=5E-1_ratio=1E-5_B=0p9_A1=%d',A1_value(i)),sprintf('_A2=%d.mat',A2_value(j))),A1_value(i),A2_value(j));
 end
 
end


