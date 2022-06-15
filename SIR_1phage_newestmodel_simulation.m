clear all;
close all;


A1_value=[100];
for i=1:length(A1_value)
    SIR_1phage_newestmodel_function(strcat(sprintf('Smesh_L=20_init=2p5E-3_beta=0p6_A1=%d.mat',A1_value(i))),A1_value(i));
end

