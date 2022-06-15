clear all;
close all;

TimeT=1*60:10:24*60;
phage_ratio=1E-5;
for m=1:length(TimeT)

filename=strcat('Smesh_L=20_init=5E-1_ratio=1E-5_B=0p9_A1=40_A2=100','_',int2str(TimeT(m)));
load(strcat(filename,'.mat')); 
Total_cell=Cell_den_S+Cell_den_I1+Cell_den_R1+Cell_den_I2+Cell_den_R2;
infected_cell_1=Cell_den_I1+Cell_den_R1;
infected_cell_2=Cell_den_I2+Cell_den_R2; 
Phage_T=Phag1+Phag2;
[Nx,Ny]=size(Total_cell);
Nx0=round(Nx/2); 
Ny0=round(Ny/2);
for i=1:Ny0;
  [~,FW_max_P_loc]=max(Phage_T(:,i));
    if ~isnan(FW_max_P_loc)
    T_left_line_phage_1(i)=Phag1(FW_max_P_loc,i);
    T_left_line_phage_2(i)=Phag2(FW_max_P_loc,i);
    T_left_line_PS_foldchange_1(i)=Phag2(FW_max_P_loc,i)/(Phag1(FW_max_P_loc,i)+Phag2(FW_max_P_loc,i));
  else
      T_left_line_phage_1(i)=nan;
      T_left_line_phage_2(i)=nan;
      T_left_line_PS_foldchange_1(i)=nan;
    end
end
  T_left_line_PS_foldchange_1t=T_left_line_PS_foldchange_1(~isnan(T_left_line_PS_foldchange_1));
  Vline_PS_fold_change(m)=max(T_left_line_PS_foldchange_1t(1:20))/phage_ratio;
   %% 
   [FW_val,FW_loc]=findpeaks(-1*Total_cell(:,Ny0));    
    central_line_phage1_t1=Phag1(Nx0:end,Ny0);
    central_line_phage2_t1=Phag2(Nx0:end,Ny0);
    central_line_phage2_ratio_t1=central_line_phage2_t1./(central_line_phage1_t1+central_line_phage2_t1);
    central_line_phage2_ratio_t2=central_line_phage2_ratio_t1(~isnan(central_line_phage2_ratio_t1));
    central_line_phage2_ratio(m)=max(central_line_phage2_ratio_t2(end-20:end))/phage_ratio; 
end  
%% 
figure;plot(Vline_PS_fold_change,'linewidth',4);hold on;
plot(central_line_phage2_ratio,'linewidth',4);hold off;
legend('Vline','central');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
title('foldchange of PS');