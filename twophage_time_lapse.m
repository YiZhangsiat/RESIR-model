clear all;
close all;

TimeT=10*60;
for m=1:length(TimeT)
filename=strcat('Smesh_L=20_init=5E-1_ratio=1E-5_B=0p9_A1=40_A2=100','_',int2str(TimeT(m)));
load(strcat(filename,'.mat'));
infection_X_range=round(LXg0/dx):1:round(LX/dx);
Total_cell=Cell_den_S+Cell_den_I1+Cell_den_R1+Cell_den_I2+Cell_den_R2;
infected_cell_1=Cell_den_I1+Cell_den_R1;
infected_cell_2=Cell_den_I2+Cell_den_R2; 
Phage_T=Phag1+Phag2;
[Nx,Ny]=size(Total_cell);
Nx0=round(Nx/2); 
Ny0=round(Ny/2);
%% 
for i=1:length(infection_X_range)
      [~,FW_Cell_Infect_half_loc]=findpeaks(Phage_T(infection_X_range(i),:)); 
      if ~isnan(FW_Cell_Infect_half_loc)
        half_loc=FW_Cell_Infect_half_loc(1);       
    left_loc(i)=FW_Cell_Infect_half_loc(1);
    left_position(i)=sqrt((infection_X_range(i)-infection_X_range(1))^2+(half_loc-left_loc(1))^2)*dx/1000;  
    T_left_line_S(i)=sum(Cell_den_S(infection_X_range(i),half_loc-1:half_loc+1))/3;
    T_left_line_infected_cell_1(i)=sum(infected_cell_1(infection_X_range(i),half_loc-1:half_loc+1))/3;
    T_left_line_infected_cell_2(i)=sum(infected_cell_2(infection_X_range(i),half_loc-1:half_loc+1))/3;
    T_left_line_phage_1(i)=sum(Phag1(infection_X_range(i),half_loc-1:half_loc+1))/3;
    T_left_line_phage_2(i)=sum(Phag2(infection_X_range(i),half_loc-1:half_loc+1))/3;
      else
    left_loc(i)=nan;
    left_position(i)=nan;
    T_left_line_S(i)=nan;
    T_left_line_infected_cell_1(i)=nan;
    T_left_line_infected_cell_2(i)=nan;
    T_left_line_phage_1(i)=nan;
    T_left_line_phage_2(i)=nan;
      end     
end
 
end
%% 
    V_line_end=**;
    left_loc(V_line_end:end)=nan;
    T_left_line_S(V_line_end:end)=0;
    T_left_line_infected_cell_1(V_line_end:end)=0;
    T_left_line_infected_cell_2(V_line_end:end)=0;
    T_left_line_phage_1(V_line_end:end)=0;
    T_left_line_phage_2(V_line_end:end)=0;    
  left_line_S=smooth(T_left_line_S',2,'sgolay',0);
  left_line_infected_cell_1=smooth(T_left_line_infected_cell_1',2,'sgolay',0);
  left_line_infected_cell_2=smooth(T_left_line_infected_cell_2',2,'sgolay',0);
  left_line_phage_1=smooth(T_left_line_phage_1',2,'sgolay',0);
  left_line_phage_2=smooth(T_left_line_phage_2',2,'sgolay',0);
  
    resmooth_value=**;   
    left_line_infected_cell_1_V2=smooth(left_line_infected_cell_1(1:resmooth_value),10,'sgolay',0);
    left_line_infected_cell_1(1:resmooth_value)=left_line_infected_cell_1_V2;
    left_line_infected_cell_2_V2=smooth(left_line_infected_cell_2(1:resmooth_value),10,'sgolay',0);
    left_line_infected_cell_2(1:resmooth_value)=left_line_infected_cell_2_V2;
    left_line_S_V2=smooth(left_line_S(1:resmooth_value),10,'sgolay',0);
    left_line_S(1:resmooth_value)=left_line_S_V2;   
    left_line_total_cell=left_line_S+left_line_infected_cell_1+left_line_infected_cell_2;  
    left_line_S(left_line_S<0)=0;
    left_line_infected_cell_1(left_line_infected_cell_1<0)=0;
    left_line_infected_cell_2(left_line_infected_cell_2<0)=0;
    left_line_total_cell(left_line_total_cell<0)=0;
    left_line_phage_1(left_line_phage_1<0)=0;
    left_line_phage_2(left_line_phage_2<0)=0;   
%% 
    central_position=0:dx/1000:(Nx-Nx0)*dx/1000;
    central_line_S=Cell_den_S(Nx0:end,Ny0);
    central_line_infected_cell_1=infected_cell_1(Nx0:end,Ny0);
    central_line_infected_cell_2=infected_cell_2(Nx0:end,Ny0);
    central_line_phage_1=Phag1(Nx0:end,Ny0);
    central_line_phage_2=Phag2(Nx0:end,Ny0);
    central_line_total_cell=central_line_S+central_line_infected_cell_1+central_line_infected_cell_2;      
%% 
figure;plot(left_line_S,'linewidth',4);hold on;
plot(left_line_infected_cell_1,'linewidth',4);hold on;
plot(left_line_infected_cell_2,'linewidth',4);hold off;
legend('cell S','infected_cell 1','infected_cell 2');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
title('left line');
%% 
figure;
plot(central_line_S,'linewidth',4);hold on;
plot(central_line_infected_cell_1,'linewidth',4);hold on;
plot(central_line_infected_cell_2,'linewidth',4);hold off;
legend('cell S','infected_cell 1','infected_cell 2');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
title('central line');