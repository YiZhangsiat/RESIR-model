clear all;
close all;



TimeT=2;

for m=1:length(TimeT)

filename=strcat('Smesh_L=20_init=2p5E-3_beta=0p6_A1=100','_',int2str(TimeT(m)));

load(strcat(filename,'.mat'));
infection_X_range=round(LXg0/dx):1:round(LX/dx);

Total_cell=Cell_den_S+Cell_den_I1+Cell_den_R1;
infected_cell=Cell_den_I1;
recovered_cell=Cell_den_R1; 
Phage_T=Phag1;

[Nx,Ny]=size(Total_cell);
Nx0=round(Nx/2); 
Ny0=round(Ny/2);
%% 
for i=1:length(infection_X_range)
      [FW_peak_value_1,FW_peak_loc_1]=findpeaks(Phage_T(infection_X_range(i),:)); 
      if ~isnan(FW_peak_loc_1)
       Phag_Maxloc_1=FW_peak_loc_1(1);       
       [~,Phag_Maxloc_1]  =min(abs(Phage_T(infection_X_range(i),1:Phag_Maxloc_1)-2*FW_peak_value_1(1)/3));
                
    left_loc(i)=Phag_Maxloc_1(1);
    left_position(i)=i*dx/1000;  
    T_left_line_S(i)=Cell_den_S(infection_X_range(i),Phag_Maxloc_1(1));
    T_left_line_infected_cell(i)=infected_cell(infection_X_range(i),Phag_Maxloc_1(1));
    T_left_line_recovered_cell(i)=recovered_cell(infection_X_range(i),Phag_Maxloc_1(1));
    T_left_line_phage(i)=Phag1(infection_X_range(i),Phag_Maxloc_1(1));
      else
    left_loc(i)=nan;
    left_position(i)=nan;
    T_left_line_S(i)=nan;
    T_left_line_infected_cell(i)=nan;
    T_left_line_recovered_cell(i)=nan;
    T_left_line_phage(i)=nan;
      end     
end
  
end
%%      
  left_line_S=smooth(T_left_line_S',5,'sgolay',0);
  left_line_infected_cell=smooth(T_left_line_infected_cell',5,'sgolay',0);
  left_line_recovered_cell=smooth(T_left_line_recovered_cell',5,'sgolay',0);
  left_line_phage=smooth(T_left_line_phage',5,'sgolay',0);
  
    resmooth_value=**;  
    left_line_S_V2=smooth(left_line_S(1:resmooth_value),20,'sgolay',0);
    left_line_S(1:resmooth_value)=left_line_S_V2;
    left_line_infected_cell_V2=smooth(left_line_infected_cell(1:resmooth_value),10,'sgolay',0);
    left_line_infected_cell(1:resmooth_value)=left_line_infected_cell_V2;
    left_line_recovered_cell_V2=smooth(left_line_recovered_cell(1:resmooth_value+10),15,'sgolay',0);
    left_line_recovered_cell(1:resmooth_value+10)=left_line_recovered_cell_V2;
  
    left_line_total_cell=left_line_S+left_line_infected_cell+left_line_recovered_cell;
  
    left_line_S(left_line_S<0)=0;
    left_line_infected_cell(left_line_infected_cell<0)=0;
    left_line_recovered_cell(left_line_recovered_cell<0)=0;
    left_line_total_cell(left_line_total_cell<0)=0;
    left_line_phage(left_line_phage<0)=0; 
%% 
    central_position=0:dx/1000:(Nx-Nx0)*dx/1000;
    central_line_S=Cell_den_S(Nx0:end,Ny0);
    central_line_infected_cell=infected_cell(Nx0:end,Ny0);
    central_line_recovered_cell=recovered_cell(Nx0:end,Ny0);
    central_line_phage=Phag1(Nx0:end,Ny0);
    central_line_total_cell=central_line_S+central_line_infected_cell+central_line_recovered_cell;
    
%% 
figure;plot(left_line_S,'linewidth',4);hold on;
plot(left_line_infected_cell,'linewidth',4);hold on;
plot(left_line_recovered_cell,'linewidth',4);hold off;
legend('cell S','infected cell ','recovered cell ');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
title('left line');
%%
figure;
plot(central_line_total_cell,'linewidth',4);hold on;
plot(central_line_S,'linewidth',4);hold on;
plot(central_line_infected_cell,'linewidth',4);hold on;
plot(central_line_recovered_cell,'linewidth',4);hold off;
legend('total cell','cell S','infected cell ','recovered cell ');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
title('central line');