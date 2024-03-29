clear all;
close all;
    TimeT=1*60:10:22*60;
%    TimeT=22*60;
for m=1:length(TimeT)
filename=strcat('Smesh_L=20_init=5E-1_ratio=1E-5_B=0p9_A1=40_A2=100','_',int2str(TimeT(m)));
load(strcat(filename,'.mat'));
%%
Total_cell=Cell_den_S+Cell_den_I1+Cell_den_R1+Cell_den_I2+Cell_den_R2;
Infected_cell_1=Cell_den_I1+Cell_den_R1;
Infected_cell_2=Cell_den_I2+Cell_den_R2;
[Nx,Ny]=size(Total_cell);
Nx0=round(Nx/2); 
Ny0=round(Ny/2);
[FW_Nx0_v,FW_Nx0_l]=findpeaks(Total_cell(Nx0,:));
FW_R_2=round((Ny0-FW_Nx0_l(1))/2);
RE_radius=2*FW_R_2;
RE_perimeter=pi*RE_radius;
[FW_R2_v,FW_R2_l]=findpeaks(Total_cell(Nx0+FW_R_2,:)); 
Y_start=FW_R2_l(1);
Y_end=FW_R2_l(end);
%% 
horiz_X=Nx0:Nx0+FW_R_2-2;
for i=1:length(horiz_X)
    [FW_horiz_X_v,FW_horiz_X_l]=findpeaks(Total_cell(horiz_X(i),:));      
     Horiz_X_position(i)=horiz_X(i);
     FirstMatrix_Y_position(i)=FW_horiz_X_l(1);
     LastMatrix_Y_position(i)=FW_horiz_X_l(end);
     FirstMatrix_Radial_length=sqrt((Horiz_X_position(i)-Nx0)^2+(FirstMatrix_Y_position(i)-Ny0)^2);
     LastMatrix_Radial_length=sqrt((Horiz_X_position(i)-Nx0)^2+(LastMatrix_Y_position(i)-Ny0)^2);   
     First_degree(i)=acos((Ny0-FirstMatrix_Y_position(i))/FirstMatrix_Radial_length);
     Last_degree(i)=acos((Ny0-LastMatrix_Y_position(i))/LastMatrix_Radial_length);
     First_cell_S(i)=Cell_den_S(Horiz_X_position(i),FirstMatrix_Y_position(i));
     First_Infected_cell_1(i)=Infected_cell_1(Horiz_X_position(i),FirstMatrix_Y_position(i));
     First_Infected_cell_2(i)=Infected_cell_2(Horiz_X_position(i),FirstMatrix_Y_position(i));
     First_Phage_1(i)=Phag1(Horiz_X_position(i),FirstMatrix_Y_position(i));
     First_Phage_2(i)=Phag2(Horiz_X_position(i),FirstMatrix_Y_position(i));    
     Last_cell_S(i)=Cell_den_S(Horiz_X_position(i),LastMatrix_Y_position(i));
     Last_Infected_cell_1(i)=Infected_cell_1(Horiz_X_position(i),LastMatrix_Y_position(i));
     Last_Infected_cell_2(i)=Infected_cell_2(Horiz_X_position(i),LastMatrix_Y_position(i));
     Last_Phage_1(i)=Phag1(Horiz_X_position(i),LastMatrix_Y_position(i));
     Last_Phage_2(i)=Phag2(Horiz_X_position(i),LastMatrix_Y_position(i)); 
end
%% 
Vertical_Y=Y_start:Y_end;
for j=1:length(Vertical_Y)  
[FW_Vertial_v,FW_Vertical_l]=findpeaks(Total_cell(Nx0:end,Vertical_Y(j)));  
     Vertial_X_position(j)=FW_Vertical_l(end)+Nx0-1;
     Vertial_Y_position(j)=Vertical_Y(j);  
     MiddleMatrix_Radial_length=sqrt((Vertial_X_position(j)-Nx0)^2+(Vertial_Y_position(j)-Ny0)^2);
     Middle_degree(j)=acos((Ny0-Vertial_Y_position(j))/MiddleMatrix_Radial_length);  
     Middle_cell_S(j)=Cell_den_S(Vertial_X_position(j),Vertial_Y_position(j));
     Middle_Infected_cell_1(j)=Infected_cell_1(Vertial_X_position(j),Vertial_Y_position(j));
     Middle_Infected_cell_2(j)=Infected_cell_2(Vertial_X_position(j),Vertial_Y_position(j));
     Middle_Phage_1(j)=Phag1(Vertial_X_position(j),Vertial_Y_position(j));
     Middle_Phage_2(j)=Phag2(Vertial_X_position(j),Vertial_Y_position(j));
end
Cartesian_X_position=[Horiz_X_position,Vertial_X_position,flip(Horiz_X_position)]';
Cartesian_Y_position=[FirstMatrix_Y_position,Vertial_Y_position,flip(LastMatrix_Y_position)]';
Polar_degree_t=[First_degree,Middle_degree,flip(Last_degree)]';
FW_cell_ST=[First_cell_S,Middle_cell_S,flip(Last_cell_S)]';
FW_Infected_cell_1T=[First_Infected_cell_1,Middle_Infected_cell_1,flip(Last_Infected_cell_1)]';
FW_Infected_cell_2T=[First_Infected_cell_2,Middle_Infected_cell_2,flip(Last_Infected_cell_2)]';
FW_Phage_1T=[First_Phage_1,Middle_Phage_1,flip(Last_Phage_1)]';
FW_Phage_2T=[First_Phage_2,Middle_Phage_2,flip(Last_Phage_2)]';
Polar_degree=Polar_degree_t(1:2:end);
FW_cell_S=FW_cell_ST(1:2:end);
FW_Infected_cell_1=FW_Infected_cell_1T(1:2:end);
FW_Infected_cell_2=FW_Infected_cell_2T(1:2:end);
FW_Phage_1=FW_Phage_1T(1:2:end);
FW_Phage_2=FW_Phage_2T(1:2:end);
Polar_degree_interp = 0:pi/900:pi;
interp1_FW_cell_S=smooth(interp1(Polar_degree,FW_cell_S,Polar_degree_interp),20);     
interp1_Infected_cell_1=smooth(interp1(Polar_degree,FW_Infected_cell_1,Polar_degree_interp),20);    
interp1_Infected_cell_2=smooth(interp1(Polar_degree,FW_Infected_cell_2,Polar_degree_interp),20);    
interp1_Phage_1=smooth(interp1(Polar_degree,FW_Phage_1,Polar_degree_interp),20);
interp1_Phage_2=smooth(interp1(Polar_degree,FW_Phage_2,Polar_degree_interp),20); 
FW_cell_S_smooth=smooth(interp1_FW_cell_S,20,'sgolay',0);
FW_infected_cell_1_smooth=smooth(interp1_Infected_cell_1,20,'sgolay',0);
FW_infected_cell_2_smooth=smooth(interp1_Infected_cell_2,20,'sgolay',0);
FW_phage_1_smooth=smooth(interp1_Phage_1,20,'sgolay',0);
FW_phage_2_smooth=smooth(interp1_Phage_2,20,'sgolay',0);
expandingfront_cell_S_data(:,m)=FW_cell_S_smooth;
expandingfront_infected_1_data(:,m)=FW_infected_cell_1_smooth;
expandingfront_infected_2_data(:,m)=FW_infected_cell_2_smooth;
expandingfront_phage_1_data(:,m)=FW_phage_1_smooth;
expandingfront_phage_2_data(:,m)=FW_phage_2_smooth;
end
 %% 
[FW_Nx,FW_Ny]=size(expandingfront_cell_S_data);
for i=1:FW_Nx    
    FW_cell_1_data_smooth(i,:)=smooth(expandingfront_infected_1_data(i,:),10,'sgolay',0);
    FW_cell_2_data_smooth(i,:)=smooth(expandingfront_infected_2_data(i,:),10,'sgolay',0);
    FW_phage_1_data_smooth(i,:)=smooth(expandingfront_phage_1_data(i,:),10,'sgolay',0);
    FW_phage_2_data_smooth(i,:)=smooth(expandingfront_phage_2_data(i,:),10,'sgolay',0);
end
FW_cell_1_data_smooth(FW_cell_1_data_smooth<1E-5)=nan;
FW_cell_2_data_smooth(FW_cell_2_data_smooth<1E-5)=nan;
FW_phage_1_data_smooth(FW_phage_1_data_smooth<1E-5)=nan;
FW_phage_2_data_smooth(FW_phage_2_data_smooth<1E-5)=nan;