
function []=SIR_1phage_newestmodel_function(name,A1)
%% define parameters
dt=1;                  %time step in unit of seconds
dx=100;                  %spacial grids in unit of um
dy=100;
Time=60*60*24;         %simulated time;
T=round(Time/dt);       %time steps
LX=200000;             % range of plate
LY=200000;
LX0=LX/2; 
LY0=LY/2; 
Rb1=2000;
Rb2=1000;
LXg0=LX/2+10000; 
LYg0=LY/2; 
Rg1=2000; 
Rg2=1000;

K1=3.5;                 %receptor sensing range in unit of uM 
K2=1000;                %receptor sensing range in unit of uM
Sk=1;                   % attractant consumption half threshold in unit of uM
Nk=0.4;                 % nutrient consumption half threshold 

Yn=0.064;               

Dn=800;                 % diffusion of small molecules %A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry,Volume 166, Issue 2, 1 November 1987, Pages 335-341 
Da=800;                 % diffusion of small molecules %A method for the determination of diffusion coefficients for small molecules in aqueous solution, Analytical Biochemistry,Volume 166, Issue 2, 1 November 1987, Pages 335-341
G0=9*0.84/60;           % basal consumption of attractant 9uM/min for OD 1 cells.
Nutr_init=30;           % initial nutrient concentration in unit of mM
Attr_init=60;           % initial attractant concentration in unit of uM
bacter_init=0.2;
phage_init=2.5E-3;


Diff=60;
chi=410;    
lamnp0=log(2)/20/60;
kap0=0.9;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
P0=0.01;


% A1=8;
eta1=0.1;
beta1=0.6;
theta1=1E-4;


%% initial conditions
Cell_den_S=zeros(round(LX/dx),round(LY/dy));
Cell_den_I1=zeros(round(LX/dx),round(LY/dy));
Cell_den_R1=zeros(round(LX/dx),round(LY/dy));

Nutr=ones(round(LX/dx),round(LY/dy))*Nutr_init;
Attr=ones(round(LX/dx),round(LY/dy))*Attr_init;
Phag1=zeros(round(LX/dx),round(LY/dy));

tic

for x=1:round(LX/dx)
    for y=1:round(LY/dy)
        RR1=Rb1^2;
        RR2=Rb2^2;
        r2=((x-1)*dx-LX0)^2+((y-1)*dy-LY0)^2;
        if r2<=RR1
            bb(x,y)=bacter_init*exp(-r2/RR2);
        else
            bb(x,y)=0;
        end
    end
end    
Cell_den_S=bb;



for x=1:round(LX/dx)
    for y=1:round(LY/dy)
        RRg1=Rg1^2;
        RRg2=Rg2^2;
        r2=((x-1)*dx-LXg0)^2+((y-1)*dy-LYg0)^2;
        if r2<=RRg1
            pp(x,y)=phage_init*exp(-r2/RRg2);
        else
            pp(x,y)=0;
        end
    end
end  
Phag1=pp;


%%
Cell_den_S_temp_a=Cell_den_S;Cell_den_S_temp_b=Cell_den_S;Cell_den_S_temp_c=Cell_den_S;Cell_den_S_temp_d=Cell_den_S;
Cell_den_I1_temp_a=Cell_den_I1;Cell_den_I1_temp_b=Cell_den_I1;Cell_den_I1_temp_c=Cell_den_I1;Cell_den_I1_temp_d=Cell_den_I1;
Cell_den_R1_temp_a=Cell_den_R1;Cell_den_R1_temp_b=Cell_den_R1;Cell_den_R1_temp_c=Cell_den_R1;Cell_den_R1_temp_d=Cell_den_R1;

Attr_temp_a=Attr;Attr_temp_b=Attr;Attr_temp_c=Attr;Attr_temp_d=Attr;
Nutr_temp_a=Nutr;Nutr_temp_b=Nutr;Nutr_temp_c=Nutr;Nutr_temp_d=Nutr;
Phag1_temp_a=Phag1;Phag1_temp_b=Phag1;Phag1_temp_c=Phag1;Phag1_temp_d=Phag1;

g1_S=Cell_den_S;g2_S=Cell_den_S;
g1_I1=Cell_den_I1;g2_I1=Cell_den_I1;
g1_R1=Cell_den_R1;g2_R1=Cell_den_R1;

%%
for t=1:T
    Cell_den_S_temp=Cell_den_S;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;
    Cell_den_I1_temp=Cell_den_I1;
    Cell_den_I1_temp(Cell_den_I1_temp<=1E-10)=0;
    Cell_den_R1_temp=Cell_den_R1;
    Cell_den_R1_temp(Cell_den_R1_temp<=1E-10)=0; 
    Attr_temp=Attr;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag1_temp=Phag1; 
    Phag1_temp(Phag1_temp<=1E-10)=0;
   
    lambda_np_temp=lamnp0.*Nutr_temp./(Nutr_temp+Nk);
    kappa1_temp=kap0.*Phag1_temp./(Phag1_temp+P0);  
    Gamma_temp=G0./(1+Sk./Attr_temp); 
    S_growth_temp=(1-kappa1_temp.*Phag1_temp).*lambda_np_temp.*Cell_den_S_temp;
    I1_growth_temp=eta1*lambda_np_temp.*Cell_den_I1_temp+kappa1_temp.*Phag1_temp.*lambda_np_temp.*Cell_den_S_temp-theta1*Cell_den_I1_temp;
    R1_growth_temp=beta1*lambda_np_temp.*Cell_den_R1_temp+theta1*Cell_den_I1_temp; 
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp);
    phage1_product_temp=A1*((1-eta1)*lambda_np_temp.*Cell_den_I1_temp+(1-beta1)*lambda_np_temp.*Cell_den_R1_temp);

    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1);    
    g1_I1(2:end-1,:)=Cell_den_I1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I1(:,2:end-1)=Cell_den_I1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I1(1,:)=-g1_I1(2,:);g1_I1(end,:)=-g1_I1(end-1,:);
    g2_I1(:,1)=-g2_I1(:,2);g2_I1(:,end)=-g2_I1(:,end-1);      
    g1_R1(2:end-1,:)=Cell_den_R1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R1(:,2:end-1)=Cell_den_R1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R1(1,:)=-g1_R1(2,:);g1_R1(end,:)=-g1_R1(end-1,:);
    g2_R1(:,1)=-g2_R1(:,2);g2_R1(:,end)=-g2_R1(:,end-1); 
 
    Cell_den_S_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I1_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_I1_temp(1:end-2,2:end-1)+Cell_den_I1_temp(3:end,2:end-1)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I1_temp(2:end-1,1:end-2)+Cell_den_I1_temp(2:end-1,3:end)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dy^2)-(g1_I1(3:end,2:end-1)-g1_I1(1:end-2,2:end-1))/2/dx-(g2_I1(2:end-1,3:end)-g2_I1(2:end-1,1:end-2))/2/dy+I1_growth_temp(2:end-1,2:end-1);
    Cell_den_R1_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_R1_temp(1:end-2,2:end-1)+Cell_den_R1_temp(3:end,2:end-1)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R1_temp(2:end-1,1:end-2)+Cell_den_R1_temp(2:end-1,3:end)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dy^2)-(g1_R1(3:end,2:end-1)-g1_R1(1:end-2,2:end-1))/2/dx-(g2_R1(2:end-1,3:end)-g2_R1(2:end-1,1:end-2))/2/dy+R1_growth_temp(2:end-1,2:end-1);   
    Attr_temp_a(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_a(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);  
    Phag1_temp_a(2:end-1,2:end-1)=phage1_product_temp(2:end-1,2:end-1);

    Cell_den_S_temp_a(1,:)= Cell_den_S_temp_a(2,:);  Cell_den_S_temp_a(:,1)= Cell_den_S_temp_a(:,2); 
    Cell_den_S_temp_a(end,:)= Cell_den_S_temp_a(end-1,:);   Cell_den_S_temp_a(:,end)= Cell_den_S_temp_a(:,end-1);    
    Cell_den_I1_temp_a(1,:)= Cell_den_I1_temp_a(2,:);  Cell_den_I1_temp_a(:,1)= Cell_den_I1_temp_a(:,2); 
    Cell_den_I1_temp_a(end,:)= Cell_den_I1_temp_a(end-1,:);   Cell_den_I1_temp_a(:,end)= Cell_den_I1_temp_a(:,end-1);    
    Cell_den_R1_temp_a(1,:)= Cell_den_R1_temp_a(2,:);  Cell_den_R1_temp_a(:,1)= Cell_den_R1_temp_a(:,2); 
    Cell_den_R1_temp_a(end,:)= Cell_den_R1_temp_a(end-1,:);   Cell_den_R1_temp_a(:,end)= Cell_den_R1_temp_a(:,end-1);    
    Attr_temp_a(1,:)=Attr_temp_a(2,:); Attr_temp_a(:,1)=Attr_temp_a(:,2); 
    Attr_temp_a(end,:)=Attr_temp_a(end-1,:);  Attr_temp_a(:,end)=Attr_temp_a(:,end-1); 
    Nutr_temp_a(1,:)=Nutr_temp_a(2,:); Nutr_temp_a(:,1)=Nutr_temp_a(:,2); 
    Nutr_temp_a(end,:)=Nutr_temp_a(end-1,:);  Nutr_temp_a(:,end)=Nutr_temp_a(:,end-1);   
    Phag1_temp_a(1,:)=Phag1_temp_a(2,:); Phag1_temp_a(:,1)=Phag1_temp_a(:,2); 
    Phag1_temp_a(end,:)=Phag1_temp_a(end-1,:);  Phag1_temp_a(:,end)=Phag1_temp_a(:,end-1);
    %%  
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_a*dt/2;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;  
    Cell_den_I1_temp=Cell_den_I1+Cell_den_I1_temp_a*dt/2;
    Cell_den_I1_temp(Cell_den_I1_temp<=1E-10)=0;    
    Cell_den_R1_temp=Cell_den_R1+Cell_den_R1_temp_a*dt/2;
    Cell_den_R1_temp(Cell_den_R1_temp<=1E-10)=0;      
    Attr_temp=Attr+Attr_temp_a*dt/2;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_a*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;
    
    Phag1_temp=Phag1+Phag1_temp_a*dt/2; 
    Phag1_temp(Phag1_temp<=1E-10)=0;

    lambda_np_temp=lamnp0.*Nutr_temp./(Nutr_temp+Nk);
    kappa1_temp=kap0.*Phag1_temp./(Phag1_temp+P0);   
    Gamma_temp=G0./(1+Sk./Attr_temp);
        
    S_growth_temp=(1-kappa1_temp.*Phag1_temp).*lambda_np_temp.*Cell_den_S_temp;
    I1_growth_temp=eta1*lambda_np_temp.*Cell_den_I1_temp+kappa1_temp.*Phag1_temp.*lambda_np_temp.*Cell_den_S_temp-theta1*Cell_den_I1_temp;
    R1_growth_temp=beta1*lambda_np_temp.*Cell_den_R1_temp+theta1*Cell_den_I1_temp;   
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp);
    phage1_product_temp=A1*((1-eta1)*lambda_np_temp.*Cell_den_I1_temp+(1-beta1)*lambda_np_temp.*Cell_den_R1_temp);

    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I1(2:end-1,:)=Cell_den_I1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I1(:,2:end-1)=Cell_den_I1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I1(1,:)=-g1_I1(2,:);g1_I1(end,:)=-g1_I1(end-1,:);
    g2_I1(:,1)=-g2_I1(:,2);g2_I1(:,end)=-g2_I1(:,end-1); 
      
    g1_R1(2:end-1,:)=Cell_den_R1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R1(:,2:end-1)=Cell_den_R1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R1(1,:)=-g1_R1(2,:);g1_R1(end,:)=-g1_R1(end-1,:);
    g2_R1(:,1)=-g2_R1(:,2);g2_R1(:,end)=-g2_R1(:,end-1); 
 
    Cell_den_S_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I1_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_I1_temp(1:end-2,2:end-1)+Cell_den_I1_temp(3:end,2:end-1)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I1_temp(2:end-1,1:end-2)+Cell_den_I1_temp(2:end-1,3:end)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dy^2)-(g1_I1(3:end,2:end-1)-g1_I1(1:end-2,2:end-1))/2/dx-(g2_I1(2:end-1,3:end)-g2_I1(2:end-1,1:end-2))/2/dy+I1_growth_temp(2:end-1,2:end-1);
    Cell_den_R1_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_R1_temp(1:end-2,2:end-1)+Cell_den_R1_temp(3:end,2:end-1)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R1_temp(2:end-1,1:end-2)+Cell_den_R1_temp(2:end-1,3:end)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dy^2)-(g1_R1(3:end,2:end-1)-g1_R1(1:end-2,2:end-1))/2/dx-(g2_R1(2:end-1,3:end)-g2_R1(2:end-1,1:end-2))/2/dy+R1_growth_temp(2:end-1,2:end-1);   
    Attr_temp_b(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_b(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);  
    Phag1_temp_b(2:end-1,2:end-1)=phage1_product_temp(2:end-1,2:end-1); 
    
    Cell_den_S_temp_b(1,:)= Cell_den_S_temp_b(2,:);  Cell_den_S_temp_b(:,1)= Cell_den_S_temp_b(:,2); 
    Cell_den_S_temp_b(end,:)= Cell_den_S_temp_b(end-1,:);   Cell_den_S_temp_b(:,end)= Cell_den_S_temp_b(:,end-1); 
    Cell_den_I1_temp_b(1,:)= Cell_den_I1_temp_b(2,:);  Cell_den_I1_temp_b(:,1)= Cell_den_I1_temp_b(:,2); 
    Cell_den_I1_temp_b(end,:)= Cell_den_I1_temp_b(end-1,:);   Cell_den_I1_temp_b(:,end)= Cell_den_I1_temp_b(:,end-1); 
    Cell_den_R1_temp_b(1,:)= Cell_den_R1_temp_b(2,:);  Cell_den_R1_temp_b(:,1)= Cell_den_R1_temp_b(:,2); 
    Cell_den_R1_temp_b(end,:)= Cell_den_R1_temp_b(end-1,:);   Cell_den_R1_temp_b(:,end)= Cell_den_R1_temp_b(:,end-1); 
    Attr_temp_b(1,:)=Attr_temp_b(2,:); Attr_temp_b(:,1)=Attr_temp_b(:,2); 
    Attr_temp_b(end,:)=Attr_temp_b(end-1,:);  Attr_temp_b(:,end)=Attr_temp_b(:,end-1); 
    Nutr_temp_b(1,:)=Nutr_temp_b(2,:); Nutr_temp_b(:,1)=Nutr_temp_b(:,2); 
    Nutr_temp_b(end,:)=Nutr_temp_b(end-1,:);  Nutr_temp_b(:,end)=Nutr_temp_b(:,end-1); 
    Phag1_temp_b(1,:)=Phag1_temp_b(2,:); Phag1_temp_b(:,1)=Phag1_temp_b(:,2); 
    Phag1_temp_b(end,:)=Phag1_temp_b(end-1,:);  Phag1_temp_b(:,end)=Phag1_temp_b(:,end-1); 
    %%
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_b*dt/2;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;   
    Cell_den_I1_temp=Cell_den_I1+Cell_den_I1_temp_b*dt/2;
    Cell_den_I1_temp(Cell_den_I1_temp<=1E-10)=0;   
    Cell_den_R1_temp=Cell_den_R1+Cell_den_R1_temp_b*dt/2;
    Cell_den_R1_temp(Cell_den_R1_temp<=1E-10)=0;    
    Attr_temp=Attr+Attr_temp_b*dt/2;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_b*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;  
    Phag1_temp=Phag1+Phag1_temp_b*dt/2; 
    Phag1_temp(Phag1_temp<=1E-10)=0;
  
    lambda_np_temp=lamnp0.*Nutr_temp./(Nutr_temp+Nk);
    kappa1_temp=kap0.*Phag1_temp./(Phag1_temp+P0);   
    Gamma_temp=G0./(1+Sk./Attr_temp);  
    S_growth_temp=(1-kappa1_temp.*Phag1_temp).*lambda_np_temp.*Cell_den_S_temp;
    I1_growth_temp=eta1*lambda_np_temp.*Cell_den_I1_temp+kappa1_temp.*Phag1_temp.*lambda_np_temp.*Cell_den_S_temp-theta1*Cell_den_I1_temp;
    R1_growth_temp=beta1*lambda_np_temp.*Cell_den_R1_temp+theta1*Cell_den_I1_temp;
  
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp);
    phage1_product_temp=A1*((1-eta1)*lambda_np_temp.*Cell_den_I1_temp+(1-beta1)*lambda_np_temp.*Cell_den_R1_temp);

    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I1(2:end-1,:)=Cell_den_I1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I1(:,2:end-1)=Cell_den_I1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I1(1,:)=-g1_I1(2,:);g1_I1(end,:)=-g1_I1(end-1,:);
    g2_I1(:,1)=-g2_I1(:,2);g2_I1(:,end)=-g2_I1(:,end-1); 
      
    g1_R1(2:end-1,:)=Cell_den_R1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R1(:,2:end-1)=Cell_den_R1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R1(1,:)=-g1_R1(2,:);g1_R1(end,:)=-g1_R1(end-1,:);
    g2_R1(:,1)=-g2_R1(:,2);g2_R1(:,end)=-g2_R1(:,end-1); 
 
    Cell_den_S_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I1_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_I1_temp(1:end-2,2:end-1)+Cell_den_I1_temp(3:end,2:end-1)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I1_temp(2:end-1,1:end-2)+Cell_den_I1_temp(2:end-1,3:end)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dy^2)-(g1_I1(3:end,2:end-1)-g1_I1(1:end-2,2:end-1))/2/dx-(g2_I1(2:end-1,3:end)-g2_I1(2:end-1,1:end-2))/2/dy+I1_growth_temp(2:end-1,2:end-1);
    Cell_den_R1_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_R1_temp(1:end-2,2:end-1)+Cell_den_R1_temp(3:end,2:end-1)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R1_temp(2:end-1,1:end-2)+Cell_den_R1_temp(2:end-1,3:end)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dy^2)-(g1_R1(3:end,2:end-1)-g1_R1(1:end-2,2:end-1))/2/dx-(g2_R1(2:end-1,3:end)-g2_R1(2:end-1,1:end-2))/2/dy+R1_growth_temp(2:end-1,2:end-1);   
    Attr_temp_c(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_c(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);  
    Phag1_temp_c(2:end-1,2:end-1)=phage1_product_temp(2:end-1,2:end-1); 
        
    Cell_den_S_temp_c(1,:)= Cell_den_S_temp_c(2,:);  Cell_den_S_temp_c(:,1)= Cell_den_S_temp_c(:,2); 
    Cell_den_S_temp_c(end,:)= Cell_den_S_temp_c(end-1,:);   Cell_den_S_temp_c(:,end)= Cell_den_S_temp_c(:,end-1); 
    Cell_den_I1_temp_c(1,:)= Cell_den_I1_temp_c(2,:);  Cell_den_I1_temp_c(:,1)= Cell_den_I1_temp_c(:,2); 
    Cell_den_I1_temp_c(end,:)= Cell_den_I1_temp_c(end-1,:);   Cell_den_I1_temp_c(:,end)= Cell_den_I1_temp_c(:,end-1);   
    Cell_den_R1_temp_c(1,:)= Cell_den_R1_temp_c(2,:);  Cell_den_R1_temp_c(:,1)= Cell_den_R1_temp_c(:,2); 
    Cell_den_R1_temp_c(end,:)= Cell_den_R1_temp_c(end-1,:);   Cell_den_R1_temp_c(:,end)= Cell_den_R1_temp_c(:,end-1);    
    Attr_temp_c(1,:)=Attr_temp_c(2,:); Attr_temp_c(:,1)=Attr_temp_c(:,2); 
    Attr_temp_c(end,:)=Attr_temp_c(end-1,:);  Attr_temp_c(:,end)=Attr_temp_c(:,end-1); 
    Nutr_temp_c(1,:)=Nutr_temp_c(2,:); Nutr_temp_c(:,1)=Nutr_temp_c(:,2); 
    Nutr_temp_c(end,:)=Nutr_temp_c(end-1,:);  Nutr_temp_c(:,end)=Nutr_temp_c(:,end-1);   
    Phag1_temp_c(1,:)=Phag1_temp_c(2,:); Phag1_temp_c(:,1)=Phag1_temp_c(:,2); 
    Phag1_temp_c(end,:)=Phag1_temp_c(end-1,:);  Phag1_temp_c(:,end)=Phag1_temp_c(:,end-1);
    %%
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_c*dt;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;   
    Cell_den_I1_temp=Cell_den_I1+Cell_den_I1_temp_c*dt;
    Cell_den_I1_temp(Cell_den_I1_temp<=1E-10)=0; 
    Cell_den_R1_temp=Cell_den_R1+Cell_den_R1_temp_c*dt;
    Cell_den_R1_temp(Cell_den_R1_temp<=1E-10)=0; 
    Attr_temp=Attr+Attr_temp_c*dt;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_c*dt; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag1_temp=Phag1+Phag1_temp_c*dt; 
    Phag1_temp(Phag1_temp<=1E-10)=0;
  
    lambda_np_temp=lamnp0.*Nutr_temp./(Nutr_temp+Nk);
    kappa1_temp=kap0.*Phag1_temp./(Phag1_temp+P0); 
    Gamma_temp=G0./(1+Sk./Attr_temp);

    S_growth_temp=(1-kappa1_temp.*Phag1_temp).*lambda_np_temp.*Cell_den_S_temp;
    I1_growth_temp=eta1*lambda_np_temp.*Cell_den_I1_temp+kappa1_temp.*Phag1_temp.*lambda_np_temp.*Cell_den_S_temp-theta1*Cell_den_I1_temp;
    R1_growth_temp=beta1*lambda_np_temp.*Cell_den_R1_temp+theta1*Cell_den_I1_temp;

    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I1_temp+Cell_den_R1_temp);
    phage1_product_temp=A1*((1-eta1)*lambda_np_temp.*Cell_den_I1_temp+(1-beta1)*lambda_np_temp.*Cell_den_R1_temp);

    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I1(2:end-1,:)=Cell_den_I1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I1(:,2:end-1)=Cell_den_I1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I1(1,:)=-g1_I1(2,:);g1_I1(end,:)=-g1_I1(end-1,:);
    g2_I1(:,1)=-g2_I1(:,2);g2_I1(:,end)=-g2_I1(:,end-1); 
      
    g1_R1(2:end-1,:)=Cell_den_R1_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R1(:,2:end-1)=Cell_den_R1_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R1(1,:)=-g1_R1(2,:);g1_R1(end,:)=-g1_R1(end-1,:);
    g2_R1(:,1)=-g2_R1(:,2);g2_R1(:,end)=-g2_R1(:,end-1); 
 
    Cell_den_S_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I1_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_I1_temp(1:end-2,2:end-1)+Cell_den_I1_temp(3:end,2:end-1)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I1_temp(2:end-1,1:end-2)+Cell_den_I1_temp(2:end-1,3:end)-2*Cell_den_I1_temp(2:end-1,2:end-1))/dy^2)-(g1_I1(3:end,2:end-1)-g1_I1(1:end-2,2:end-1))/2/dx-(g2_I1(2:end-1,3:end)-g2_I1(2:end-1,1:end-2))/2/dy+I1_growth_temp(2:end-1,2:end-1);
    Cell_den_R1_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_R1_temp(1:end-2,2:end-1)+Cell_den_R1_temp(3:end,2:end-1)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R1_temp(2:end-1,1:end-2)+Cell_den_R1_temp(2:end-1,3:end)-2*Cell_den_R1_temp(2:end-1,2:end-1))/dy^2)-(g1_R1(3:end,2:end-1)-g1_R1(1:end-2,2:end-1))/2/dx-(g2_R1(2:end-1,3:end)-g2_R1(2:end-1,1:end-2))/2/dy+R1_growth_temp(2:end-1,2:end-1);   

    Attr_temp_d(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_d(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);  
    Phag1_temp_d(2:end-1,2:end-1)=phage1_product_temp(2:end-1,2:end-1);

    Cell_den_S_temp_d(1,:)= Cell_den_S_temp_d(2,:);  Cell_den_S_temp_d(:,1)= Cell_den_S_temp_d(:,2); 
    Cell_den_S_temp_d(end,:)= Cell_den_S_temp_d(end-1,:);   Cell_den_S_temp_d(:,end)= Cell_den_S_temp_d(:,end-1);   
    Cell_den_I1_temp_d(1,:)= Cell_den_I1_temp_d(2,:);  Cell_den_I1_temp_d(:,1)= Cell_den_I1_temp_d(:,2); 
    Cell_den_I1_temp_d(end,:)= Cell_den_I1_temp_d(end-1,:);   Cell_den_I1_temp_d(:,end)= Cell_den_I1_temp_d(:,end-1);    
    Cell_den_R1_temp_d(1,:)= Cell_den_R1_temp_d(2,:);  Cell_den_R1_temp_d(:,1)= Cell_den_R1_temp_d(:,2); 
    Cell_den_R1_temp_d(end,:)= Cell_den_R1_temp_d(end-1,:);   Cell_den_R1_temp_d(:,end)= Cell_den_R1_temp_d(:,end-1);  
    Attr_temp_d(1,:)=Attr_temp_d(2,:); Attr_temp_d(:,1)=Attr_temp_d(:,2); 
    Attr_temp_d(end,:)=Attr_temp_d(end-1,:);  Attr_temp_d(:,end)=Attr_temp_d(:,end-1); 
    Nutr_temp_d(1,:)=Nutr_temp_d(2,:); Nutr_temp_d(:,1)=Nutr_temp_d(:,2); 
    Nutr_temp_d(end,:)=Nutr_temp_d(end-1,:);  Nutr_temp_d(:,end)=Nutr_temp_d(:,end-1);   
    Phag1_temp_d(1,:)=Phag1_temp_d(2,:); Phag1_temp_d(:,1)=Phag1_temp_d(:,2); 
    Phag1_temp_d(end,:)=Phag1_temp_d(end-1,:);  Phag1_temp_d(:,end)=Phag1_temp_d(:,end-1);
    %%
    Cell_den_S=Cell_den_S+1/6*(Cell_den_S_temp_a+2*Cell_den_S_temp_b+2*Cell_den_S_temp_c+Cell_den_S_temp_d)*dt;   
    Cell_den_I1=Cell_den_I1+1/6*(Cell_den_I1_temp_a+2*Cell_den_I1_temp_b+2*Cell_den_I1_temp_c+Cell_den_I1_temp_d)*dt;   
    Cell_den_R1=Cell_den_R1+1/6*(Cell_den_R1_temp_a+2*Cell_den_R1_temp_b+2*Cell_den_R1_temp_c+Cell_den_R1_temp_d)*dt;  
    Attr=Attr+1/6*(Attr_temp_a+2*Attr_temp_b+2*Attr_temp_c+Attr_temp_d)*dt;
    Nutr=Nutr+1/6*(Nutr_temp_a+2*Nutr_temp_b+2*Nutr_temp_c+Nutr_temp_d)*dt;
    Phag1=Phag1+1/6*(Phag1_temp_a+2*Phag1_temp_b+2*Phag1_temp_c+Phag1_temp_d)*dt;
   
    Cell_den_S(1,:)=Cell_den_S(2,:); Cell_den_S(:,1)=Cell_den_S(:,2); 
    Cell_den_S(end,:)=Cell_den_S(end-1,:);  Cell_den_S(:,end)=Cell_den_S(:,end-1);   
    Cell_den_S(Cell_den_S<=1E-10)=0;   
    Cell_den_I1(1,:)=Cell_den_I1(2,:); Cell_den_I1(:,1)=Cell_den_I1(:,2); 
    Cell_den_I1(end,:)=Cell_den_I1(end-1,:);  Cell_den_I1(:,end)=Cell_den_I1(:,end-1);   
    Cell_den_I1(Cell_den_I1<=1E-10)=0;   
    Cell_den_R1(1,:)=Cell_den_R1(2,:); Cell_den_R1(:,1)=Cell_den_R1(:,2); 
    Cell_den_R1(end,:)=Cell_den_R1(end-1,:);  Cell_den_R1(:,end)=Cell_den_R1(:,end-1);   
    Cell_den_R1(Cell_den_R1<=1E-10)=0;    
    Attr(1,:)=Attr(2,:); Attr(:,1)=Attr(:,2); 
    Attr(end,:)=Attr(end-1,:);  Attr(:,end)=Attr(:,end-1);
    Attr(Attr<=0)=0;   
    Nutr(1,:)=Nutr(2,:); Nutr(:,1)=Nutr(:,2); 
    Nutr(end,:)=Nutr(end-1,:);  Nutr(:,end)=Nutr(:,end-1); 
    Nutr(Nutr<=0)=0;   
    Phag1(1,:)=Phag1(2,:); Phag1(:,1)=Phag1(:,2); 
    Phag1(end,:)=Phag1(end-1,:);  Phag1(:,end)=Phag1(:,end-1);  
    Phag1(Phag1<=1E-10)=0; 
       
    if t/3600==round(t/3600)
    save(strrep(name,'.mat',['_' num2str(t/3600) '.mat']));
    end   
    
end

toc

end