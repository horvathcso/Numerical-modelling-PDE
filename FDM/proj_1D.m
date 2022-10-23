%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                      %%%
%%% Advection diffusion eq.              %%%
%%%                                      %%%
%%% Solved with SBP-Projection           %%%
%%%                                      %%%
%%% Solution to part d) in the second    %%%
%%% problemsolving in the course         %%%
%%% Scientific Computing for PDE         %%%
%%%                                      %%%
%%% Date: Sep 14 2022                    %%%
%%% Author: Ken Mattsson                 %%%
%%%                                      %%%
%%%                                      %%%
%%%                                      %%%
%%%  u_t =au_x+bu_xx,  -1<= x <=1        %%%
%%%  u=0            , x=-1               %%%
%%%  au+2bu_x=0     , x=1                %%% 
%%%                                      %%%
%%%                                      %%%
%%% Initial data                         %%%
%%% u(x,0)=exp(-(6*x)^2)                 %%%
%%%                                      %%%
%%%                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


disp(' '); 
filnamn = input('The name of the avi-file: ','s');
disp(' '); 

order=0;
while (order ~=2 &&  order ~=4 && order ~=6)
order = input('Order of accuracy explicit (2) (4) (6) : ');
end

 m=0;
 while ( m <=15 )
    m = input('How many grid-intervals (>15) : ');
 end

 t_1=0;
 while ( t_1 <=0 )
 t_1 = input('End time (T): ');
 end

  BC=0;
  while (BC ~=1 &&  BC ~=2 && BC ~= 3 )
      disp(' ');
      disp('3 different types of boundary conditions:')
      disp('-------------------------------------------------------------------------');
      disp('(1) Dirichlet ')
      disp('(2) Neuman ')
      disp('(3) System settings. Right now the absorbing boundary conditions ')
      disp('-------------------------------------------------------------------------');
      disp(' ');
  BC = input('Type of boundry condition (1) (2) (3): ');
  end
   
    
    if BC == 3
        %Set BC
        al = 1;
        bl = 0;
        cl = -c;
        ar = 1;
        br = 0;
        cr = c;
    elseif BC == 2
        al = 0;
        bl = 0;
        cl = 1  ;
        ar = 0;
        br = 0;
        cr = 1;
    else
        al = 0;
        bl = 1;
        cl = 0;
        ar = 0;
        br = 1;
        cr = 0;
    end


    


%equation data
c = 1;
x_l=-1;x_r=1;               % The boundaries of the domain
h=(x_r-x_l)/(m-1);          % Grid step


k=0;
k=input('Give k value. If you give k=<=0, than we use k=0.1*h');
if k<=0
    k=0.1*h;
end


e1=[1,0]; e2=[0,1]; I2=eye(2); zer=zeros(m);

    max_itter=floor(t_1/k);

    % Construct SBP op
    if order == 2
        SBP2;
    elseif order == 4
        SBP4;                   
    elseif order == 6
        SBP6;
    end

    % Construct the projection operator
    H_I=kron(I2,HI);   
    L=[al*kron(e2,e_1')+bl*kron(e1,e_1')+cl*kron(e1,d_1);al*kron(e2,e_m')+br*kron(e1,e_m')+cr*kron(e1,d_m)];  % Boundary operator
    I2m=eye(2*m);
    Im=eye(m);
    P=I2m-H_I*L'*((L*H_I*L')\L);
    Q=[zer,Im;c*c*D2,zer];
    %Solution matrix 
    A=P*(Q)*P;
    
    V=zeros(m,1);           % Numerical solution

    x=linspace(x_l,x_r,m)';	% discrete x values
    
    %%% Initialize

    t=0.0;

    V=exp(-(x/0.2).^2);
    V=[V;zeros(m,1)];


    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)])
    vidObj = VideoWriter(filnamn);
    open(vidObj);


felet=zeros(max_itter+1,1);	      % Error vector in time
felet(1)=0;



   n_step=10;

    %%%%%%%%%%%%% RK4 time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for nr_itter=1:max_itter

        w1=A*V;
        w2=A*(V+k/2*w1);
        w3=A*(V+k/2*w2);
        w4=A*(V+k*w3);

        V=V+k/6*(w1+2*w2+2*w3+w4);

        t=t+k;


        %l2
        if BC == 2
            exact= 1/2*exp(-((x+2-t)/0.2).^2)+1/2*exp(-((x-2+t)/0.2).^2);
            felet(nr_itter+1)=sqrt(h)*norm(V(1:m)-exact);
        elseif BC==1
           exact= -1/2*exp(-((x+2-t)/0.2).^2)-1/2*exp(-((x-2+t)/0.2).^2);
            felet(nr_itter+1)=sqrt(h)*norm(V(1:m)-exact);
        end
       
       if mod(nr_itter,n_step)==0
            figure(2);  
            %plot(x,V(1:m),x,exact,'--','LineWidth',1);
            plot(x,V(1:m),'LineWidth',1);
            title(['Numerical solution at t = ',num2str(t)]);
            xlabel('x');ylabel('v');
            %legend('v', 'exact');
 
            ax = gca;          % current axes
            ax.FontSize = 16;
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
            
       end
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 tiden=linspace(0,t_1,max_itter+1);
 close(vidObj);

 if BC==1 || BC==2
disp(['The l_2-error is: ', num2str(felet(nr_itter+1))])
 end

figure(2);
plot(x,V(1:m),'LineWidth',1);
title(['Numerical solution at t = ',num2str(t_1)]);

xlabel('x');ylabel('v');
ax = gca;          % current axes
ax.FontSize = 16;
