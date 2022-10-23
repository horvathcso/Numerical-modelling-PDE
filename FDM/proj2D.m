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
%%%                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

k=0;
 while ( k <=0 )
 k = input('Time step (k): ');
 end


  BC=0;
  while (BC ~=1 &&  BC ~=2)
      disp(' ');
      disp('2 different types of boundary conditions:')
      disp('-------------------------------------------------------------------------');
      disp('(1) Dirichlet ')
      disp('(2) Neuman ')
      disp('-------------------------------------------------------------------------');
      disp(' ');
  BC = input('Type of boundry condition (1) (2) : ');
  end

   if BC == 2
        al = 0;
        bl = 0;
        cl = 1  ;
        ar = 0;
        br = 0;
        cr = 1;
   elseif BC==1
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


%k=0.005;


    max_itter=floor(t_1/k);

    % Construct SBP op
    if order == 2
        SBP2;
    elseif order == 4
        SBP4;                   
    elseif order == 6
        SBP6;
    end

    % Construct the projection operato   
    L=[bl*e_1'+cl*d_1;br*e_m'+cr*d_m];  % Boundary operator
    Im=eye(m);
    P=Im-HI*L'*((L*HI*L')\L);
   
    %Solution matrix 
    A=sparse(P*(c*c*D2)*P);
    Im=sparse(eye(m));
    D= sparse(kron(A,Im)+kron(Im,A));
    zer=sparse(zeros(m*m));
    Imm=sparse(eye(m*m));
    Q=sparse([zer,Imm;D,zer]);

            
    x=linspace(x_l,x_r,m)';	% discrete x values
    y=linspace(x_l,x_r,m)';
    %%% Initialize

    t=0.0;

   

    V=sparse(zeros(m*m,1));

    for i = 1:m
        v=sparse(zeros(m,1));  
        v=sparse(exp(-100*(x.^2+(x_l+(i-1)*h).^2)));
        V((i-1)*m+1:i*m) = sparse(v) ;
    end
    v=sparse(zeros(2*m*m,1));
    v(1:m*m)=V;
    V=v;

    [X, Y ] = meshgrid(x, y);
    U=reshape(V(1:m*m),m,m);
    surf(X,Y,U);
    title(['Numerical solution at t = ',num2str(t)]);
    xlabel('x');
ylabel('y');
colormap(jet(256))
shading interp
material dull
lighting phong
view(-70,20)
axis equal;
axis([-1 1 -1 1 -1 1]);
camlight left;
caxis([-0.1,0.1])
ax = gca;          % current axes
ax.FontSize = 16;

   
    figure(1);
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2 scrsz(4) scrsz(3)/2 scrsz(4)])
    vidObj = VideoWriter(filnamn);
    open(vidObj);

   n_step=10;

    %%%%%%%%%%%%% RK4 time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for nr_itter=1:max_itter

        w1=Q*V;
        w2=Q*(V+k/2*w1);
        w3=Q*(V+k/2*w2);
        w4=Q*(V+k*w3);
        
        V=V+k/6*(w1+2*w2+2*w3+w4);
        
        t=t+k;

       if mod(nr_itter,n_step)==0
           figure(1);   
           [X, Y ] = meshgrid(x, y);
                U=reshape(V(1:m*m),m,m);
                surf(X,Y,U);
                title(['Numerical solution at t = ',num2str(t)]);
               
xlabel('x');
ylabel('y');
colormap(jet(256))
shading interp
material dull
lighting phong
view(-70,20)
axis equal;
axis([-1 1 -1 1 -1 1]);
camlight left;
caxis([-0.1,0.1])
ax = gca;          % current axes
ax.FontSize = 16;
currFrame = getframe;
writeVideo(vidObj,currFrame);

       end
       
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %tiden=linspace(0,t_1,max_itter+1);
 close(vidObj);

figure(2);
[X, Y ] = meshgrid(x, y);
U=reshape(V(1:m*m),m,m);
surf(X,Y,U);
    