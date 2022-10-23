function u=myproject(n)

a=-1;
b=1;
N=n;
%N=641;
h=(b-a)/(N-1);
c=2;
e=0.1;
T=0.4;
M=Ma(N,h);
A=Aa(N);
S=e*Sa(N,h);
M_=inv(M);
A_=M_*A;
S_=M_*S;

    %%%%%%%%%%%%% RK4 time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x= a:h:b;
x=x';
u= c-tanh((x+0.5)/(2*e));
t=0;
k=h^2;
%disp(A);

    for nr_itter=1:(T/k)
        u1=func(u,A_,S_);
        u2=func(u+k/2*u1,A_,S_);
        u3=func(u+k/2*u2,A_,S_);
        u4=func(u+k*u3,A_,S_);
        u=(u+k/6*(u1+2*u2+2*u3+u4));
        u(1)=u_ex(a,t,c,e);
        
        u(N)=u_ex(b,t,c,e);
        t=t+k;
   
            if mod(nr_itter,100)==0
            figure(1);  
            clf('reset');
            
            plot(x,u)
            %plot(x,u,x,u_ex(x,t,c,e));
           
            title(['Numerical solution at t = ',num2str(t)]);
            xlabel('x');ylabel('u');
            legend('v');
            end
    end
    disp(u);
end

function M=Ma(N,h)
a=2*h/3;
b=h/6;
c=h/6;
M = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
M(1,1)=h/3;
M(N,N)=h/3;
end

function A=Aa(N)
a=0;
b=0.5 ;
c=-0.5;
A = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
A(1,1)=-0.5;
A(N,N)=0.5;
end

function S=Sa(N,h)
a=2/h ;
b=-1/h ;
c=-1/h ;
S = diag(a*ones(1,N)) + diag(b*ones(1,N-1),1) + diag(c*ones(1,N-1),-1);
S(1,1)=1/h;
S(N,N)=1/h;
end

function u=func(u,A_,S_)
u = -A_*(u.^2*0.5)-S_*u;
end

function u=u_ex(x,t,c,e)
u=c-tanh((x+0.5-c*t)/(2*e));
end