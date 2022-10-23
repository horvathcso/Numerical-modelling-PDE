function myproject()

a=0;
b=2*pi;
N=201;
h=(b-a)/(N-1);
c=2;
e=0.1;
T=2;
M=Ma(N,h);
A=Aa(N);
S=e*Sa(N,h);
M_=inv(M);
A_=M_*A;
S_=M_*S;

    %%%%%%%%%%%%% RK4 time integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x= a:h:b;
x=x';
t=0;
k=0.2*h^2;

i=1;
for e = [ 1, 0.1, 0.001, 0]
t=0;
M=Ma(N,h);
A=Aa(N);
S=e*Sa(N,h);
M_=inv(M);
A_=M_*A;
S_=M_*S;
u=sin(x);
    for nr_itter=1:(T/k)
        u1=func(u,A_,S_);
        u2=func(u+k/2*u1,A_,S_);
        u3=func(u+k/2*u2,A_,S_);
        u4=func(u+k*u3,A_,S_);
        u=(u+k/6*(u1+2*u2+2*u3+u4));
        u(1)=0;
        u(N)=0;
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
    figure(2)
    subplot(2,2,i)
    plot(x,u)
    title("eps= " +num2str(e) )
    i=i+1;
end
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