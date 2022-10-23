clear();
i=1;
h=[];
e=[];
figure(2)

for N = [ 41, 81, 161, 321, 641]
    u=myproject(N);
    x=-1:(2/(N-1)):1;
    ex=2-tanh((x+0.5-2*0.4)/(2*0.1));
    err=sum(ex-u,'all')*2/(N-1);
    h(i)=(2/(N-1));
    e(i)=err;
    figure(2);
    subplot(2,3,i);
    plot(x,u )
    disp(err)
    title(" N= "+num2str(N)+", error: "+num2str(err))
    i = i+1;
end

x=-1:0.0001:1;
u=2-tanh((x+0.5-2*0.4)/(2*0.1));
subplot(2,3,i);
plot(x,u )
title("Exact solution")

figure(3);
plot(h,e)
title("Error as function of mesh-size")