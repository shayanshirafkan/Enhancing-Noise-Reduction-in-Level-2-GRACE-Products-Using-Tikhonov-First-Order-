clc
clear all
close all
%% read first data
name='ITSG_Lmax96_7512_DOS.ts'
[y,yh,xh,time,A,m,n ] = read( name );
Qy=eye(m,m);
ys1=y;
%% functional of regularization
L=zeros(n-1,n);
for i=1:n-1;
    L(i,i)=-1;
    L(i,i+1)=1;
end
[U,sm,X,V] = cgsvd(A,L);
%% VCE
Rx=L'*L;
Ry = eye(length(y));
s2_y =1; s2_x = 1;
[x1,s2_y,s2_x,lambda_VCE ] = VCE(y,A,Ry,Rx,s2_y,s2_x );
y1=A*x1;
%% GCV
[reg_min,G,reg_param] = gcv(U,sm,y,'Tikh');
alpha=reg_min;
x2=(A'*inv(Qy)*A+alpha*L'*L)\(A'*inv(Qy)*y);
y2=A*x2;
%% L-Curve
[reg_corner,rho,eta,reg_param] = l_curve(U,sm,y,'Tikh');
x3=(A'*inv(Qy)*A+reg_corner*L'*L)\(A'*inv(Qy)*y);
y3=A*x3;
%% Plot
figure;
plot(time+2003,y,'b','markersize',20,'linewidth',2)
grid on
hold on
plot(time+2003,y1,'r','markersize',20,'linewidth',2)
hold on
plot(time+2003,y2,'m','markersize',20,'linewidth',2)
hold on
plot(time+2003,y3,'k','markersize',20,'linewidth',2)
xlabel('Time [year]')
ylabel('Equivalent Water Height [meter]')
saveas(gcf,'im1.png')
