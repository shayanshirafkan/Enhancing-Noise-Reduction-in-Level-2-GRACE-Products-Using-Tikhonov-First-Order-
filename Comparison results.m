clc
clear all
close all
%% read first data
name='ITSG_Lmax60_7512_DOS.ts'
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
%% GCV
[reg_min,G,reg_param] = gcv(U,sm,y,'Tikh');
alpha=reg_min;
x1=(A'*inv(Qy)*A+alpha*L'*L)\(A'*inv(Qy)*y);
y11=A*x1;
%% read second data
name='ITSG_Lmax96_7512_DOS.ts'
[y,yh,xh,time,A,m,n ] = read( name );
Qy=eye(m,m);
ys2=y;
%% functional of regularization
L=zeros(n-1,n);
for i=1:n-1;
    L(i,i)=-1;
    L(i,i+1)=1;
end
[U,sm,X,V] = cgsvd(A,L);
%% GCV
[reg_min,G,reg_param] = gcv(U,sm,y,'Tikh');
alpha=reg_min;
x2=(A'*inv(Qy)*A+alpha*L'*L)\(A'*inv(Qy)*y);
y22=A*x2;
%% Plot differnce of first and its reg-data
figure;
plot(time+2003,ys1-y11,'b','markersize',20,'linewidth',2)
grid on
xlabel('Time [year]')
ylabel('Difference [meter]')
% saveas(gcf,'im3.png')
%% Plot differnce of second and its reg-data
figure;
plot(time+2003,ys2-y22,'b','markersize',20,'linewidth',2)
grid on
xlabel('Time [year]')
ylabel('Difference [meter]')
% saveas(gcf,'im4.png')

%% Plot differnce of first and second reg-data 
figure;
plot(time+2003,y11,'b','markersize',20,'linewidth',2)
grid on
hold on
plot(time+2003,y22,'r','markersize',20,'linewidth',2)
xlabel('Time [year]')
ylabel('Equivalent Water Height [meter]')
% saveas(gcf,'im5.png')


