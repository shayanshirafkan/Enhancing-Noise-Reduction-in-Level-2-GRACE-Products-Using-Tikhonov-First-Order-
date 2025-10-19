function [ y,yh,xh,time,A,m,n ] = read( name )
[filename pathname]=uigetfile(name);
fileID = fopen([pathname filename]);
C = textscan(fileID,'%f %f %f ');
fclose(fileID);
celldisp(C)
a=[C(1,1),C(1,2),C(1,3)];
a1=a{1, 1};  
a2=a{1, 2};
y=a{1, 3};
u=[];
[m,k]=size(a1);
time = juliandate(a1,a2,1);
time=time-time(1)+1;
time=time/365.25;
A=eye(m,m);
n=m;
xh=inv(A'*A)*A'*y;
yh=A*xh;
end

