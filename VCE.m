function [ x,s2_y,s2_x,lambda_VCE ] = VCE( y,A,Ry,Rx,s2_y,s2_x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[m,n]=size(A);
lambda_VCE = s2_y / s2_x;
tl=0;
tll=0;
P=inv(Ry);
 while norm(tl-s2_y) && norm(tll-s2_x) >5*10^-3
Ny=(1/s2_y)*(A'*P*A);
Nx=(1/s2_x)*Rx;
N=Nx+Ny;
x = (1/s2_y)*inv(N)*A'*P*y;
tao_y = trace(Ny * inv(N));
tao_x = trace(Nx * inv(N));
e=y - A*x;
tl=s2_y;
tll=s2_x;
% [tl,tll]
s2_y =(e'*P*e)/(m - tao_y)
s2_x =(x'*Rx*x)/(n - tao_x);
lambda_VCE = s2_y / s2_x;
 end
end

