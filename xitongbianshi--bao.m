clc;
clear all;
A=[0.96,0.5,0.27,0.28;-0.125,0.96,-0.08,-0.07;0,0,0.85,0.97;0,0,0,0.99];
C=[0,2,0,0];
B=sqrt(15)*[0.5,0,0,1]';
I=[0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
for t=1:1:1000
v1(:,:,t)=B*B'*randn(4,4);
v2(t)=sqrt(1000)*randn;
end
x(:,1)=sqrt(0.5)*randn(4,1);
P(:,:,1)=0.5*eye(4);
P1(:,:,1)=0.5*eye(4);
x1(:,1)=[0,0,0,0]';
x2(:,1)=[0,0,0,0]';
x3(:,1)=[0,0,0,0]';
x4(:,1)=[0,0,0,0]';
%---------------------------------------%
sys=ss(A,[B B],C,0,-1,'inputname',{'u' 'w'},'outputname','y');
q=15;
r=1000;
[kalmf,L,P,K2]=kalman(sys,q,r);
%--------------------------------------%Steady-state Kalman ?lter (F1)
sys2=ss(A,B,C,0,-1);
yy4=lsim(sys2,v2);
a=A;
b=[B B B 0*B];
c=[C; C; C];
d=[0 0 0 0;0 0 0 0 ;0 0 0 1];
kalmf=kalmf(1,:);
P2=ss(a,b,c,d,-1,'inputname',{'u' 'w' 'v' 's'},'outputname',{'y' 'yv' 'yvx'});
sys2=parallel(P2,kalmf,1,1,[ ],[ ]);
SimModel=feedback(sys2,1,4,2,1);
SimModel=SimModel([1 3],[1 2 3]);
t1=[1:1000]';
u=sin(t1/5);
n=length(u);
w=sqrt(q)*rand(n,1);
v=sqrt(1)*randn(n,1);
[out,xx]=lsim(SimModel,[w,v,u]);
y3=out(:,1);
ye=out(:,2);
yv=y3+v2';
measerr=x(:,1);
measerr2=x(:,1);
%---------------------------------------%
for t=1:1:1000
    x(:,t+1)=A*x(:,t)+v1(:,t);
    y(t)=C*x(:,t)+v2(t);
%--------------------------------------%Dynamic Kalman 1-step predictor (K)&Dynamic Kalman fliter
    K(:,:,t)=A*P(:,:,t)*C'*pinv(C*P(:,:,t)*C'+v2(t));
    P(:,:,t+1)=A*P(:,:,t)*A'+B*B'-K(:,:,t)*(C*P(:,:,t)*C'+v2(t))*K(:,:,t)';
    y1(t)=C*x1(:,t);
    e(t)=y(t)-y1(t);
    x1(:,t+1)=A*x1(:,t)+K(:,:,t)*e(t);
    
%--------------------------------------%Steady-state Kalman 1-step predictor (K2)
    y2(t)=C*x2(:,t);
    e(t)=y(t)-y2(t);
    x2(:,t+1)=A*x2(:,t)+K2*e(t);
    % x3(:,t+1)=x3(:,t)+K3(:,:,t)*e(t);

%-------------------------------------(Optional) Dynamic Kalman 1-step predictor in predictor/corrector form (Kpc)
       K0(:,:,t)=P1(:,:,t)*C'*pinv(C*P1(:,:,t)*C'+v2(t));
       P0(:,:,t)=(I-K0(:,:,t)*C)*P1(:,:,t);
       P1(:,:,t+1)=A*P0(:,:,t)*A'+v1(:,:,t);
       y4(t)=C*x4(:,t);
       e2(t)=y(t)-y4(t);
       x4(:,t+1)=A*(x4(:,t)+K0(:,:,t)*e2(t));
       measerr=measerr+x(:,t+1)-x1(:,t+1);
       measerr2=measerr+x(:,t+1)-x2(:,t+1);
       measerr4=measerr+x(:,t+1)-x4(:,t+1);
end
%----------------------------------------RMSE-----------------------------------------------------------------------
s1=sum(measerr);
l1=length(measerr);
measerrcov=s1/sqrt(l1);

s2=sum(measerr2);
l2=length(measerr2);
meas2errcov=s2/sqrt(l2);

s4=sum(measerr4);
l4=length(measerr4);
meas4errcov=s4/sqrt(l4);
