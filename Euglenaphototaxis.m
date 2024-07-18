close all;
clc;
clear;

rng('shuffle');
count=0;

ti = 0.00;      % start time
tf = 10;        % end time
dt = 0.001;     % time step
tn=ti:dt:tf;

% Kepsilon=0.00;
alpha=3*pi/2;   % alpha


% model=1:light-dependent angular turn; 
% model=2:delay in photoresponse; 
% model 3:photoresponse inversion
model=1; 

omega=2*pi;
U=1;

Ix=-0.05;       % light signal
Iy=0;
Iz=0;

K=0;

Ka=0*omega;
Kb=0;
Kc=0.0*omega;
Kd=0.2*omega;

Aa0=-pi/2;
Ab=0;
Ac=0;
Ad=0;

gamma=0;

theta0=2*pi*rand;
phi0=-pi/2+pi*rand;
psi0=2*pi*rand;


px0=cos(theta0)*cos(phi0);
py0=sin(theta0)*cos(phi0);
pz0=sin(phi0);
pvec0=[px0;py0;pz0];

e_angle1=0*pi/180;
e_angle2=0*pi/180;  
eye_angle1=0*pi/180;    
eye_angle2=0*pi/180;


bevec=[sin(e_angle1);-cos(e_angle1)*sin(e_angle2);cos(e_angle1)*cos(e_angle2)];
beye=-[sin(eye_angle1);-cos(eye_angle1)*sin(eye_angle2);cos(eye_angle1)*cos(eye_angle2)];
bhead=[1;0;0];

bpvec=[1;0;0];
bdvec=[0;1;0];

% mapping the body frame to initial frame
if (theta0==0 && phi0==0) 
    pvec1=bpvec;
    evec1=bevec;
    dvec1=bdvec;
    eye1=beye;
    head1=bhead;
else
    cvec=[0;-pz0;py0];
    cvec=cvec/sqrt(sum(cvec.^2));
    angle=acos(sum(pvec0.*bpvec));
    angle=-angle;
    
    qangle=[sin(angle/2)*cvec(1); sin(angle/2)*cvec(2); sin(angle/2)*cvec(3); cos(angle/2)];
    Qangle=[qangle(1)^2-qangle(2)^2-qangle(3)^2+qangle(4)^2 2*(qangle(1)*qangle(2)+qangle(3)*qangle(4)) 2*(qangle(1)*qangle(3)-qangle(2)*qangle(4));...
        2*(qangle(1)*qangle(2)-qangle(3)*qangle(4)) -qangle(1)^2+qangle(2)^2-qangle(3)^2+qangle(4)^2 2*(qangle(2)*qangle(3)+qangle(1)*qangle(4));...
        2*(qangle(1)*qangle(3)+qangle(2)*qangle(4)) 2*(qangle(2)*qangle(3)-qangle(1)*qangle(4)) -qangle(1)^2-qangle(2)^2+qangle(3)^2+qangle(4)^2];

    pvec1=Qangle*bpvec;
    evec1=Qangle*bevec;
    dvec1=Qangle*bdvec;
    eye1=Qangle*beye;
    head1=Qangle*bhead;
end

qpsi0=[sin(-psi0/2)*px0; sin(-psi0/2)*py0; sin(-psi0/2)*pz0; cos(-psi0/2)];
Qpsi0=[qpsi0(1)^2-qpsi0(2)^2-qpsi0(3)^2+qpsi0(4)^2 2*(qpsi0(1)*qpsi0(2)+qpsi0(3)*qpsi0(4)) 2*(qpsi0(1)*qpsi0(3)-qpsi0(2)*qpsi0(4));...
    2*(qpsi0(1)*qpsi0(2)-qpsi0(3)*qpsi0(4)) -qpsi0(1)^2+qpsi0(2)^2-qpsi0(3)^2+qpsi0(4)^2 2*(qpsi0(2)*qpsi0(3)+qpsi0(1)*qpsi0(4));...
    2*(qpsi0(1)*qpsi0(3)+qpsi0(2)*qpsi0(4)) 2*(qpsi0(2)*qpsi0(3)-qpsi0(1)*qpsi0(4)) -qpsi0(1)^2-qpsi0(2)^2+qpsi0(3)^2+qpsi0(4)^2];

evec0=Qpsi0*evec1;
dvec0=Qpsi0*dvec1;
eye0=Qpsi0*eye1;
head0=Qpsi0*head1;

head=head0;
pvec=pvec0;
evec=evec0;
dvec=dvec0;
eye=eye0;
alpha0=0;
qalpha0=[sin(alpha0/2)*head(1); sin(alpha0/2)*head(2); sin(alpha0/2)*head(3); cos(alpha0/2)];

Qalpha0=[qalpha0(1)^2-qalpha0(2)^2-qalpha0(3)^2+qalpha0(4)^2 2*(qalpha0(1)*qalpha0(2)+qalpha0(3)*qalpha0(4)) 2*(qalpha0(1)*qalpha0(3)-qalpha0(2)*qalpha0(4));...
            2*(qalpha0(1)*qalpha0(2)-qalpha0(3)*qalpha0(4)) -qalpha0(1)^2+qalpha0(2)^2-qalpha0(3)^2+qalpha0(4)^2 2*(qalpha0(2)*qalpha0(3)+qalpha0(1)*qalpha0(4));...
            2*(qalpha0(1)*qalpha0(3)+qalpha0(2)*qalpha0(4)) 2*(qalpha0(2)*qalpha0(3)-qalpha0(1)*qalpha0(4)) -qalpha0(1)^2-qalpha0(2)^2+qalpha0(3)^2+qalpha0(4)^2];
T0=-Qalpha0*evec;
xi=0;
x=[xi; 0; 0];
xn=zeros(3,length(tn));
pn=zeros(3,length(tn));
en=zeros(3,length(tn));
eyen=zeros(3,length(tn));
dn=zeros(3,length(tn));
dthetan=zeros(1,length(tn)-1);

Kn=zeros(3,length(tn));
Kvalue=zeros(1,length(tn)-1);


I1matrix=zeros(1,length(tn)-1);
Smatrix=zeros(1,length(tn)-1);

xn(:,1)=x;
pn(:,1)=pvec0;
en(:,1)=evec0;
eyen(:,1)=eye0;
dn(:,1)=dvec0;
head=head0;

dpsi=-omega*dt;
dpsi0=dpsi;

% diffusion coefficient
kb=1.380649*10^(-23); %m^2 kg/s^2.K
T=298;
eta=8.90*10^(-4); % kg/m.s
a=50*10^(-6); % m
Dt=kb*T/(6*pi*eta*a);
Dr=kb*T/(8*pi*eta*a^3);   % small value compared to flagella rotational noise, can be neglected
Dt=Dt/a^2; % to scale length of Euglena to 1

Dt=2*10^(-6)/50^2;
Dr=1.5*10^(-6);

Iepsilon=0.025; % light signal noise level
epsilon=0.02;   % rotational noise due to active reorientation


%%

inttime=0;

delay=480;      % delay time step for model 2
for iter=1:length(tn)-1
    x=x+U*pvec*dt+randn*sqrt(2*Dt*dt);

    I1noisevec=Iepsilon*randn(1,3);  
    I1=-((Ix+I1noisevec(1))*evec(1)+(Iy+I1noisevec(2))*evec(2)+(Iz+I1noisevec(3))*evec(3));
    S=I1;
    Smatrix(iter)=S;

    if (Ix==0 && Iy==0 && Iz==0)
            Ic=0;      
    else
        if Ix~=0
            Ic=Ix;
        end
        if Iy~=0
            Ic=Iy;
        end 
        if Iz~=0
            Ic=Iz;
        end
        Aa=Aa0;
    end
    if model==1
        I1_r=Ka*abs(Ic)+Kd*I1*heaviside(I1);
        I1matrix(:,iter)=I1_r; 
    end
%%%%%%%%%%%%%%%        
	if model==2 
        I1_r=Ka*abs(Ic)+Kd*I1*heaviside(I1);
        I1matrix(:,iter)=I1_r;

        if iter<delay+1
            I1_r=0;
        else
            I1_r=I1matrix(iter-delay);
        end
    end
%%%%%%%%%%%%%%   
	if model==3
         kick=0.1*rand;     % random switching in beat
         I1_r=Ka*abs(Ic)+Kd*kick*heaviside(-I1);
         I1matrix(:,iter)=I1_r;
    end
%%%%%%%%%%%%%%   
    if Ic==0
        I1_r=0;
    end
    K=I1_r+Kc;
    
    
    dtheta=K*dt;
    qalpha=[sin(alpha/2)*head(1); sin(alpha/2)*head(2); sin(alpha/2)*head(3); cos(alpha/2)];

        Qalpha=[qalpha(1)^2-qalpha(2)^2-qalpha(3)^2+qalpha(4)^2 2*(qalpha(1)*qalpha(2)+qalpha(3)*qalpha(4)) 2*(qalpha(1)*qalpha(3)-qalpha(2)*qalpha(4));...
            2*(qalpha(1)*qalpha(2)-qalpha(3)*qalpha(4)) -qalpha(1)^2+qalpha(2)^2-qalpha(3)^2+qalpha(4)^2 2*(qalpha(2)*qalpha(3)+qalpha(1)*qalpha(4));...
            2*(qalpha(1)*qalpha(3)+qalpha(2)*qalpha(4)) 2*(qalpha(2)*qalpha(3)-qalpha(1)*qalpha(4)) -qalpha(1)^2-qalpha(2)^2+qalpha(3)^2+qalpha(4)^2];
        T=-Qalpha*evec;

    qt=[T(1)*sin(-dtheta/2); T(2)*sin(-dtheta/2); T(3)*sin(-dtheta/2); cos(-dtheta/2)];
    Qt=[qt(1)^2-qt(2)^2-qt(3)^2+qt(4)^2 2*(qt(1)*qt(2)+qt(3)*qt(4)) 2*(qt(1)*qt(3)-qt(2)*qt(4));...
        2*(qt(1)*qt(2)-qt(3)*qt(4)) -qt(1)^2+qt(2)^2-qt(3)^2+qt(4)^2 2*(qt(2)*qt(3)+qt(1)*qt(4));...
        2*(qt(1)*qt(3)+qt(2)*qt(4)) 2*(qt(2)*qt(3)-qt(1)*qt(4)) -qt(1)^2-qt(2)^2+qt(3)^2+qt(4)^2];
    
    pvec=Qt*pvec;
    evec=Qt*evec;
    eye=Qt*eye;
    head=Qt*head;
    dvec=Qt*dvec;
    
    if epsilon>0
        depsilon_1=randn*sqrt(2*epsilon*dt)+randn*sqrt(2*Dr*dt);
        depsilon_2=randn*sqrt(2*epsilon*dt)+randn*sqrt(2*Dr*dt);
        depsilon_3=randn*sqrt(2*epsilon*dt)+randn*sqrt(2*Dr*dt);
    
        depsilon=sqrt(depsilon_1^2+depsilon_2^2+depsilon_3^2);
        Ep=[depsilon_1,depsilon_2,depsilon_3]/depsilon;
    else
        Ep=[1;0;0];
        depsilon=0;
    end
    
    qn=[Ep(1)*sin(-depsilon/2); Ep(2)*sin(-depsilon/2); Ep(3)*sin(-depsilon/2); cos(-depsilon/2)];
    Qn=[qn(1)^2-qn(2)^2-qn(3)^2+qn(4)^2 2*(qn(1)*qn(2)+qn(3)*qn(4)) 2*(qn(1)*qn(3)-qn(2)*qn(4));...
        2*(qn(1)*qn(2)-qn(3)*qn(4)) -qn(1)^2+qn(2)^2-qn(3)^2+qn(4)^2 2*(qn(2)*qn(3)+qn(1)*qn(4));...
        2*(qn(1)*qn(3)+qn(2)*qn(4)) 2*(qn(2)*qn(3)-qn(1)*qn(4)) -qn(1)^2-qn(2)^2+qn(3)^2+qn(4)^2];
    
    
    pvec=Qn*pvec;
    evec=Qn*evec;
    eye=Qn*eye;
    head=Qn*head;
    dvec=Qn*dvec;
    
    
    qp=[pvec(1)*sin(-dpsi/2);pvec(2)*sin(-dpsi/2);pvec(3)*sin(-dpsi/2);cos(-dpsi/2)];
    Qp=[qp(1)^2-qp(2)^2-qp(3)^2+qp(4)^2 2*(qp(1)*qp(2)+qp(3)*qp(4)) 2*(qp(1)*qp(3)-qp(2)*qp(4));...
        2*(qp(1)*qp(2)-qp(3)*qp(4)) -qp(1)^2+qp(2)^2-qp(3)^2+qp(4)^2 2*(qp(2)*qp(3)+qp(1)*qp(4));...
        2*(qp(1)*qp(3)+qp(2)*qp(4)) 2*(qp(2)*qp(3)-qp(1)*qp(4)) -qp(1)^2-qp(2)^2+qp(3)^2+qp(4)^2];
    
    pvec=Qp*pvec;
    evec=Qp*evec;
    eye=Qp*eye;
    head=Qp*head;
    dvec=Qp*dvec;
    

    xn(:,iter+1)=x;
    pn(:,iter+1)=pvec;
    en(:,iter+1)=evec;
%    eyen(:,iter+1)=eye;
%     Kn(:,iter+1)=T;
%     Kvalue(iter)=K;
%     dn(:,iter+1)=dvec;
%     dthetan(:,iter)=(dtheta)/dt;

end

%%
framenovec=[1 2001 4001 6001 8001 10001];

bodylength=1;
thickness=bodylength/5;

for k=1:length(framenovec)
set(gcf,'color','white');
hold on;
box on;
frameno=framenovec(k);
xc=xn(1,frameno);
yc=xn(2,frameno);
zc=xn(3,frameno);


% 
px_rot=pn(1,frameno);
py_rot=pn(2,frameno);
pz_rot=pn(3,frameno);
bodyrot=[0;-pz_rot;py_rot];
bodyrot=bodyrot/sqrt(sum(bodyrot.^2));
bodyangle=acos(sum(pn(1,frameno).*bpvec));

if (px_rot==1) 
    bodyrot=[1 0 0];
    bodyangle=0;
end

ex_rot=eyen(1,frameno);
ey_rot=eyen(2,frameno);
ez_rot=eyen(3,frameno);
eyerot=[ex_rot;-ey_rot;0];
eyerot=eyerot/sqrt(sum(eyerot.^2));
eyeangle=acos(sum(eyen(1,frameno).*bevec));
Euglenaplot(xc,yc,zc,bodylength,thickness,bodyrot,bodyangle,eyerot,eyeangle,0.6,1/0.6)

[x2, y2, z2] = ellipsoid(xc+pn(1,frameno)*bodylength*0.3+eyen(1,frameno)*thickness*0.35,yc+pn(2,frameno)*bodylength*0.3+eyen(2,frameno)*thickness*0.35,zc+pn(3,frameno)*bodylength*0.3+eyen(3,frameno)*thickness*0.35,thickness*0.25,thickness*0.25,thickness*0.25,10);
E = surfl(x2, y2, z2);
set(E,'FaceColor',[1 0 0],'FaceAlpha',1,'EdgeColor','none'); 
axis equal;
end
hold on;
plot3(xn(1,1:frameno),xn(2,1:frameno),xn(3,1:frameno),'b-','linewidth',2);
set(gca,'FontSize',16);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
ylim([-5 5]);
zlim([-5 5]);
xlim([-10 10]);
view([45 30]);