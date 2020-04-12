clc;
clear all;
close all;
N=512;
wo=N/2;
d1=1*10^-6;
NoiseLevel=5;
lambda=637*10^-9;
x=(-N/2:N/2-1)*d1;
y=(-N/2:N/2-1)*d1;
[X,Y]=meshgrid(x,y);
umax=1/(2*d1); 
df=1/(N*d1);
u=(-N/2:N/2-1)*df;
v=(-N/2:N/2-1)*df;
[U,V]=meshgrid(u,v);
RR=sqrt(U.^2 + V.^2);
% I = imread('testpat.tiff');
% I = I(257:768,257:768);
% I = double(I);
I=imread('Lena.png');
I=rgb2gray(I);
I=double(I);
figure
imagesc(x*10^3,y*10^3,I);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['Input Image']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"1.png");

figure
imagesc(log(abs(fftshift(fft2(I)))));
colorbar
colormap(gray);
axis on
title(['Object Fourier transform']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"2.png");

scale=800000;
R=sqrt(X.^2 +Y.^2);
psf=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2; 

figure
imagesc(x*10^3,y*10^3,psf);
colorbar
colormap(jet);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['PSF']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"3.png");

Otf=fftshift(fft2(psf));
Otf= real(Otf./max(max(abs(Otf))));
Otf=abs(Otf);

figure
imagesc(u,v,(Otf));
colorbar
colormap(jet);
xlabel('u') ;
ylabel('v');
axis on
title(['OTF']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"4.png");

figure
mesh(U,V,(Otf));
colorbar
colormap(jet);
xlabel('u') ;
ylabel('v');
axis on
title(['OTF']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"5.png");

[Im,IA1,IA2,IA3,IB1,IB2,IB3,IC1,IC2,IC3,ILt1p1]= SIMimages(I,Otf,x,y,X,Y,u,v,U,V,N,NoiseLevel);

fftemp = fftshift(fft2(ILt1p1));
figure
hold on
plot(u,abs(Otf(wo+1,:)),'--','LineWidth',2,'MarkerSize',6)
plot(u,0.5*abs(fftemp(wo+1,:))./max(max(abs(fftemp))),'r--','LineWidth',2,'MarkerSize',6)
legend('OTFo','illumination spectrum')
grid on
box on
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"6.png");

figure
imagesc(x*10^3,y*10^3,Im);
colorbar
colormap(gray);
xlabel('x(mm)') ;
ylabel('y(mm)');
axis on
title(['Image with uinform illumination']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"7.png");

figure
mesh(log(abs(fftshift(fft2(Im)))));
colorbar
colormap(gray);
axis on
title(['Fourier transform of Image with uinform illumination']);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,"8.png");

% 
% figure
% imagesc(x*10^3,y*10^3,IA1);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IA1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"9.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IA1)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IA1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"10.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IA2);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IA2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"11.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IA2)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IA2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"12.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IA3);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IA3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"13.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IA3)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IA3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"14.png");
% 
% 
% figure
% imagesc(x*10^3,y*10^3,IB1);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IB1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"15.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IB1)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IB1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"16.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IB2);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IB2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"17.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IB2)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IB2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"18.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IB3);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IB3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"19.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IB3)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IB3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"20.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IC1);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IC1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"21.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IC1)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IC1']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"22.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IC2);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IC2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"23.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IC2)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IC2']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"24.png");
% 
% figure
% imagesc(x*10^3,y*10^3,IC3);
% colorbar
% colormap(gray);
% xlabel('x(mm)') ;
% ylabel('y(mm)');
% axis on
% title(['IC3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"25.png");
% 
% figure
% imagesc(log(abs(fftshift(fft2(IC3)))));
% colorbar
% colormap(gray);
% axis on
% title(['Fourier transform of IC3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"26.png");

[Kotf,cindex] = OTFcutoff(Otf);
psfe=psfforedgetaper(psf);
% figure
% imagesc(psfe);
% colorbar
% colormap(jet);
% axis on
% title(['PSF for edge tapering']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"27.png");

Kotf=abs(u(cindex));

[fvecA]= fvector_optimizer(IA1,IA2,IA3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y);

[fvecB]= fvector_optimizer(IB1,IB2,IB3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y);

[fvecC]= fvector_optimizer(IC1,IC2,IC3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y);
fvecC=-fvecC;
[phaseA]= phase_optimizer(IA1,IA2,IA3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y,fvecA);
[phaseB]= phase_optimizer(IB1,IB2,IB3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y,fvecB);
[phaseC]= phase_optimizer(IC1,IC2,IC3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y,fvecC);

[AF1,AF2,AF3]=separate_noisyfcomponenets(IA1,IA2,IA3,fvecA,phaseA,psfe);

% figure
% imagesc(log(abs(AF1)));
% title('S(k)H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(AF2)));
% title('S(k-p_{\theta_1})H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(AF3)));
% title('S(k+p_{\theta_1})H(k)');
% colormap(gray)

[BF1,BF2,BF3]=separate_noisyfcomponenets(IB1,IB2,IB3,fvecB,phaseB,psfe);

% figure
% imagesc(log(abs(BF1)));
% title('S(k)H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(BF2)));
% title('S(k-p_{\theta_2})H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(BF3)));
% title('S(k+p_{\theta_2})H(k)');
% colormap(gray)


[CF1,CF2,CF3]=separate_noisyfcomponenets(IC1,IC2,IC3,fvecC,phaseC,psfe);

% figure
% imagesc(log(abs(CF1)));
% title('S(k)H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(CF2)));
% title('S(k-p_{\theta_3})H(k)');
% colormap(gray)
% 
% figure
% imagesc(log(abs(CF3)));
% title('S(k+p_{\theta_3})H(k)');
% colormap(gray)

NoisySkHk=(AF1+BF1+CF1)./3;

figure
imagesc(log(abs(NoisySkHk)));
title('Noisy S(k)H(k)');
colormap(gray)

[A,alpha]=Object_Para_Estimation(NoisySkHk,Otf,Kotf,U,V);

K=sqrt(U.^2+V.^2);
Aobj = A;
Bobj = -alpha;
OBJo = Aobj*(K.^-Bobj);
SIGp = OBJo.*Otf;
pp = 3;
figure;
hold on
plot([0:512-1]-256,log(abs(NoisySkHk(256+pp,:))),'k--','LineWidth',3,'MarkerSize',6)
plot([0:512-1]-256,log(abs(SIGp(256+pp,:))),'mo-','LineWidth',2,'MarkerSize',6)
legend('acutal signal power','avg. signal power')
grid on
box on

function [Obj1,Obj2]=Object_Para_Estimation(NoisySkHk,Otf,Kotf,U,V)
K=sqrt(U.^2+V.^2);
CC=(K>0.3*Kotf).*(K<0.4*Kotf);
NSK=NoisySkHk.*CC;
A = sum(sum(abs(NSK)))./sum(sum(CC));
alpha = -0.5;
OBJparam = [A alpha];
Objparaopt=@(OBJparam)Objoptfunc(OBJparam,NoisySkHk,Kotf,U,V,K,Otf,CC)
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
[Objpara,fval]=fminsearch(Objparaopt,OBJparam,options);
Obj1=Objpara(1);
Obj2= Objpara(2);
end
function C= Objoptfunc(OBJparam,NoisySkHk,Kotf,U,V,K,Otf,CC)
K(257,257)=1;
A= OBJparam(1,1);
alpha= OBJparam(1,2);
Objamp=A.*(K.^(alpha));
Signalamp=Objamp.*(abs(Otf));
NOisef=1.5*Kotf;
ZZ=K>NOisef;
Noisespectrum=NoisySkHk.*ZZ;
NoisePower = sum(sum(abs(Noisespectrum).^2))./sum(sum(ZZ));
Noisefreesignalamp=(abs(NoisySkHk).^2)-NoisePower;
Noisefreesignalamp=sqrt(Noisefreesignalamp);
Error=(Noisefreesignalamp)-Signalamp;
Zloop = (K<0.75*Kotf).*((K>0.25*Kotf));
invK=1./K;
C=sum(sum(((Error.^2.*invK).*Zloop)));
end
function [FC1,FC2,FC3]= separate_noisyfcomponenets(I1,I2,I3,fvector,phase,psfe)
I1 = edgetaper(I1,psfe);
I2 = edgetaper(I2,psfe);
I3 = edgetaper(I3,psfe);

FI1=fftshift(fft2(I1));
FI2=fftshift(fft2(I2));
FI3=fftshift(fft2(I3));

m=1;
phi1=phase(1,1);
phi2=phase(1,2);
phi3=phase(1,3);

M11= 1;
M12=-(m/2)*exp(-1j.*phi1);
M13=(-m/2)*exp(1j.*phi1);

M21= 1;
M22=-(m/2)*exp(-1j.*phi2);
M23=-(m/2)*exp(1j.*phi2);

M31= 1;
M32=-(m/2)*exp(-1j.*phi3);
M33=-(m/2)*exp(1j.*phi3);

M= [M11 M12 M13; M21 M22 M23; M31 M32 M33;];

Minv=inv(M);
FC1=Minv(1,1).*FI1 + Minv(1,2).*FI2 + Minv(1,3).*FI3;
FC2=Minv(2,1).*FI1 + Minv(2,2).*FI2 + Minv(2,3).*FI3;
FC3=Minv(3,1).*FI1 + Minv(3,2).*FI2 + Minv(3,3).*FI3;
end
function [C] =phaseoptfunc(phiest,II,fvec,X,Y)
C1=II;
C2=-cos((2*pi*(X.*fvec(1)+ Y.*fvec(2)))-phiest);
C=-(sum(sum(C1.*C2)));
end
function [phase]= phase_optimizer(I1,I2,I3,Otf,psf,x,y,Kotf,U,V,u,~,X,Y,fvec)
N = size(I1,1);
wo = N/2;
h = 30;
psfe = psf(wo-h+1:wo+h,wo-h+1:wo+h);
I1_et = edgetaper(I1,psfe); 
I2_et = edgetaper(I2,psfe);
I3_et = edgetaper(I3,psfe);

ph = @(phiest)phaseoptfunc(phiest,I1_et,fvec,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
phiest = 0;
[phase1,fval] = fminsearch(ph,phiest,options);

ph = @(phiest)phaseoptfunc(phiest,I2_et,fvec,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
phiest = 2*pi/3;
[phase2,fval] = fminsearch(ph,phiest,options);

ph = @(phiest)phaseoptfunc(phiest,I3_et,fvec,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
phiest = 4*pi/3;
[phase3,fval] = fminsearch(ph,phiest,options);
phase=[phase1 phase2 phase3];
end
function [fvector,Ix,Iy] = fvectorapproximate(FiSMap,Kotf,U,V,u,v)
% figure
% imagesc(log(abs(FiSMap)));
% colormap(gray)

R=sqrt(U.^2+V.^2);
Z0 = R> round(0.5*Kotf);
Z1 = U > 0;
FiSMap=FiSMap.*Z0.*Z1;
% figure
% imagesc(log(abs(FiSMap)));
% colormap(gray)

dumX = max( FiSMap,[],2 );
[dummx Ix] = max(dumX);
dumY = max( FiSMap,[],1);
[~, Iy] = max(dumY);

fvector = [U(Ix,Iy) V(Ix,Iy)];

end
function [fvec]= fvector_optimizer(I1,I2,I3,Otf,psf,x,y,Kotf,U,V,u,v,X,Y)

N = size(I1,1);
wo = N/2;
h = 30;
psfe = psf(wo-h+1:wo+h,wo-h+1:wo+h);
I1_et = edgetaper(I1,psfe);
I2_et = edgetaper(I2,psfe);
I3_et = edgetaper(I3,psfe);
FI1_et=fftshift(fft2(I1_et));
FI2_et=fftshift(fft2(I2_et));
FI3_et=fftshift(fft2(I3_et));
[fvector1,Ix,Iy] = fvectorapproximate(FI1_et,Kotf,U,V,u,v);
OPT = 1;
fv = @(k2fa0)frequencyoptfunc(k2fa0,Otf,OPT,I1_et,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');%,'PlotFcns',@optimplotfval);
k2fa0 = fvector1;
[fvec1,fval] = fminsearch(fv,k2fa0,options);

[fvector2,Ix,Iy] = fvectorapproximate(FI2_et,Kotf,U,V,u,v);
OPT = 1;
fv = @(k2fa0)frequencyoptfunc(k2fa0,Otf,OPT,I2_et,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');%,'PlotFcns',@optimplotfval);
k2fa0 = fvector2;
[fvec2,~] = fminsearch(fv,k2fa0,options);

[fvector3,Ix,~] = fvectorapproximate(FI3_et,Kotf,U,V,u,v);

OPT = 1;
fv = @(k2fa0)frequencyoptfunc(k2fa0,Otf,OPT,I3_et,X,Y);
options = optimset('LargeScale','off','Algorithm','active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');%,'PlotFcns',@optimplotfval);
k2fa0 = fvector3;
[fvec3,fval] = fminsearch(fv,k2fa0,options);

fvec = (fvec1 + fvec2 + fvec3)./3;

end
function [C] = frequencyoptfunc(fm,Otf,OPT,I,X,Y)

Inew=I.*exp(-1j*2*pi*(fm(1).*X + fm(2).*Y));
Inew=fftshift(fft2(Inew));
C1=fftshift(fft2(I)).*conj(Otf);
C2=conj(((Inew)).*Otf);
C=(sum(sum(C1.*C2)));
C = C./sum(sum( C2.*conj(C2) ));
C=abs(C);
if (OPT > 0)
    C = -abs(C);
else
    C = C;    
    
end
end
function [Kotf,cindex] = OTFcutoff(Otf)
w = size(Otf,1);
wo = w/2;
OTFarray=abs(Otf(wo+1,:));
v=0.001;
[Kotf,cindex] = (min(abs(OTFarray - v)));
end
function [Im,IA1,IA2,IA3,IB1,IB2,IB3,IC1,IC2,IC3,ILt1p1]= SIMimages(I,Otf,x,y,X,Y,u,v,U,V,N,NoiseLevel)

fm=175781.25;
m=0.8;
a=0.5;
theta1=0;
theta2=pi/3;
theta3=2*pi/3;
phi1=0;
phi2=2*pi/3;
phi3=4*pi/3;

randomphaseerror=(0.5-rand(9,1))*pi/18;

ILt1p1= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta1)+ Y.*fm.*sin(theta1))) - (phi1+randomphaseerror(1,1)));
ILt1p2= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta1)+ Y.*fm.*sin(theta1))) - (phi2+randomphaseerror(2,1)));
ILt1p3= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta1)+ Y.*fm.*sin(theta1))) - (phi3+randomphaseerror(3,1)));
ILt2p1= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta2)+ Y.*fm.*sin(theta2))) - (phi1+randomphaseerror(4,1)));
ILt2p2= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta2)+ Y.*fm.*sin(theta2))) - (phi2+randomphaseerror(5,1)));
ILt2p3= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta2)+ Y.*fm.*sin(theta2))) - (phi3+randomphaseerror(6,1)));
ILt3p1= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta3)+ Y.*fm.*sin(theta3))) - (phi1+randomphaseerror(7,1)));
ILt3p2= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta3)+ Y.*fm.*sin(theta3))) - (phi2+randomphaseerror(8,1)));
ILt3p3= a - 0.5*m*cos((2*pi*(X.*fm.*cos(theta3)+ Y.*fm.*sin(theta3))) - (phi3+randomphaseerror(9,1)));
% 
% figure
% imagesc(x,y,ILt1p1);
% colorbar
% colormap(gray);
% axis on
% title(['ILt1p1']);
% 
% figure
% imagesc(x,y,ILt1p2);
% colorbar
% colormap(gray);
% axis on
% title(['ILt1p2']);
% 
% figure
% imagesc(x,y,ILt1p3);
% colorbar
% colormap(gray);
% axis on
% title(['ILt1p3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"28.png");

% figure
% imagesc(x,y,ILt2p1);
% colorbar
% colormap(gray);
% axis on
% title(['ILt2p1']);
% 
% figure
% imagesc(x,y,ILt2p2);
% colorbar
% colormap(gray);
% axis on
% title(['ILt2p2']);

% figure
% imagesc(x,y,ILt2p3);
% colorbar
% colormap(gray);
% axis on
% title(['ILt2p3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"29.png");
% 
% 
% figure
% imagesc(x,y,ILt3p1);
% colorbar
% colormap(gray);
% axis on
% title(['ILt3p1']);
% 
% figure
% imagesc(x,y,ILt3p2);
% colorbar
% colormap(gray);
% axis on
% title(['ILt3p2']);
% 
% figure
% imagesc(x,y,ILt3p3);
% colorbar
% colormap(gray);
% axis on
% title(['ILt3p3']);
% set(gcf, 'PaperPositionMode', 'auto');
% saveas(gcf,"30.png");

Im=    ifft2(ifftshift(Otf.*fftshift(fft2(I))));
IA1 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt1p1)))));
IA2 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt1p2)))));
IA3 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt1p3)))));
IB1 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt2p1)))));
IB2 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt2p2)))));
IB3 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt2p3)))));
IC1 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt3p1)))));
IC2 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt3p2)))));
IC3 = (ifft2(ifftshift(Otf.*fftshift(fft2(I.*ILt3p3)))));
 
Noisepercentage=NoiseLevel/100;
SNR_image=1/Noisepercentage;
nn=1;

w=N;
nA1 = random('norm', 0, Noisepercentage*std2(IA1), w , w);
nA2 = random('norm', 0, Noisepercentage*std2(IA2), w , w);
nA3 = random('norm', 0, Noisepercentage*std2(IA3), w , w);
nB1 = random('norm', 0, Noisepercentage*std2(IB1), w , w);
nB2 = random('norm', 0, Noisepercentage*std2(IB2), w , w);
nB3 = random('norm', 0, Noisepercentage*std2(IB3), w , w);
nC1 = random('norm', 0, Noisepercentage*std2(IC1), w , w);
nC2 = random('norm', 0, Noisepercentage*std2(IC2), w , w);
nC3 = random('norm', 0, Noisepercentage*std2(IC3), w , w);
nIm = random('norm', 0, Noisepercentage*std2(Im), w , w);

Im=Im+nn.*nIm;
IA1=IA1+nn.*nA1;
IA2=IA2+nn.*nA2;
IA3=IA3+nn.*nA3;
IB1=IB1+nn.*nB1;
IB2=IB2+nn.*nB2;
IB3=IB3+nn.*nB3;
IC1=IC1+nn.*nC1;
IC2=IC2+nn.*nC2;
IC3=IC3+nn.*nC3;


end
function psfe=psfforedgetaper(psf)
N = size(psf,1);
wo = N/2;
h = 30;
psfe = psf(wo-h+1:wo+h,wo-h+1:wo+h);
end
