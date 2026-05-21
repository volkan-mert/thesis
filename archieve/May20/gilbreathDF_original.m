%This M-file will calculate the closed-loop describing function
%for the aircraft model example from AIAA-95-3204-CP
clear all
global w R K magGc magGac phGc phGac dtr qco z
%System input frequency range, amplitude, rate limit
w=logspace(-1,2,100); qco=1.1; R=60;

%System numbers
K=13.68;
Gc=tf(5.21*conv([1 -57.36],conv([1 4.26],[1 .55])),conv([1 2*.442*22.85 22.85^2],conv([1 0],[1 1.16])));
Gac=tf(-10.524*conv([1 1.562],conv([1 .038],[1 0])),conv([1 2*.212*.088 .088^2],conv([1 3.75],[1 -1.44])));
num=K*Gc;
den=zpk([],[],1)+series(Gc,Gac);
pcl=num/den;
dtr=pi/180;

%Find linear response amplitude and phase for initial guess
[magpcl,phi2(1)]=bode(pcl,w(1));
deltao(1)=qco*magpcl; phi2(1);

for z=1:length(w)
    [magGc,phGc]=bode(Gc,w(z));
    [magGac,phGac]=bode(Gac,w(z));

    xo=[deltao(z) phi2(z)]; %initial guess for optimization algorithm

    [a,fval(z)]=fminsearch('eqs',xo,optimset('Display','off'));

    deltait=a(1);
    phi2it=a(2);

    [magNit,phNit]=dfunction(w(z),R,deltait);
    A=deltait*magGac*magNit/(K*qco);
    phi=phi2it+phGac+phNit*180/pi;

    %Calculate Open-loop describing function for Nichols plot
    NoMag(z)=20*log10(A/sqrt(1-2*A*cos(dtr*phi)+A^2));
    NoPh(z)=dtr*phi-atan2(-A*sin(dtr*phi),1-A*cos(dtr*phi));

    %Give new initial guess if function cannot converge to zero
    %Disply onset frequency where discontinuity starts to occur
    if fval(z)>.0001
        deltao(z+1)=24; phi2(z+1)=-237;onset=w(z)
    else
        %Else, use previous solution for next starting point
        deltao(z+1)=a(1); phi2(z+1)=a(2);
    end
end

[mlin,plin]=bode(Gc*Gac,w);
figure(1) %Plot results on Nichols Chart
plot(plin(1,:)-360,20*log10(mlin(1,:)),'-o',NoPh/dtr,NoMag,'-*')
axis([-300,-50,-20,12])

figure(2) %Plot function evaluation vs. frequency to see convergence
semilogx(w,fval)


function f=eqs(x)
global w R dtr magGc magGac phGc phGac K qco z
deltao=x(1); phi2=x(2);

f=(deltao/magGc*cos(dtr*(phi2-phGc))+deltao*magGac*dfuncmag(w(z),R,deltao)...
    *cos(dtr*(phi2+phGac)+dfuncphase(w(z),R,deltao))-K*qco)^2.0...
    +(deltao/magGc*sin(dtr*(phi2-phGc))+deltao*magGac*dfuncmag(w(z),R,deltao)...
    *sin(dtr*(phi2+phGac)+dfuncphase(w(z),R,deltao)))^2.0;

%Function file that determines magnitude and phase of describing function given the input frequency
%rate limit and input amplitude. Includes cubic spline interpolation coefficients.
function [magN,phN]=dfunction(freq,rate,inamp)
%Describing Function of a rate limiting element
x=freq*inamp/rate;
if x<1, magN=1;phN=0;
elseif x<1.862
    magN=polyval([.2908 -1.4396 1.9232 .223],x);%From cubic spline interpolation
    phN=polyval([.528 -2.6213 3.5056 -1.4171],x);% From cubic spline interpolation
else
    magN=4/x/pi;
    phN=-acos(pi/x/2);
end