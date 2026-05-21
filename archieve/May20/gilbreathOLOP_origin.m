%*****************************************olop.m***************************
%This M-file is used for application of the Open-loop Onset Point (OLOP) to an aircraft system
%with a rate limiting element in it.  Two simulink models are needed for the analysis.  The first is
%the linear model of the aircraft system with the pilot output as the input and theta as the output
%to include any feedback gains.  The linear model should be named 'namelin.mdl'.  Once the pilot 
%model 
%is chosen the second model is developed called 'nameolop.mdl'.  It has the pilot model in the loop
%and has the output of the rate limiter as its input and the output is the input into the rate limiter.
%**************************************************************************
clear
format short
lin_abcd %Initialization File to get aircraft model
w=logspace(-2,2,1000);
%Define the rate limiter value and the pilot input amplitude
R=input(' Enter Vector of Rate Limits to evaluate[default = 60]: ');
if(isempty(R))
 R=60;
end

Per=input(' Enter Stick Input Amplitude in percentage (range 0-100) [default = 100% ]: ');
if(isempty(Per))
 Per=100;
end

Po=Per/100*3.6
if Po<.8 %Takes into account non-linear stik gradient
 Po=2.5*Po;
else Po=2.5*(2.86*Po-1.48);end

tag=input(' Choose desired pilot model (1=Low Gain 2=Med Gain 3=High gain 4=Neal-Smith)  [default Neal-Smith]:');
if(isempty(tag))
 tag=4;
elseif tag==1; xover=-90; %Set crossover angle for gain pilot model
elseif tag==2; xover=-110;
elseif tag==3; xover=-130;
end

Kp=-.05*2.5; Tp1=.06; Tp2=.01; %Gain, lead, and lag time constants for Neal-Smith Pilot Model
%Use LINMOD to find appropriate transfer function for Yc=theta/pilot
[A,B,C,D]=LINMOD('nt33v2Dlin');
Yc=ss(A,B,C(4,:),D(4,:)); %Make sure appropriate row is selected for theta output

%Find pure gain pilot model
Kg=.1*2.5:.01:10;
for i=1:length(Kg)
 [gm,pm,wcg,wcp]=margin(Kg(i)*Yc);
 if tag==4; break, end
 if -180+pm > xover
  xover_angle=-180+pm,
  Pilot_gain=Kg(i)
  break
 end
end

figure(1),ngridneal,nichols(Kg(i)*Yc,w)

%Find Neal-Smith Modified Pilot model
bw=3.5; %Set bandwidth requirement
Yp=Kp*tf([5 1],[1 0])*tf([Tp1 1],[Tp2 1]); set(Yp,'InputDelay',.25);
Ypp=pade(Yp,2);%Need pade approximation for linmod to get OLOP transfer function
figure(2),ngridneal,nichols(Yp*Yc,{.1,3*bw})

%Find closed-loop response amplitude from stick input to RLE input
P_RLE=ss(A,B,C(6,:),D(6,:));
[magP,phP]=bode(P_RLE,w);

for n=1:length(R)
%Calculate Rate Limit Line for frequency plot
[magR,phR]=bode(tf(R(n),[1 0]),w);
%Plot closed-loop amplitude and rate limit line
figure(3),semilogx(w,20*log10(Po*magP(1,:)),w,20*log10(magR(1,:)),'k:')
hold on,axis([.01,10,10,60])
title('Determine Closed-loop Onset Frequency')
%Calculate omega onset for closed-loop system
iter=100;
for i=1:length(w)
 dif(i)=abs(magR(i)-Po*magP(i));
 if dif(i)<= iter
  iter=dif(i); k=i;
 end
end
w_onset=w(k),Rate=R(n)

%Calculate OLOP on Nichols chart

if tag==4
Pilot=Ypp; %Neal-Smith Pilot Model with Pade approximation for time delay
else
Pilot=tf(Pilot_gain,1); %Simple Gain Pilot Model
end

[AA,BB,CC,DD]=linmod('nt33v2dolop');
OLOP=ss(AA,BB,-CC,-DD);
[magN,phN]=nichols(OLOP,w);
r=150; %Plot range
v1=[-60 -90 -100 -120 -140 -160 -180];v2=[13.5 7.5 5.5 2.5 1.1 0 0];%OLOP stability boundary
%Plot OLOP points for each rate value
if n==1
figure(4),plot(phN(k),20*log10(magN(k)),'o','MarkerFaceColor','k','MarkerSize',7)
elseif n==2
figure(4),plot(phN(k),20*log10(magN(k)),'d','MarkerFaceColor','k','MarkerSize',7)
elseif n==3
figure(4),plot(phN(k),20*log10(magN(k)),'o','MarkerSize',7)
elseif n==4
figure(4),plot(phN(k),20*log10(magN(k)),'p','MarkerFaceColor','k','MarkerSize',7)
elseif n==5
figure(4),plot(phN(k),20*log10(magN(k)),'s','MarkerFaceColor','k','MarkerSize',7)
elseif n==6
figure(4),plot(phN(k),20*log10(magN(k)),'d','MarkerSize',7)
else
figure(4),plot(phN(k),20*log10(magN(k)),'s','MarkerSize',7)
end
hold on
%Plot mag and phase on Nichols chart
plot(phN(1,k-r:r+k),20*log10(magN(1,k-r:r+k)),...
 [-360 0],[0 0],'-k',[-180 -180],[-50 50],'-k')
axis([-200,-60,-10,20]);
plot(v1,v2,'-k','LineWidth',1.5)%Plot OLOP stability boundary
end