%%
%Conversor DC-DC
Vin = 1.2;
Vout = 0.6;
Verror=0.99*Vout;
Rl = 1000;
C1 = 0.5*10^-9;

Lmin=120*10^-9;


Fclk=Verror/(2*Rl*C1*(Vin-2*Verror))
Tclk=1/Fclk;

Pout=Vout^2/Rl

rendimento = abs(-Vin*2*C1+4*C1)/abs(Vin*(-Vin*C1+2*C1))

 
%syms tau
%vpasolve(Vout-(Vout-0)*exp(-Tclk/(2*tau))==0.01,tau)

tau=-Tclk/(2*log(0.01))

RonTotal=tau/C1

Wp=Lmin/(RonTotal*145*10^-6*(Vin-0.28)) %Wp=205um

Wn=Lmin/(RonTotal*545*10^-6*(0.6-0.3))  %Wn=167um

%Wp1=Wp2=205um;L1=L2=120nm
%Wn1=Wn2=167um;L1=L2=120nm

%%
%Comparador
%%
%Parte 2

Vout2=0.39;

Tclk2=(12*C1*Rl*Vin-36*C1*Rl*Vout2)/8*Vout2

Fclk2=1/Tclk2

Ctotal=C1*C1/(C1+C1)
tau=-Tclk2/(2*log(0.01))

RonTotal2=tau/Ctotal

Ron1=RonTotal/3;

Wp1=Lmin/(Ron1*145*10^-6*(Vin-0.28))

Wn1=Lmin/(Ron1*545*10^-6*(0.6-0.3))

Ron2=0.5*RonTotal

Wp2=Lmin/(Ron2*145*10^-6*(Vin-0.28))

Wn2=Lmin/(Ron2*545*10^-6*(0.6-0.3))
