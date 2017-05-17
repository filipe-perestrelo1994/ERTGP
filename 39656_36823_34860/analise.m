%%
Vdd = 0.7;
P = 4*10^(-6);
I_total = P/Vdd;

Id_1 = 1.4*10^-6; %Corrente que é injectada no par diferencia

Id_2 = 2*10^-6; %Corrente que é injectada no andar de saída

Iref = 1.4*10^-6; %Corrente de referencia


Lmin = 120*10^-9; %L mínimo

%%Vai buscar os ficheiros NMOS e PMOS e passá-los para variáveis
NMOS = csvread('NMOS(data).csv');
Vdc_NMOS = NMOS(:,1);
gmOverId_NMOS = NMOS(:,2);
gm_NMOS = NMOS(:,3);
id_NMOS = NMOS(:,4);

PMOS = csvread('PMOS(data).csv');
Vdc_PMOS = PMOS(:,1);
gmOverId_PMOS = PMOS(:,2);
gm_PMOS = PMOS(:,3);
id_PMOS = PMOS(:,4);

%% Gráficos dos I0 de NMOS e PMOS

%%%NMOS
y_NMOS = gmOverId_NMOS;
x_NMOS = id_NMOS;

Id_max_NMOS = max(y_NMOS);

recta1_NMOS = Id_max_NMOS*ones(size(x_NMOS));

declive_NMOS = (log10(y_NMOS(30))-log10(y_NMOS(36)))/(log10(x_NMOS(30))-log10(x_NMOS(36)))
b_NMOS = log10(y_NMOS(40))- declive_NMOS*log10(x_NMOS(40));

recta2_NMOS = 10.^(declive_NMOS*log10(x_NMOS)+b_NMOS);

figure(1)

loglog(x_NMOS,y_NMOS)
grid on
hold on
loglog(x_NMOS,recta1_NMOS)

loglog(x_NMOS,recta2_NMOS)
hold off



I0_NMOS = 5.847*10^(-7); %Por observação


%%%PMOS
x_PMOS = -1.*id_PMOS;
y_PMOS = gmOverId_PMOS;

Id_max_PMOS = max(y_PMOS);
recta1_PMOS = Id_max_PMOS*ones(size(x_PMOS));

declive_PMOS = (log10(y_PMOS(20))-log10(y_PMOS(19)))/(log10(x_PMOS(20))-log10(x_PMOS(19)))
b_PMOS = log10(y_PMOS(20))- declive_PMOS*log10(x_PMOS(20));

recta2_PMOS = 10.^(declive_PMOS*log10(x_PMOS)+b_PMOS);

figure(2)
loglog(x_PMOS,y_PMOS)
grid on
hold on
loglog(x_PMOS,recta1_PMOS)

loglog(x_PMOS,recta2_PMOS)
hold off


I0_PMOS =1.399*10^(-7); %Por observação 

%%
%Cálculo dos W's
%Par diferencial




L = 7*Lmin;
IC=0.3;
W_12 = (Id_1/2)/((I0_PMOS*IC)/L) %(W_12 > 2*Lmin)

L = 5*Lmin;
IC = 0.3;
W_34 = (Id_1/2)/((I0_NMOS*IC)/L) %(W_34 > 2*Lmin)


L = 10*Lmin;
IC = 10
W_5 = Id_1/((I0_PMOS*IC)/L) %(W_5 > 9*Lmin)

%Segundo andar

L = 10*Lmin;
IC = 6.75;
W_6 = Id_2/((I0_PMOS*IC)/L) %(W_6 > 2*Lmin)


L = 10*Lmin;
IC = 1;
W_7 = Id_2/((I0_NMOS*IC)/L) %

%%
% Tirar os Vearly dos NMOS para se calcular os rds e gds
NMOS_Early = csvread('NMOS_Vearly.csv');
NMOS_EARLY_ID = csvread('NMOS_ID.csv');

Vds_NMOS = NMOS_Early(:,1);

Vearly_NMOS_240n = NMOS_Early(:,2);
Id_NMOS_240n = NMOS_EARLY_ID(:,2);


Vearly_NMOS_360n = NMOS_Early(:,3);
Id_NMOS_360n = NMOS_EARLY_ID(:,3);


Vearly_NMOS_480n = NMOS_Early(:,4);
Id_NMOS_480n = NMOS_EARLY_ID(:,4);



Vearly_NMOS_600n = NMOS_Early(:,5);
Id_NMOS_600n = NMOS_EARLY_ID(:,5);


Vearly_NMOS_720n = NMOS_Early(:,6);
Id_NMOS_720n = NMOS_EARLY_ID(:,6);


Vearly_NMOS_840n = NMOS_Early(:,7);
Id_NMOS_840n = NMOS_EARLY_ID(:,7);


Vearly_NMOS_960n = NMOS_Early(:,8);
Id_NMOS_960n = NMOS_EARLY_ID(:,8);


Vearly_NMOS_1u80n = NMOS_Early(:,9);
Id_NMOS_1u80n = NMOS_EARLY_ID(:,9);


Vearly_NMOS_1u200n = NMOS_Early(:,10);
Id_NMOS_1u200n = NMOS_EARLY_ID(:,10);


Vearly_NMOS_1u320n = NMOS_Early(:,11);
Id_NMOS_1u320n = NMOS_EARLY_ID(:,11);


Vearly_NMOS_1u440n = NMOS_Early(:,12);
Id_NMOS_1u440n = NMOS_EARLY_ID(:,12);


Vearly_NMOS_1u560n = NMOS_Early(:,13);
Id_NMOS_1u560n = NMOS_EARLY_ID(:,13);


Vearly_NMOS_1u680n = NMOS_Early(:,14);
Id_NMOS_1u680n = NMOS_EARLY_ID(:,14);


Vearly_NMOS_1u800n = NMOS_Early(:,15);
Id_NMOS_1u800n = NMOS_EARLY_ID(:,15);


Vearly_NMOS_1u920n = NMOS_Early(:,16);
Id_NMOS_1u920n = NMOS_EARLY_ID(:,16);

Vearly_NMOS_2u40n = NMOS_Early(:,17);
Id_NMOS_2u40n = NMOS_EARLY_ID(:,17);

Vearly_NMOS_2u160n = NMOS_Early(:,18);
Id_NMOS_2u160n = NMOS_EARLY_ID(:,18);

Vearly_NMOS_2u280n = NMOS_Early(:,19);
Id_NMOS_2u280n = NMOS_EARLY_ID(:,19);

Vearly_NMOS_2u400n = NMOS_Early(:,20);
Id_NMOS_2u400n = NMOS_EARLY_ID(:,20);

Vearly_NMOS_2u520n = NMOS_Early(:,21);
Id_NMOS_2u520n = NMOS_EARLY_ID(:,21);

Vearly_NMOS_2u640n = NMOS_Early(:,22);
Id_NMOS_2u640n = NMOS_EARLY_ID(:,22);

Vearly_NMOS_2u760n = NMOS_Early(:,23);
Id_NMOS_2u760n = NMOS_EARLY_ID(:,23);

Vearly_NMOS_2u880n = NMOS_Early(:,24);
Id_NMOS_2u880n = NMOS_EARLY_ID(:,24);

Vearly_NMOS_3u = NMOS_Early(:,25);
Id_NMOS_3u = NMOS_EARLY_ID(:,25);

Vearly_NMOS_3u120n = NMOS_Early(:,26);
Id_NMOS_3u120n = NMOS_EARLY_ID(:,26);

Vearly_NMOS_3u240n = NMOS_Early(:,27);
Id_NMOS_3u240n = NMOS_EARLY_ID(:,27);

Vearly_NMOS_3u360n = NMOS_Early(:,28);
Id_NMOS_3u360n = NMOS_EARLY_ID(:,28);

Vearly_NMOS_3u480n = NMOS_Early(:,29);
Id_NMOS_3u480n = NMOS_EARLY_ID(:,29);

Vearly_NMOS_3u600n = NMOS_Early(:,30);
Id_NMOS_3u600n = NMOS_EARLY_ID(:,30);

Vearly_NMOS_3u720n = NMOS_Early(:,31);
Id_NMOS_3u720n = NMOS_EARLY_ID(:,31);

Vearly_NMOS_3u840n = NMOS_Early(:,32);
Id_NMOS_3u840n = NMOS_EARLY_ID(:,32);

Vearly_NMOS_4u = NMOS_Early(:,33);
Id_NMOS_4u = NMOS_EARLY_ID(:,33);




%(Transistor_34)
figure(3)
IC_NMOS = Id_NMOS_600n/I0_NMOS;
grid on
loglog(IC_NMOS,Vearly_NMOS_600n)
title('Transistor 3,4')


figure(4)
IC_NMOS = Id_NMOS_1u200n/I0_NMOS;
grid on
loglog(IC_NMOS,Vearly_NMOS_1u200n)
title('Transistor 7')

% loglog(Id_NMOS_240n,Vearly_NMOS_240n)
% loglog(Id_NMOS_360n,Vearly_NMOS_360n)
% loglog(Id_NMOS_480n,Vearly_NMOS_480n)
% loglog(Id_NMOS_600n,Vearly_NMOS_600n)
% loglog(Id_NMOS_720n,Vearly_NMOS_720n)
% loglog(Id_NMOS_840n,Vearly_NMOS_840n)
% loglog(Id_NMOS_960n,Vearly_NMOS_960n)
% loglog(Id_NMOS_1u80n,Vearly_NMOS_1u80n)
% loglog(Id_NMOS_1u200n,Vearly_NMOS_1u200n)
% loglog(Id_NMOS_1u320n,Vearly_NMOS_1u320n)
% loglog(Id_NMOS_1u440n,Vearly_NMOS_1u440n)
% loglog(Id_NMOS_1u560n,Vearly_NMOS_1u560n)
% loglog(Id_NMOS_1u680n,Vearly_NMOS_1u680n)
% loglog(Id_NMOS_1u800n,Vearly_NMOS_1u800n)
% loglog(Id_NMOS_2u,Vearly_NMOS_2u)
%hold off





%%
%Tirar os Veary dos PMOS para se calcular os rds e gds
PMOS_EARLY_ID = csvread('PMOS_ID.csv');
PMOS_Early = csvread('PMOS_Vearly.csv');
Vds_PMOS = PMOS_Early(:,1);

Vearly_PMOS_240n = -1.*PMOS_Early(:,2);
Id_PMOS_240n = PMOS_EARLY_ID(:,2);


Vearly_PMOS_360n = -1.*PMOS_Early(:,3);
Id_PMOS_360n = PMOS_EARLY_ID(:,3);


Vearly_PMOS_480n = -1.*PMOS_Early(:,4);
Id_PMOS_480n = PMOS_EARLY_ID(:,4);



Vearly_PMOS_600n = -1.*PMOS_Early(:,5);
Id_PMOS_600n = PMOS_EARLY_ID(:,5);


Vearly_PMOS_720n = -1.*PMOS_Early(:,6);
Id_PMOS_720n = PMOS_EARLY_ID(:,6);


Vearly_PMOS_840n = -1.*PMOS_Early(:,7);
Id_PMOS_840n = PMOS_EARLY_ID(:,7);


Vearly_PMOS_960n = -1.*PMOS_Early(:,8);
Id_PMOS_960n = PMOS_EARLY_ID(:,8);


Vearly_PMOS_1u80n = -1.*PMOS_Early(:,9);
Id_PMOS_1u80n = PMOS_EARLY_ID(:,9);


Vearly_PMOS_1u200n = -1.*PMOS_Early(:,10);
Id_PMOS_1u200n = PMOS_EARLY_ID(:,10);


Vearly_PMOS_1u320n = -1.*PMOS_Early(:,11);
Id_PMOS_1u320n = PMOS_EARLY_ID(:,11);


Vearly_PMOS_1u440n = -1.*PMOS_Early(:,12);
Id_PMOS_1u440n = PMOS_EARLY_ID(:,12);


Vearly_PMOS_1u560n = -1.*PMOS_Early(:,13);
Id_PMOS_1u560n = PMOS_EARLY_ID(:,13);


Vearly_PMOS_1u680n = -1.*PMOS_Early(:,14);
Id_PMOS_1u680n = PMOS_EARLY_ID(:,14);


Vearly_PMOS_1u800n = -1.*PMOS_Early(:,15);
Id_PMOS_1u800n = PMOS_EARLY_ID(:,15);

Vearly_PMOS_1u920n = -1.*PMOS_Early(:,16);
Id_PMOS_1u920n = PMOS_EARLY_ID(:,16);

Vearly_PMOS_2u40n = -1.*PMOS_Early(:,17);
Id_PMOS_2u40n = PMOS_EARLY_ID(:,17);

Vearly_PMOS_2u160n = -1.*PMOS_Early(:,18);
Id_PMOS_2u160n = PMOS_EARLY_ID(:,18);

Vearly_PMOS_2u280n = -1.*PMOS_Early(:,19);
Id_PMOS_2u280n = PMOS_EARLY_ID(:,19);

Vearly_PMOS_2u400n = -1.*PMOS_Early(:,20);
Id_PMOS_2u400n = PMOS_EARLY_ID(:,20);

Vearly_PMOS_2u520n = -1.*PMOS_Early(:,21);
Id_PMOS_2u520n = PMOS_EARLY_ID(:,21);

Vearly_PMOS_2u640n = -1.*PMOS_Early(:,22);
Id_PMOS_2u640n = PMOS_EARLY_ID(:,22);

Vearly_PMOS_2u760n = -1.*PMOS_Early(:,23);
Id_PMOS_2u760n = PMOS_EARLY_ID(:,23);

Vearly_PMOS_2u880n = -1.*PMOS_Early(:,24);
Id_PMOS_2u880n = PMOS_EARLY_ID(:,24);

Vearly_PMOS_3u = -1.*PMOS_Early(:,25);
Id_PMOS_3u = PMOS_EARLY_ID(:,25);

Vearly_PMOS_3u120n = -1.*PMOS_Early(:,26);
Id_PMOS_3u120n = PMOS_EARLY_ID(:,26);

Vearly_PMOS_3u240n = -1.*PMOS_Early(:,27);
Id_PMOS_3u240n = PMOS_EARLY_ID(:,27);

Vearly_PMOS_3u360n = -1.*PMOS_Early(:,28);
Id_PMOS_3u360n = PMOS_EARLY_ID(:,28);

Vearly_PMOS_3u480n = -1.*PMOS_Early(:,29);
Id_PMOS_3u480n = PMOS_EARLY_ID(:,29);

Vearly_PMOS_3u600n = -1.*PMOS_Early(:,30);
Id_PMOS_3u600n = PMOS_EARLY_ID(:,30);

Vearly_PMOS_3u720n = -1.*PMOS_Early(:,31);
Id_PMOS_3u720n = PMOS_EARLY_ID(:,31);

Vearly_PMOS_3u840n = -1.*PMOS_Early(:,32);
Id_PMOS_3u840n = PMOS_EARLY_ID(:,32);

Vearly_PMOS_4u = -1.*PMOS_Early(:,33);
Id_PMOS_4u = PMOS_EARLY_ID(:,33);


%Transistor_12
figure(5)
%hold on
IC_PMOS = Id_PMOS_840n/I0_PMOS;
grid on
loglog(IC_PMOS,Vearly_PMOS_840n)
title('Transistor 1,2')

%Transistor_5
figure(6)
IC_PMOS = Id_PMOS_1u200n/I0_PMOS;
grid on
loglog(IC_PMOS,Vearly_PMOS_1u200n)
title('Transistor 5')

figure(7)
IC_PMOS = Id_PMOS_1u200n/I0_PMOS;
grid on
loglog(IC_PMOS,Vearly_PMOS_1u200n)
title('Transistor 6')

% loglog(Id_PMOS_240n,Vearly_PMOS_240n)
% loglog(Id_PMOS_360n,Vearly_PMOS_360n)
% loglog(Id_PMOS_480n,Vearly_PMOS_480n)
% loglog(Id_PMOS_600n,Vearly_PMOS_600n)
% loglog(Id_PMOS_720n,Vearly_PMOS_720n)
% loglog(Id_PMOS_840n,Vearly_PMOS_840n)
% loglog(Id_PMOS_960n,Vearly_PMOS_960n)
% loglog(Id_PMOS_1u80n,Vearly_PMOS_1u80n)
% loglog(Id_PMOS_1u200n,Vearly_PMOS_1u200n)
% loglog(Id_PMOS_1u320n,Vearly_PMOS_1u320n)
% loglog(Id_PMOS_1u440n,Vearly_PMOS_1u440n)
% loglog(Id_PMOS_1u560n,Vearly_PMOS_1u560n)
% loglog(Id_PMOS_1u680n,Vearly_PMOS_1u680n)
% loglog(Id_PMOS_1u800n,Vearly_PMOS_1u800n)
% loglog(Id_PMOS_2u,Vearly_PMOS_2u)
%hold off


%%
%Análise de pequenos sinais para se calcular os ganhos
GBW = 20*10^3;

%%%Cálculo da tensão térmica
Ut = 25.8*10^(-3); %V
%%%Arbitração de n
n=1.35;
%%%IC utilizado no transistor 1
IC=0.3;
Vdsat1 = 2*Ut*sqrt(IC+0.5)+3*Ut; %Super importante!!!!!!!!!

gmp1 = (Id_1/2)/(n*Ut*sqrt(IC+0.25)+0.5)
gds_34 = (Id_1/2)/4.301 %Por observação
rds_34 = 1/gds_34;
gds_12 = (Id_1/2)/16.75%Por observação
rds_12 = 1/gds_12;
gds_5 = Id_1/22.51%Por observação
rds_5 = 1/gds_5;

%gds_34=320*10^(-9);

rout_1 = (rds_12*rds_34)/(rds_34+rds_12)

Av1 = gmp1*rout_1

%Par diferencial
%%%Rout1 = rds4//(rdsPF*(gmp1/(1/gds2))

Rout1 = (rds_34*(rds_5*(gmp1/gds_12)))/(rds_34+(rds_5*(gmp1/gds_12)))


Gm_eff = gmp1
Av_CM = (gmp1*Rout1)/(1+gmp1*(1/gds_5))

CMRR = 20*log10(Av1/Av_CM)


Cc = gmp1/(2*pi*GBW)


%M7 -> L=1.2u; W=7.7u; IC=1
%M6 -> L=1.2u; W=11u; IC=6.75
IC = 1;

Vdsat = 2*Ut*sqrt(IC+0.5)+3*Ut

gmp7 = Id_2/(n*Ut*sqrt(IC+0.25)+0.5);
gds7=Id_2/8.16; %por observação
rds7=1/gds7;
gds6=Id_2/24.74; %por observação
rds6=1/gds6

rout2=(rds6*rds7)/(rds6+rds7);

Av2=gmp7*rout2


Av_total=Av1*Av2

Ganho_dB=20*log10(Av_total)
