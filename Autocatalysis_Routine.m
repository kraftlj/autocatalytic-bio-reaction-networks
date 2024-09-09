% Autocatalysis routine, Lewis J. Kraft, Whitesides Research Group, Harvard University.
clear all;
close all;

simTime=2000;

%Initiate variables
X0 = zeros(1,9);
% Initialize Indices
AlaSEt =1;
CSH =2;
Cmpd5=3;
EtSH=4;
CSSC=5;
Cmpd6=6;
Cmpd7=7;
Mal=8;
AAm=9;

%Initiate initial concentrations
X0(AlaSEt)=0.046;%AlaSEt initial concentration in reactor
X0(CSH)=0;%CSH initial concentration in reactor
X0(Cmpd5)=0;%Cmpd5 initial concentration in ractor
X0(EtSH)=0;%EtSH initial concentration in reactor
X0(CSSC)=0.046;%CSSC initial concentration in reactor
X0(Cmpd6)=0;%Cmpd6 initial concentration in reactor
X0(Cmpd7)=0;%Cmpd7 initial concentration in reactor
X0(Mal)=0;%Malemide initial concentration in reactor
X0(AAm)=0;%Acrylamide initial concentration in reactor

%Initialize indices for reactions
AAm_rxn=1;
Mal_rxn=2;
Disulf_ex=3;
Ligation=4;
Hydrolysis=5;

%Initiate rate constants
k(AAm_rxn)= 0.014; % EtSH + AAm -> inhibited %Cmpd5 + AAm -> inhibited %CSH + AAm -> inhibited
k(Mal_rxn)= 150; % EtSH + Mal -> inhibited %Cmpd5 + Mal -> inhibited %CSH + Mal -> inhibited
k(Disulf_ex) = 0.65; %CSSC + EtSH <-> CSH + Cmpd6 %Cmpd5 + Cmpd6 <-> Cmpd7 + EtSH %CSSC + Cmpd5 <-> CSH + Cmpd7
k(Ligation) = 0.411; %AlaSEt+CSH -> Cmpd5 + EtSH
k(Hydrolysis) = 7e-6; % AlaSEt -> ROH + EtSH

%Initialize indices for flow/volume
FvV=1;
inputTE=2;
inputHSR=3;
inputDCys=4;
inputMal=5;
inputAAm=6;

%Initiate constants
A(FvV)=0;%flow rate (total of two input streams with equal flow rates)
A(inputTE)=0;%Initial concentration of AlaSEt in input stream
A(inputHSR)=0;%Initial concentration of EtSH in input stream 
A(inputDCys)=0;%Initial concentration of DCYS in input stream 
A(inputMal)=0;%Initial concentration of Malemide in input stream
A(inputAAm)=0;%Initial concentration of AAm in input stream

chemsyst=@(t,X)OscillationEquationsFlow(t,X,k,A);
tspan=linspace(0,simTime,300);
[T,X]=ode45(chemsyst,tspan,X0);
NHCO = (X(:,Cmpd5)+X(:,Cmpd7))*1e3; % Convert to mM


plot(T,NHCO);
grid on;
xlabel('Time (s)');
ylabel('[NHCO]    (mM)');
