% Bistability routine, Lewis J. Kraft, Whitesides Research Group, Harvard University.
clear all;
close all;

simTime=5*60*60;

%Initiate variables
X0 = zeros(1,8);
% Initialize Indices
AlaSEt =1;
CSH =2;
Cmpd5=3;
EtSH=4;
CSSC=5;
Cmpd6=6;
Cmpd7=7;
Mal=8;

%Initiate initial concentrations
X0(AlaSEt)=0.046;%AlaSEt initial concentration in reactor
X0(CSH)=0;%CSH initial concentration in reactor
X0(Cmpd5)=0;%Cmpd5 initial concentration in ractor
X0(EtSH)=0;%EtSH initial concentration in reactor
X0(CSSC)=0.092;%CSSC initial concentration in reactor
X0(Cmpd6)=0;%Cmpd6 initial concentration in reactor
X0(Cmpd7)=0;%Cmpd7 initial concentration in reactor
X0(Mal)=3.47e-3;%Malemide initial concentration in reactor

%Initialize indices for reactions
Mal_rxn=2;
Disulf_ex=3;
Ligation=4;
Hydrolysis=5;

%Initiate rate constants
k(1)= 0; % This was for the acrylamide in the oscillatory model
k(Mal_rxn)= 150; % EtSH + Mal -> inhibited %Cmpd5 + Mal -> inhibited %CSH + Mal -> inhibited
k(Disulf_ex) = 0.65; %CSSC + EtSH <-> CSH + Cmpd6 %Cmpd5 + Cmpd6 <-> Cmpd7 + EtSH %CSSC + Cmpd5 <-> CSH + Cmpd7
k(Ligation) = 0.411; %AlaSEt+CSH -> Cmpd5 + EtSH
k(Hydrolysis) = 9.26e-6; % AlaSEt -> ROH + EtSH

%Initialize indices for flow/volume
FvV=1;
inputTE=2;
inputHSR=3;
inputDCys=4;
inputMal=5;

%Initiate constants
flows=[linspace(0.01,0.0001,20), linspace(0.0001,0.01,20)];
num=length(flows);
A(inputTE)=0.046;%Initial concentration of AlaSEt in input stream
A(inputHSR)=0;%Initial concentration of EtSH in input stream
A(inputDCys)=0.092;%Initial concentration of DCYS in input stream
A(inputMal)=0.00347;%Initial concentration of Malemide in input stream
Thiols=zeros(num,1);
A(FvV)=flows(1);
[Thiols(1),X]=FlowSolver(X0,k,A,simTime);
for j=2:num
    X0=X(end,:);
    A(FvV)=flows(j);%flow rate (total of two input streams with equal flow rates)
    [Thiols(j),X]=FlowSolver(X0,k,A,simTime);
end

figure
plot(flows,Thiols.*1000,'ko')
ylim([-5 50])
xlim([0 1.1e-2])
xlabel('Flow/Volume (s^{-1})')
ylabel('[Thiols]    (mM)')
set(gca,'FontSize',6)
set(gcf,'Position',[935 667 238 195])
