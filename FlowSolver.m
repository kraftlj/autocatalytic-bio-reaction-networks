function [Thiols,X]=FlowSolver(X0,k,A,simTime)
%Oct 4, 2014

% Initialize Indices

AlaSEt =1;
CSH =2;
Cmpd5=3;
EtSH=4;
CSSC=5;
Cmpd6=6;
Cmpd7=7;
Mal=8;

chemsyst=@(t,X)BistabilityEquationsFlow(t,X,k,A);
tspan=linspace(0,simTime,300);
[T,X]=ode45(chemsyst,tspan,X0);
Thiols=X(end,CSH)+X(end,Cmpd5)+X(end,EtSH);
end

