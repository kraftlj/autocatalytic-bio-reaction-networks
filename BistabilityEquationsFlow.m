function dX=BistabilityEquationsFlow(t,X,k,A)

% Initialize Indices
AlaSEt =1;
CSH =2;
Cmpd5=3;
EtSH=4;
CSSC=5;
Cmpd6=6;
Cmpd7=7;
Mal=8;
%Initialize indices for reactions
Mal_rxn=2;
Disulf_ex=3;
Ligation=4;
Hydrolysis=5;
%Initialize indices for flow/volume
FvV=1;
inputTE=2;
inputHSR=3;
inputCSSC=4;
inputMal=5;

dX=zeros(8,1);

dX(AlaSEt)=-k(Ligation)*X(AlaSEt)*X(CSH) -k(Hydrolysis)*X(AlaSEt) -A(FvV)*X(AlaSEt)...
        +A(inputTE)*A(FvV);
    
dX(CSH)=-k(Disulf_ex)*X(CSH)*X(Cmpd6) -k(Ligation)*X(AlaSEt)*X(CSH) -k(Mal_rxn)*X(CSH)*X(Mal) -k(Disulf_ex)*X(Cmpd7)*X(CSH) -A(FvV)*X(CSH)...
        +k(Disulf_ex)*X(CSSC)*X(EtSH) +k(Disulf_ex)*X(CSSC)*X(Cmpd5);
    
dX(Cmpd5)=-A(FvV)*X(Cmpd5) -k(Disulf_ex)*X(Cmpd5)*X(Cmpd6) -k(Mal_rxn)*X(Cmpd5)*X(Mal) -k(Disulf_ex)*X(Cmpd5)*X(CSSC)...
        +k(Disulf_ex)*X(Cmpd7)*X(EtSH) +k(Disulf_ex)*X(Cmpd7)*X(CSH) +k(Ligation)*X(AlaSEt)*X(CSH);

dX(EtSH)=-k(Disulf_ex)*X(CSSC)*X(EtSH) -k(Disulf_ex)*X(Cmpd7)*X(EtSH) -k(Mal_rxn)*X(EtSH)*X(Mal) -A(FvV)*X(EtSH)...
        +k(Disulf_ex)*X(CSH)*X(Cmpd6) +k(Ligation)*X(AlaSEt)*X(CSH) +A(FvV)*A(inputHSR) +k(Disulf_ex)*X(Cmpd5)*X(Cmpd6) +k(Hydrolysis)*X(AlaSEt);

dX(CSSC)=-k(Disulf_ex)*X(CSSC)*X(EtSH) -k(Disulf_ex)*X(Cmpd5)*X(CSSC) -A(FvV)*X(CSSC)...
        +k(Disulf_ex)*X(CSH)*X(Cmpd6) +k(Disulf_ex)*X(CSH)*X(Cmpd7) +A(FvV)*A(inputCSSC);

dX(Cmpd6)=-k(Disulf_ex)*X(CSH)*X(Cmpd6) -k(Disulf_ex)*X(Cmpd5)*X(Cmpd6) -A(FvV)*X(Cmpd6)...
        +k(Disulf_ex)*X(CSSC)*X(EtSH) +k(Disulf_ex)*X(Cmpd7)*X(EtSH);

dX(Cmpd7)=-A(FvV)*X(Cmpd7) -k(Disulf_ex)*X(Cmpd7)*X(EtSH) -k(Disulf_ex)*X(Cmpd7)*X(CSH)...
        +k(Disulf_ex)*X(Cmpd5)*X(Cmpd6) +k(Disulf_ex)*X(Cmpd5)*X(CSSC);

dX(Mal)=-A(FvV)*X(Mal) -k(Mal_rxn)*X(CSH)*X(Mal) -k(Mal_rxn)*X(EtSH)*X(Mal) -k(Mal_rxn)*X(Cmpd5)*X(Mal)...
        +A(FvV)*A(inputMal);
