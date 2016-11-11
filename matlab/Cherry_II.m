% Cherry

cherry = [0.5*10^(-17) 32.78245 ; 0.5*10^(-16) 29.09606 ; 0.5*10^(-15) 23.545156 ; 0.5*10^(-14) 20.3741995 ; 0.5*10^(-13) 17.39074 ; 0.5*10^(-12) 14.1602545] ;

% % CHERRY
% [cubicCoefR,statsR,ctrR] = polyfit(log(cherry(:,1)),cherry(:,2),2) ;
% cubicFitR = polyval(cubicCoefR,log(cherry(:,1)),[],ctrR);
% figure; plot(log(cherry(:,1)),cubicFitR,'r-',log(cherry(:,1)),cherry(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[1 0 0],'LineWidth',2)
% xlabel('massa [g]'); ylabel('CT');

R_RY41 = [29.594206+A(1) 29.409613+A(1) 29.350616+A(1) 22.242758+A(2) 22.17024+A(2) 22.222555+A(2) 22.255655+A(3) 22.19835+A(3) 22.178944+A(3)] ;
R_RY41PRE = [28.994684+A_pre(1) 29.246027+A_pre(1) 29.09437+A_pre(1) 23.424353+A_pre(2) 23.52039+A_pre(2) 23.5354+A_pre(2) 29.786695+A_pre(3) 29.93465+A_pre(3) 29.664057+A_pre(3)] ;

[cubicCoefR] = polyfit(log(cherry(:,1)),cherry(:,2),2) ;

RY41 = R_RY41 ;
RY41PRE = R_RY41PRE ; 

cubicCoef = cubicCoefR ;

for i=1:9
    pippo(i,:) = roots(cubicCoef - [0 0 RY41(i)]) ;
    pippo_PRE(i,:) = roots(cubicCoef - [0 0 RY41PRE(i)]) ;
    pluto(i) = find(pippo(i,:)>-45 & pippo(i,:)<-28) ;
    pluto_PRE(i) = find(pippo_PRE(i,:)>-45 & pippo_PRE(i,:)<-28) ; 
    massa_prePCR(i) = exp(pippo(i,pluto(i))) ;
    massa_prePCR_PRE(i) = exp(pippo_PRE(i,pluto_PRE(i))) ;    
end

% Per avere la massa iniziale, contenuta in 500 nanogrammi di RNA, moltiplico
% per le diluizioni del cDNA (1/50) -> *50 e poi 1/2 -> *2

for i=1:9 
    massa_iniziale(i) = massa_prePCR(i)*2*50 ;
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*50 ;
end

PM = 139260 ; % PM cherry

for i=1:9
    massa_cDNA_totale(i) = massa_iniziale(i)*(massa_RNA_totale(i)/500)*correzione_rna_perso(i) ; % massa_iniziale*nanogrammi_totali
    massa_cDNA_totale_PRE(i) = massa_iniziale_PRE(i)*(massa_RNA_totale_PRE(i)/500)*correzione_rna_perso_PRE(i) ; % massa_iniziale*nanogrammi_totali

    numero_moli(i) = massa_cDNA_totale(i)/PM;
    numero_moli_PRE(i) = massa_cDNA_totale_PRE(i)/PM ;
    
    numero_molecole_totale(i) = numero_moli(i)*(6.022*10^23) ;
    numero_molecole_totale_PRE(i) = numero_moli_PRE(i)*(6.022*10^23) ;
    
    numero_molecole_per_cellula(i) = numero_molecole_totale(i)/numero_cellule(i) ;
    numero_molecole_per_cellula_PRE(i) = numero_molecole_totale_PRE(i)/numero_cellule_PRE(i) ;
end

molecole_cherry = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

CHERRY_II = [mean(molecole_cherry(1,1:3)) mean(molecole_cherry(1,4:6)) mean(molecole_cherry(1,7:9)) ; mean(molecole_cherry(2,1:3)) mean(molecole_cherry(2,4:6)) mean(molecole_cherry(2,7:9))] ;
CHERRY_std_II = [std(molecole_cherry(1,1:3)) std(molecole_cherry(1,4:6)) std(molecole_cherry(1,7:9)) ; std(molecole_cherry(2,1:3)) std(molecole_cherry(2,4:6)) std(molecole_cherry(2,7:9))] ;

