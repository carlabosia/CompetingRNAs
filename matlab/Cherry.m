% Cherry

cherry = [0.5*10^(-17) 32.026333 ; 0.5*10^(-16) 28.73617 ; 0.5*10^(-15) 23.307676 ; 0.5*10^(-14) 20.34642 ; 0.5*10^(-13) 17.143194 ; 0.5*10^(-12) 13.296725] ;

% % CHERRY
% [cubicCoefR,statsR,ctrR] = polyfit(log(cherry(:,1)),cherry(:,2),2) ;
% cubicFitR = polyval(cubicCoefR,log(cherry(:,1)),[],ctrR);
% figure; plot(log(cherry(:,1)),cubicFitR,'r-',log(cherry(:,1)),cherry(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[1 0 0],'LineWidth',2)
% hold on
% xlabel('massa [g]'); ylabel('CT');

R_RY41 = [19.34519+A(1) 19.155756+A(1) 19.199903+A(1) 20.402935+A(2) 20.361588+A(2) 20.470684+A(2) 18.871069+A(3) 18.670658+A(3) 18.690641+A(3)] ;
R_RY41PRE = [20.454283+A_pre(1) 20.264828+A_pre(1) 20.341682+A_pre(1) 23.731686+A_pre(2) 23.599752+A_pre(2) 23.801167+A_pre(2) 23.31532+A_pre(3) 23.136175+A_pre(3) 23.105368+A_pre(3)] ;

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

% Per avere la massa iniziale, contenuta in 1000 nanogrammi di RNA, moltiplico
% per le diluizioni del cDNA (1/50) -> *50 e poi 1/2 -> *2

for i=1:9  % conto per cherry e yfp
    massa_iniziale(i) = massa_prePCR(i)*2*50 ;
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*50 ;
end

PM = 139260 ; % PM cherry

for i=1:9
    massa_cDNA_totale(i) = massa_iniziale(i)*(massa_RNA_totale(i)/1000)*correzione_rna_perso(i) ; % massa_iniziale*nanogrammi_totali
    massa_cDNA_totale_PRE(i) = massa_iniziale_PRE(i)*(massa_RNA_totale_PRE(i)/1000)*correzione_rna_perso_PRE(i) ; % massa_iniziale*nanogrammi_totali

    numero_moli(i) = massa_cDNA_totale(i)/PM;
    numero_moli_PRE(i) = massa_cDNA_totale_PRE(i)/PM ;
    
    numero_molecole_totale(i) = numero_moli(i)*(6.022*10^23) ;
    numero_molecole_totale_PRE(i) = numero_moli_PRE(i)*(6.022*10^23) ;
    
    numero_molecole_per_cellula(i) = numero_molecole_totale(i)/numero_cellule(i) ;
    numero_molecole_per_cellula_PRE(i) = numero_molecole_totale_PRE(i)/numero_cellule_PRE(i) ;
end

molecole_cherry = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

CHERRY = [mean(molecole_cherry(1,1:3)) mean(molecole_cherry(1,4:6)) mean(molecole_cherry(1,7:9)) ; mean(molecole_cherry(2,1:3)) mean(molecole_cherry(2,4:6)) mean(molecole_cherry(2,7:9))] ;
CHERRY_std = [std(molecole_cherry(1,1:3)) std(molecole_cherry(1,4:6)) std(molecole_cherry(1,7:9)) ; std(molecole_cherry(2,1:3)) std(molecole_cherry(2,4:6)) std(molecole_cherry(2,7:9))] ;

