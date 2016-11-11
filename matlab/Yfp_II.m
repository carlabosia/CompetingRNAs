% Yfp

yfp = [0.5*10^(-17) 32.694695 ; 0.5*10^(-16) 29.017807 ; 0.5*10^(-15) 23.467371 ; 0.5*10^(-14) 19.518457 ; 0.5*10^(-13) 16.612638 ; 0.5*10^(-12) 12.30845] ;

% % YFP
% [cubicCoefY,statsY,ctrY] = polyfit(log(yfp(:,1)),yfp(:,2),2) ;
% cubicFitY = polyval(cubicCoefY,log(yfp(:,1)),[],ctrY);
% figure; plot(log(yfp(:,1)),cubicFitY,'r-',log(yfp(:,1)),yfp(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[1 1 0],'LineWidth',2)
% xlabel('massa [g]'); ylabel('CT');

Y_RY41 = [25.731197+A(1) 25.846992+A(1) 25.902403+A(1) 21.842686+A(2) 22.412416+A(2) 21.818607+A(2) 22.653858+A(3) 22.690063+A(3) 22.401299+A(3)] ;
Y_RY41PRE = [25.598999+A_pre(1) 25.446535+A_pre(1) 25.791355+A_pre(1) 21.458939+A_pre(2) 21.473625+A_pre(2) 21.481153+A_pre(2) 27.419674+A_pre(3) 27.412123+A_pre(3) 27.601395+A_pre(3)] ;

[cubicCoefY] = polyfit(log(yfp(:,1)),yfp(:,2),2) ;

RY41 = Y_RY41 ;
RY41PRE = Y_RY41PRE ; 

cubicCoef = cubicCoefY ;

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

for i=1:9  % conto per cherry e yfp
    massa_iniziale(i) = massa_prePCR(i)*2*50 ;
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*50 ;
end

PM = 93060 ; % PM yfp

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

molecole_yfp = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

YFP_II = [mean(molecole_yfp(1,1:3)) mean(molecole_yfp(1,4:6)) mean(molecole_yfp(1,7:9)) ; mean(molecole_yfp(2,1:3)) mean(molecole_yfp(2,4:6)) mean(molecole_yfp(2,7:9))] ;
YFP_std_II = [std(molecole_yfp(1,1:3)) std(molecole_yfp(1,4:6)) std(molecole_yfp(1,7:9)) ; std(molecole_yfp(2,1:3)) std(molecole_yfp(2,4:6)) std(molecole_yfp(2,7:9))] ;

