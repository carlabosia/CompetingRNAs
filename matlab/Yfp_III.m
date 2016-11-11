% Yfp

yfp = [0.5*10^(-17) 33.756916 ; 0.5*10^(-16) 30.369036 ; 0.5*10^(-15) 23.974335 ; 0.5*10^(-14) 19.974829 ; 0.5*10^(-13) 16.946165 ; 0.5*10^(-12) 13.366947] ;

% % YFP
% [cubicCoefY,statsY,ctrY] = polyfit(log(yfp(:,1)),yfp(:,2),2) ;
% cubicFitY = polyval(cubicCoefY,log(yfp(:,1)),[],ctrY);
% figure; plot(log(yfp(:,1)),cubicFitY,'r-',log(yfp(:,1)),yfp(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[1 1 0],'LineWidth',2)
% xlabel('massa [g]'); ylabel('CT');

Y_RY41 = [20.805595+A(1) 20.660673+A(1) 20.8026+A(1) 20.41863+A(2) 20.389132+A(2) 20.341545+A(2) 34.547253+A(3) 33.992214+A(3) 34.042892+A(3)] ;
Y_RY41PRE = [21.929329+A_pre(1) 21.593794+A_pre(1) 21.804052+A_pre(1) 22.792189+A_pre(2) 21.867678+A_pre(2) 21.654968+A_pre(2) 23.654741+A_pre(3) 23.329863+A_pre(3) 23.460886+A_pre(3)] ;

[cubicCoefY] = polyfit(log(yfp(:,1)),yfp(:,2),2) ;

RY41 = Y_RY41 ;
RY41PRE = Y_RY41PRE ; 

cubicCoef = cubicCoefY ;

clear pippo pippo_NEG pippo_PRE

for i=1:9
    pippo(i,:) = roots(cubicCoef - [0 0 RY41(i)]) ;
    pippo_PRE(i,:) = roots(cubicCoef - [0 0 RY41PRE(i)]) ;
    pluto(i) = find(pippo(i,:)>-45 & pippo(i,:)<-26) ;
    pluto_PRE(i) = find(pippo_PRE(i,:)>-45 & pippo_PRE(i,:)<-26) ; 
    massa_prePCR(i) = exp(pippo(i,pluto(i))) ;
    massa_prePCR_PRE(i) = exp(pippo_PRE(i,pluto_PRE(i))) ;    
end

% Per avere la massa iniziale, contenuta in 250 nanogrammi di RNA, moltiplico
% per le diluizioni del cDNA (1/50) -> *50 e poi 1/2 -> *2

for i=1:9  % conto per cherry e yfp
    massa_iniziale(i) = massa_prePCR(i)*2*50 ;
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*50 ;
end

PM = 93060 ; % PM yfp

for i=1:9
    massa_cDNA_totale(i) = massa_iniziale(i)*(massa_RNA_totale(i)/250)*correzione_rna_perso(i) ; % massa_iniziale*nanogrammi_totali
    massa_cDNA_totale_PRE(i) = massa_iniziale_PRE(i)*(massa_RNA_totale_PRE(i)/250)*correzione_rna_perso_PRE(i) ; % massa_iniziale*nanogrammi_totali

    numero_moli(i) = massa_cDNA_totale(i)/PM;
    numero_moli_PRE(i) = massa_cDNA_totale_PRE(i)/PM ;
    
    numero_molecole_totale(i) = numero_moli(i)*(6.022*10^23) ;
    numero_molecole_totale_PRE(i) = numero_moli_PRE(i)*(6.022*10^23) ;
    
    numero_molecole_per_cellula(i) = numero_molecole_totale(i)/numero_cellule(i) ;
    numero_molecole_per_cellula_PRE(i) = numero_molecole_totale_PRE(i)/numero_cellule_PRE(i) ;
end

molecole_yfp = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

YFP_III = [mean(molecole_yfp(1,1:3)) mean(molecole_yfp(1,4:6)) mean(molecole_yfp(1,7:9)) ; mean(molecole_yfp(2,1:3)) mean(molecole_yfp(2,4:6)) mean(molecole_yfp(2,7:9))] ;
YFP_std_III = [std(molecole_yfp(1,1:3)) std(molecole_yfp(1,4:6)) std(molecole_yfp(1,7:9)) ; std(molecole_yfp(2,1:3)) std(molecole_yfp(2,4:6)) std(molecole_yfp(2,7:9))] ;

