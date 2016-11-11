% Cerulean

cerulean = [0.5*10^(-17) 30.183733 ; 0.5*10^(-16) 29.010271 ; 0.5*10^(-15) 24.996275 ; 0.5*10^(-14) 21.274153 ; 0.5*10^(-13) 18.035173 ; 0.5*10^(-12) 13.210649] ;

% %CERULEAN
% [cubicCoefC,statsC,ctrC] = polyfit(log(cerulean(:,1)),cerulean(:,2),2) ;
% cubicFitC = polyval(cubicCoefC,log(cerulean(:,1)),[],ctrC);
% figure; plot(log(cerulean(:,1)),cubicFitC,'r-',log(cerulean(:,1)),cerulean(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[0 1 1],'LineWidth',2)
% xlabel('massa [g]'); ylabel('CT');

C_RY41 = [18.227121+A(1) 18.119673+A(1) 18.198444+A(1) 20.823072+A(2) 20.595587+A(2) 20.539675+A(2) 20.895905+A(3) 20.42778+A(3) 20.463964+A(3)] ;
C_RY41PRE = [18.736418+A_pre(1) 18.656206+A_pre(1) 18.632595+A_pre(1) 19.197348+A_pre(2) 19.014267+A_pre(2) 19.096895+A_pre(2) 21.485899+A_pre(3) 21.176214+A_pre(3) 20.989914+A_pre(3)] ;

[cubicCoefC] = polyfit(log(cerulean(:,1)),cerulean(:,2),2) ;

RY41 = C_RY41 ;
RY41PRE = C_RY41PRE ; 

cubicCoef = cubicCoefC ;

for i=1:9
    pippo(i,:) = roots(cubicCoef - [0 0 RY41(i)]) ;
    pippo_PRE(i,:) = roots(cubicCoef - [0 0 RY41PRE(i)]) ;
    pluto(i) = find(pippo(i,:)>-45 & pippo(i,:)<-28) ;
    pluto_PRE(i) = find(pippo_PRE(i,:)>-45 & pippo_PRE(i,:)<-28) ; 
    massa_prePCR(i) = exp(pippo(i,pluto(i))) ;
    massa_prePCR_PRE(i) = exp(pippo_PRE(i,pluto_PRE(i))) ;    
end

% Per avere la massa iniziale, contenuta in 1000 nanogrammi di RNA, moltiplico
% per le diluizioni del cDNA (1/5) -> *5 e poi 1/2 -> *2

for i=1:9  % conto per cerulean e orange
    massa_iniziale(i) = massa_prePCR(i)*2*5 ; 
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*5 ;
end

PM = 118800 ; % PM cerulean

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

molecole_cerulean = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

CERULEAN = [mean(molecole_cerulean(1,1:3)) mean(molecole_cerulean(1,4:6)) mean(molecole_cerulean(1,7:9)) ; mean(molecole_cerulean(2,1:3)) mean(molecole_cerulean(2,4:6)) mean(molecole_cerulean(2,7:9))] ;
CERULEAN_std = [std(molecole_cerulean(1,1:3)) std(molecole_cerulean(1,4:6)) std(molecole_cerulean(1,7:9)) ;  std(molecole_cerulean(2,1:3)) std(molecole_cerulean(2,4:6)) std(molecole_cerulean(2,7:9))] ;

