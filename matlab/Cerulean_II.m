% Cerulean

cerulean = [0.5*10^(-17) 31.697641 ; 0.5*10^(-16) 30.010946 ; 0.5*10^(-15) 27.515776 ; 0.5*10^(-14) 21.487556 ; 0.5*10^(-13) 17.69022 ; 0.5*10^(-12) 13.988193] ;

%  %CERULEAN
%  [cubicCoefC,statsC,ctrC] = polyfit(log(cerulean(:,1)),cerulean(:,2),2) ;
%  cubicFitC = polyval(cubicCoefC,log(cerulean(:,1)),[],ctrC);
%  figure; plot(log(cerulean(:,1)),cubicFitC,'r-',log(cerulean(:,1)),cerulean(:,2),'<k', 'MarkerSize',7,'MarkerFaceColor',[0 1 1],'LineWidth',2)
%  xlabel('massa [g]'); ylabel('CT');

C_RY41 = [27.963337+A(1) 27.602585+A(1) 27.602585+A(1) 21.52707+A(2) 22.195679+A(2) 21.431292+A(2) 22.18691+A(3) 22.072508+A(3) 22.330772+A(3)] ;
C_RY41PRE = [25.171381+A_pre(1) 25.042969+A_pre(1) 25.126356+A_pre(1) 20.906298+A_pre(2) 20.74485+A_pre(2) 20.808542+A_pre(2) 26.968414+A_pre(3) 26.680132+A_pre(3) 27.016623+A_pre(3)] ;

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

% Per avere la massa iniziale, contenuta in 500 nanogrammi di RNA, moltiplico
% per le diluizioni del cDNA (1/5) -> *5 e poi 1/2 -> *2

for i=1:9  % conto per cerulean e orange
    massa_iniziale(i) = massa_prePCR(i)*2*5 ; 
    massa_iniziale_PRE(i) = massa_prePCR_PRE(i)*2*5 ;
end

PM = 118800 ; % PM cerulean

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

molecole_cerulean = [numero_molecole_per_cellula ; ...
    numero_molecole_per_cellula_PRE] ;

CERULEAN_II = [mean(molecole_cerulean(1,1:3)) mean(molecole_cerulean(1,4:6)) mean(molecole_cerulean(1,7:9)) ; mean(molecole_cerulean(2,1:3)) mean(molecole_cerulean(2,4:6)) mean(molecole_cerulean(2,7:9))] ;
CERULEAN_std_II = [std(molecole_cerulean(1,1:3)) std(molecole_cerulean(1,4:6)) std(molecole_cerulean(1,7:9)) ; std(molecole_cerulean(2,1:3)) std(molecole_cerulean(2,4:6)) std(molecole_cerulean(2,7:9))] ;

