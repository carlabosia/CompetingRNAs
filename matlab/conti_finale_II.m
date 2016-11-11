
% sample order: HIGH - MEDIUM - LOW

s18 = [28.060469 18.046976 15.104772] ; 
s18_pre = [23.175064 17.841488 20.385714];

S18 = [28.060469 18.046976 15.104772 23.175064 17.841488 20.385714] ;

minnie = min(S18) ;

for i=1:3
    A(i) = minnie-s18(i) ;
    A_pre(i) = minnie-s18_pre(i) ;  
end

massa_RNA_totale = [54.3*30 54.3*30 54.3*30 79*30 79*30 79*30 110.6*30 110.6*30 110.6*30] ; % massa contenuta in RNA totale estratto (nanogrammi)
massa_RNA_totale_PRE = [16.3*30 16.3*30 16.3*30 27*30 27*30 27*30 20.5*30 20.5*30 20.5*30] ;

% l'ordine del numero di cellule è HIGH - MEDIUM - LOW

numero_cellule = [197791 197791 197791 369284 369284 369284 556726 556726 556726] ;
numero_cellule_PRE = [74389 74389 74389 166527 166527 166527 156415 156415 156415] ;

correzione = [54.3/197791 79/369284 110.6/556726 25/63185 31.2/111932 15.3/100336 16.3/74389 27/166527 20.5/156415] ;

max_corr = max(correzione) ;

correzione_rna_perso = [max_corr/(54.3/197791) max_corr/(54.3/197791) max_corr/(54.3/197791) ...
    max_corr/(79/369284) max_corr/(79/369284) max_corr/(79/369284) ...
    max_corr/(110.6/556726) max_corr/(110.6/556726) max_corr/(110.6/556726)] ;
correzione_rna_perso_PRE = [max_corr/(16.3/74389) max_corr/(16.3/74389) max_corr/(16.3/74389) ...
    max_corr/(27/166527) max_corr/(27/166527) max_corr/(27/166527) ...
    max_corr/(20.5/156415) max_corr/(20.5/156415) max_corr/(20.5/156415)] ;

Cerulean_II
Cherry_II
Yfp_II

