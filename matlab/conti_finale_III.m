
% sample order: HIGH - MEDIUM - LOW

s18 = [21.749464 19.586191 32.943645] ; 
s18_pre = [25.784464 16.353994 18.693254];

S18 = [21.749464 19.586191 32.943645 25.784464 16.353994 18.693254] ;

minnie = min(S18) ;

clear A A_neg A_pre

for i=1:3
    A(i) = minnie-s18(i) ;
    A_pre(i) = minnie-s18_pre(i) ;  
end

massa_RNA_totale = [36.3*30 36.3*30 36.3*30 57*30 57*30 57*30 37*30 37*30 37*30] ; % massa contenuta in RNA totale estratto (nanogrammi)
massa_RNA_totale_PRE = [25.5*30 25.5*30 25.5*30 61.4*30 61.4*30 61.4*30 51.5*30 51.5*30 51.5*30] ;

% l'ordine del numero di cellule è HIGH - MEDIUM - LOW

numero_cellule = [126127 126127 126127 258461 258461 258461 488849 488849 488849] ;
numero_cellule_PRE = [83369 83369 83369 211679 211679 211679 300243 300243 300243] ;
 
correzione = [36.3/126127 57/258461 37/488849 33.7/157105 54.7/288058 95.1/403472 61.4/211679 51.5/300243] ;

max_corr = max(correzione) ;

correzione_rna_perso = [max_corr/(36.3/126127) max_corr/(36.3/126127) max_corr/(36.3/126127) ...
    max_corr/(57/258461) max_corr/(57/258461) max_corr/(57/258461) ...
    max_corr/(37/488849) max_corr/(37/488849) max_corr/(37/488849)] ;
correzione_rna_perso_PRE = [max_corr/(25.5/83369) max_corr/(25.5/83369) max_corr/(25.5/83369) ...
    max_corr/(61.4/211679) max_corr/(61.4/211679) max_corr/(61.4/211679) ...
    max_corr/(51.5/300243) max_corr/(51.5/300243) max_corr/(51.5/300243)] ;

Cerulean_III
Cherry_III
Yfp_III
