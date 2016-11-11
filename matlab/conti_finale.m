
% sample order: HIGH - MEDIUM - LOW

s18 = [14.0590 14.6532 12.7058] ;
s18_pre = [15.7591 14.6695 13.6044];

minnie = min(s18) ;

for i=1:3
    A(i) = minnie-s18(i) ;
    A_pre(i) = minnie-s18_pre(i) ;  
end

massa_RNA_totale = [181*30 181*30 181*30 214*30 214*30 214*30 207*30 207*30 207*30] ; % massa contenuta in RNA totale estratto (nanogrammi)
massa_RNA_totale_PRE = [119*30 119*30 119*30 163*30 163*30 163*30 135*30 135*30 135*30] ;

% l'ordine del numero di cellule è HIGH - MEDIUM - LOW

numero_cellule = [701405 701405 701405 736181 736181 736181 906714 906714 906714] ;
numero_cellule_PRE = [400972 400972 400972 579752 579752 579752 528356 528356 528356] ;

correzione = [181/701405 214/736181 207/906714 150.5/500263 151/505278 190.5/790984 119/400972 163/579752 135/528356] ;

max_corr = max(correzione) ;

correzione_rna_perso = [max_corr/(181/701405) max_corr/(181/701405) max_corr/(181/701405) ...
    max_corr/(214/736181) max_corr/(214/736181) max_corr/(214/736181) ...
    max_corr/(207/906714) max_corr/(207/906714) max_corr/(207/906714)] ;
correzione_rna_perso_PRE = [max_corr/(119/400972) max_corr/(119/400972) max_corr/(119/400972) ...
    max_corr/(163/579752) max_corr/(163/579752) max_corr/(163/579752) ...
    max_corr/(135/528356) max_corr/(135/528356) max_corr/(135/528356)] ;

Cerulean
Cherry
Yfp


