
miR = [22.661 22.497 22.224 ; 21.079 20.989 20.939 ; 23.507 23.329 23.377 ; 23.082 23.028 23.023 ; 22.624 22.529 22.419 ; ...
    21.389 21.396 21.404 ; 21.422 21.426 21.254 ; 20.072 20.136 20.022 ; 19.000 18.963 18.774 ; 18.188 18.183 18.232 ; 17.371 17.207 17.191 ; ...
    16.340 16.182 15.161] ;

vuoto = [24.126 24.059 24.137] ;
pre_miR = [16.872 16.722 16.789] ;

u44 = [21.500 21.317 21.454 ; 21.371 21.304 21.371 ; 21.675 21.547 21.494 ; 22.005 21.532 22.028 ; ...
    21.507 21.410 21.399 ; 21.199 21.157 21.080 ; 21.764 21.676 21.659 ; 21.705 21.572 22.001 ; 21.590 21.717 21.409 ; ...
    21.643 21.600 21.569 ; 21.584 21.555 21.575 ; 21.858 21.627 21.698] ;

vuoto_u44 = [21.783 21.418 21.405] ;
pre_miR_u44 = [21.420 21.700 21.593] ;

for i=1:12
    miR_mean(i) = mean(miR(i,:)) ;
    u44_mean(i) = mean(u44(i,:)) ;
end

molecole_per_cell = [78 156 312 625 1250 2500 5000 10000 20000 40000 80000 160000] ;

max_u44 = max(u44_mean) ;
u44_corretto(1:12) = max_u44 - u44_mean(1:12) ;

for i=1:12
    miR_corretto(i,1:3) = miR(i,1:3) + u44_corretto(i) ;
end

miR_mean_corretto = miR_mean + u44_corretto ;

for i=1:12
    miR_corretto_mean(i) = mean(miR_corretto(i,:)) ;
end

fit_cT = polyfit(log(molecole_per_cell(5:end)),miR_corretto_mean(5:end),1) ; 
x1 = linspace(7.2,12.5);
y1 = polyval(fit_cT,x1);

fit_cT2 = polyfit(log(molecole_per_cell([1 4 5])),miR_corretto_mean([1 4 5]),1) ; 
x2 = linspace(4,7.25);
y2 = polyval(fit_cT2,x2);


figure; 
plot(log(molecole_per_cell([1 4:end])), miR_corretto_mean([1 4:end]),'ob','linewidth',3)
hold on
plot(x1,y1,'r--','linewidth',2)
plot(x2,y2,'r--','linewidth',2)
hold off
xlim([4 13])
ylim([15 25])

vuoto_u44_mean = mean(vuoto_u44) ;
vuoto_u44_corretto = max_u44 - vuoto_u44_mean ;

vuoto_corretto(1:3) = vuoto(1:3) + vuoto_u44_corretto ;
vuoto_mean = mean(vuoto_corretto) ;


