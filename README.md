# CompetingRNAs

data and scripts


This repository contains the Matlab scripts and the RealTime raw data
presented in the paper [Quantitative study of crossregulation, noise
and synchronization between microRNA targets in single
cells](https://arxiv.org/abs/1503.06696) by Carla Bosia, Francesco
Sgrò, Laura Conti, Carlo Baldassi, Federica Cavallo, Ferdinando Di
Cunto, Emilia Turco, Andrea Pagnani, Riccardo Zecchina.

Installation
------------

To create a local copy of this repository please

   `git clone https://github.com/carlabosia/CompetingRNAs`

Matlab script
-------------

The matlab scripts can be run without external parameters and will
return the mean molecule amount for the exogenous transcripts corresponding to the 
CTs present in the RealTime data directory in this repository for three different 
experiments:

* "conti_finale.m" calls Yfp.m, Cherry.m and Cerulean.m (I experiment)
* "conti_finale_II.m" calls Yfp_II.m, Cherry_II.m and Cerulean_II.m (II experiment)
* "conti_finale_III.m" calls Yfp_III.m, Cherry_III.m and Cerulean_III.m (III experiment)

Finally:

* "real_time_mirna.m" evaluates the number of endogenous miR-20a from
  the CTs from the corresponding real time experiment.


RealTime data
-------------
The RealTime data are in following four files: 

* "24_07_14_quantificazione_campioni_sortati_excel.xls" -> I experiment
* "quantificazione_fluorofori_021014.xls" -> II experiment
* "quantificazione_fluorofori_161014.txt" -> III experiment
* "080716_analisi_miR20a_t.xls" -> RealTime for the quantification of miR-20a






