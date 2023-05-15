# BNI_N_Cycling
 
## This repository contains soil datasets associated with Schaedel et al. (in review) : Nitrogen cycling functional gene abundance and potential activity in novel perennial forage cropping systems in Rwanda
### Data and processing scripts for conventional and high-throuphut qPCR can be found at https://github.com/schaedem/NiCE_Chip_Paper

The goal of this project was to:

1. Test BNI under field conditions
2. Determine the impact of a legume intercrop on BNI and soil Fertility
3. Identify genetic drivers of NP and DEA across time

The data in this repository derives from bulk soil sampling from replicated trials in Rwanda (Nyagatare, Rubona).
Field trials were sampled for a total of 6 timepoints:
  T1 (forage anthesis, dry season) : September 2020
  T2 (forage regrowth, dry season) : October 2020
  T3 (forage anthesis, early rainy season) : November 2020
  T4 (forage regrowth, early rainy season) : December 2020
  T5 (forage anthesis, rainy season) : January 2021
  T6 (forage regrowth, rainy season) : February/March 2021

T1 samples have ID numbers in the 200s.
T2 samples have ID numbers in the 300s.
T3 samples have ID numbers in the 400s.
T4 samples have ID numbers in the 500s.
T5 samples have ID numbers in the 600s.
T6 samples have ID numbers in the 700s.

Treatments (7), 4 replicates:
Maize monoculture, maize + Desmodium distortum intercrop
Napier grass monoculture, Napier + D. distortum intercrop
Brachiaria cv. Mulato II monoculture, Brachiaria + D. distortum intercrop
D. distortum monocrop

Laboratory procedures: the same set of procedures was collected on all soil samples.
pH, GWC
POXC
Mineral N and PMN
DEA
NP
Shimadzu (NPOC, TON)
HT-qPCR (NiCE Chip)

## Datasets

_16S_data.csv_ : 16S rRNA gene abundance data from conventional qPCR

_final_nice_std_data.csv_: standardized high-throughput qPCR N cycle gene abundances used in analysis

_assay_list_2.csv_: metadata for high-throughput qPCR

_dea_final.csv_: denitrification enzyme activity dataset

_np_final.csv_: nitrification potential dataset

_poxc_final_snc.csv_: permanganate-oxidizable carbon dataset

_SNC_metadata_2020.csv_: metadata for soil samples

_SNC_NPOC_TN_infal.csv_: non-purgable organic carbon and dissolved organic nitrogen dataset

_NH4_final.csv_: ammonium dataset

_NO3_final.csv_: nitrate dataset

## Scripts
