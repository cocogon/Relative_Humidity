# Relative_Humidity
Data and scripts for the paper 'Relative humidity predicts day-to-day variations in COVID-19 cases in the city of Buenos Aires'
https://www.medrxiv.org/content/10.1101/2021.01.29.21250789v2

This project is licensed under the terms of the CC BY-NC-ND 4.0 license.


The AllData.db file contains all data in Matlab format '-mat'. These data are public and can be obtained in other formats from the sources cited in the publication. The file contains a single structure 'data' with fields 'time' for the date in MATLAB date format and 'covid' for the number of cases reported that day. The other fields correspond to daily meteorological observations and are named as in the abbreviations of the publication.

The Scripts.m file contains MATLAB code for all major figures in the publication. As the script runs, it loads the data.
