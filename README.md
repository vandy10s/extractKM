# extractKM
The algorithm by Guyot et al.(2012) derives data from KM curve that are close to the original individual patient time-to-event data (IPD). The functions in this package map from digitised curves back to KM data using the algorithm in an easier and more convenient way. There are 3 functions which require different information types or conditions. To install the package using R, please run the commands below. 

```
library(remotes)
remotes::install_github("vandy10s/extractKM")
```

### 1. km.table
`km.table` function derives or extracts IPD when only number at risk table information are available. To check the exmaple, run the following command. 

```
example(km.table)
```


### 2. km.event
`km.event` function derives or extracts IPD when only total number of events information are available. To check the exmaple, run the following command. 

```
example(km.event)
```


### 3. km.tab.event
`km.tab.event` function derives or extracts IPD when both number at risk table and total number of events information are available. To check the exmaple, run the following command. 

```
example(km.tab.event)
```
