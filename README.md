
# conisi

<!-- badges: start -->
<!-- badges: end -->

The goal of conisi is to let users simulate the Covid-19 pandemic.

## Installation

You can install the released version of conisi from github with:

``` r
devtools::install_github(
  "wimmyteam/conisi",
  ref="master"
  )
```

## Example

How to run a single simulation:

``` r

  file_path <- system.file("test-data", "manually_edited_parameters.csv", package="conisi")
  # You would normally provide your own parameter file
  
  parameters <- conisi::load_parameters(file_path, sim_name = "sa_37", start_date = lubridate::ymd("2020-03-10"))
  
  pop <- 1000000 #The size of the population to simulate
  
  days <- 100 #The number of days to simulate
  
  mod_result <- COVIDmodel_run_and_mutate(par_table, pop, days) # Runs the simulation


```

How to run many experiments:

``` r

  parameters <- get_par_table() # Get the parameters somewhere, needs "experiment" column.
  # Need to update the conisi::load_parameters function to allow including many experiments in a single file.
  
  pop <- 1000000 #The size of the population to simulate
  
  days <- 100 #The number of days to simulate
  
  
  availablecores <- 2
  doParallel::registerDoParallel(cores = availablecores)
  
  ## Alternatively to run on azure:
  # library(doAzureParallel)
  # setCredentials("../credentials.json")
  # cluster <- makeCluster("cluster.json")
  # registerDoAzureParallel(cluster)
  # getDoParWorkers()
  
  mod_result <- COVIDmodel_run_and_mutate_many(xparm_table, pop, 100)

  ## If using azure, remember to turn off cluster
  # stopCluster(cluster)
```



