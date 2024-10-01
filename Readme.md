# SEITR Network Analysis

This project implements a SEITR (Susceptible, Exposed, Infected, Treated, Recovered) model for network analysis using R. The model simulates the spread of a disease through a network of individuals, allowing for various network types and parameters.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
  - [SEITR](#seitrs)
  - [SEITR_network](#seitr_network)
  - [compare_experiment_sets](#compare_experiment_sets)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

To use this project, you need to have R installed on your system. You can install the required packages using the following commands:

```r
install.packages(c("igraph", "deSolve", "ggplot2", "dplyr"))
```

### Installing the Package from GitHub
You can install the SEITR Network Analysis package directly from GitHub using the devtools package. First, make sure you have devtools installed:
```r
install.packages("devtools")
```
Then, use the ```install_github``` function to install the package:
```r
devtools::install_github("skaraoglu/SEITRNet")
```
### Loading the Package
After installing the package, you can load it using the library function:
```r
library(SEITRNet)
```

## Usage

### SEITR Model
The SEITR function defines the SEITR model for network analysis.
### SEITR Network Analysis
The SEITR_network function performs SEITR network analysis.
### Compare Experiment Sets
The compare_experiment_sets function compares the results of multiple experiment sets.

## Examples
Here are some examples of how to use the functions in this project:

### Example 1: SEITR Network Analysis
```r
# Perform SEITR network analysis
ws_p.1_k20 <- SEITR_network("WS", n_par1=0.1, n_par2=20, num_exp = 3)
ws_p.3_k20 <- SEITR_network("WS", n_par1=0.3, n_par2=20, num_exp = 3)
ws_p.5_k20 <- SEITR_network("WS", n_par1=0.5, n_par2=20, num_exp = 3)
ws_p.9_k20 <- SEITR_network("WS", n_par1=0.9, n_par2=20, num_exp = 3)
```

### Example 2: Compare Experiment Sets
```r
# Compare the results of multiple experiment sets
compare_experiment_sets(list(ws_p.1_k20, ws_p.3_k20, ws_p.5_k20, ws_p.9_k20))
```

### Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an issue if you have any suggestions or improvements.

### License
This project is licensed under the GPL-3 License. See the LICENSE file for details.

