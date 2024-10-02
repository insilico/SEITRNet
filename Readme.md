# SEITR Network Analysis

This project implements a SEITR (Susceptible, Exposed, Infected, Treated, Recovered) model for network analysis using R. The model simulates the spread of a disease through a network of individuals, allowing for various network types and parameters.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [SEITR_network](#seitr_network)
    - [Parameters](#parameters)
  - [compare_experiment_sets](#compare_experiment_sets)
- [Examples](#examples)
- [Well-Posed SEITR Mathematical Model](#seitrmodel)
- [Algorithm](#algorithm)
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

### SEITR_network
The SEITR_network function performs SEITR network analysis.

#### Parameters

- n: (int) Number of nodes in the network.
- network_type: (str) The type of network to create. This can be one of the following:

  - ER: Erdős-Renyi random graph
  - BA: Barabasi-Albert scale-free network
  - WS: Watts-Strogatz small-world network
  - LN: Lattice network
  - RR: Random regular network

- n_par1: (float) First network parameter. Assigned for "p" argument in Erdős-Renyi and Watts-Strogatz graphs, assigned for "k" argument for Lattice and Random regular, n * n_par1 is used as "m" for Barabási-Albert networks.

- n_par2: (float) Second network parameter that is only needed for Watts-Strogatz networks to assign "k" argument.

- Lambda: (float) Birth rate.

- alpha1: (float) Treatment rate
- alpha2: (float) Recovery rate from treatment
- delta_I: (float) Death rate from infection
- delta_T: (float) Death rate from treatment
- mu: (float) Natural death rate
- beta1: (float) The infection rate parameter.
- beta2: (float) The exposure rate parameter.
- beta3: (float) The recovery rate parameter.
- initial_statuses: (int) The initial status of the nodes in the network. Where each status can be one of the following:

    S: Susceptible, E: Exposed, I: Infected, T: Treatment, R: Recovered

- t: (int) Time period.
- num_exp: (int) The number of experiments to run.
- verbose: (bool) Verbose output.

```r
SEITR_network <- function(network_type="ER", n=100, n_par1=.9, n_par2=10, Lambda=1.1, beta1=.8, beta2=.18, beta3=.02, alpha1=.1, alpha2=.055, delta_I=.03, delta_T=.03, mu=.01, S=85, E=5, I=10, Tt=0, R=0, N=100, t=100, num_exp = 10, verbose = F, state = NULL, parameters = NULL) {
  # Function implementation
}
```
### Compare Experiment Sets
The compare_experiment_sets function compares the results of multiple experiment sets. Currently inactive.

## Examples
Here are some examples of how to use the functions in this project:

### Example 1: SEITR Network Analysis
```r
# Perform SEITR network analysis
ws_p.1_k20 <- SEITR_network("WS", n_par1=0.1, n_par2=20, num_exp = 3)
er_p.2 <- SEITR_network("ER", n_par1=0.2, num_exp = 3)
ba_m.75 <- SEITR_network("BA", n_par1=0.75, num_exp = 3)
```

### Example 2: Compare Experiment Sets
```r
# Compare the results of multiple experiment sets
compare_experiment_sets(list(ws_p.1_k20, er_p.2, ba_m.75))
```

## Well-Posed SEITR Mathematical Model

Here we present the mathematical model used to analyze the disease transmission dynamics. This model is an extension of the classic Susceptible-Infectious-Recovered model.

### Model Formulation
The model is an extension of the classic Susceptible-Infectious-Recovered (SIR) model with Exposed and Recovered statuses. The model consists of a system of non-linear ordinary differential equations (ODEs) expressing the transmission dynamics:
- Susceptible (S): These are individuals who are susceptible to the disease. The rate of change of S is given by the equation: $$\frac{dS}{dt}= \Lambda-\frac{\beta_1 SI}{N}-\eta S$$
Here, $\Lambda$ is the birth rate, $\beta_1$ is the transmission rate of the disease, I is the number of infectious individuals, N is the total population, and $\eta$ is the natural death rate.
- Exposed (E): These are individuals who have been infected but are not yet infectious. The rate of change of E is given by the equation: $$\frac{dE}{dt}= \frac{\beta_1SI}{N}-(\gamma_2+\eta) E$$
Here, $\gamma_2$ is the rate at which exposed individuals become infectious.
- Infectious (I): These are individuals who are infectious. The rate of change of I is given by the equation: $$\frac{dI}{dt}=\gamma_2 E-(\gamma_3+\eta+d_I+\kappa_1)I$$
Here, $\gamma_3$ is the recovery rate, $dI$ is the disease-induced death rate, and $\kappa_1$ is the rate at which infectious individuals are treated.
- Treated (T): These are individuals who have been infected and are receiving treatment. The rate of change of T is given by the equation: $$\frac{dT}{dt}=\kappa_1 I-(\eta+d_T+\kappa_2)T$$
Here, $dT$ is the treatment-induced death rate, and $\kappa_2$ is the rate at which treated individuals recover.
- Recovered (R): These are individuals who have recovered from the disease. The rate of change of R is given by the equation: $$\frac{dR}{dt}= \gamma_3 I+\kappa_2 T-\eta R$$

with non-negative initial conditions;
$$S(0)=S_0,E(0)=E_0,I(0)=I_0,T(0)=T_0,R(0)=\mathcal{R}_0.$$

### Equilibrium Points

The model also includes two equilibrium points: the disease-free equilibrium and the endemic equilibrium. The disease-free equilibrium represents the state where the disease has been eradicated, while the endemic equilibrium represents the state where the disease persists in the population. The reproduction number $\mathcal{R}_0$ is a key quantity that determines the number of infections produced by an infectious individual. If $\mathcal{R}_0<1$, the disease-free equilibrium is stable, meaning the disease will die out. If $\mathcal{R}_0>1$, the endemic equilibrium is stable, meaning the disease will persist in the population.

Disease-Free Equilibrium: $\zeta =\left(\frac{\Lambda}{\nu} , 0,0,0,0\right)$.

Endemic Equilibrium: $\Gamma^{1}=(S^{1},E^{1},I^{1},T^{1},R^{1})$, where $S^{1}$, $E^{1}$, $I^{1}$, $T^{1}$, and $R^{1}$ are given by the following equations:
$$S^{1}=\frac{\Lambda(-k_1 k_2 k_3+k_3\delta I \beta_2+\delta_T \alpha_1 \beta_2)} {(-k_3 \beta_1+k_3\delta_I+\delta_T \alpha_1)\beta_2 \mu}$$
$$E^{1}=\frac{-\Lambda\beta_1 \beta_2 k_1 k_3+ k_1^2 k_2 k_3 \Lambda}{-k_1 k_2 k_3 \beta_1 \beta_2 + k_1 k_2 k_3\delta_I \beta_2+k_1 k_2 \delta_T \alpha_1 \beta_2}$$
$$I^{1}=\frac{-\Lambda \beta_1 \beta_2 k_3+ k_1 k_2 k_3 \Lambda}{-k_1 k_2 k_3 \beta_1 + k_1 k_2 k_3\delta_I +k_1 k_2 \delta_T \alpha_1}$$
$$T^{1}=\frac{\alpha(-\Lambda \beta_1 \beta_2 k_3+ k_1 k_2 k_3\Lambda)}{k_3(-k_1 k_2 k_3 \beta_1 +k_1 k_2 k_3\delta_I +k_1 k_2 \delta_T \alpha_1)}$$
$$R^{1}=\frac{-\Lambda \beta_1 \beta_2 \beta_3 k_3+ k_1 k_2 k_3\beta_3 \Lambda-\Lambda \beta_1 \beta_2 \alpha_1 \alpha_2+k_1 k_2\alpha_1 \alpha_2 \Lambda}{\mu(-k_1 k_2 k_3 \beta_1 + k_1 k_2 k_3 \delta_I+ k_1 k_2 \alpha_1 \delta_T)}$$
$$N^{1}=\frac{\Lambda(-k_1 k_2 k_3+k_3\delta_I \beta_2+\delta_T\alpha_1 \beta_2)}{\mu \beta_2(-k_3 \beta_1 +k_3 \delta_I+\delta_T\alpha_1)}$$
### Reproduction Number and Stability
The reproduction number $\mathcal{R}_0$ is a mathematical quantity that determines the number of infections produced by an infectious individual. For our model, it is given by:

$$\mathcal{R}_0=\frac{\beta_1 \beta_2 }{(\beta_3+\mu+\delta_I+\alpha_1)(\beta_2+\mu)}$$

The model is locally asymptotically stable (LAS) and globally asymptotically stable (GAS) at the disease-free equilibrium point $\zeta$ when $\mathcal{R}_0 < 1$ and unstable when $\mathcal{R}_0 > 1$, provided $\mathcal{R}_0 \ne 1$. The endemic equilibrium point $\Gamma^{1}$ of the model is LAS and GAS stable, provided $\mathcal{R}_0>1$.

## Algorithm: SEITR Network Analysis Algorithm

**Input**: Network Type `ER`, n=100, n_par1=.9, n_par2=10, Lambda=1.1, beta1=.8, beta2=.18, beta3=.02, alpha1=.1, alpha2=.055, delta_I=.03, delta_T=.03, mu=.01,
                          S=85, E=5, I=10, Tt=0, R=0, N=100, t=100, num_exp = 10, verbose = F

**Output**: avgs 'Average results of the num_exp experiments"

1. **For each experiment in num_exp experiments:**
    1. Create network 'g' with given parameters
    2. Assign statuses to nodes randomly with given parameters
    3. **For each time step t in the simulation:**
        1. Save the current status of each node.
        2. **For each type of removal (Infection, Treatment, natural death):**
            1. Calculate the number of nodes to remove `nodes_to_remove_count`
            2. Separate `nodes_to_remove_count` into a floor value `floor_value` and a fractional part `fractional_part`
            3. **if** `floor_value` > 0 **then**
                1. Remove `floor_value` number of nodes
            4. **if** `fractional_part` > 0 **then**
                1. Remove an additional node with a probability equal to `fractional_part`
        3. Calculate the number of nodes to add `nodes_to_add_count`
        4. Separate `nodes_to_add_count` into a floor value `floor_value` and a fractional part `fractional_part`
        5. **if** `floor_value` > 0 **then**
            1. Add `floor_value` number of nodes 
            2. Add connections to the new nodes according to Network Type
        6. **if** `fractional_part` > 0 **then**
            1. Add an additional node with a probability equal to `fractional_part`
            2. Add connections to the new nodes according to Network Type
        7. **for** each node in the graph **do**
            1. Get the current status of the node.
            2. Generate a random number `rand`
            3. If the node is Susceptible and `rand` < `beta1` * `I` / `N`, change the status to Exposed
            4. If the node is Exposed and `rand` < `beta2`, change the status to Infected
            5. If the node is Infected:
                1. Generate another random number `rand2`
                2. If `rand` < `beta3`, change the status to Recovered
                3. If `rand2` < `alpha1`, change the status to Treatment
            6. If the node is under Treatment and `rand` < `alpha2`, change the status to Recovered
        8. Count the number of nodes in each status and store these counts in `S_count`, `E_count`, `I_count`, `Tt_count`, `R_count`, `N_count`
        9. Calculate various network metrics and store these metrics in `degree_dist`, `clustering_coeff`, `avg_path_length`, `largest_comp_size`
        10. Store the experiment result
    4. Calculate and store the averages of num_exp experiments
2. Return averages


## Contributing
Contributions are welcome! Please feel free to submit a Pull Request or open an issue if you have any suggestions or improvements.

## License
This project is licensed under the GPL-3 License. See the LICENSE file for details.

