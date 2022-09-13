# psydiff

**Numerical Integration of Panel Stochastic Differential Equations**

**WARNING**: This package is under development and should **not** be used in any real analyses. **Please use dynr instead** (https://github.com/mhunter1/dynr). psydiff can be used to get some ideas of how to implement stochastic differential equations in R using the square root unscented Kalman filter procedure developed by Sarkka (2007).

psydiff is an R package which implements numerical integration of panel stochastic differential equations. The square root unscented Kalman filter procedure developed by Sarkka (2007) is implemented.

The package uses Rcpp and the interface to the Boost odeint package for numerical integration. It also offers a simple approach to estimate parameters person-specific or group specific.

Sarkka, S. (2007). On Unscented Kalman Filtering for State Estimation of Continuous-Time Nonlinear Systems. IEEE Transactions on Automatic Control, 52(9), 1631–1641. https://doi.org/10.1109/TAC.2007.904453


# Installation

If you want to install psydiff from GitHub, use the following commands in R:

    if(!require(devtools))install.packages("devtools")

    devtools::install_github("jhorzek/psydiff")

# Getting Started

A good place to start is the help page of the main function for setting up models:

    ?psydiff::newPsydiff
