
# script to set up the parameters used in simulations
# Things to track or calculate:

# psi - prob a site is occupied in year 1
# r - prob that if a site is occupied in year 1, that it is occupied with many bats

# gamma1 - prob an unoccupied site in the previous year gets colonized
# gamma2 - prob that if an unoccupied site in the previous year gets colonized, it is with many bats

# phi1 - prob a site that was occupied with a few bats the previous year remains occupied with a few bats
# G - prob a site with few bats the previous year sees growth in the population to many bats

# phi2 - prob a site with many bats last year persists with many bats
# D - prob a site with many bats last year declines to have few bats


calc_psi <- function(coeffs = NULL, celldata = NULL){
  
  lpsi <- 0
  psi <- exp(lpsi)/(1+exp(lsi))
  
}
