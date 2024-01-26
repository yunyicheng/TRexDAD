# Purpose: Run TRexDAD shiny app
# Author: Yunyi Cheng
# Date: 01.25.2024
# Version: 1.0.0
# Bugs and Issues: N/A

#' Run Shiny App For Package TRexDAD
#'
#' A function that launches the shiny app for this package.
#' The shiny app offers users an interactive interface that allows them to
#' put in their gene of interest for assembly design.
#' Supplementary results, including the final positions, optimal tiles, 
#' gene codons, oligo cost, overhang_fidelity chart, 
#' are also provided to the user.
#'
#' @return It opens the window for the shiny app.
#'
#' @examples
#' \dontrun{
#' run_TRexDAD()
#' }
#'
#' @references
#' #' @references 
#' Vladimir Potapov, Jennifer L. Ong, Rebecca B. Kucera, Bradley W. Langhorst, 
#'      Katharina Bilotti, John M. Pryor, Eric J. Cantor, Barry Canton, 
#'      Thomas F. Knight, Thomas C. Jr. Evans, and Gregory J. S. Lohman. 
#'      Comprehensive Profiling of Four Base Overhang Ligation Fidelity by 
#'      T4 DNA Ligase and Application to DNA Assembly. 
#'      ACS Synthetic Biology, 7(11):2665–2674, 2018.
#' Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, et al. 
#'      Enabling one-pot Golden Gate assemblies of unprecedented 
#'      complexity using data-optimized assembly design. 
#'      PLOS ONE 15(9): e0238592, 2020.
#' Alberts B., Johnson A., Lewis J., et al. 
#'      "Molecular Biology of the Cell." 6th edition. 
#'      Garland Science, 2014.
#'
#' @export
#' @importFrom shiny runApp
run_TRexDAD <- function() {
    dir <- system.file("shiny_script", package = "TRexDAD")
    shiny::runApp(dir, display.mode = "normal")
    return()
}