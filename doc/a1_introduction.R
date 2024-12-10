## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options(crayon.enabled = F)

## ----setup, warning=FALSE-----------------------------------------------------
library(tickTack)
require(dplyr)

## ----echo=TRUE----------------------------------------------------------------
tickTack::pcawg_example$mutations %>% head()

## ----echo=TRUE----------------------------------------------------------------
tickTack::pcawg_example$cna %>% head()

## ----echo=TRUE----------------------------------------------------------------
tickTack::pcawg_example$metadata

## ----fig.width=9.5, fig.height=3----------------------------------------------
tickTack::chr_coordinates_hg19

tickTack::chr_coordinates_GRCh38

