# Spatial-effects-in-parasite-produced-marine-diseases

This repository contains the code used to produce the results published in 

## Abstract
Marine infectious diseases are more prevalent in recent times due to climate change and other anthropogenic pressures, posing a substantial threat to marine ecosystems and the conservation of their biodiversity. An important subset of marine organisms are sessile, for which the most common mechanism for disease transmission is direct contact with waterborne parasites. Only recently, some deterministic compartmental models have been proposed to describe this kind of epidemics, being these models based on non-spatial descriptions where space is homogenised and parasite mobility is not explicitly accounted for. However, in realistic situations, epidemic transmission is conditioned by the spatial distribution of hosts and the parasites mobility patterns. Thus, the interplay between these factors  is expected to have a crucial effect in the evolution of the epidemic, so calling for a explicit description of space. In this work we develop a spatially-explicit individual-based model to study disease transmission by waterborne parasites in sessile marine populations. We investigate the impact of spatial disease transmission, performing extensive numerical simulations and analytical approximations. Specifically, the effects of parasite mobility into the epidemic threshold and the temporal evolution of the epidemic are assessed. We show that larger values of pathogen mobility have two main implications: more severe epidemics, as the number of infections increases, and shorter time-scales to extinction. Moreover, an analytical expression for the basic reproduction number of the spatial model is derived as function of the non-spatial counterpart, which characterises a transition between a disease-free and a propagation phase, in which the disease propagates over a large fraction of the system. This allows to determine a phase diagram for the epidemic model as function of the parasite mobility and the basic reproduction number of the non-spatial model.

## Requirements

- Julia v 1.6 or higher with the following libraries
  
  - StatsBase
  - Plots (for animations, but not necessary to run the simulations)
  
## Tutorial
