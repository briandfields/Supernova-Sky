# Supernova-Sky

Codes for calcuating the sky distributions of supernovae of different types, brighter than a given apparent magnitude for an observer at Earth (or other specified Galactic location).

## SNvisMap_int.py
Calculation via integration along an array of sightlines.  

The integration depth (radius upper limit) for each sightline is found via a root-finding approach, using Newton's method.

## mc_map.py
Monte Carlo calculation.  Based on Tanner's Monte Carlo code for Milky Way supernova distributions.  

Initial version includes changes from BDF, including:
* diagnostics of number of supernovae in different quadrants, to test for symmetries
* plots tweaked to show magnitude cutoff, and in other smaller ways
* now also makes a scatter plot showing visible and invisible supernovae
