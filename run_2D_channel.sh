#!/bin/bash

# Define the filename variable
FILENAME="2D_channel.jl"
SUCSESS="Skript ferdig! :) \n\nHilsen\n$(hostname)"
FAIL="Oi, nå har det skjedd noe galt. Skript feila :( \nDet går bra, dette fikser du! \n\nHilsen\n$(hostname)"
julia = "/itf-fi-ml/home/alsjur/.julia/juliaup/julia-1.10.1+0.x64.linux.gnu/bin/julia"
# multithreading 
#export JULIA_NUM_THREADS=2

nice $julia --project=. channel/$FILENAME && echo -e "$SUCSESS" | mail -s "$FILENAME" alsjur@uio.no || echo -e "$FAIL" | mail -s "$FILENAME" alsjur@uio.no 
