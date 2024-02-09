#!/bin/bash

# Define the filename variable
FILENAME="2D_BLOM_channel.jl"
SUCSESS="Skript ferdig! :) \n\nHilsen\n$(hostname)"
FAIL="Oi, her har det skjedd noe galt. Skript feila :( \nDet går bra, dette fikser du! \n\nHilsen\n$(hostname)"

nice julia --project=. channel/$FILENAME && echo -e "$SUCSESS" | mail -s "$FILENAME" alsjur@uio.no  || echo -e "$FAIL" | mail -s "$FILENAME" alsjur@uio.no
