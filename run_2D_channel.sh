#!/bin/bash

# Define the filename variable
FILENAME="2D_channel.jl"
SUCSESS="Skript ferdig! :) \n\nHilsen\n$(hostname)"
FAIL="Oi, nå har det skjedd noe galt. Skript feila :( \nDet går bra, dette fikser du! \n\nHilsen\n$(hostname)"

nice julia --project=. channel/$FILENAME && echo -e "$SUCSESS" | mail -s "$FILENAME" alsjur@uio.no alpsjur@runbox.com || echo -e "$FAIL" | mail -s "$FILENAME" alsjur@uio.no alpsjur@runbox.com
