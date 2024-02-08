#!/bin/bash

nice julia --project=. channel/2D_BLOM_channel.jl && echo -e "Skript ferdig! :) \n\nHilsen\n$(hostname)" | mail -s "2D_BLOM_channel.jl" alsjur@uio.no  || echo -e "Oi, her har det skjedd noe galt. Skript feila :( \n\nHilsen\n$(hostname)" | mail -s "2D_BLOM_channel.jl" alsjur@uio.no
