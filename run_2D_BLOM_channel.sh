#!/bin/bash

nice julia --project=. channel/2D_BLOM_channel.jl && echo "Skript ferdig på ml3! :)" | mail -s "2D kanal" alsjur@uio.no  || echo "Skript feila på ml3 :(" | mail -s "2D kanal" alsjur@uio.no
