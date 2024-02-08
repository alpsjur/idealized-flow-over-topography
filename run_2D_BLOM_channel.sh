#!/bin/bash

nice julia --project=. channel/2D_BLOM_channel.jl && echo "Skript ferdig på ml3! :)" | sendmail alsjur@uio.no "2D kanal" || echo "Skript feila på ml3 :(" | sendmail alsjur@uio.no "2D kanal"