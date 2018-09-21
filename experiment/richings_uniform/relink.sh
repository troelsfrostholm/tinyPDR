#!/bin/bash

rm reactions_verbatim.dat info.log output.dat plots
ln -s $1/reactions_verbatim.dat .
ln -s $1/info.log .
ln -s $1/output.dat .
ln -s plots-$1 plots
