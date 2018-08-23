#!/bin/bash

rm reactions_verbatim.dat info.log output.dat
ln -s $1/reactions_verbatim.dat .
ln -s $1/info.log .
ln -s $1/output.dat .