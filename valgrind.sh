#!/usr/bin/env bash

set -eu

INNAME=test.hep.gz
OUTNAME=valgrind_delphes_output.root
CARD=cards/delphes_card_ATLAS.tcl

function check-input() {
    if [[ ! -f $INNAME ]] ; then
	echo "no input file ${INNAME}... quitting" >&2
	return 1
    fi
}

check-input $INNAME
check-input $CARD

gunzip $INNAME -c > test.hep
valgrind ./DelphesSTDHEP $CARD $OUTNAME test.hep