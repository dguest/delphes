#!/usr/bin/env bash

set -eu

OUTNAME=delphes_output.root
INNAME=test.hep.gz
CARD=cards/delphes_card_ATLAS.tcl

function check-input() {
    if [[ ! -f $INNAME ]] ; then
	echo "no input file ${INNAME}... quitting" >&2
	return 1
    fi
}

check-input $INNAME
check-input $CARD

rm -f $OUTNAME
gunzip $INNAME -c | ./DelphesSTDHEP $CARD $OUTNAME -
