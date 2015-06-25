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
rm -f $OUTNAME

gunzip $INNAME -c >| test.hep
SUPP=--suppressions=${ROOTSYS}/etc/valgrind-root.supp
valgrind $SUPP ./DelphesSTDHEP $CARD $OUTNAME test.hep 2>&1 | tee valgrind.log
