#!/usr/bin/env bash

set -eu

OUTDIR=output
INNAME=${1-test.hep.gz}
OUTNAME=${2-delphes_output.root}
CARD=cards/delphes_tracksmear.tcl

function check-input() {
    if [[ ! -f $INNAME ]] ; then
	echo "no input file ${INNAME}... quitting" >&2
	return 1
    fi
}

check-input $INNAME
check-input $CARD

mkdir -p $OUTDIR
OUTPATH=${OUTDIR}/${OUTNAME}
if [[ -f $OUTPATH ]] ; then
    mv $OUTPATH ${OUTDIR}/$(date +%F-%R)-$OUTNAME
fi

gunzip $INNAME -c | ./DelphesSTDHEP $CARD $OUTPATH -
