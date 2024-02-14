#!/bin/bash

if [ "$1" == "path" ]
then
    echo "Testing stxtyper command in your \$PATH"
    which stxtyper
    STXTYPER=stxtyper
else
    echo "Testing ./stxtyper"
    echo "  To test stxtyper in your path run 'test_stxtyper.sh path'"
    STXTYPER=./stxtyper
fi

if ! $STXTYPER -n test/basic.fa > test/basic.got
then
    echo "not ok: $STXTYPER returned a non-zero exit value indicating a failure of the software"
    echo "#  $STXTYPER -n test/basic.fa > test/basic.got"
    exit 1
else
    if ! diff -q test/basic.expected test/basic.got
    then
        echo "not ok: $STXTYPER returned output different from expected"
        echo "#  $STXTYPER -n test/basic.fa > test/basic.got"
        echo "# diff test/basic.expected test/basic.got"
        diff test/basic.expected test/basic.got
        echo "#  To approve run:"
        echo "#     mv test/basic.got test/basic.expected "
        exit 1
    else
        echo "ok: test/basic.fa"
    fi
fi

if ! $STXTYPER -n test/synthetics.fa > test/synthetics.got
then
    echo "not ok: $STXTYPER returned a non-zero exit value indicating a failure of the software"
    echo "#  $STXTYPER -n test/synthetics.fa > test/synthetics.got"
    exit 1
else
    if ! diff -q test/synthetics.expected test/synthetics.got
    then
        echo "not ok: $STXTYPER returned output different from expected"
        echo "#  $STXTYPER -n test/synthetics.fa > test/synthetics.got"
        echo "# diff test/synthetics.expected test/synthetics.got"
        diff test/synthetics.expected test/synthetics.got
        echo "#  To approve run:"
        echo "#     mv test/synthetics.got test/synthetics.expected "
        exit 1
    else
        echo "ok: test/synthetics.fa"
    fi
fi


if ! $STXTYPER -n test/virulence_ecoli.fa > test/virulence_ecoli.got
then
    echo "not ok: $STXTYPER returned a non-zero exit value indicating a failure of the software"
    echo "#  $STXTYPER -n test/virulence_ecoli.fa > test/virulence_ecoli.got"
    exit 1
else
    if ! diff -q test/virulence_ecoli.expected test/virulence_ecoli.got
    then
        echo "not ok: $STXTYPER returned output different from expected"
        echo "#  $STXTYPER -n test/virulence_ecoli.fa > test/virulence_ecoli.got"
        echo "# diff test/virulence_ecoli.expected test/virulence_ecoli.got"
        diff test/virulence_ecoli.expected test/virulence_ecoli.got
        echo "#  To approve run:"
        echo "#     mv test/virulence_ecoli.got test/virulence_ecoli.expected "
        exit 1
    else
        echo "ok: test/virulence_ecoli.fa"
    fi
fi
echo "ok: all tests passed"
