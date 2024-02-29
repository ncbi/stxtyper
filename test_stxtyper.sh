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

FAILURES=0
TESTS=4

function test_input_file {
    local test_base=$1
    if ! $STXTYPER -n "test/$test_base.fa" > "test/$test_base.got"
    then
        echo "not ok: $STXTYPER returned a non-zero exit value indicating a failure of the software"
        echo "#  $STXTYPER -n test/$test_base.fa > test/$test_base.got"
        return 1
    else
        if ! diff -q "test/$test_base.expected" "test/$test_base.got"
        then
            echo "not ok: $STXTYPER returned output different from expected"
            echo "#  $STXTYPER -n test/$test_base.fa > test/$test_base.got"
            echo "# diff test/$test_base.expected test/$test_base.got"
            diff "test/$test_base.expected" "test/$test_base.got"
            echo "#  To approve run:"
            echo "#     mv test/$test_base.got test/$test_base.expected "
            return 1
        else
            echo "ok: test/$test_base.fa"
            return 0
        fi
    fi
}

test_input_file 'basic'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'synthetics'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'virulence_ecoli'
FAILURES=$(( $? + $FAILURES ))

test_input_file 'cases'
FAILURES=$(( $? + $FAILURES ))

echo "Done."
echo ""
if [ "$FAILURES" -gt 0 ]
then
    PASSED=$(( $TESTS - $FAILURES ))
    echo "not ok overall: $FAILURES out of $TESTS stxtyper tests failed"
    exit 1
else
    echo "ok: all $TESTS stxtyper tests passed"
fi
