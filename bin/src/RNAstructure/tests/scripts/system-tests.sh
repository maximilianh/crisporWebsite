verifyInputFiles -t $SINGLESEQ && echo "GOOD" || echo "BAD"

verifyInputFiles -t banana.txt && echo "BAD" || echo "GOOD"

(verifyInputFiles banana.txt; echo "BAD") && echo "BAD" || echo "GOOD"
