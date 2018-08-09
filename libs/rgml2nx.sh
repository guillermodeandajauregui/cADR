#!/bin/bash
tail -n+3 $1 > $2
sed -i 'N;s/\s+\[/ \[/g;P;D' $2
sed -i 's/name/label/g' $2
