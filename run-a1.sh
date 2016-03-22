#!/bin/sh

# param[0]: resolution
# param[1]: unit length of the routing data
# param[2]: defect density
# param[3]: clustering parameter
# param[4]: wire failure rate
params="30 8.8e-7 0.004e+6 2 0.000010"
#params="30 0.2058333e-6 0.004e+6 2 1e-4"

pushd a1
echo $params | ../part adaptec1.capo70.3d.35.50.90.gr a1.out
popd
