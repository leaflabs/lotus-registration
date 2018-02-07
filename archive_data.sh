#!/bin/bash

declare -a arr=(
"bead/20170915/registration"
"bead/20171011/registration"
"bead/20180109/registration"
"fish/20171202/interpolate"
"fish/20171221/interpolate"
"worm/20170720/interpolate"
"worm/20171012/interpolate"
"worm/20171114/registration"
"worm/20171117/registration"
"worm/20171215/interpolate"
"worm/20180103/registration"
)

LPATH=/om/user/jkinney/DLFM

for i in "${arr[@]}"; do
    echo "##"
    echo "## $i"
    echo "##"

    DEST="$LPATH/$i"
    if [ -d "$DEST" ]; then
        echo ""
        echo $DEST
        echo 'Folder does exist.'
	today=`date +%Y%m%d_%H%M%S` # or whatever pattern you desire
	NEW='_'
	cmd="mv $DEST $DEST$NEW$today"
	echo $cmd
	mv $DEST $DEST$NEW$today
    fi
done
