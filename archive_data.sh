#!/bin/bash

declare -a arr=(
"bead/20180328/registration"
)

LPATH=/om/user/jkinney/DLFM

for i in "${arr[@]}"; do
    echo "##"
    echo "## $i"
    echo "##"

    ORIGIN="$LPATH/$i"

    if [ -d "$ORIGIN" ]; then

        echo ""
        echo $ORIGIN
        echo 'Folder exists.'
	TODAY=`date +%Y%m%d_%H%M%S` # or whatever pattern you desire
	DEST="${ORIGIN}_$TODAY"
	CMD="mv $ORIGIN $DEST"
	echo $CMD
	eval $CMD

	# find any *null.mat from $DEST
	CMD='find '$DEST' -name '\"'*null.mat'\"' -exec basename {} .mat \;'
	echo $CMD
	NULL=$(eval $CMD)
	echo NULL=\'$NULL\'
	# if any found
	if [[ ! -z $NULL ]]; then
		CMD="mkdir $ORIGIN"
		echo $CMD
		eval $CMD
		# for each file
		for F in $NULL; do
			CMD='cp '$DEST'/'$F'.mat '$ORIGIN
			echo $CMD
			eval $CMD
			echo ""
		done
	fi
    fi
done
