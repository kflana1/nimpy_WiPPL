#!/bin/bash
# filename: nimupdate.sh
# description: contains nimupdate function that should be sourced
# via ~/.bashrc checks if pwd is descendant of nimrod source
# directory and sets the binary executables in ~/bin appropriately

nimupdate ()
{
    read -r -d '' help <<- EOM
--------------------------------------------------------------------
Usage: nimupdate [--help]
--------------------------------------------------------------------
    This script will help you keep track of which nimrod install you
    are working with. To use, cd to a working nimrod install
    directory and run this. It will automatically update the softlinks
    in ~/bin for all the nimrod executables (nimrod, nimset, nimfl,
    stitch, nimplot, nimeq, sol, xlog, and fluxgrid) to this nimrod
    install.
    There is an 'are you sure' check, so don't sweat it.

    --help: displays this helpful description

    Now go update like nobody has updated before!
-------------------------------------------------------------------
EOM

    if [ "$1" == "--help" ]; then
        echo "$help"; return
    fi

    updatelinks ()
    {
        echo "NIMSRCDIR is now set to ${NIMSRCDIR}"
        ln -sfn ${NIMSRCDIR}/nimset/stitch ${HOME}/bin/stitch
        ln -sfn ${NIMSRCDIR}/rundir/nimplot ${HOME}/bin/nimplot 
        ln -sfn ${NIMSRCDIR}/rundir/nimrod ${HOME}/bin/nimrod
        ln -sfn ${NIMSRCDIR}/rundir/nimeq ${HOME}/bin/nimeq 
        ln -sfn ${NIMSRCDIR}/rundir/nimset ${HOME}/bin/nimset 
        ln -sfn ${NIMSRCDIR}/rundir/fluxgrid ${HOME}/bin/fluxgrid
        ln -sfn ${NIMSRCDIR}/rundir/nimfl ${HOME}/bin/nimfl 
        ln -sfn ${NIMSRCDIR}/rundir/sol ${HOME}/bin/sol
        ln -sfn ${NIMSRCDIR}/rundir/xlog ${HOME}/bin/xlog 
		echo "Printing contents of ${HOME}/bin..."
		echo "-----------------------------------"
        ls -lah ${HOME}/bin
    }

    createmenu ()
    {
        PS3="Pick a source: "
        select option in "$@"; do
            if [ 1 -le "$REPLY" ] && [ "$REPLY" -le "$#" ]; then
                break;
            else
                echo "Incorrect Input: Select a number 1-$#"
            fi
        done
    }                                                                                 
    echo ""
	echo "Current NIMSRCDIR: ${NIMSRCDIR}"
    export NIMRESULTPATHS=$(echo $(find ${NIMTOP} -type d -print | grep results) | sed 's/ /:/g')
    export NIMSRCPATHS=$(echo $(find ${NIMTOP} -type d -print | grep .*/rundir | sed 's/\/rundir.*//g') | sed 's/ /:/g')
    IFS=':' read -ra srcdirs <<< "$NIMSRCPATHS"
    result=""
    for src in "${srcdirs[@]}"; do
        if [[ -n "$(find "$src" -type d -wholename "$PWD")" ]];then
            result="$src"
        fi
    done

    if [[ -n "$result" ]]; then
        read -rp "Are you sure you want to set NIMSRCDIR to $result? (yY/nN): "
        case $REPLY in
            n|N ) echo "Quitting nimupdate, no change to NIMSRCDIR: $NIMSRCDIR";;
            y|Y ) export NIMSRCDIR="$result"
                  sed -i "s@NIMSRCDIR=.*@NIMSRCDIR=$NIMSRCDIR@" ${HOME}/.bashrc
                  updatelinks;;
            * ) echo "Quitting nimupdate, no change to NIMSRCDIR: $NIMSRCDIR";;
        esac
    else
        echo "$PWD is not a child of a nimrod source directory,"
        echo "Please select one of the following valid sources..."
        createmenu "${srcdirs[@]}" "Quit"
        if [ "$option" == "Quit" ]; then
            echo "Quitting nimupdate, no change to NIMSRCDIR: $NIMSRCDIR"
        else
            export NIMSRCDIR="$option"
            sed -i "s@NIMSRCDIR=.*@NIMSRCDIR=$NIMSRCDIR@" ${HOME}/.bashrc
            updatelinks
        fi
    fi

    return
}
