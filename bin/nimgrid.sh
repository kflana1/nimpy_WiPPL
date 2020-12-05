#!/bin/bash
#-----------------------------------------------------------------
#Script: nimgrid 
#Usage: nimgrid [--help]
#-----------------------------------------------------------------
nimgrid ()
{

    read -r -d '' help <<- EOM
--------------------------------------------------------------------
Usage: nimgrid [--help]
--------------------------------------------------------------------
    A convient way to make the WiPAl grid nimrod.in files and stitch
    them. Takes a nimrod.in file and appends everything below grid
    to nimrod.inb, nimrod.int, nimrod.inlo, nimrod.inlc, nimrod.inuo
    and nimrod.inuc then stitches them together for you making a
    ready to go WiPAL nimrod simulation.

    --help: displays this helpful description

    Hooray WiPAL!
-------------------------------------------------------------------
EOM

    if [ "$1" == "--help" ]; then
        echo "$help"; return
    fi

    if [ ! -f nimrod.in ]; then 
        echo "No nimrod.in exists, try again."
        return
    fi

    createmenu ()
    {
		echo ""
        echo "List of grids in grid source directory"
        echo "Please select a grid to use..."
        PS3="Pick an option: "
        export COLUMNS=1
        select option in "$@"; do
            if [ 1 -le "$REPLY" ] && [ "$REPLY" -le "$#" ]; then
                break;
            else
                echo "Incorrect Input: Select a number 1-$#"
            fi
        done
    }
    
    grids=($(ls $NIMGRIDDIR))
    for i in ${!grids[@]}; do
        if [ "${grids[$i]}" == "README" ]; then
            unset grids[$i]
        fi
    done

    case ${#grids[@]} in
        0 ) echo "No grids available, exiting"
            return;;
        1 ) grid=${grids[0]};;
        * ) createmenu "${grids[@]}" "Quit"
            grid=$option
            if [ ! "$grid" == "Quit" ];then
                grid=${grids[$(($REPLY))]};
            else
                echo "Ok, exiting"
                return;
            fi;;
    esac

    griddir=$NIMGRIDDIR/$grid

    if [ -f $griddir/stitin ]; then 
        cp $griddir/stitin .
    fi
    
    cp $griddir/nimrod.in* .
    regions=($(ls $griddir/nimrod.in* | xargs -n 1 basename))
    for i in "${regions[@]}"
    do
        sed '/&grid_input/,/\//d' nimrod.in >> $i
    done
    
    cores=$(cat $griddir/cores) 
    echo ""
    echo "Excellent choice!"
    echo "$grid requires $cores cores."

    echo ""
    echo "Would you like to perform stitch?"
    read -rp "Enter (yY/nN): " choice
    case $choice in 
          y|Y ) stitch < stitin;;
          n|N ) echo "ok, exiting";;
            * ) echo "invalid choice, exiting";;
    esac

    return
}
