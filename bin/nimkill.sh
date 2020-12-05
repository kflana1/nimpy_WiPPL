#!/bin/bash
#-----------------------------------------------------------------
#Script: nimkill
#Usage: nimkill [--help]
#-----------------------------------------------------------------
nimkill ()
{

    read -r -d '' help <<- EOM
--------------------------------------------------------------------
Usage: nimkill [--help]
--------------------------------------------------------------------
    A nimrod process exterminator. When run, will look for any
    nimrod processes that are running and ask you which to kill. If
    multiple processes are running they will be listed via the
    directory they are in.
    Don't worry, we'll ask if you are sure before killing anything!

    --help: displays this helpful description

    Happy Exterminating!
-------------------------------------------------------------------
EOM

    if [ "$1" == "--help" ]; then
        echo "$help"; return
    fi

    createmenu ()
    {
		echo ""
        echo "List of currently running nimrod working directories"
        echo "Please select nimrod run to kill or exit program..."
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
    
    nimdirs=($(pwdx $(pgrep nimrod) 2>/dev/null | sort -uk 2 | sed 's/.*: //g'))
    nimpids=($(pwdx $(pgrep nimrod) 2>/dev/null | sort -uk 2 | sed 's/:.*//g'))

    case ${#nimdirs[@]} in
        0 ) echo "No current nimrod processes, exiting"
            return;;
        1 ) nimdir=${nimdirs[0]}
            nimpid=${nimpids[0]};;
        * ) createmenu "${nimdirs[@]}" "Quit"
            nimdir=$option
            if [ ! "$nimdir" == "Quit" ];then 
                nimpid=${nimpids[$(($REPLY - 1))]};
            else
                echo "No processes have been killed, exiting"
                return;
            fi;;
    esac

    echo ""
    echo "Are you sure you want to kill PID#: $nimpid ?"
    echo "nimrod working directory: $nimdir"
    read -rp "Enter (yY/nN): " choice
    case $choice in 
          y|Y ) kill $nimpid && echo "killed PID#: $nimpid";;
          n|N ) echo "ok, exiting";;
            * ) echo "invalid choice, exiting";;
    esac

    return
}
