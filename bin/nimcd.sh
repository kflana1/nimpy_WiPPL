#!/bin/bash
#-----------------------------------------------------------------
#Script: nimcd
#Usage: nimcd [--help]
#-----------------------------------------------------------------
nimcd ()
{
    read -r -d '' help <<- EOM
--------------------------------------------------------------------
Usage: nimcd [--all|-a] [--help]
--------------------------------------------------------------------
    This script will aide you in getting into a nimrod results
    directory. If only one process is running, it will take you to
    that particular directory. If more than one process is running
    it'll let you choose. There is also a flag to show all nimrod
    related results and source directories

    --all|-a: displays all nimrod results and source directories
    --help: displays this helpful description

    cd, schm-cd am I right?
-------------------------------------------------------------------
EOM

    createmenu ()
    {
        echo "Please select nimrod directory or exit program"
        echo "(... corresponds to NIMTOP: $NIMTOP)"
        echo "----------------------------------------------"
        local ii=0
		[ -n "$1" ] && declare -a act=("${!1}")
        [ -n "$2" ] && declare -a inact=("${!2}")
        [ -n "$3" ] && declare -a src=("${!3}")
        allopts=( "${act[@]}" "${inact[@]}" "${src[@]}" "$4" )

        if [ -n "$act" ]; then
            echo "List of active nimrod directories:"
            for adir in "${act[@]}"; do
                let "ii++"
                echo "$ii ) ...${adir#$NIMTOP}"
            done
		else
			echo "List of active nimrod directories:"
			echo
        fi

        if [ -n "$inact" ]; then
            echo "List of inactive nimrod results directories:"
            for idir in "${inact[@]}"; do
                let "ii++"
                echo "$ii ) ...${idir#$NIMTOP}"
            done
        fi

        if [ -n "$src" ]; then
			echo "List of nimrod source directories:"
			for sdir in "${src[@]}"; do
				let "ii++"
				echo "$ii ) ...${sdir#$NIMTOP}"
			done
		fi

		let "ii++"
		echo "$ii ) $4"
		echo

        while true; do
            read -rp "Please select a directory: "
            if [ 1 -le "$REPLY" ] && [ "$REPLY" -le "${#allopts[@]}" ]; then
                break;
            else
                echo "Incorrect Input: Select a number 1-${#allopts[@]}"
            fi
        done
    }                                                                                 

	echo
    usage_line="Usage: nimcd [--all|-a] [--help]"
	list_all=
    while [ -n "$1" ]; do
        case "$1" in
            --all | -a ) list_all="yes"; shift;;
            --help ) echo "$help"; return;;
            * ) printf "Invalid option, $usage_line \n\n"; return;;
        esac
    done

    activedirs=($(pwdx $(pgrep nimrod) 2>/dev/null | sort -uk 2 | sed 's/.*: //g'))
    IFS=':' read -ra resultdirs <<< "$NIMRESULTPATHS"
    IFS=':' read -ra srcdirs <<< "$NIMSRCPATHS"
	# Move to active directory if there is only one and --all is not set.
	([ -z "$list_all" ] && [ ${#activedirs[@]} -eq 1 ]) && cd "${activedirs[0]}" && echo "You are now in $PWD" && return
	# Otherwise, find inactive directories
    inactivedirs=()
    for i in "${resultdirs[@]}"; do
        skip=
        for j in "${activedirs[@]}"; do
            [[ $i == $j ]] && { skip=1; break; }
        done
        [[ -n $skip ]] || inactivedirs+=("$i")
    done

    alldirs=( "${activedirs[@]}" "${inactivedirs[@]}" "${srcdirs[@]}")
    if [ "${#alldirs[@]}" -ne "0" ]; then
        if [ -z "$list_all" ]; then
			if [ "${#activedirs[@]}" -ne "0" ]; then
				createmenu activedirs[@] "" srcdirs[@] "Quit"
			else
				echo "No active nimrod directories found... defaulting to all"
				createmenu activedirs[@] inactivedirs[@] srcdirs[@] "Quit"
			fi
        else
            createmenu activedirs[@] inactivedirs[@] srcdirs[@] "Quit"
        fi

        option="${allopts[$(($REPLY - 1))]}"
        if [ ! "$option" == "Quit" ]; then
            cd "$option"
            echo "You are now in $PWD"
            return
        else
            echo "Exiting without changing directories"
            return
        fi

    else
        echo "No valid nimrod directories could be found, check installation"
        return
    fi
}
