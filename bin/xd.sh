#!/bin/bash
#-----------------------------------------------------------------
#Script: xd
#Usage: xd <file>
#-----------------------------------------------------------------
# Input Args #
xd () 
{
    read -r -d '' help <<- EOM
--------------------------------------------------------------------
Usage: xd filename [--edit|-e] [--help]
--------------------------------------------------------------------
    A handy little xdraw wrapper.

    filename: either a .bin file in current directory or any valid
                xdraw option name: grid, dis, en, fl, jpar, polfl,
                prof, pss, stab, t, xt, y, xy, xyeq, con, coneq, or
                vec
    --edit|-e: flag allows you to edit xdraw config file before
                plotting using your editor of choice
    --save|-s: flag will save the draw.in file in your directory
                this way you can have a custom draw.in that is 
                always called in this results directory.  
    --help: displays this helpful description

    Enjoy!
-------------------------------------------------------------------
EOM
    edit_file=
    in_file=
    save_file=
    if [ "$#" -eq "0" ]; then
        echo "$help"; return
    else
        while [ -n "$1" ]; do
            case "$1" in
                --edit | -e ) edit_file="yes"; shift;;
                --save | -s ) save_file="yes"; shift;;
                --help ) echo "$help"; return;;
                * ) filename="$1"; shift;;
            esac
        done
    fi

    filelabel=$(basename "$filename" .bin)
    case "$filelabel" in
        grid ) xdoption="grid";;
        discharge | dis ) xdoption="dis";;
        energy | en ) xdoption="en";;
        flsurf | fl ) xdoption="fl";;
        parallel_current | jpar ) xdoption="jpar";;
        polflux | polfl ) xdoption="polfl";;
        fluxgrid | prof ) xdoption="prof";;
        nimfl | pss ) xdoption="pss";;
        stability | stab ) xdoption="stab";;
        nimhist | t ) xdoption="t";;
        xt_slice | xt ) xdoption="xt";;
        yt_slice | y ) xdoption="y";;
        xy_slice ) read -p "xy or xyeq? :" xdoption;;
        contour ) read -p "con, coneq, or vec? :" xdoption;;
        * ) xdoption="$filelabel";;
    esac

    # Make draw file ##
    new_draw_file=draw$xdoption.in
    # 
    if [ -f "$new_draw_file" ]; then
        in_file="yes"
    else
        cp $NIMSRCDIR/draw/$new_draw_file .
    fi
    #------------------------------------------------------------------
    # Edit if requested
    #------------------------------------------------------------------
    if [ -n "$edit_file" ]; then
        $EDITOR $new_draw_file
    fi
    #------------------------------------------------------------------
    # Run xdraw
    #------------------------------------------------------------------
    xdraw $xdoption
    #------------------------------------------------------------------
    # Remove draw file
    #------------------------------------------------------------------
    if [ -z "$in_file" ]; then
        if [ -z "$save_file" ]; then
            rm $new_draw_file
        fi
    fi
    return
}
