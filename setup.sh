#!/bin/bash
# Install script for nimpy_wipal "package"
# To run this script and source changes, run the following command
# ./setup.sh && source ~/.bash_profile
#
# Primary objective is to set up access to python modules 
# for parsing nimrod binary/vtk files and plotting fields
# from simulations. (Not really taken care of yet)
#
# Secondary objective is setting environment variables and
# building/sourcing shell functions nimcd,nimupdate,nimkill,
# and xd for easy navigation between builds, result directories
# and active nimrod processes.
#
###############################################################

# Start with some dependency checks
# 1.) Make sure user is running setup.sh from nimpy_wipal root
# 2.) Define some important variables
# 3.) Check if ~/.bashrc exists. If no, create it.
# 4.) Check if ~/.bash_profile exists. If no, create it.
# 5.) Check if ~/bin exists. If no, create it.
if [ "$(basename $(pwd))" != "nimpy_WiPAL" ]; then
	echo "ERROR: $(pwd) is not the root directory for the nimpy_WiPAL package"
	exit 1
fi	

NIMPYHOME=$PWD
BASH_RC=${HOME}/.bashrc
BASH_PROFILE=${HOME}/.bash_profile
HEADER="## Added by nimpy_wipal ##"
TAILER="## End of nimpy_wipal adds ##"

# Check for ~/.bashrc, create if necessary
if [ ! -f $BASH_RC ]; then
	echo "Creating $BASH_RC..."
	read -r -d '' brctxt<<-EOF
	# .bashrc
	
	# Source global definitions
	if [ -f /etc/bashrc ]; then
	\t. /etc/bashrc
	fi
	
	# User defined functions and variables
	EOF
	printf "$brctxt\n\n" > $BASH_RC
fi    

# Check for ~/.bash_profile and make sure it sources ~/.bashrc
if [ ! -f $BASH_PROFILE ]; then
	echo "Creating $BASH_PROFILE..."
	read -r -d '' bprotxt<<-EOF
	# .bash_profile

	# Get the aliases and functions
	if [ -f ~/.bashrc ]; then
	\t. ~/.bashrc
	fi

	# User specific environment and startup programs
	EOF
	printf "$bprotxt\n\n" > $BASH_PROFILE
else
	grep -E -q "source|\. ~/\.bashrc|${HOME}/\.bashrc" $BASH_PROFILE
	if [ $? -eq 0 ]; then
		:
	else
		echo "Creating backup of ~/.bash_profile before editing..."
		cp $BASH_PROFILE ${BASH_PROFILE}-nimpy.bak
		echo "Prepending sourcing of ~/.bashrc to ~/.bash_profile..."
		read -r -d '' brcsrctxt<<-EOF
		# .bash_profile

		# Get the aliases and functions
		if [ -f ~/.bashrc ]; then
		\t. ~/.bashrc
		fi

		# User specific environment and startup programs
		EOF
		bpro_old=$(cat $BASH_PROFILE) 
		printf "$brcsrctxt\n\n" > $BASH_PROFILE && printf "$bpro_old" >> $BASH_PROFILE
	fi

fi    

# Check for ~/bin directory
if [ ! -d "${HOME}/bin" ]; then
    echo "Creating ${HOME}/bin directory..." 
    mkdir ${HOME}/bin
fi

# Begin interactive prompt for install
# 1.) Get top level nimrod directory from user
# 2.) Find a valid nimrod build to initialize the install
# 3.) Get path to xdraw executable automatically or from user
# 4.) If EDITOR isn't defined get it from the user

echo ""
read -e -p "Please specify nimrod top level directory: " NIMTOP
NIMTOP="${NIMTOP/#\~/$HOME}"
NIMTOP="${NIMTOP%/}"
NIMSRCPATHS=$(echo $(find $NIMTOP -type d -print | grep .*/rundir | sed 's/\/rundir.*//g') | sed 's/ /:/g')
IFS=':' read -ra NIMDIRS <<< "$NIMSRCPATHS"
if [ "${#NIMDIRS[@]}" -eq "0" ]; then
    echo "ERROR: Unable to find a valid nimrod source directory under $NIMTOP"
    exit 1
else
    echo ""
    options=( "${NIMDIRS[@]}" "Abort install" )
    echo "Please select one of the following nimrod source directories"
    echo "to initialize NIMSRCDIR"
    echo "------------------------------------------------------------"
    select opt in "${options[@]}"; do
        if [ 1 -le "$REPLY" ] && [ "$REPLY" -le "${#options[@]}" ]; then
            break;
        else
            echo "Incorrect Input: Select a number 1-${#options[@]}"
        fi
    done
fi

if [ "$opt" == "Abort install" ];then
    echo "Aborting install..."
    exit 1
else
    NIMSRCDIR="$opt"
fi

EXISTXDRAW=$(which xdraw 2>/dev/null)
if [ -x "$EXISTXDRAW" ];then
	echo "xdraw executable found at $EXISTXDRAW"
else
	echo "Could not find xdraw executable"
	read -e -p "Please provide the path to the xdraw executable: " XDRAWPATH
	XDRAWPATH="${XDRAWPATH/#\~/$HOME}"
	if [ -f "$XDRAWPATH" -a -x "$XDRAWPATH" -a "$(basename $XDRAWPATH)" == "xdraw" ];then
		:
	else
		echo "ERROR: $XDRAWPATH is not an appropriate xdraw executable, aborting install..."
		exit 1
	fi
fi

if [ -z "$EDITOR" ]; then
    read -e -p "Please provide your text editor of choice: " TEXTEDITOR
fi

python_mod_path=$NIMPYHOME

# Begin install
# 1.) create backup ~/.bashrc file before editing
# 2.) uninstall previous version of nimpy_wipal
# 3.) append necessary environment variables/aliases to ~/.bashrc
# 4.) create nimrod and xdraw softlinks in ~/bin

# Create backup file for ~/.bashrc before making changes
cp $BASH_RC $BASH_RC-nimpy.bak

# Delete previous install if possible
sed -i /"$HEADER"/,/"$TAILER"/d "$BASH_RC"

# Append necessary commands to ~/.bashrc
cat <<EOF >> $BASH_RC

$HEADER
export NIMPYHOME=$NIMPYHOME
export NIMTOP=$NIMTOP
export NIMSRCDIR=$NIMSRCDIR
export NIMGRIDDIR=$NIMSRCDIR/results/grids
export NIMRESULTPATHS=\$(echo \$(find ${NIMTOP} -type d -print | grep results) | sed 's/ /:/g')
export NIMSRCPATHS=\$(echo \$(find ${NIMTOP} -type d -print | grep .*/rundir | sed 's/\/rundir.*//g') | sed 's/ /:/g')
export DRAWDIR=${NIMPYHOME}/draw
source ${NIMPYHOME}/bin/nimcd.sh
source ${NIMPYHOME}/bin/nimupdate.sh
source ${NIMPYHOME}/bin/nimkill.sh
source ${NIMPYHOME}/bin/xd.sh
source ${NIMPYHOME}/bin/nimgrid.sh
alias cdnimpy="cd ${NIMPYHOME}"
EOF
grep -q "^export PYTHONPATH=.*$python_mod_path" $BASH_RC
if grep -q "^export PYTHONPATH=.*" $BASH_RC; then
	if grep -q "^export PYTHONPATH=.*$python_mod_path" $BASH_RC; then
		:
	else
		echo "export PYTHONPATH=\$PYTHONPATH:$python_mod_path" >> $BASH_RC
	fi
else
	echo "export PYTHONPATH=$python_mod_path" >> $BASH_RC
fi
([ -n "$EDITOR" ] && grep -q "^export EDITOR=" $BASH_RC) || (echo "export EDITOR=$TEXTEDITOR" >> $BASH_RC && echo "exported EDITOR choice to $BASHRC")
grep -q "^export PATH=.*\(\${HOME}\|${HOME}\)/bin" $BASH_RC || (echo "export PATH=${HOME}/bin:\$PATH" >> $BASH_RC && echo "Prepended ${HOME}/bin to $BASH_RC")
echo "$TAILER" >> $BASH_RC

## Create softlinks to nimrod and xdraw executables in ~/bin
ln -sfn ${NIMSRCDIR}/nimset/stitch ${HOME}/bin/stitch
ln -sfn ${NIMSRCDIR}/rundir/nimplot ${HOME}/bin/nimplot 
ln -sfn ${NIMSRCDIR}/rundir/nimrod ${HOME}/bin/nimrod
ln -sfn ${NIMSRCDIR}/rundir/nimeq ${HOME}/bin/nimeq 
ln -sfn ${NIMSRCDIR}/rundir/nimset ${HOME}/bin/nimset 
ln -sfn ${NIMSRCDIR}/rundir/fluxgrid ${HOME}/bin/fluxgrid
ln -sfn ${NIMSRCDIR}/rundir/nimfl ${HOME}/bin/nimfl 
ln -sfn ${NIMSRCDIR}/rundir/sol ${HOME}/bin/sol
ln -sfn ${NIMSRCDIR}/rundir/xlog ${HOME}/bin/xlog 
ln -sfn ${NIMPYHOME}/bin/nimh5.sh ${HOME}/bin/nimh5
ln -sfn ${NIMPYHOME}/bin/nimen.sh ${HOME}/bin/nimen
[ -x "$EXISTXDRAW" ] || ln -sfn ${XDRAWPATH} ${HOME}/bin/xdraw
echo
echo "Printing contents of ${HOME}/bin..."
echo "-----------------------------------"
ls -lah --color ${HOME}/bin
echo
echo "nimpy_wipal install complete! Checkout $BASH_RC to see changes." 
echo "Please run 'source ~/.bash_profile' for these changes to take effect"
exit 0
