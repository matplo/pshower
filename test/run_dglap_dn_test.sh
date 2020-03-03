#!/bin/bash

function thisdir()
{
        SOURCE="${BASH_SOURCE[0]}"
        while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
          DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
          SOURCE="$(readlink "$SOURCE")"
          [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
        done
        DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
        echo ${DIR}
}
THISD=$(thisdir)
if [ -e ${THISD}/../scripts/util.sh ]; then
	source ${THISD}/../scripts/util.sh 

	need_help=$(get_opt "help" $@)
	if [ ! -z ${need_help} ]; then
		echo "$0 [--help] [--nev=<integer>]"
		exit 0
	fi

	separator "dglap: ${BASH_SOURCE}"

	echo_warning "using modules from ${THISD}/../modules/dglap/python"
	module use ${THISD}/../modules/dglap/python
	echo_warning "loading dglap module"

	module load dglap
	module list

	nevents=$(get_opt "nev" $@)
	if [ -z ${nevents} ] || [ "xyes" == "x$nevents" ]; then
		nevents=10000
	fi

	echo_info "running for ${nevents} events..."

	example_command="dglap_dn_exe 101 ${nevents} 250 1.5708 1 1"

	separator "dglap: ${example_command}"

	${example_command}

	separator "done."
fi
