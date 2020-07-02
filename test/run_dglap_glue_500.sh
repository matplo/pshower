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
	separator "dglap: ${BASH_SOURCE}"

	echo_warning "using modules from ${THISD}/../modules/dglap/python"
	module use ${THISD}/../modules/dglap/python
	echo_warning "loading dglap module"

	module load dglap
	module list

	example_command="python3 ${THISD}/run_dglap_dn_test.py --initfile Initialize_Parton_Shower_glue_500.txt -Q 500 --flavor 0 --gg-only $@"

	separator "dglap: ${example_command}"

	${example_command}

	separator "done."
fi
