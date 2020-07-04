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
if [ -e ${THISD}/util.sh ]; then
	source ${THISD}/util.sh 
	separator "dglap: ${BASH_SOURCE}"

	echo_warning "using modules from ${THISD}/../modules/dglap/python"
	module use ${THISD}/../modules/dglap/python
	echo_warning "loading dglap module"

	module load dglap
	module list

	separator "done."
fi
