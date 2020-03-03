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
export -f thisdir

THISD=$(thisdir)
source ${THISD}/util.sh

function make_module_package()
{
	dirinst=${1}
	module_name=$(basename ${dirinst})
	package_name=$(basename ${dirinst})
	[ ! -z ${2} ] && package_name=${2}
	[ ! -z ${3} ] && package_version=${3}

	[ "x$(get_opt "python2" $@)" == "xyes" ] && X_USER_PYTHON_VERSION=python2
	[ "x$(get_opt "python3" $@)" == "xyes" ] && X_USER_PYTHON_VERSION=python3
	[ -z ${X_USER_PYTHON_VERSION} ] && X_USER_PYTHON_VERSION=python

	X_PYTHON_EXECUTABLE=$(which ${X_USER_PYTHON_VERSION})
	X_PYTHON_CONFIG_EXECUTABLE=$(which ${X_USER_PYTHON_VERSION}-config)
	# if [ -f "${X_PYTHON_EXECUTABLE}" ] && [ -f "${X_PYTHON_CONFIG_EXECUTABLE}" ]; then

	if [ -d ${dirinst} ]; then
		mkdir -p ${THISD}/../modules/${package_name}
		modulefiledir=$(abspath_python_expand ${THISD}/../modules/${package_name})
		[ ! -z ${X_USER_PYTHON_VERSION} ] && modulefiledir=${modulefiledir}/${X_USER_PYTHON_VERSION}
		[ ! -z ${package_name} ] && modulefiledir=${modulefiledir}/${package_name}
		modulefile="${modulefiledir}/${module_name}"
		[ ! -z ${package_version} ] && modulefile="${modulefiledir}/${package_version}"
		mkdir -p ${modulefiledir}
		separator "making ${package_name} module ${modulefile}"
		[ -f ${modulefile} ] && warning "removing ${modulefile}" && rm -f ${modulefile}

		X_PYTHON_VERSION=$(${X_PYTHON_EXECUTABLE} --version 2>&1 | cut -f 2 -d' ' | cut -f 1-2 -d.)
		X_PYTHON_BIN_DIR=$(dirname ${X_PYTHON_EXECUTABLE})
		X_PYTHON_INCLUDE_DIR=$(${X_PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
		X_PYTHON_LIBDIR=$(${X_PYTHON_EXECUTABLE} -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
		X_PYTHON_NUMPY_INCLUDE_DIR=$(${X_PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())")
		if [ ! -d "${X_PYTHON_NUMPY_INCLUDE_DIR}" ]; then
			error "missing numpy and/or headers "
			error "we are strongly relying on numpy - numpy AND headers must be installed/accessible - anything below does not matter..."
			error "try: pip install numpy [--user]"
			return 0
		fi
		X_PYTHON_LIBS=$(${X_PYTHON_CONFIG_EXECUTABLE} --libs)
		X_PYTHON_LIBS_LINK="-L${X_PYTHON_LIBDIR} ${X_PYTHON_LIBS}"
		X_PYTHON_CONFIG_LDFLAGS=$(${X_PYTHON_CONFIG_EXECUTABLE} --ldflags)
		X_PYTHON_CONFIG_INCLUDES=$(${X_PYTHON_CONFIG_EXECUTABLE} --includes)
		X_PYTHON_SETUP=TRUE

		bin_path="${dirinst}/bin"
		lib_path="${dirinst}/lib"
		lib64_path="${dirinst}/lib64"
		python_path="${dirinst}/lib/python${X_PYTHON_VERSION}/site-packages"
		python_path64="${dirinst}/lib64/python${X_PYTHON_VERSION}/site-packages"

		setenv_module ${modulefile} ${package_name}DIR ${dirinst}
		setenv_module ${modulefile} ${package_name}_DIR ${dirinst}
		setenv_module ${modulefile} ${package_name}_ROOT ${dirinst}

		setenv_module ${modulefile} ${package_name}_INCLUDE_DIR ${dirinst}/include

		[ $(os_linux) ] && add_path_module ${modulefile} PATH ${bin_path}
		[ $(os_darwin) ] && add_path_module ${modulefile} PATH ${bin_path}

		for sp in ${lib_path} ${lib64_path} ${python_path} ${python_path64}
		do
			[ $(os_linux) ] && add_path_module ${modulefile} LD_LIBRARY_PATH ${sp}
			[ $(os_darwin) ] && add_path_module ${modulefile} DYLD_LIBRARY_PATH ${sp}
		done

		for sp in ${python_path} ${python_path64} ${lib_path} ${lib64_path}
		do
			[ $(os_linux) ] &&  add_path_module ${modulefile} PYTHONPATH ${sp}
			[ $(os_darwin) ] &&  add_path_module ${modulefile} PYTHONPATH ${sp}
		done

		if [ -d "${bin_path}" ]; then
			add_path_module ${modulefile} PATH ${bin_path}
		fi

		if [ ! -z ${X_PYTHON_MODULE_LOADED} ]; then
			echo "prereq ${X_PYTHON_MODULE_LOADED}" >> ${modulefile}
		fi

	else
		error "${dirinst} does not exists - no module generation"
	fi
}
export -f make_module_package

packagedir=$(get_opt "dir" $@)
packagename=$(get_opt "name" $@)
packageversion=$(get_opt "version" $@)
if [ -d "${packagedir}" ]; then
	separator "make_modules.sh :: package module"
	make_module_package ${packagedir} ${packagename} ${packageversion}
	separator "make_modules.sh - done"
fi
