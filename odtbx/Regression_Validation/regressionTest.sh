#!/bin/bash
# Driver script for the ODTBX regression testing.  This script updates
# the git working copies, builds with Apache Maven, then runs the
# MATLAB ODTBX regressionTesting.m in batch mode.  The results are 
# emailed out.  Log files are placed in a logs directory under the pwd
# when this script is run.  If nothing changed in the develop branch,
# then the regression test is not run.
#
# Arguments: 
# [-f]       Force-run the regression test, even if nothing changed in
#            the develop branch
# [-e email] The recipient's email address (default odtbx@lists.nasa.gov)
# [-s smtp]  The smtp server to use (default mailhost.gsfc.nasa.gov)
# [-m cmd]   Run the specified MATLAB command verbatim (default regressionTesting)
#
# ODTBX: Orbit Determination Toolbox
# 
# Copyright (c) 2003-2011 United States Government as represented by the
# administrator of the National Aeronautics and Space Administration. All
# Other Rights Reserved.
# 
# This file is distributed "as is", without any warranty, as part of the
# ODTBX. ODTBX is free software; you can redistribute it and/or modify it
# under the terms of the NASA Open Source Agreement, version 1.3 or later.
# 
# You should have received a copy of the NASA Open Source Agreement along
# with this program (in a file named License.txt); if not, write to the 
# NASA Goddard Space Flight Center at opensource@gsfc.nasa.gov.

# Sends email using MATLAB (since that supports SMTP servers)
function send_email
{
  email=$1
  subject=$2
  message=$3
  server=$4
  attachfile=$5

  matlab -nodisplay -r "\
    setpref('Internet','SMTP_Server','$server');\
    setpref('Internet','E_mail','odtbx-builder@emergentspace.com');\
    sendmail('$email','$subject','$message','$attachfile');\
    quit;"
}

# Add paths to built-in commands (e.g. 'date' and 'tee')
PATH=/usr/bin:/bin:$HOME/bin

export LANG=en_US.UTF-8

# Set paths to needed files/folders
ODTBX="$HOME/projects/odtbxsync-git/odtbx"
VENDOR="$HOME/projects/odtbxsync-git/vendor"
LOG_PATH="$HOME/log"
INFO_LOG_FILE="$LOG_PATH/supporting_info.txt"

echo \e "\nODTBX Regression Test script running on `date`" | tee $INFO_LOG_FILE
echo "ODTBX Directory = $ODTBX"
echo "ODTBX Vendor Directory = $VENDOR"

# Set default variable values
EMAIL="odtbx@lists.nasa.gov"
SMTP_SERVER="mailhost.gsfc.nasa.gov"
FORCE=false
MATLABCMD="regressionTesting('$LOG_PATH','$EMAIL','$SMTP_SERVER','$INFO_LOG_FILE')"

# Parse input arguments
USAGE="Usage: regressionTest.sh [-f] [-h] [-e emailaddr] [-s smtpserver] [-m matlabcmd]"
while getopts fhe:s:m: OPT; do
	case "$OPT" in
	f)
		FORCE=true
		;;
	h)
		echo $USAGE
		exit 0
		;;
	e)
		EMAIL=$OPTARG
		;;
	s)
		SMTP_SERVER=$OPTARG
		;;
	m)
		MATLABCMD="$OPTARG"
		;;
	\?)
		# getopts issues an error message
		echo $USAGE | tee -a $INFO_LOG_FILE
		exit 1
		;;
	esac
done
shift `expr $OPTIND - 1`

# Make sure there are no other arguments
if [ $# != 0 ]; then
	echo "Invalid input arguments to regressionTest.sh"
	echo $USAGE | tee -a $INFO_LOG_FILE
	exit 1
fi

# rmathur
# The process of merging public and internal Git repositories is documented
# in ~/projects/README-GIT.txt .
git_subj="ODTBX Git failure"
git_msg="There was an error while merging Git repositories. See the attached log."

# Switch to the master branch, so that changes to master can be 
# properly merged
echo -e "\n***Switching to master branch." | tee -a $INFO_LOG_FILE
git checkout master >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git checkout master
if [ $gitstatus != 0 ]; then
	msg="git checkout master failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus
fi

# Merge changes from internal master branch
echo -e "\n***Pulling internal master branch." | tee -a $INFO_LOG_FILE
git pull internal master >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git pull
if [ $gitstatus != 0 ]; then
	msg="git pull internal master failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus
fi

# Merge changes from public master branch
echo -e "\n***Pulling public master branch." | tee -a $INFO_LOG_FILE
git pull public master >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git pull
if [ $gitstatus != 0 ]; then
	msg="git pull public master failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus
fi

# Push merged master branch to both repos
echo -e "\n***Pushing master branch." | tee -a $INFO_LOG_FILE
git push internal master >> $INFO_LOG_FILE 2>&1
git push public master >> $INFO_LOG_FILE 2>&1

# Switch to the develop branch, so that changes to develop can be 
# properly merged
echo -e "\n***Switching to develop branch." | tee -a $INFO_LOG_FILE
git checkout develop >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git checkout develop
if [ $gitstatus != 0 ]; then
	msg="git checkout develop failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus
fi

# Check develop commit before merging internal & public repos
prevcommit=`git log -1 --pretty=format:%H`
echo -e "\nBefore merge, develop commit = $prevcommit" | tee -a $INFO_LOG_FILE

# Merge changes from internal develop branch
echo -e "\n***Pulling internal develop branch." | tee -a $INFO_LOG_FILE
git pull internal develop >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git pull
if [ $gitstatus != 0 ]; then
	msg="git pull internal develop failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus $gitstatus
fi

# Merge changes from public develop branch
echo -e "\n***Pulling public develop branch." | tee -a $INFO_LOG_FILE
git pull public develop >> $INFO_LOG_FILE 2>&1
gitstatus=$? # Get error code of git pull
if [ $gitstatus != 0 ]; then
	msg="git pull public develop failed! Git status = $gitstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "$git_subj" "$git_msg" $SMTP_SERVER "$INFO_LOG_FILE"
	exit $gitstatus
fi

# Check develop commit after merging internal & public repos
nextcommit=`git log -1 --pretty=format:%H`
echo -e "\nAfter merge, develop commit = $nextcommit" | tee -a $INFO_LOG_FILE

# If the develop branch commit hash didn't change after merging, then
# nothing was done during the merge, so there's nothing to test.
# But if the "-f" flag was given, then always run the regression test.
if [ "$prevcommit" == "$nextcommit" ] && ! $FORCE ; then
	echo "Branch develop has not changed. Nothing to test." >> $INFO_LOG_FILE
	exit 0
fi

# Push merged develop branch to both repos
echo -e "\n***Pushing develop branch." | tee -a $INFO_LOG_FILE
git push internal develop >> $INFO_LOG_FILE 2>&1
git push public develop >> $INFO_LOG_FILE 2>&1

# Build JAT using Apache Maven
echo -e "\nBuilding JAT with Apache Maven." | tee -a $INFO_LOG_FILE
mvn -f "$VENDOR/Jat/maven/pom.xml" clean compile >> $INFO_LOG_FILE 2>&1
mvnstatus=$? # Get error code of mvn compile
if [ $mvnstatus != 0 ]; then
	msg="mvn clean compile failed! Maven status = $mvnstatus"
	echo $msg >> $INFO_LOG_FILE
	send_email $EMAIL "ODTBX Maven failure" "There was an error while building JAT with Maven. See the attached log." $SMTP_SERVER "$INFO_LOG_FILE"
	exit $mvnstatus
fi

# Run Matlab regression test suite
echo -e "\nRunning ODTBX regression tests in Matlab." | tee -a $INFO_LOG_FILE
matlab -nodisplay -r "\
        basePath = '$ODTBX';\
        addpath(basePath);\
        startup();\
	$MATLABCMD;\
        quit;"
