#!/bin/bash
# Driver script for the ODTBX regression testing.  This script updates
# the git working copies, builds with Apache Maven, then runs the
# MATLAB ODTBX regressionTesting.m in batch mode.  The results are 
# emailed out.  Log files are placed in a logs directory under the pwd
# when this script is run.
#
# Arguments: 
# (optional) The recipient's email address, if not provided then
#            the DEFAULT_EMAIL below is used
# (optional) The email server to use, if not provided then
#            the DEFAULT_SMTP_SERVER below is used
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

# Set these to enable email, if no email is desired, then set each to
# ''.
# keb - send mail to the nasa mail list
DEFAULT_EMAIL='odtbx@lists.nasa.gov'
DEFAULT_SMTP_SERVER='mailhost.gsfc.nasa.gov'

export LANG=en_US.UTF-8

# Set outgoing email address if specified
if [ -n "$1" ] ; then
        EMAIL=$1
else
        EMAIL=$DEFAULT_EMAIL
fi

# Set outgoing email server if specified
if [ -n "$2" ] ; then
        SMTP_SERVER=$2
else
        SMTP_SERVER=$DEFAULT_SMTP_SERVER
fi

# Add paths to built-in commands (e.g. 'date' and 'tee')
PATH=/usr/bin:/bin:$HOME/bin

ODTBX="$HOME/projects/odtbxsync-git/odtbx"
echo "ODTBX Directory = $ODTBX"

VENDOR="$HOME/projects/odtbxsync-git/vendor"
echo "ODTBX Vendor Directory = $VENDOR"

LOG_PATH="$HOME/log"
INFO_LOG_FILE="$LOG_PATH/supporting_info.txt"

echo "ODTBX Regression Test update,build,run script running at `date`" | tee $INFO_LOG_FILE

# rmathur
# The process of merging public and internal Git repositories is documented
# in ~/projects/README-GIT.txt .

# Get changes to internal master branch
echo -e "\n***Pulling internal master branch." | tee -a $INFO_LOG_FILE
git pull internal master >> $INFO_LOG_FILE 2>&1
giterr=${PIPESTATUS[0]} # Get error code of first command above (git pull)
if [ "$giterr" != 0 ]; then
	echo "git pull internal master failed!" >> $INFO_LOG_FILE
	exit
fi

# Get changes to public master branch
echo -e "\n***Pulling public master branch." | tee -a $INFO_LOG_FILE
git pull public master >> $INFO_LOG_FILE 2>&1
giterr=${PIPESTATUS[0]} # Get error code of first command above (git pull)
if [ "$giterr" != 0 ]; then
	echo "git pull public master failed!" >> $INFO_LOG_FILE
	exit
fi

# Push merged master branch to both repos
echo -e "\n***Pushing master branch." | tee -a $INFO_LOG_FILE
git push internal master >> $INFO_LOG_FILE 2>&1
git push public master >> $INFO_LOG_FILE 2>&1

# Get changes to internal develop branch
echo -e "\n***Pulling internal develop branch." | tee -a $INFO_LOG_FILE
git pull internal develop >> $INFO_LOG_FILE 2>&1
giterr=${PIPESTATUS[0]} # Get error code of first command above (git pull)
if [ "$giterr" != 0 ]; then
	echo "git pull internal develop failed!" >> $INFO_LOG_FILE
	exit
fi

# Get changes to public develop branch
echo -e "\n***Pulling public develop branch." | tee -a $INFO_LOG_FILE
git pull public develop >> $INFO_LOG_FILE 2>&1
giterr=${PIPESTATUS[0]} # Get error code of first command above (git pull)
if [ "$giterr" != 0 ]; then
	echo "git pull public develop failed!" >> $INFO_LOG_FILE
	exit
fi

# Push merged develop branch to both repos
echo -e "\n***Pushing develop branch." | tee -a $INFO_LOG_FILE
git push internal develop >> $INFO_LOG_FILE 2>&1
git push public develop >> $INFO_LOG_FILE 2>&1

# Build JAT using Apache Maven
echo -e "\nBuilding JAT with Apache Maven." | tee -a $INFO_LOG_FILE
mvn -f "$VENDOR/Jat/maven/pom.xml" clean compile >> $INFO_LOG_FILE 2>&1

echo -e "\nRunning ODTBX regression tests in Matlab." | tee -a $INFO_LOG_FILE
matlab -nodisplay -r "\
        basePath = '$ODTBX';\
        addpath(basePath);\
        startup();\
        regressionTesting('$LOG_PATH','$EMAIL','$SMTP_SERVER','$INFO_LOG_FILE');\
        quit;\
"

