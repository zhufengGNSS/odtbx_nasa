#!/bin/bash
# Driver script for the ODTBX regression testing.  This script updates
# the svn working copies, builds with Apache Maven, then runs the
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
#DEFAULT_EMAIL='allen.brown@emergentspace.com'
#DEFAULT_SMTP_SERVER='mailhost.gsfc.nasa.gov'

# aet - send email to an alias on cain
#DEFAULT_EMAIL=odtbx
#DEFAULT_SMTP_SERVER=localhost

# keb - send mail to the nasa mail list
DEFAULT_EMAIL='odtbx@lists.nasa.gov'
DEFAULT_SMTP_SERVER='mailhost.gsfc.nasa.gov'

export LANG=en_US.UTF-8

if [ -n "$1" ] ; then
        EMAIL=$1
else
        EMAIL=$DEFAULT_EMAIL
fi

if [ -n "$2" ] ; then
        SMTP_SERVER=$2
else
        SMTP_SERVER=$DEFAULT_SMTP_SERVER
fi

PATH=/usr/bin:/bin:$HOME/bin

# I'm sure there's a better way to go up 2 directories...
DIR=`dirname $0`
DIR=`dirname $DIR`
DIR=`dirname $DIR`

LOG_PATH=$DIR/../log
INFO_LOG_FILE=$LOG_PATH/supporting_info.txt

echo "ODTBX Regression Test update,build,run script running at `date`" > $INFO_LOG_FILE
echo "ODTBX Regression Test update,build,run script running at `date`"

# The --config-dir is set for cypher:
echo "Updating source trees."
svn update $DIR/odtbx
svn update $DIR/mice
svn update $DIR/jat
svn update $DIR/jat-lib
svn update $DIR/gmat

# record the svn info:
echo "Gathering information on working copies for log files."
echo "Subversion working copy information:" >> $INFO_LOG_FILE
echo -e "\n" >> $INFO_LOG_FILE

svn info $DIR/odtbx >> $INFO_LOG_FILE
svn info $DIR/mice >> $INFO_LOG_FILE
svn info $DIR/jat >> $INFO_LOG_FILE
svn info $DIR/jat-lib >> $INFO_LOG_FILE
svn info $DIR/gmat >> $INFO_LOG_FILE

echo -e "\n" >> $INFO_LOG_FILE
echo "Maven build results:" >> $INFO_LOG_FILE
echo "Building checkout with Maven."
mvn -f $DIR/jat/maven/pom.xml clean compile >> $INFO_LOG_FILE

echo "Running Matlab."
matlab -nodisplay -r "\
        basePath = '$DIR/odtbx';\
        addpath(basePath);\
        startup();\
        regressionTesting('$LOG_PATH','$EMAIL','$SMTP_SERVER','$INFO_LOG_FILE');\
        quit;\
"

