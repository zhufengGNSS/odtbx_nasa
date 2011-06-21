#!/bin/csh

############################################################################
# Find the path to java.  First look in JAT_JAVA, then JAVA_HOME, then in PATH
if ((! $?JAT_JAVA) && ($?JAVA_HOME)) then
  if (-r ${JAVA_HOME}/bin/java) then
    setenv JAT_JAVA ${JAVA_HOME}
  endif
endif

if (! $?JAT_JAVA) then
  foreach dir ($path)
    if (-r $dir/java) then
      setenv JAT_JAVA $dir/..
      break
    endif
  end
endif

if (! $?JAT_JAVA) then
    echo -n "Java could not be found.  Please specify "
    echo "location in environment variable JAT_JAVA."
    exit 1
else if (! -r ${JAT_JAVA}/bin/java) then
    echo -n "JAT_JAVA environment variable is set, but does not "
    echo -n "seem to point to Java installation.  JAT_JAVA/bin/java"
    echo " not found"
    exit 1
endif

#############################################################
# Determine the JAT root directory.  We look first in the
# JAT_HOME environment variable then in the current directory
if (! $?JAT_HOME) then
  if (-r ./jat.sh) then
    setenv JAT_HOME .
  endif
endif

if (! $?JAT_HOME) then
    echo -n "JAT could not be found.  Please run from JAT's root "
    echo "directory or specify location in environment variable JAT_HOME."
    exit 1
else if !(-r ${JAT_HOME}/jat.sh) then
    echo -n "JAT_HOME environment variable is set, but does not "
    echo -n "seem to point to JAT system.  JAT_HOME/jat.sh "
    echo "not found."
    exit 1
endif

########################################################################
# Setup the classpath - the JAT source tree and everything in external_jar

setenv JAT_CLASSPATH ${JAT_HOME}

# Add all jars in external_jar/ directory to classpath.
foreach lib (`find ${JAT_HOME}/external_jar -name '*.jar'`)
  setenv JAT_CLASSPATH ${JAT_CLASSPATH}:$lib
end

###################################################################
# Now figure out the Java class to execute.  Defaults to CEVSim.

if ($#argv == 0) then
  setenv _MAIN_CLASS jat.sim.CEVSim
else
  setenv _MAIN_CLASS $1
endif

#########################################################################
# Now run the command

${JAT_JAVA}/bin/java -classpath ${JAT_CLASSPATH} $_MAIN_CLASS $argv[2-$#argv]

unsetenv _MAIN_CLASS
