#!/bin/sh

mvn install:install-file -Dfile=lib/BigWig.jar -DgroupId=org.broad -DartifactId=bigwig -Dversion=0.1 -Dpackaging=jar
