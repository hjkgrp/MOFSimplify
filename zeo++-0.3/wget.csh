#!/bin/csh

set num_args = 1

if( $#argv != $num_args ) then

        echo "This script downloads net files (cgd format) from the RCSR database, saving them in a new directory called nets."
        echo "Requires $num_args arguments: net_name_string. Please try again."
        exit

endif

set net = $1

wget http://rcsr.anu.edu.au/webfiles/systreinput/$net.cgd

mkdir nets
mv $net.cgd nets/

