
# Start of OS-dependent functions
proc remove_file {fileName} {
    exec rm $fileName
}

proc copy_file {oldName newName} {
    exec cp $oldName $newName
}

proc run_network_program {list_of_args} {
    eval exec ./network $list_of_args
}
# End of OS-dependent functions




proc readFileName {} {
    global filename
    puts "Please enter the filename of your input file" 
    gets stdin filename
    return ""
}

proc checkInputFile {} {
    global filename
    global ending
    if {[string match "*.cuc" $filename]} {
	set ending ".cuc"
    } elseif {[string match "*.cssr" $filename]} {
	set ending ".cssr"
    } elseif {[string match "*.cif" $filename]} {
	set ending ".cif"
    } else {
	puts "ERROR: INPUT FILE MUST BE OF TYPE .cuc or .cssr \nExiting"
	exit
    }
}


readFileName
checkInputFile
set inputFile "ZeoVisInput$ending"
copy_file $filename $inputFile
puts "CREATING zeovis INPUT FILE"

#run_network_program [list "-r" "example.rad" "-zvis" "ZeoVisInput.zvis" $inputFile]
run_network_program [list "-r" "-zvis" "ZeoVisInput.zvis" $inputFile]
#run_network_program [list "-zvis" "ZeoVisInput.zvis" $inputFile]
puts "INPUT FILE SUCCESSFULLY CREATED"

source "ZeoVisInput.zvis"
remove_file "ZeoVisInput.zvis"

source "ZeoVisCommands.tcl"
set num_segments 0
set num_features 0
set num_cages    0
initialize
puts "ENJOY ZEOVIS!" 
