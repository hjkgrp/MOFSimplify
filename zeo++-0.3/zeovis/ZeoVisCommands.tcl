set axs_node_color   green
set inaxs_node_color red
set color_list      {blue orange orange  lime    cyan  purple  orange red}
set element_types   {atom node  vorcell unitcell face  channel vornet cage segment feature}
set colorable_types {atom node  vorcell unitcell face  channel vornet cage segment feature}
set sizeable_types  {atom node}
set available_colors [colorinfo colors]
set draw_count 0

proc initialize {} {
    global list_of_actions
    global element_types colorable_types sizeable_types
    global available_colors
    global color_list
    global draw_count
    
    set draw_count 0

    set index [lsearch $available_colors black]
    if {$index != -1} {
	set available_colors [lreplace $available_colors $index $index]
    }

    foreach type $element_types {
	global "num_$type\s"
	global "draw_$type\_log"
	global "draw_$type\_count"
 
	set num_elements [subst "\$num_$type\s"]
	for {set i 0} {$i < $num_elements} {incr i} {
	    set "draw_$type\_log($i)" [dict create]
	}
    }
    
    set colorIndex 0
    foreach type $colorable_types {
	global "num_$type\s"
	global "old_$type\_colors"
	global "$type\Colors"
	set num_elements [subst "\$num_$type\s"]
	set newColor [lindex $color_list $colorIndex]
	
	puts "COLORED $type $newColor"

	for {set i 0} {$i < $num_elements} {incr i} {
	    set "$type\Colors($i)" $newColor
	    set "old_$type\_colors($i)" {}
	}
	incr colorIndex
    }

    foreach type $sizeable_types {
	global "num_$type\s"
	global "old_$type\_radii"
	set num_elements [subst "\$num_$type\s"]
	for {set i 0} {$i < $num_elements} {incr i} {
	    set "old_$type\_radii($i)" {}
	}
    }
    
    set list_of_actions {}
    set probe_radius -1
    return ""
}

proc clear_env {} {
    reset size   atom    all
    reset color  atom    all
    reset size   node    all
    reset color  node    all
    reset color  vorcell all
    reset color  face    all
    initialize
    draw  delete         all
    return "ENVIRONMENT CLEARED"
} 

proc translate_point {pt da db dc} {
    if {[llength $pt] != 3} {
	return "ERROR: A POINT MUST BE A LIST CONTAINING THREE NUMBERS.";
    }
    global uc_a_vector uc_b_vector uc_c_vector 
    set newX [expr [lindex $pt 0] + $da * [lindex $uc_a_vector 0] + $db * [lindex $uc_b_vector 0] + $dc * [lindex $uc_c_vector 0] ] 
    set newY [expr [lindex $pt 1] + $da * [lindex $uc_a_vector 1] + $db * [lindex $uc_b_vector 1] + $dc * [lindex $uc_c_vector 1] ] 
    set newZ [expr [lindex $pt 2] + $da * [lindex $uc_a_vector 2] + $db * [lindex $uc_b_vector 2] + $dc * [lindex $uc_c_vector 2] ] 
    return [list $newX $newY $newZ]
}

proc translate_command {comm da db dc} {
    set result ""
    foreach element $comm {
	if {[llength $element] == 3} {
	    set result [concat $result  "{" [translate_point $element $da $db $dc] "}"]
	} else {
	    set result [concat $result $element]
	}
    }
    set result [concat $result ""]
    return $result;
}

proc show_trans {da db dc type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    set suffix "s"
    helperInvoker "show_trans_helper" $type $indices "draw_" "$type$suffix DRAWN WITH A TRANSLATION OF ($da $db $dc)" [list $da $db $dc]
}

proc show_trans_helper {type index extra_args} {
    global "$type\s"
    set da [lindex $extra_args 0]
    set db [lindex $extra_args 1]
    set dc [lindex $extra_args 2]

    set item_name "$type\s($index)"
    set old_record [subst "\$$item_name"]
    
    set result {}
    foreach line $old_record {
	lappend result [translate_command $line $da $db $dc]
    }

    set $item_name $result
    show_helper $type $index $extra_args
    set $item_name $old_record
    return "$type #$index DRAWN WITH A TRANSLATION OF ($da $db $dc)"
}

proc checkElementTypes {type list_of_types} {
    set numTypes [llength $list_of_types]
    for {set i 0} {$i < $numTypes} {incr i} {
	if {[string compare [lindex $list_of_types $i] $type] == 0} {
	    return "APPROPRIATE TYPE SPECIFIED"
	}
    }
    error "ERROR: This procedure only accepts these types: $list_of_types";
}

proc checkElementIndex {type index} {
    global "num_$type\s"
    set num_elements [subst "\$num_$type\s"] 
    if {$num_elements == 0} {
	error "ERROR: THERE AREN'T ANY $type\s"
    } elseif {$index >= $num_elements} {
	error "ERROR: THERE ARE ONLY $num_elements $type\s INDEXED FROM 0 TO [expr $num_elements - 1]"
    }
    return "OK"
}

proc is_list {element} {
    return [expr [string compare [lindex $element 0] $element] != 0]
}

proc checkListIndices {lst type} {
    for {set i 0} {$i < [llength $lst]} {incr i} {
	checkElementIndex $type [lindex $lst $i]
    }
}


proc undo {} {
    global list_of_actions 
    if {[llength $list_of_actions] == 0} {
	return "NOTHING TO BE UNDONE" 
    }
    set last_action [lindex $list_of_actions [expr [llength $list_of_actions] - 1]]
    set list_of_actions [lreplace $list_of_actions end end]
    set numPairs [expr [llength $last_action] / 2] 
    for {set i 0} {$i < $numPairs} {incr i} {
	set element_name  [lindex $last_action [expr 2 * $i]]
	set indices [lindex $last_action [expr 1 + [expr 2 * $i]]]
	
	if {[string first "color_" $element_name] != -1} {
	    set type [string range $element_name 6 end];
	    revert color $type $indices
	} elseif {[string first "size_" $element_name] != -1} {
	    set type [string range $element_name 5 end];
	    revert size $type $indices
	} elseif {[string first "hide_" $element_name] != -1} {
	    set type [string range $element_name 5 end];
	    if {[is_list $indices]} {
		set disp [lindex $indices end]
		set real_indices [lreplace $indices end end]
	    } else {
		set real_indices [string range $indices 0 [expr [string first "_" $indices] - 1]]
		set disp [string range $indices [expr [string first "_" $indices] + 1] end]
	    }
	    show_trans [expr [string index $disp 0] - 4] [expr [string index $disp 1] - 4]  [expr [string index $disp 2] - 4] $type $real_indices
	    set list_of_actions [lreplace $list_of_actions end end]
	} elseif {[string first "draw_" $element_name] != -1} {
	    set type [string range $element_name 5 end];
	    if {[is_list $indices]} {
		set disp [lindex $indices end]
		set real_indices [lreplace $indices end end]
	    } else {
		set real_indices [string range $indices 0 [expr [string first "_" $indices] - 1]]
		set disp [string range $indices [expr [string first "_" $indices] + 1] end]
	    }
	    hide_trans [expr [string index $disp 0] - 4] [expr [string index $disp 1] - 4] [expr [string index $disp 2] - 4] $type $real_indices 
	    set list_of_actions [lreplace $list_of_actions end end]
	} elseif {[string first "sample_sa" $element_name] != -1} { 
	    set first_index [lindex $indices 0]
	    set last_index  [lindex $indices end]
	    for {set i $first_index} {$i <= $last_index} {incr i} {
		draw delete $i
	    }
	} elseif {[string first "sample_ray_atom" $element_name] != -1} { 
	    set first_index [lindex $indices 0]
	    set last_index  [lindex $indices end]
	    for {set i $first_index} {$i <= $last_index} {incr i} {
		draw delete $i
	    }
        } elseif {[string first "draw_ray_file" $element_name] != -1} {
            set first_index [lindex $indices 0]
            set last_index  [lindex $indices end]
            for {set i $first_index} {$i <= $last_index} {incr i} {
                draw delete $i
            }
	} elseif {[string first "sample_ray_sphere" $element_name] != -1} { 
	    set first_index [lindex $indices 0]
	    set last_index  [lindex $indices end]
	    for {set i $first_index} {$i <= $last_index} {incr i} {
		draw delete $i
	    }
	} elseif {[string first "sample_ray_andrew_sphere" $element_name] != -1} { 
	    set first_index [lindex $indices 0]
	    set last_index  [lindex $indices end]
	    for {set i $first_index} {$i <= $last_index} {incr i} {
		draw delete $i
	    }
	} elseif {[string first "sample_ray_andrew_atom" $element_name] != -1} { 
	    set first_index [lindex $indices 0]
	    set last_index  [lindex $indices end]
	    for {set i $first_index} {$i <= $last_index} {incr i} {
		draw delete $i
	    }
	} else {
	    puts "LAST ACTION $element_name HAS NOTHING TO BE UNDONE"
	}
    }
    return UNDONE
}

proc do_all {proc_name log_prefix type msg_name {extra_args {}}}  {
    global "num_$type\s" list_of_actions
    set num_elements [subst "\$num_$type\s"]
    set actionCount 0
    for {set i 0} {$i < $num_elements} {incr i} {
	set command [list $proc_name $type $i $extra_args]
	eval $command
	incr actionCount
    }

    if {[string compare $log_prefix "revert_"] != 0} {
	global list_of_actions
	if {[is_collapseable $proc_name]} {
	    if {[has_translation $proc_name]} {
		set vec [expr [lindex $extra_args 0] + 4][expr [lindex $extra_args 1] + 4][expr [lindex $extra_args 2] + 4]
		set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
		set actions [list $log_prefix$type "all_$vec"]
		lappend list_of_actions $actions
	    } else {
		set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
		lappend list_of_actions [list "$log_prefix$type" "all"]
	    }
	} else {
	    set new_actions [lrange $list_of_actions [expr [llength $list_of_actions] - $actionCount]  end]
	    set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
	    lappend list_of_actions [join $new_actions]
	}	
    }
    return "$msg_name"
}


proc do_list {proc_name lst log_prefix type msg_name {extra_args {}}}  {
    global "num_$type\s"
    set num_elements [llength $lst]
    set actionCount 0
    for {set i 0} {$i < $num_elements} {incr i} {
	set command [list $proc_name $type [lindex $lst $i] $extra_args]
	eval $command
	incr actionCount
    }
    if {[string compare $log_prefix "revert_"] != 0} {
	global list_of_actions
	if {[is_collapseable $proc_name]} {
	    if {[has_translation $proc_name]} {
		set vec [expr [lindex $extra_args 0] + 4][expr [lindex $extra_args 1] + 4][expr [lindex $extra_args 2] + 4]
		set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
		set new_indices $lst
		lappend new_indices $vec
		lappend list_of_actions [list $log_prefix$type $new_indices]
	    } else {
		set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
		lappend list_of_actions [list "$log_prefix$type" $lst]
	    }
	} else {
	    set new_actions [lrange $list_of_actions [expr [llength $list_of_actions] - $actionCount]  end]
	    set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $actionCount + 1]]]
	    lappend list_of_actions [join $new_actions]
	}
    }
    return "$msg_name"
}


proc helperInvoker {helperProc type indices log_message msg extra_args} {
    if {[is_list $indices]} {
	checkListIndices $indices $type
	do_list $helperProc $indices $log_message $type $msg $extra_args 
    } elseif {[string compare $indices "all"] == 0} {
	do_all $helperProc $log_message $type $msg $extra_args
    } else {
	checkElementIndex $type $indices
	$helperProc $type $indices $extra_args
    }
}

proc show {type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    helperInvoker "show_helper" $type $indices "draw_" "$type\s DRAWN" {0 0 0}
}

proc show_helper {type index extra_args} {
    global list_of_actions draw_count
    global "draw_$type\_log" "$type\s" "$type\Colors" nodeColors vorcellColors
    global atomRadii nodeRadii sphere_resolution
    
    set vec [expr [lindex $extra_args 0] + 4][expr [lindex $extra_args 1] + 4][expr [lindex $extra_args 2] + 4]
    if {[dict exists [subst "\$draw\_$type\_log($index)"] $vec]} {
	return -1
    }

    incr "draw_$type\_count($index)" 

    set suffix "s($index)"
    set element [subst "\$$type$suffix"]	
    set draw_id_list {}
    
    foreach line $element {
	set command [subst $line]
	set id [draw $command] 
	lappend draw_id_list $id
	set draw_count $id
    }
    
    lappend list_of_actions [list "draw_$type" "$index\_$vec"]
    dict set [subst "draw_$type\_log($index)"] $vec $draw_id_list
    return "$type #$index DRAWN"
}

proc hide_trans {da db dc type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    helperInvoker "hide_trans_helper" $type $indices "hide_" "$type\s HIDDEN" [list $da $db $dc]
}

proc hide_trans_helper {type index extra_args} {
    global "draw_$type\_log"
    global list_of_actions
    set vec [expr [lindex $extra_args 0] + 4][expr [lindex $extra_args 0] + 4][expr [lindex $extra_args 0] + 4]
    if {[dict exists [subst "\$draw\_$type\_log($index)"] $vec]} {
	set values [dict get [subst "\$draw\_$type\_log($index)"] $vec]
	set [subst "draw\_$type\_log($index)"] [dict remove [subst "\$draw\_$type\_log($index)"] $vec]
	
	for {set i 0} {$i < [llength $values]} {incr i} {
	    foreach drawID $values {draw delete $drawID}
	}
	lappend list_of_actions [list "hide_$type" "$index\_$vec"]
	return 1;
    } else {
	return 0;
    }
}

proc hide {type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    hide_trans 0 0 0 $type $indices
}


proc hide_all_images {type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    helperInvoker "hide_all_images_helper" $type $indices "hideall_" "$type\s HIDDEN" {}
}

proc hide_all_images_helper {type index extra_args} {
    global "draw_$type\_log"
    global list_of_actions
    set logname "draw_$type\_log($index)"
    
    if {[dict size [subst "\$draw\_$type\_log($index)"]] == 0} {
	return -1
    }
    
    set keys   [dict keys   [subst "\$draw\_$type\_log($index)"]]
    set values [dict values [subst "\$draw\_$type\_log($index)"]]
    
    set actions ""
    for {set i 0} {$i < [llength $keys]} {incr i} {
	foreach drawID [lindex $values $i] {draw delete $drawID}
	set key [lindex $keys $i]
	append actions "hide_$type $index\_$key "
    }

    set "draw_$type\_log($index)" [dict create]
    lappend list_of_actions $actions
    return "$type #$index HIDDEN"
}


proc color {newColor type {indices 0}} {
    global colorable_types available_colors
    checkElementTypes $type $colorable_types
    checkElementTypes $newColor $available_colors
    helperInvoker "color_helper" $type $indices "color_" "$type\s COLORED $newColor" [list $newColor]

    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc color_helper {type index extra_args} {
    global "$type\Colors" "old_$type\_colors"
    global list_of_actions
    set newColor [lindex $extra_args 0]
    lappend "old_$type\_colors($index)" [subst "\$$type\Colors($index)"] 
    lappend list_of_actions [list "color_$type" $index]
    set "$type\Colors($index)" $newColor
    update_env $type $index
    if {[string compare $type "unitcell"] == 0} {
	return "COLOR OF $type CHANGED TO $newColor"
    } else {
	return "COLOR OF $type #$index CHANGED TO $newColor" 
    }
}

proc adjust_size {radiusChange type {indices 0}} {
    global sizeable_types
    checkElementTypes $type $sizeable_types
    helperInvoker "adjust_size_helper" $type $indices "size_" "$type\s RADIUS CHANGED BY $radiusChange" [list $radiusChange]

    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc adjust_size_helper {type index extra_args} {
    global "$type\Radii" "old_$type\_radii"
    set radius_change [lindex $extra_args 0]
    set old_radius [subst "\$$type\Radii($index)"]
    set new_radius [expr $radius_change + $old_radius]
    size_helper $type $index [list $new_radius]
}

proc all_updater {type} {
    global "num_$type\s"
    set num_elements [subst "\$num_$type\s"]
    for {set i 0} {$i < $num_elements} {incr i} {
	update_env $type $i
    }
}

proc update_env {type index} {
    global list_of_actions
    global "draw_$type\_log"
    if {[expr [dict size [subst "\$draw_$type\_log($index)"] ] > 0]} {
	set num_actions [llength $list_of_actions]
	set keys [dict keys [subst "\$draw\_$type\_log($index)"]]
	hide_all_images $type $index
	foreach key $keys {show_trans [expr [string index $key 0] - 4] [expr [string index $key 1] - 4] [expr [string index $key 2] - 4] $type $index}
	set list_of_actions [lreplace $list_of_actions $num_actions end] 
    }
}

proc size {newRadius type {indices 0}} {
    global sizeable_types
    checkElementTypes $type $sizeable_types
    helperInvoker "size_helper" $type $indices "size_" "$type\s RADIUS SET TO $newRadius" [list $newRadius]

    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc size_helper {type index extra_args} {
    set newRadius [lindex $extra_args 0]
    global "$type\Radii" "old_$type\_radii"
    global list_of_actions
    lappend "old_$type\_radii($index)" [subst "\$$type\Radii($index)"]
    lappend list_of_actions [list "size_$type" $index]
    set "$type\Radii($index)" $newRadius
    update_env $type $index
    return "RADIUS OF $type #$index CHANGED TO $newRadius"
}

proc scale {factor type {indices 0}} {
    global sizeable_types
    checkElementTypes $type $sizeable_types
    helperInvoker "scale_helper" $type $indices "size_" "$type\s SCALED BY A FACTOR OF $factor" [list $factor]

    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc scale_helper {type index extra_args} {
    global "$type\Radii"
    set old_radius [subst "\$$type\Radii($index)"]
    set factor [lindex $extra_args 0]
    set new_radius [expr $factor * $old_radius]
    size_helper $type $index [list $new_radius]
}

proc reset {characteristic type {indices 0}} { 
    checkElementTypes $characteristic {color size}
    if {[string compare $characteristic "color"] == 0} {
	global colorable_types
	checkElementTypes $type $colorable_types
	set logname "color_"
    } elseif {[string compare $characteristic "size"] == 0} {
	global sizeable_types
	checkElementTypes $type $sizeable_types
	set logname "size_"
    }
    helperInvoker "reset_helper" $type $indices $logname "$characteristic OF $type\s RESET" [list $characteristic]
    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc reset_helper {type index extra_args} {
    set characteristic [lindex $extra_args 0]
    if {[string compare $characteristic "size"] == 0} {
	global "old_$type\_radii"
	if {[llength [subst "\$old_$type\_radii($index)"]] == 0} {
	    return "$characteristic OF $type #$index HAS NOT BEEN CHANGED FROM ITS ORIGINAL VALUE"
	}
	size_helper $type $index [list [lindex [subst "\$old_$type\_radii($index)"] 0]]
    } elseif {[string compare $characteristic "color"] == 0} {
	global "old_$type\_colors"
	if {[llength [subst "\$old_$type\_colors($index)"]] == 0} {
	    return "$characteristic OF $type #$index HAS NOT BEEN CHANGED FROM ITS ORIGINAL VALUE"
	}
	color_helper $type $index [list [lindex [subst "\$old_$type\_colors($index)"] 0]]
    }
    update_env $type $index
    return "$characteristic OF $type #$index RESET"
}

proc revert {characteristic type {indices 0}} {
    checkElementTypes $characteristic {color size}
    if {[string compare $characteristic "color"] == 0} {
	global colorable_types
	checkElementTypes $type $colorable_types
    } elseif {[string compare $characteristic "size"] == 0} {
	global sizeable_types
	checkElementTypes $type $sizeable_types
    }
    helperInvoker "revert_helper" $type $indices "revert_" "$characteristic OF $type\s RESET" [list $characteristic]
    if {[expr [string compare $type node] == 0]} {
	all_updater vorcell
	all_updater vornet
	all_updater channel
    }
}

proc revert_helper {type index extra_args} {
    set characteristic [lindex $extra_args 0]
    if {[string compare $characteristic "size"] == 0} {
	global "$type\Radii"
	global "old_$type\_radii"
	set suffix1 "_radii($index)"
	set suffix2 "Radii($index)"
    } elseif {[string compare $characteristic "color"] == 0} {
	global "$type\Colors"
	global "old_$type\_colors"
	set suffix1 "_colors($index)"
	set suffix2 "Colors($index)"
    }
    global list_of_actions
    if {[llength [subst "\$old_$type$suffix1"]] >= 1} {
	set "$type$suffix2" [lindex [subst "\$old_$type$suffix1"] end]
	set "old_$type$suffix1" [lreplace [subst "\$old_$type$suffix1"] end end]
	update_env $type $index
	return "$characteristic OF $type #$index REVERTED"
    } else {
	return "$characteristic OF $type #$index HAS NO PREVIOUS VALUE"
    }
}

proc color_code_nodes {probeRad} {
    global num_nodes
    global axs_node_color
    global inaxs_node_color 
    global inputFile
    run_network_program [list "-r" "-axs" $probeRad "ZeoVisInput.axs" $inputFile]
    set fp [open "ZeoVisInput.axs" r]
    set file_data [read $fp]
    set nodeID 0

    set new_list_actions {}
    foreach line $file_data {
	if {[string compare $line "true"] == 0} {
	    color $axs_node_color node $nodeID 
	} else {
	    color $inaxs_node_color node $nodeID
	}
	incr nodeID
    }

    remove_file "ZeoVisInput.axs"
    
    global list_of_actions
    set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - [expr $num_nodes + 1]]]
    lappend list_of_actions {color_node all}
    return "NODES COLORED ACCORDING TO ACCESSIBILITY. \nACCESSIBLE NODES ARE $axs_node_color WHILE INACCESSIBLE NODES ARE $inaxs_node_color "
}

proc sample_ray_atom {chanRad probeRad  numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zray_atom" $chanRad $probeRad $numSamples "ZeoVisInput.zray_atom" $inputFile]
    set fp [open "ZeoVisInput.zray_atom"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_ray_atom" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_atom"
}

proc draw_ray_file {} {
    global inputFile num_atoms list_of_actions draw_count
#    run_network_program [list "-r" "-zray_atom" $chanRad $probeRad $numSamples "ZeoVisInput.zray_atom" $inputFile]
    set fp [open "ZeoVisInput.zray_file"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
        draw $line
        set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
        lappend list_of_actions [list "sample_ray_atom" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_file"
}


proc sample_ray_sphere {chanRad probeRad numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zray_sphere" $chanRad $probeRad $numSamples "ZeoVisInput.zray_sphere" $inputFile]
    set fp [open "ZeoVisInput.zray_sphere"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_ray_sphere" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_sphere"
}

proc sample_ray_node {chanRad probeRad numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zray_node" $chanRad $probeRad $numSamples "ZeoVisInput.zray_node" $inputFile]
    set fp [open "ZeoVisInput.zray_node"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_ray_node" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_node"
}


proc sample_ray_andrew_sphere {chanRad probeRad  numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zray_andrew_sphere" $chanRad $probeRad $numSamples "ZeoVisInput.zray_andrew_sphere" $inputFile]
    set fp [open "ZeoVisInput.zray_andrew_sphere"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_ray_andrew_sphere" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_andrew_sphere"
}

proc sample_ray_andrew_atom {chanRad probeRad  numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zray_andrew_atom" $chanRad $probeRad $numSamples "ZeoVisInput.zray_andrew_atom" $inputFile]
    set fp [open "ZeoVisInput.zray_andrew_atom"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_ray_andrew_atom" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zray_andrew_atom"
}

proc sample_surface_area {probeRad numSamples} {
    global inputFile num_atoms list_of_actions draw_count
    run_network_program [list "-r" "-zsa" $probeRad $probeRad $numSamples "ZeoVisInput.zsa" $inputFile]
    set fp [open "ZeoVisInput.zsa"]
    set file_data [read $fp]
    set first_index $draw_count

    foreach line $file_data {
	draw $line
	set draw_count [expr $draw_count + 1]
    }
    set last_index [expr $draw_count - 1]

    if {$last_index >= $first_index} {
	lappend list_of_actions [list "sample_sa" [list $first_index "to" $last_index]]
    }
    remove_file "ZeoVisInput.zsa"
}

proc id_channels {probeRad} {
    global channels num_channels channelColors draw_channel_log old_channel_colors
    global inputFile probe_radius
    global list_of_actions
    set probe_radius $probeRad
    
    if {$num_channels != 0} {
	hide_all_images channel all
	set list_of_actions [lreplace $list_of_actions end end]
	unset channels
	unset draw_channel_log
	unset channelColors
	unset old_channel_colors
    }

    set list_of_actions [lsearch -regexp -all -inline -index 0 -not $list_of_actions ".*channel"]

    run_network_program [list "-r" "-zchan" $probeRad "ZeoVisInput.zchan" $inputFile]
    source "ZeoVisInput.zchan"
    
    remove_file "ZeoVisInput.zchan"
    if {$num_channels > 0} {
	puts "$num_channels CHANNELS IDENTIFIED INDEXED 0 THROUGH [expr $num_channels - 1]"
    } else {
	puts "NO CHANNELS WERE IDENTIFIED. THE PROBE SIZE IS TOO LARGE FOR TRAVERSAL. \nTRY A SMALLER PROBE SIZE"
    }
    
    global color_list
    global element_types
    set chanColor [lindex $color_list [lsearch $element_types channel]]

    for {set i 0} {$i < $num_channels} {incr i} {
	set channelColors($i)  $chanColor
	set old_channel_colors($i) {}
	set draw_channel_log($i) [dict create]
    }
}

proc helperInvoker {helperProc type indices log_message msg extra_args} {
    if {[is_list $indices]} {
	checkListIndices $indices $type
	do_list $helperProc $indices $log_message $type $msg $extra_args 
    } elseif {[string compare $indices "all"] == 0} {
	do_all $helperProc $log_message $type $msg $extra_args
    } else {
	checkElementIndex $type $indices
	$helperProc $type $indices $extra_args
    }
}


proc id_segments {} {
    global num_segments segmentColors old_segment_colors segments draw_segment_log
    global list_of_actions
     if {$num_segments != 0} {
	 hide_all_images segment all
	 set list_of_actions [lreplace $list_of_actions end end]
	 unset segments
	 unset segmentColors
	 unset old_segment_colors
	 unset  draw_segment_log
    }

    set list_of_actions [lsearch -regexp -all -inline -index 0 -not $list_of_actions ".*segment"]
    set num_segments 0
    global num_channels segments probe_radius inputFile segmentColors available_colors 
    run_network_program [list "-r" "-zseg" $probe_radius "ZeoVisInput.zseg" $inputFile]
    source "ZeoVisInput.zseg"
    remove_file "ZeoVisInput.zseg"

    source "seginfo.data"
    remove_file "seginfo.data"
   
    for {set i 0} {$i < $num_channels} {incr i} {
	set old_seg_index $num_segments
	set num_segments [expr $num_segments + [lindex $seg_counts $i]]
	set new_seg_index [expr $num_segments - 1]
	puts "CHANNEL #$i\'s SEGMENTS INDEXED FROM $old_seg_index TO $new_seg_index"

	for {set j $old_seg_index} {$j <= $new_seg_index} {incr j} {
	    set draw_segment_log($j) [dict create]
	    set segmentColors($j) [lindex $available_colors [expr $j % [llength $available_colors]]]
	    set old_segment_colors($j) {}
	}
    }
    return "ALL SEGMENTS IDENTIFIED"
}


proc id_features {} {
    global num_features featureColors old_feature_colors features draw_feature_log
    global list_of_actions
     if {$num_features != 0} {
	 hide_all_images feature all
	 set list_of_actions [lreplace $list_of_actions end end]
	 unset features
	 unset featureColors
	 unset old_feature_colors
	 unset draw_feature_log
    }

    set list_of_actions [lsearch -regexp -all -inline -index 0 -not $list_of_actions ".*feature"]

    set num_features 0
    global num_channels probe_radius inputFile available_colors
    run_network_program [list "-r" "-zfeat" $probe_radius "ZeoVisInput.zfeat" $inputFile]
    source "ZeoVisInput.zfeat"
    remove_file "ZeoVisInput.zfeat"

    source "featinfo.data"
    remove_file "featinfo.data"

    for {set i 0} {$i < $num_channels} {incr i} {
	set old_index $num_features
	set num_features [expr $num_features + [lindex $feat_counts $i]]
	set new_index [expr $num_features - 1]
	puts "CHANNEL #$i\'s FEATURES INDEXED FROM $old_index TO $new_index"


	for {set j $old_index} {$j <= $new_index} {incr j} {
	    set draw_feature_log($j) [dict create]
	    set featureColors($j) [lindex $available_colors [expr $j % [llength $available_colors]]]
	    set old_feature_colors($j) {}
	}
    }
    return "ALL FEATURES IDENTIFIED"
}

proc id_cages {probeRad} {
    global num_cages cages inputFile 
    global old_cage_colors cageColors available_colors draw_cage_log
    global list_of_actions
    if {$num_cages != 0} {
	hide_all_images cage all
	set list_of_actions [lreplace $list_of_actions end end]
	unset cages
	unset draw_cage_log
	unset cageColors
	unset old_cage_colors
    }

    set list_of_actions [lsearch -regexp -all -inline -index 0 -not $list_of_actions ".*cage"]
    set num_cages 0
    run_network_program [list "-r" "-zcage" $probeRad "ZeoVisInput.zcage" $inputFile]
    source "ZeoVisInput.zcage"
    remove_file "ZeoVisInput.zcage"

    for {set i 0} {$i < $num_cages} {incr i} {
	set draw_cage_log($i) [dict create]
	set cageColors($i) [lindex $available_colors [expr $i % [llength $available_colors]]]
	set old_cage_colors($i) {}
    }

    if {$num_cages == 0} {
	return "NO CAGES WERE IDENTIFIED"
    } else {
	set max_index [expr $num_cages - 1]
	return "CAGES INDENTIFIED AND INDEXED FROM 0 TO $max_index"
    }
}

proc load_new_file {fileName incRadius} {
    global num_vorcells num_faces num_vornets num_nodes num_atoms num_unitcells num_channels uc_a_vector uc_b_vector uc_c_vector 
    global sphere_resolution atoms atomRadii nodes nodeRadii vornets vorcells faces num_faces num_channels num_features num_segments num_cages 
    clear_env
    run_network_program [list "-r" "-zvor" $incRadius "ZeoVisInput.zor" $fileName]
    source "ZeoVisInput.zvor"
    remove_file "ZeoVisInput.zvor"
}

proc help {{procName "NONE"}} {
    set types {}
    set characteristics {}
    set special ""
    set indexMsg "VALUES FOR indices: a number, a list such as {1 5 7} or all"
    global sizeable_types colorable_types element_types

    if {[string compare $procName "NONE"] == 0} {
	puts "Please specify about which function you'd like more information:"
	list_commands
	return
    } elseif {[string compare $procName "show"] == 0} {
	set types $element_types
	set description "DISPLAYS THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES IN THE (0,0,0) UNITCELL"
    } elseif {[string compare $procName "show_trans"] == 0} {
	set types $element_types
	set description "DISPLAYS THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES, SHIFTED IN UNIT CELL COORDINATES BY (da, db, dc)"
    } elseif {[string compare $procName "show_range"] == 0} {
	set types $element_types
	set description "DISPLAYS THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES OVER THE RANGE OF UNIT CELLS SPECIFIED.  "
	set special "VALUES FOR minCoord AND maxCoord: a 3-element list such as {1 0 1}"
    } elseif {[string compare $procName "hide"] == 0} {
	set types $element_types
	set description "HIDES THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES IN THE (0,0,0) UNITCELL"
    } elseif {[string compare $procName "hide_trans"] == 0} {
	set types $element_types
	set description "HIDES THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES IN THE UNITCELL (da, db, dc)"
    } elseif {[string compare $procName "hide_all_images"] == 0} {
	set types $element_types
	set description "HIDES THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES IN ALL UNITCELLS"
    } elseif {[string compare $procName "color"] == 0} {
	set types $colorable_types
	set description "CHANGES THE COLOR OF THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES"
    } elseif {[string compare $procName "size"] == 0} {
	set types $sizeable_types
	set description "CHANGES THE RADIUS OF THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES TO THE PROVIDED VALUE"
    } elseif {[string compare $procName "adjust_size"] == 0} {
	set types $sizeable_types
	set description "CHANGES THE RADIUS OF THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES BY THE PROVIDED AMOUNT"
    } elseif {[string compare $procName "scale"] == 0} {
	set types $sizeable_types
	set description "SCALES THE RADIUS OF THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES BY THE PROVIDED FACTOR"
    } elseif {[string compare $procName "reset"] == 0} {
	set characteristics {color size}
	set special "TYPES ARE SPECIFIED BY THE CORRESPONDING TYPES IN THE COLOR/SIZE COMMANDS"
	set description "RESETS THE SIZE/COLOR OF THE ELEMENT(s) REFERRED TO BY THE TYPE AND INDEX/INDICES TO THEIR ORIGINAL VALUES"
    } elseif {[string compare $procName "undo"] == 0} {
	set description "UNDOES THE LAST ACTION IF POSSIBLE"
    } elseif {[string compare $procName "clear_env"] == 0} {
	set description "REMOVES ALL DISPLAYED ELEMENTS. RESETS SIZE AND COLORS OF ALL ELEMENTS TO THEIR ORIGINAL VALUES. CAN NOT BE UNDONE"
    }  elseif {[string compare $procName "color_code_nodes"] == 0} {
	global axs_node_color
	global inaxs_node_color 
	set description "COLORS THE NODES ACCORDING TO ACCESSIBILITY FOR THE PROBE SIZE. \nACCESSIBLE AND INACCESSIBLE NODES ARE $axs_node_color AND $inaxs_node_color RESPECTIVELY."
    } elseif {[string compare $procName "sample_surface_area"] == 0} {
	set description "SAMPLES THE ACCESSIBLE SURFACE AREA FOR THE PROVIDED PROBE.\n POINTS WITHIN ATOMS ARE NOT SHOWN. ACCESSIBLE AND INACESSIBLE POINTS ARE GREEN AND RED RESPECTIVELY."
	puts "**CAUTION: PROGRAM MAY HALT IF PROBE SIZE IS TOO SMALL OR NUMBER OF SAMPLES IS TOO LARGE.**"
    } elseif {[string compare $procName "id_channels"] == 0} {
	set description "IDENTIFIES THE CHANNELS ACCESSIBLE FOR PROBE OF RADIUS probeRad, INDEXES THEM FROM 0 TO N-1,\nAND ADDS THEM TO THE ENVIRONMENT. DESTROYS ANY PREVIOUSLY IDENTIFIED CHANNELS"
    } elseif {[string compare $procName "id_segments"] == 0}  {
	set description "IDENTIFIES THE SEGMENTS PRESENT IN THE MOST RECENTLY IDENTIFIED CHANNELS AND DESTROYS ANY PREVIOUSLY IDENTIFIED SEGMENTS"
    } elseif {[string compare $procName "id_features"] == 0}  {
	set description "IDENTIFIES THE FEATURES PRESENT IN THE MOST RECENTLY IDENTIFIED CHANNELS AND DESTROYS AND PREVIOUSLY IDENTIFIED FEATURES"
    } elseif {[string compare $procName "sample_ray_atom"] == 0}  {
	set description "Samples the accessible volume of a void netork randomly shooting lines within the cell until it hits an atom or maximum ray distance"
    } elseif {[string compare $procName "draw_ray_file"] == 0}  {
        set description "Draws rays from a precomputed file (ZeoVisInput.zray_file) "
    } elseif {[string compare $procName "sample_ray_sphere"] == 0}  {
	set description "Sample a random point within vornode radius's and shoot the line until it hits a point outside of vorcell"
    } elseif {[string compare $procName "sample_ray_andrew_sphere"] == 0}  {
	set description "Samples a random point within a unit cell and shoots in random orientation until it reaches max ray distance. Only segments that go through vornodes are recorded"
    } elseif {[string compare $procName "sample_ray_andrew_atom"] == 0}  {
	set description "Samples a random point within a unit cell and shoots in random orientation until it reaches max ray distance. Only segments that go through accessible volume are recorded"
    } elseif {[string compare $procName "help"] == 0} {
	puts "YOU MUST FEEL REALLY CLEVER"
	after 2000
	puts "BUT YOU'RE NOT"
	return
    } else {
	error "ERROR: COMMAND $procName DOES NOT EXIT"
	return
    }

    set numArgs [llength [info args $procName]]
    puts "FUNCTION NAME: $procName"
    puts "DESCRIPTION: $description"
    puts "# OF ARGUMENTS: $numArgs"
    if {$numArgs != 0} {
	puts "ARGUMENTS: [info args $procName]"
    }
    if {[llength $types] != 0} {
	puts "VALUES FOR type: $types"
    }
    if {[llength $characteristics] != 0} {
	puts "VALUES FOR characteristic: $characteristics"
    }
    if {[lsearch [info args $procName] indices] != -1} {
	puts $indexMsg
    }
    if {[string compare $special ""] != 0} {
	puts $special
    }
}

proc max {x y} {
    if {[expr $x < $y] == 1} {
	return $y
    } else {
	return $x
    }
}

proc min {x y} {
    if {[expr $x > $y] == 1} {
	return $y
    } else {
	return $x
    }
}


proc show_range {minCoord maxCoord type {indices 0}} {
    global element_types
    checkElementTypes $type $element_types
    if { [expr [is_list $minCoord] + [is_list $maxCoord] + [expr [llength $minCoord] == 3] + [expr [llength $maxCoord] == 3] ] == 4} {
	show_range_helper $type $indices [list $minCoord $maxCoord]
	return "RANGE OF $type DRAWN"
    } else {
	error "ERROR: minCoord AND maxCoord MUST BE 3-ELEMENT LISTS"
    }
}

proc show_range_helper {type index extra_args} {
    set coord1 [lindex $extra_args 0]; set coord2 [lindex $extra_args 1];
    set minX [lindex $coord1 0];  set minY [lindex $coord1 1];  set minZ [lindex $coord1 2];
    set maxX [lindex $coord2 0];  set maxY [lindex $coord2 1];  set maxZ [lindex $coord2 2];
    
    set minX [max $minX -4]
    set minY [max $minY -4]
    set minZ [max $minZ -4]

    set maxX [min $maxX 4]
    set maxY [min $maxY 4]
    set maxZ [min $maxZ 4]

    set showCount 0
    for {set i $minX} {$i <= $maxX} {incr i} {
	for {set j $minY} {$j <= $maxY} {incr j} {
	    for {set k $minZ} {$k <= $maxZ} {incr k} {
		show_trans $i $j $k $type $index
		incr showCount
	    }
	}
    }
    global list_of_actions
    set new_list_actions [join [lrange $list_of_actions  [expr [llength $list_of_actions] - $showCount] end]]
    set list_of_actions [lrange $list_of_actions 0 [expr [llength $list_of_actions] - $showCount - 1]]
    lappend list_of_actions $new_list_actions
}

proc list_commands {} {
    set commands {show show_trans show_range hide hide_trans hide_all_images color size adjust_size scale reset undo clear_env color_code_nodes sample_surface_area sample_ray_atom draw_ray_file sample_ray_sphere sample_ray_andrew_atom sample_ray_andrew_sphere id_channels id_segments id_features}
    foreach command $commands {puts $command}
    return "COMMANDS SHOWN"
}

proc is_collapseable {procName} {
    set collapseable_commands {show show_trans show_range hide hide_trans color size adjust_size scale reset color_code_nodes}
    set index [string last "_" $procName]
    if {[expr $index == -1]} {
	set origName $procName
    } else {
	set origName [string range $procName 0 [expr $index - 1]]
    }
    return [expr [lsearch -exact $collapseable_commands $origName] != -1];
}

proc has_translation {procName} {
    set index [string first "_" $procName]
    if {[expr $index == -1]} {
	set origName $procName
    } else {
	set origName [string range $procName 0 [expr $index - 1]]
    }
    return [expr [lsearch -exact {show hide} $origName] != -1];
}
