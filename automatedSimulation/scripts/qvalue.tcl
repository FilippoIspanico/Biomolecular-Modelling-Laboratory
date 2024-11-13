source scripts/VMDextensions.tcl

## The atom selection
set protein [atomselect top "protein"]


set end_frame [lindex $argv 0]

## Use the first trajectory frame as a reference

$protein frame 0
set ref [ prepareNativeContacts 7 $protein ]

set nnc [ llength $ref ]
puts "There are $nnc contacts in the native state"




$protein frame $end_frame
set nc [measureNativeContacts $ref 7 $protein] 
puts "There are $nc contacts in the last frame"

set qnc [ expr 100.0 * $nc / $nnc ]


## Write the end contacts to a file
set output_file "log/native_contacts.txt"
set outfile [open $output_file w]
puts $outfile "$nc, $qnc"  ;# Write the number of contacts at the top

close $outfile
exit
