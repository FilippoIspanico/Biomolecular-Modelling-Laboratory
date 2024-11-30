source scripts/VMDextensions.tcl

## The atom selection
set protein [atomselect top "protein and name CA"]


set end_frame [lindex $argv 0]

## Use the first trajectory frame as a reference

$protein frame 0
set ref [ prepareNativeContacts 7 $protein ]

set nnc [ llength $ref ]
puts "There are $nnc contacts in the native state"




$protein frame $end_frame
set nc [measureNativeContacts $ref 7 $protein] 
puts "There are $nc contacts in the last frame"

# set qnc [ expr 100.0 * $nc / $nnc ]

set outfile [open "log/native_contact.csv" w]

puts $outfile "frame,nc,qnc"

              # Now, for each frame,
forFrames fno $protein {
               # compute number of native contacts,
         set nc [measureNativeContacts $ref 7 $protein]
               # their fraction,
         set qnc [ expr 100.0 * $nc / $nnc ]
               # and print both.
         puts $outfile "$fno, $nc, $qnc"
 }


close $outfile
exit
