package require psfgen
package require solvate
package require autoionize

set input_protein_name [lindex $argv 0]

puts "input pdb: $input_protein_name" 


# Load the input protein pdb 
mol new $input_protein_name
puts "loaded input protein successfully"



######## Removing Water ###########

set a [atomselect top protein]
$a writepdb proteins/protein.pdb


# Close the input protein
mol delete top


######## Generating PSF ###########

mol new proteins/protein.pdb
puts "new protein loaded"



topology top_all27_prot_lipid.inp

pdbalias residue HIS HSE ;# maybe we can remove this?
pdbalias atom ILE CD1 CD

segment U {pdb proteins/protein.pdb} 
coordpdb proteins/protein.pdb U
guesscoord ;# we guess the coordines of missing atoms (mainly H, but also some >
writepdb proteins/system.pdb
writepsf proteins/system.psf

mol delete top

######## Solvating ###########

puts "Solvating the protein"
mol new proteins/system.psf 
mol addfile proteins/system.pdb
solvate proteins/system.psf proteins/system.pdb -t 5 -o proteins/solvated
mol delete top
puts "Solvating completed" 



######## Autoionizing ###########

autoionize -psf proteins/solvated.psf -pdb proteins/solvated.pdb -neutralize -o proteins/ionized 



######## Writing Box Measures ###########
mol new proteins/ionized.psf 
mol addfile proteins/ionized.pdb


# Get the selection of all atoms
set all [atomselect top all]

# Get the minimum and maximum coordinates
set minmax [measure minmax $all]
set min_x [lindex $minmax 0 0]
set min_y [lindex $minmax 0 1]
set min_z [lindex $minmax 0 2]
set max_x [lindex $minmax 1 0]
set max_y [lindex $minmax 1 1]
set max_z [lindex $minmax 1 2]

# Get the center coordinates
set center [measure center $all]
set center_x [lindex $center 0]
set center_y [lindex $center 1]
set center_z [lindex $center 2]

# Write the output to a CSV file
set output_file "log/box_measure.csv"
set fp [open $output_file "w"]
puts $fp "min_x,min_y,min_z,max_x,max_y,max_z,center_x,center_y,center_z"
puts $fp "$min_x,$min_y,$min_z,$max_x,$max_y,$max_z,$center_x,$center_y,$center_z"
close $fp

puts "Output written to $output_file"

exit 
