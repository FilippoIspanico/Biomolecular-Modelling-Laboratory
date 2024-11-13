mol new proteins/ionized.psf 
mol addfile simulation_output/result.dcd

set outfile [open results/rmsf.csv w]
set ca [atomselect top "name CA"]
set numca [$ca num]
set rmsf [measure rmsf $ca]

puts $outfile "idx,rmsf"

for {set i 0} {$i < $numca} {incr i} {
    puts $outfile "[expr {$i+1}],[lindex $rmsf $i]"
}
close $outfile
exit