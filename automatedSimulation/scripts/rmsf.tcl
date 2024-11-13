set psf_file [lindex $argv 0]
set dcd_file [lindex $argv 1]

puts "input psf_file: $psf_file"
puts "input Dcd_file: $dcd_file"

mol new $psf_file
mol addfile $dcd_file




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