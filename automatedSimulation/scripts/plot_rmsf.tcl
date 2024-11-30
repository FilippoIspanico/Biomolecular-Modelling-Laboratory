set rmsf_file "results/rmsf.csv"

set sel [atomselect top "protein"]

# Open and read RMSF data
set data [open $rmsf_file r]
set rmsf_values [read $data]
close $data

# Split the RMSF data by lines
set rmsf_lines [split $rmsf_values "\n"]

# Skip the first line (header)
set rmsf_lines [lrange $rmsf_lines 1 end]

# Assign RMSF values to the beta column
foreach line $rmsf_lines {
    if {[string length $line] > 0} {
        set fields [split $line ","]
        set resid [lindex $fields 0]
        set rmsf [lindex $fields 1]

        # Select residues by resid and assign RMSF value to beta column
        set resid_sel [atomselect $molid "resid $resid"]
        $resid_sel set beta $rmsf
        $resid_sel delete
    }
}

