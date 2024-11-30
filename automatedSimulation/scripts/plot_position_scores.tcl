set position_file "results/positions_scores.csv"

set sel [atomselect top "protein"]

# Open and read position data
set data [open $position_file r]
set position_values [read $data]
close $data

# Split the position data by lines
set position_lines [split $position_values "\n"]

# Skip the first line (header)
set position_lines [lrange $position_lines 1 end]

# Assign scores values to the beta column
foreach line $position_lines {
    if {[string length $line] > 0} {
        set fields [split $line ","]
        set resid [lindex $fields 0]
        set position [lindex $fields 1]

        # Select residues by resid and assign position score value to beta column
        set resid_sel [atomselect $molid "resid $resid"]
        $resid_sel set beta $position
        $resid_sel delete
    }
}

