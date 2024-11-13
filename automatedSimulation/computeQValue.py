import subprocess
import pandas as pd
import math
import sys
import time
import yaml



def computeQValue(frame_arg, config):
    
    print("\n####################### Q-VALUE #######################")
    print(f"Frame at which we are computing Q-value: {frame_arg}")

    print("Computing Q-value...")
    f = open("log/Qvalue.log", "w")


    command = [
        "vmd", "proteins/ionized.psf", "simulation_output/result.dcd",
        "-dispdev", "text",
        "-e", "scripts/qvalue.tcl",
        "-args", str(frame_arg)
    ] 


    try:
        # Run the VMD command
        subprocess.run(command, check=True,  stdout=f)
        
        fin = open('log/native_contacts.txt', 'r')


        # Reading results from tcl process: qvalue.tcl        
        for line in fin:
            qvalue = line.split(',')[0]
            fraction = line.split(',')[1]
        fin.close()

        # Writing those result in a csv file along with other information
        fout = open('results/qvalues.csv', 'a')
        fout.write(
            config["pdb"]+','  +
            str(config["namd"]["temperature"]) + ',' +
            str(int(config["namd"]["run"])*2/1e6)  + ','  +
            qvalue + ',' +
            fraction + '\n'
            )
        fout.close()


        
        print(f"Computed a Q-value    of : {qvalue}")
        print(f"Computed a Q-value(%) of : {fraction}")
        print("Log at: log/Qvalue.log and log/native_contacts.txt")
        print("Result written at: results/qvalues.csv")

        # return qvalue # if we want to use it again...
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running VMD: {e}")
   
def compute_box_size():
    df = pd.read_csv('log/box_measure.csv')
    x =  math.ceil(df["max_x"][0] - df["min_x"][0])
    y = math.ceil(df["max_y"][0] - df["min_y"][0])
    z = math.ceil(df["max_z"][0] - df["min_z"][0])
    return x, y, z

def compute_center():
    df = pd.read_csv('log/box_measure.csv')
    x = df['center_x'][0]
    y = df['center_y'][0]
    z = df['center_z'][0]
    return x, y, z

def write_conf(config):

    x, y, z = compute_box_size()

    xx, yy, zz = compute_center()


    conf = f"""
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# Minimization and Equilibration of 
# Ubiquitin in a Water Box


#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          proteins/ionized.psf
coordinates        proteins/ionized.pdb

set temperature    {config["namd"]["temperature"]}
set outputname     simulation_output/result

firsttimestep      0


#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm	    on
parameters          par_all27_prot_lipid.inp 
temperature         {config["namd"]["temperature"]}


# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.0
switching           on
switchdist          10.0
pairlistdist        14.0


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2  
stepspercycle       10


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     1     ;# damping coefficient (gamma) of 1/ps
langevinTemp        $temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens


# Periodic Boundary Conditions
cellBasisVector1     {x}    0.   0.0
cellBasisVector2     0.0  {y}   0.0
cellBasisVector3     0.0    0   {z} 
cellOrigin           {xx}  {yy}  {zz} ;#this must be the previously computed center


wrapAll             off


# PME (for full-system periodic electrostatics)
PME                 yes
PMEGridSpacing      1.0

#manual grid definition
#PMEGridSizeX        45
#PMEGridSizeY        45
#PMEGridSizeZ        48


# Constant Pressure Control (variable volume)
useGroupPressure      yes ;# needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  100.0
langevinPistonDecay   50.0
langevinPistonTemp    $temperature


# Output
outputName          $outputname

;# restartfreq         5000     ;# 5000steps = every 10ps
dcdfreq                {config["namd"]["dcdfreq"]}
;# xstFreq             250 
outputEnergies      100
outputPressure      100


#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################


#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

# Minimization
minimize            500
reinitvels          $temperature

run {config["namd"]["run"]} ;# 1000ps = 1 ns

"""

    
    with open("autoSimulation.conf","w") as f:
        f.writelines(conf)

def prepare_protein(protein_pdb: str):
    

    f = open("log/prepare_protein.log", "w")

    print("####################### Protein cleaning #######################")
    print(f"Input pdb: {protein_pdb}")
    print("Log at: log/prepare_protein.log")
    try:
        subprocess.run(
            [
                "vmd",
                "-dispdev", "text",
                "-e", "scripts/prepare_protein.tcl",
                "-args", "proteins/" + protein_pdb + ".pdb"
            ], 
            stdout=f
        )
        print("Protein cleaned successfully!")
      
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while cleaning : {e}")

    f.close()

def run(config):
    f = open("log/run.log", "w")
    print("\n####################### RUNNING MDS #######################")
    print("Configuration at: configuration.yaml")
    print("Log at: log/run.log\n")
    print("Yaml configuration for current run: ")
    print(yaml.dump(config["namd"], allow_unicode=True, default_flow_style=False))
    print("Simulating...")

    
    start_time = time.time() 

    try:
        subprocess.run(
            [
                "namd3", "+p8", 
                "autoSimulation.conf",
            ], 
            stdout=f
        )

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the simulation: {e}")
    
    end_time = time.time()

    print(f"Simulation completed in: {int((start_time - end_time)/60)} mins")

def compute_rmsf(config, plot: bool):
    
    f = open("log/rmsf.log", "w")

    print("####################### RMSF Computing #######################")

    try:
        subprocess.run(
            [
                "vmd",
                "-dispdev", "text",
                "-e", "scripts/rmsf.tcl",
                # we may want to specify which pdb we want to read... useful for past runs...
            ], 
            stdout=f
        )
        print("RMSF computed successfully! ")
        print("You can find the result at: result/rmsf.csv")

        # we may want to map this file to a csv: but we just need to modify tcl file!

        if plot:
            return "Feature not yet implemented"

      
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while cleaning : {e}")

def write_backup():
    print("Not yet implemented")


with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)



run_time = config["namd"]["run"]
dcdfreq = config["namd"]["dcdfreq"]
QValue_frame = run_time/dcdfreq # we need to check this

#protein_pdb = config["pdb"]

#write_conf(config)

#prepare_protein(protein_pdb)

#run(config)

#computeQValue(QValue_frame, config)

compute_rmsf(config, True)
# We miss  compute rmsf






