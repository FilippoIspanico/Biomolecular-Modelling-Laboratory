import subprocess
import pandas as pd
import math
import sys
import time
import yaml
import seaborn as sns
import matplotlib.pylab as plt


class MDRunner:

    def __init__(self, yaml_configuration_path: str = "config.yaml"):
        
        
        try: 
            config_file = open(yaml_configuration_path, "r")
            self.yaml_configuration = yaml.safe_load(config_file)
        
        except:
            print(f'Cannot open input file: {yaml_configuration_path}')


        run_time = self.yaml_configuration["namd"]["run"]
        dcdfreq = self.yaml_configuration["namd"]["dcdfreq"]
        self.QValue_frame = run_time/dcdfreq # we need to check this

    def prepare_protein(self):
    
        f = open("log/prepare_protein.log", "w")

        protein = self.yaml_configuration["pdb"]

        print("####################### Protein cleaning #######################")
        print(f"Input pdb: { protein }")
        print("Log at: log/prepare_protein.log")
        try:
            subprocess.run(
                [
                    "vmd",
                    "-dispdev", "text",
                    "-e", "scripts/prepare_protein.tcl",
                    "-args", "proteins/" + protein + ".pdb"
                ], 
                stdout=f
            )
            print("Protein cleaned successfully!")
        
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while cleaning : {e}")

        f.close()

    def compute_box_size(self):
        df = pd.read_csv('log/box_measure.csv')
        x =  math.ceil(df["max_x"][0] - df["min_x"][0])
        y = math.ceil(df["max_y"][0] - df["min_y"][0])
        z = math.ceil(df["max_z"][0] - df["min_z"][0])
        return x, y, z

    def compute_center(self):
        df = pd.read_csv('log/box_measure.csv')
        x = df['center_x'][0]
        y = df['center_y'][0]
        z = df['center_z'][0]
        return x, y, z

    def write_conf(self):

        x, y, z = self.compute_box_size()

        xx, yy, zz = self.compute_center()


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

    set temperature    {self.yaml_configuration["namd"]["temperature"]}
    set outputname     simulation_output/result

    firsttimestep      0


    #############################################################
    ## SIMULATION PARAMETERS                                   ##
    #############################################################

    # Input
    paraTypeCharmm	    on
    parameters          par_all27_prot_lipid.inp 
    temperature         {self.yaml_configuration["namd"]["temperature"]}


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
    dcdfreq                {self.yaml_configuration["namd"]["dcdfreq"]}
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

    run {self.yaml_configuration["namd"]["run"]} ;# 1000ps = 1 ns

    """

        
        with open("autoSimulation.conf","w") as f:
            f.writelines(conf)

    def run(self):

        f = open("log/run.log", "w")
        print("\n####################### RUNNING MDS #######################")
        print("Configuration at: configuration.yaml")
        print("Log at: log/run.log\n")
        print("Yaml configuration for current run: ")
        print(yaml.dump(self.yaml_configuration["namd"], allow_unicode=True, default_flow_style=False))
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


        if self.yaml_configuration["backup"]:
            self.write_backup()

    def computeQValue(self):
        

        frame_arg =  self.QValue_frame


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
                self.yaml_configuration["pdb"]+','  +
                str(self.yaml_configuration["namd"]["temperature"]) + ',' +
                str(int(self.yaml_configuration["namd"]["run"])*2/1e6)  + ','  +
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

    def compute_rmsf(self, plot: bool = True, pastRun : str = None):
        
        """
        This function computes the RMSF of the current run or of a past run. It can also visualize such result in a plot comparing it with the one of the original protein
        """



        if pastRun is not None:
            
            psf_file = "pastRuns/" + pastRun + ".psf"
            dcd_file = "pastRuns/" + pastRun + ".dcd"

            print(f"Computing rmsf of a past run: {pastRun}")
        
        else:
            
            psf_file = "proteins/ionized.psf"
            dcd_file = "simulation_output/result.dcd"

            print(f'Computing rmsf of the current run.')



        f = open("log/rmsf.log", "w")

        print("\n####################### RMSF Computing #######################")

        try:
            subprocess.run(
                [
                    "vmd",
                    "-dispdev", "text",
                    "-e", "scripts/rmsf.tcl",
                    "-args",  psf_file, dcd_file
                ], 
                stdout=f
            )
            print("RMSF computed successfully! ")
            print("You can find the result at: result/rmsf.csv")

            # we may want to map this file to a csv: but we just need to modify tcl file!

            if plot:
                
                df = pd.read_csv('results/rmsf.csv')
                df_reference = pd.read_csv(self.yaml_configuration['reference_rmsf'])

                run_name = (self.yaml_configuration["pdb"] + '_' 
                    +   str(int(self.yaml_configuration["namd"]["run"])/1e6*2) + 'ns_at_' 
                    +   str(self.yaml_configuration["namd"]["temperature"]) + 'K'
                    )

                current_rmsf_name = pastRun if pastRun is not None else run_name

                plt.plot( df_reference['idx'],  df_reference['rmsf'])
                plt.plot( df['idx'], df['rmsf'])

                plt.title(
                    f'rmsf of {current_rmsf_name}'
                    )

                plt.legend(labels=['reference rmsf', f'{current_rmsf_name}'])
                plt.show()




      
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while cleaning : {e}")

    def write_backup(self):
        

        run_name = (self.yaml_configuration["pdb"] + '_' 
        +   str(int(self.yaml_configuration["namd"]["run"])/1e6*2) + 'ns_at_' 
        +   str(self.yaml_configuration["namd"]["temperature"]) + 'K'
        )
        run_name = self.yaml_configuration["backup_folder"] + run_name
        print(f"Writing backup as: {run_name}")



        try:
            subprocess.run(
                [
                    "cp",
                    "proteins/ionized.psf",
                    run_name + '.psf'
                    
                ]
            )

            subprocess.run(
                [
                    "cp",
                    "simulation_output/result.dcd",
                    run_name + '.dcd'
                    
                ]
            )
            print("Protein cleaned successfully!")
        
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while cleaning : {e}")




runner = MDRunner()
# runner.prepare_protein()
# runner.write_conf()
# runner.run()
# runner.computeQValue()
runner.compute_rmsf()
