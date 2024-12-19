from src.MDRunner import MDRunner

"""
In this example program we run the simulation on two different proteins, compute the Q-Value and plot it
"""

proteins = ["2lyz", "GLY72"]  # add here any protein you want (as long as you have its pdb under the protein/ folder

runner = MDRunner()  # we initialize the runner: pass a config.yaml file in the constructor with your settings!
pastRuns = []  # we initialize an array that will contain the name of the simulation we are running, it is useful for
# comparing them later

for protein in proteins:
    runner.set_protein(protein)  # We set the protein
    runner.prepare_protein()  # We prepare it: removing water, adding solvent, runs the autoionize program, compute size of the box
    runner.write_conf()  # We write the NAMD configuration file
    runner.run()  # We launch the simulation
    runner.computeQValue()  # We compute Q-Value
    pastRuns.append(runner.get_run_name())  # We save the name of the run

runner.compareQValue(pastRuns)  # We plot the Q value for each of the past Runs
