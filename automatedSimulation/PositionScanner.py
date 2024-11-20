import yaml
import pandas as pd
import re

class PositionScanner:
        
    def __init__(self, yaml_configuration_path: str = "config.yaml"):
    
    
        try: 
            config_file = open(yaml_configuration_path, "r")
            self.yaml_configuration = yaml.safe_load(config_file)
            self.pdb: str = self.yaml_configuration["pdb"]

            print("Inizialing Position Scanner: ")
            print(f"Current pdb: {self.pdb}")

        except:
            print(f'Cannot open input file: {yaml_configuration_path}')

    

    class Mutation:
        
        """
        This class represents a mutation: an aminoacid in a specific position.
        """

        def __init__(self, position: int, amino_acid: str):
            
            if len(amino_acid) != 1:
                
                raise ValueError('Not valid amino_acid input: please use one letter code for specifing the aminoacid')


            self.position = position
            
            self.amino_acid = amino_acid
        
        def __init__(self, foldX_name: str):
            self.position, self.amino_acid  = self.toSequence(foldX_name)

        def toFoldX(self):
            """
            Return the mutation in FoldX format i.e. amino_acid + "A" + position
            """
            return self.amino_acid + "A" + str(self.position)

        def __repr__(self):
            return f"({self.position}, {self.amino_acid})"

        def toSequence(self, foldX_name):
            amino_acid = foldX_name[-1]
            position = re.findall(r'\d+', foldX_name)
            
            return position[0], amino_acid



    def foldX_scan(self, a: int, b: int):
        """
        This function uses foldex PositionScan to scan the position in [a, b] and returns the best one. 
        FoldX only work if the protein to be mutated resides in the same folder of such program. 
        For this reason in this function we will first copy the current pdb (specified in the yaml file) into such folder, 
        and then run PositionScan. We will also need to modify the configuration file. Due to the large number of files produced by FoldX, 
        each time the program is called (which is b-a) the produced files will be removed. Once the best position has been identified 
        """

        



    def parse_FXPositionScan_output(self, epsilon: float):
        """
        This method parses PositionsCan output file. It returns mutations that:
        1. have a negative energy gain. In particular energy_gain<=epsilon 
        2. have a negative energy gain but epislon<energy_gain<0 and the proposed aminoacid in that position is the most common one (for that position). 
        """



        scanning_output_name = "/PS_" + self.pdb + "_scanning_output.txt"         
        scanning_output_name = self.yaml_configuration["foldX_installation_folder"] + scanning_output_name


        # scanning_output = open(scanning_output_name)

        scanning_output = pd.read_csv(scanning_output_name, header=None, sep = "\t")
        scanning_output.columns = ["Mutation", "EnergyGain"]
        
     

        mask = scanning_output["EnergyGain"] < epsilon 
        



        for idx, foldXmutation in enumerate(scanning_output["Mutation"]):
            print(idx, foldXmutation)
            if foldXmutation[-1] == self.yaml_configuration["most_common_sequence"][idx]: 
                mask[idx] = True


        result = []
        for mutation in scanning_output[mask]["Mutation"]:
            result.append(self.Mutation(foldX_name=mutation))
        
        print(result) 

scanner = PositionScanner()
scanner.parse_FXPositionScan_output(-1)