from auxiliar_functions import coordinate_extractor, ensembl_request, get_sequence


class SearchInfo: 
    """This class stores both, the cleaned input and the newly acquired information of the transcript.
    Printing the object will return a string with all the information of the transcript and its sequence"""
    
    def __init__(self, transcript=None, coordinates="whole", flank=0, fragment=None, enzymes=None, include_introns=False, species="human"):
        
        #All available information about the transcript is stored as atributes. Even the information acquired from this class methods
        if transcript is None:
            raise ValueError("A transcript must be inputed for a SearchInfo-class object to be created. Valid examples = 'ENST00000357654' or 'ENST00000357654.9'")
        else:
            self.transcript = transcript
        self.coordinates = coordinates
        self.flank = flank
        self.fragment = fragment
        self.enzymes = enzymes
        self.include_introns = include_introns
        self.species = species
        #This class uses its own methods to further increase the completeness of its atributes
        self.chromosome, self.genomic_coordinates, self.strand, self.start_cdna, self.end_cdna, self.refseq_id = self.get_info()
        self.sequence = self.get_sequence()
        
    def __str__(self): #Printing the object will return a brief desccription of the transcript
         
        return f"The sequence inside the transcript {self.transcript} of {self.species} species, with its cDNA coordinates being {self.start_cdna}-{self.end_cdna} in chromosome {self.chromosome}, inside the strand {self.strand} and flanking regions of {self.flank} nucleotides is {self.sequence}"
    
    def get_info(self): 
        """This method uses the Ensembl REST API, alongside the inputed information, to retrieve valuable information from the transcript"""
        
        #===EXTRACTING COORDINATES====
        #The program uses the previously defined coordinate_extractor to check and separate the starting and ending coordinates from input.
        input_start, input_end, coord_type = coordinate_extractor(self.coordinates)
        
        #___API REST ENSEMBL___
        #===cDNA COORDINATES===
        if coord_type == "cDNA": #If the inputed coordinates are in cDNA form (The url will be different if the coordinate input was exonic or full transcript)
            
            url = f"https://rest.ensembl.org/map/cds/{self.transcript}/{input_start}..{input_end}?"
            headers = {"Content-Type": "application/json"}
            
            #The program uses the previously defined ensembl_requests() function to access Ensembl and retrieve the additional information
            response = ensembl_request(url, headers)
            data = response.json()
            
            #The information is retrieved from the json (which is represented as a dictionary of lists of dictionaries)
            chromosome = data["mappings"][0]["seq_region_name"]
            strand = data["mappings"][0]["strand"]
            refseq_id = None #This json does not contain the refseq_id (NM...) used to identify transcripts by NCBI
            
            #This variables will contain the starting and ending genomic coordinates from every exon in the interval
            start = []
            end = []
         
            for i in range (0, len(data["mappings"])):
                #The reason not to asume data["mappings"][i]["start"] is the smallest of the couple(start/end) and viceversa, is due to de existance of genes in the reverse strand
                #Currently, ensembl does store the smallest coordinate as the starting position of the exon, even on the reverse strand. But in case that may change, the program will be safer by not assuming it.
                start.append(min(data["mappings"][i]["start"], data["mappings"][i]["end"])) 
                end.append(max(data["mappings"][i]["start"], data["mappings"][i]["end"]))

            start = tuple(start) #This coordinates will be used in some methods ahead, so it's safer to make them immutable
            end = tuple(end)
            genomic_coordinates = {"start": start, "end": end} #Starting and ending genomic coordinates are stored in a dictionary        
            
            #The reason why the program returns the cDNA original coordinates too (in the form of input_start/end), which serve
            #no apparent use in here, is because in the nex step, with exonic coordinates or whole transcript, the program will need 
            #to acquire the cDNA coordinates (as those are not present in the input)        
            return chromosome, genomic_coordinates, strand, input_start, input_end, refseq_id 
        
        #===EXONIC AND WHOLE TRANSCRIPT===
        else:
           
            url = f"https://rest.ensembl.org/lookup/id/{self.transcript}?expand=1;utr=1;mane=1" #The url is different from the cDNA one
            headers = {"Content-Type": "application/json"}
            response = ensembl_request(url, headers)
            data= response.json()
            
            chromosome = data["Exon"][0]["seq_region_name"]
            strand = data["Exon"][0]["strand"]
            refseq_id = data["MANE"][0]["refseq_match"] #This json contains the refseq_id (NM...), used to identify transcripts by NCBI
            
            #Variable definition for further steps
            start = []
            end = []
            five_start = [] #When using exonic coordinates, UTR regions might be included
            five_end = []
            three_start = []
            three_end = []
            
            #Extraction of genomic coordinates from json    
            for i, coordinate in enumerate(data["Exon"]):
                coordinate = [data["Exon"][i]["start"], data["Exon"][i]["end"]]
                start.append(min(coordinate))
                end.append(max(coordinate))   

            #Extraction of UTR coordinates from json
            for i, utrinstance in enumerate(data.get("UTR", [])):
                coordinatesutr = [data["UTR"][i]["start"], data["UTR"][i]["end"]]
                
                if utrinstance["object_type"] == "five_prime_UTR":
                    five_start.append(min(coordinatesutr))
                    five_end.append(max(coordinatesutr))
                    
                elif utrinstance["object_type"] == "three_prime_UTR":
                    three_start.append(min(coordinatesutr))
                    three_end.append(max(coordinatesutr))
                    
            #UTR length calculus (Done in an accumulative way)
            accumulative_five = 0
            for start_utr, end_utr in zip(five_start, five_end):
                accumulative_five += max(start_utr, end_utr) - min(start_utr, end_utr) + 1
                
            accumulative_three = 0
            for start_utr, end_utr in zip(three_start, three_end):
                accumulative_three += max(start_utr, end_utr) - min(start_utr, end_utr) + 1
            
            #The program calculates the cDNA length without UTR regions
            cdna_length = -(accumulative_five + accumulative_three)
            for i,(startcoord, endcoord) in enumerate(zip(start, end)):
                cdna_length += max(startcoord, endcoord) - min(startcoord, endcoord) + 1
               
            #Only chosen exons coordinates are selected
            accumulative_length = -(accumulative_five) 
            start_cdna = None
            end_cdna = None
            no_utr_start = []
            no_utr_end = []
            
            if coord_type == "exon": #If inputed coordinates were exonic
                
                if input_start > len(start) or input_end > len(start): #If inputed exons are out of the real exon range of the transcript
                    raise ValueError(f"{self.transcript} only has len{start} exons and you introduced {input_start}_{input_end}. Please, introduce a valid number of exons")
                
                for i,(startcoord, endcoord) in enumerate(zip(start, end)): 
                    accumulative_length += max(startcoord, endcoord) - min(startcoord, endcoord) + 1
                    
                    if i+1 < input_start:
                        if accumulative_length < 0: #cDNA position is on 5'UTR
                            start_cdna = accumulative_length - 1
                        if accumulative_length >= 0:
                            start_cdna = accumulative_length + 1 #As c.0 does not exist, accumulative length needs to be corrected by + 1
                    
                    if input_start <= i+1 <= input_end:
                        no_utr_start.append(startcoord)
                        no_utr_end.append(endcoord)
                        if i+2 > input_end:
                            if start_cdna is None: #Then exon 1 has been chosen and the first conditional on the loop was skipped
                                start_cdna = -accumulative_five
                            if accumulative_length < 0: #Then end_cdna position will be in 5'UTR
                                end_cdna = accumulative_length - 1 
                            elif accumulative_length > cdna_length: #Then end_cdna position will be in 3'UTR
                                end_cdna = f"*{accumulative_length - cdna_length }"
                            elif accumulative_length <= cdna_length: #Then end_cdna position will be in the codificant region
                                end_cdna = accumulative_length  

            #The starting coordinates for each exon are included in its respective tuples and merged in a dictionary (same as with cDNA)
            no_utr_start = tuple(no_utr_start)
            no_utr_end = tuple(no_utr_end)
            genomic_coordinates = {"start": no_utr_start, "end": no_utr_end}
            
            if coord_type == "whole transcript": #If there were no inputed coordinates, the program will take the whole transcript
                start_cdna = -accumulative_five #Starting from the beggining of 5'UTR
                end_cdna = f"*{accumulative_three}" #And finishing at the end of 3'UTR
                
                start = tuple(start)
                end = tuple(end)
                genomic_coordinates = {"start": start, "end": end}
            
            return chromosome, genomic_coordinates, strand, start_cdna, end_cdna, refseq_id
    
    def get_sequence(self, genomic_coords = None, chromosome = None, strand = None, flank = None, include_introns = None, species = "human"): 
        """This method uses the previously defined get_sequence() auxiliar function to retrieve the DNA sequence, out from the information available
As stated before, get_sequence() auxiliar function uses a Ensembl REST API, as well as the genomic positions previously acquired, to extract the sequence"""

        #Inputed information is checked. If retreieved from get_info() method, everything should be alright, but in case the
        #user is using this method independently, the input is checked.
        if genomic_coords is None:
            genomic_coords = self.genomic_coordinates
        
        if chromosome is None:
            chromosome = self.chromosome
        
        if strand is None:
            strand = self.strand
        
        #The reason flank is set to None and not 0, is due to the fact that the default setting will set the flank value to self.flank. 
        #But maybe the user wants explicitly flank to be 0 for this specific method, while setting another value on class atributes. 
        #So by being None the default, if the user sets flank to 0 explicitly, self.flank will not overwrite it
        if flank is None: 
            flank = self.flank
            
        if include_introns is None: #The same as with the flank happens with setting include_introns to false.
            include_introns = self.include_introns
            
        if species == "human":
            species = self.species

        #The program uses the previously defined get_sequence() auxiliar function to retrieve the sequence from Ensembl
        return get_sequence(genomic_coords, chromosome, strand, flank , include_introns, species)
     