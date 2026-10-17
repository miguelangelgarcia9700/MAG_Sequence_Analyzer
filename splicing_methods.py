import re
import csv
from maxentpy import maxent
from auxiliar_functions import get_sequence, snp_mutator
from pathlib import Path

#This module contains two Sequence-class methods, extracted to reduce Sequence_class.py module size and make it more readable

def hexplorer_splicing_analysis(self, hex_score_matrix, reference_sequence=None, mutated_sequences=None, file_output=None, file_name=None, file_path="./results/hexplorer_splicing_analysis"):
    """Evaluate potential splicing-regulatory effects of single-nucleotide variants.

    Calculates the HEXplorer ΔHZEI score for variants in a DNA sequence.
    This ΔHZEI can either be positive or negative. The more positive, the more likely 
    is for the mutated sequence to have generated new enhancer sites around the mutation
    The more negative it is, the more likely is for the mutated sequence to have 
    generated new silencer sites around the mutation.

    If you want to study the whole inputed coordinates, you must include a flanking of 
    AT LEAST 11(RECOMENDED). If not, the first and last 11 instances will not be reliable
    Every flanking possition will be cut off the results, as its only use in here is to 
    offer enough context for the interval to be correctly analyzed.

    The analysis uses the supplied reference sequence and mutated sequences,
    either inputed in this method parameters, or present in the class atributes
    (for example, the mutated sequences are obtainde by default by interval_mutator()
    and the default reference is self.sequence)

    Optionally writes the results to a CSV file inside a "results" directory
    created locally.

    Parameters:
        hex_score_matrix: Score matrix required by the HEXplorer scoring
            function. Contains a score for each possible hexamer
        reference_sequence (str | None): Reference DNA sequence used for
            comparison. If None, the method uses self.sequence.
        mutated_sequences (dict[str, str] | None): Dictionary mapping variant
            labels to mutated sequences. If None, the method uses
            self.mutated_sequences (obtained by interval_mutator).
        file_output (str | None): Optional csv output format
            May be introduced as a string (file_output='csv')
        file_name (str | None): Set the name to the csv. If not set, default can be:
            'HEXplorer_splicing_analysis' if self.transcript is empty
            "{self.transcript}_'HEXplorer_splicing_analysis'" if self.transcript not empty
            "{self.transcript}_{self.coordinates}_'HEXplorer_splicing_analysis'" if 
            neither self.transcript nor self.coordinates are empty
        file_path(str): Route where the output file will be created.
            Default to "./results/hexplorer_splicing_analysis".  
            This are directories, the file name will be concatenated after

    Returns:
        dict [str, float]: Returns a dictionary, mapping variant identifiers to
            their ΔHZEI scores. {Variant: ΔHZEI}
        float: If only one variant was introduced, the method will just return
            its ΔHZEI. 

    Notes:
        The score estimates changes in predicted splicing-regulatory
        potential; it should not be interpreted as direct experimental
        evidence of altered splicing.
        
        Both reference_sequence and mutated_sequences must be inputed
        in the parameters for them to be used. If only one has been inputed
        the program will ignore it and use the atribute equivalent instead.
        
        The score matrix for the hexamer scores is present in a dictionary 
        form inside a different module called hex_matrix
        
        The method returns the results dictionary even when CSV output
        is requested.
    """
    if self.flank < 11: #The program will alert the user if a flanking shorter than one ha been inputed
        print("[ALERT] remember that hexplorer_splicing_analysis needs 11 nucleotides of flanking to be reliable. If you included them there should be no issue, but if you didn't, you probably should")

    csv_output = False
    if file_output is not None: #The program checks if the user wants a csv output
        
        if type(file_output) is str and file_output in ["csv"]:
            csv_output = True
            
        else:
            raise ValueError ("The only supported file output for hexplorer_splicing_analysis() method is 'csv', and must be introduced as a string. Example file_output='csv'")
        

    if file_name is None:
        file_name = "HEXplorer_splicing_analysis"
        
        if self.transcript is not None:
            file_name = f"{self.transcript}_HEXplorer_splicing_analysis"
            
            if self.coordinates is not None:
                file_name = f"{self.transcript}_{self.coordinates}_HEXplorer_splicing_analysis"
        
    elif file_name is not None and re.search(r'[\\/:*?"<>|]', file_name):
        raise ValueError("The 'file_name' parameter contains invalid characters for file paths (e.g. /, \\, :, *, ?, \", <, >, |).")

    #If the user has not inputed, either the mutated_sequences or the reference sequence in the parameters, and thus the program
    #would be facing the need to take one of those from the atributes and the other from the parameters.
    #The program will assume there is risk of both data not corresponding and will take BOTH from the atributes
    if mutated_sequences is None or reference_sequence is None:
        mutated_sequences = self.get_mutated_sequences()
        reference_sequence = self.sequence
        
    if mutated_sequences is None:
        #This should not be happening if the user is using the program correctly, as, in order to create a Sequecne-class object
        #the user must input a sequence, and thus, the program should always be able to perform an interval_mutator() method
        raise ValueError (""""No mutated sequences dictionary was included for the hexplorer_splicing_analysis(), please, either: 
    Input a valid dictionary (keys being the variant names and values being the corresponding mutated sequences) and assure the reference sequence is the wild type form of the mutated variants in the dictionary
    Perform a regular analysis by introducing the wild type sequence during the Sequence-class object creation, to enable the program to automatically perform a interval_mutator() method, generating the dictionary (recomended)
    """)
    if reference_sequence is None:
    #This should not be happening if the user is using the program correctly, as, in order to create a Sequecne-class object the user must input a sequence
        raise ValueError (""""No reference sequence was included for the hexplorer_splicing_analysis(), please, either: 
    Input a valid wild type sequence as reference, and assure the reference sequence is the wild type form of the mutated variants in the mutated_sequences dictionary
    Perform a regular analysis by introducing the wild type sequence during the Sequence-class object creation (recomended)""")
            
    if len(reference_sequence) >= 11: #A HEXplorer analysis needs at least a sequence 11 nucleotides long
        ref_seq = reference_sequence
        ref_seq_length = len(reference_sequence)
    else:
        raise ValueError ("reference_sequence must be at least 11 nucleotides long")

    if type(mutated_sequences) is str: #A single string may be analyzed, returning a single value in the terminal (if printed)
        if len(mutated_sequences) >= 11:
            name_list = []
            sequence_list = [mutated_sequences]
            csv_output = False
        else:
            raise ValueError ("mutated_sequences must be at least 11 nucleotides long")
        
    elif type(mutated_sequences) is dict: #A whole dictionary of (variants: sequences) may be analyzed, returning a dictionary and a .csv(if indicated)
        name_list = list(mutated_sequences.keys())
        sequence_list = list(mutated_sequences.values())
        for name, mut_seq in zip(name_list, sequence_list):
            if len(mut_seq) == ref_seq_length:
                continue
            else:
                raise ValueError (f"All mutated sequences must have the same length as the reference one, {name} does not")
    else:
        raise ValueError ("mutated_sequences input must be a string of a sequence or a dictionary with its keys being the sequence name and its values being the sequences(variants: sequences)")

    hexplorer_results = {} #This will be the dictionary containing the results
    ref_hzei_sum = 0

    #The program checks every hexamer containing the analyzed nucleotide in the current instance and retrieves its corresponding
    #value from hex_score_matrix. Then, the program calculates the mean of the 6 retrieved values, obtaining the positional hzei
    #Finally, every positional hzei is added, calculating the whole sequence hzei score (in this case for the reference)
    for i in range(0, ref_seq_length-10):
        hexamer_batch = {
            "hexamer1": ref_seq[0 + i : 6 + i],
            "hexamer2": ref_seq[1 + i : 7 + i],
            "hexamer3": ref_seq[2 + i : 8 + i],
            "hexamer4": ref_seq[3 + i : 9 + i],
            "hexamer5": ref_seq[4 + i : 10 + i],
            "hexamer6": ref_seq[5 + i : 11 + i]
        }
        hzei_score = 0
        for h in range(1,7):
            hzei_score += hex_score_matrix[hexamer_batch[f"hexamer{h}"]]
        hzei_score = hzei_score/6
        ref_hzei_sum += hzei_score    

    #The program executes the exact same analysis for each mutated sequence
    for name, mut_seq in zip(name_list, sequence_list):
        mut_hzei_sum = 0
        for i in range(0, ref_seq_length-10):
            hexamer_batch = {
                "hexamer1": mut_seq[0 + i : 6 + i],
                "hexamer2": mut_seq[1 + i : 7 + i],
                "hexamer3": mut_seq[2 + i : 8 + i],
                "hexamer4": mut_seq[3 + i : 9 + i],
                "hexamer5": mut_seq[4 + i : 10 + i],
                "hexamer6": mut_seq[5 + i : 11 + i]
            }
            hzei_score = 0
            for h in range(1,7):
                hzei_score += hex_score_matrix[hexamer_batch[f"hexamer{h}"]]
            hzei_score = hzei_score/6
            mut_hzei_sum += hzei_score
        
        #The reference sequence hzei score is substracted to the mutated sequence score, obtaining the resulting ΔHZEI of the mutation
        delta_score = mut_hzei_sum - ref_hzei_sum
        
        if not name_list: #If the input was a string, the program returns the single ΔHZEI value
            return round(delta_score, 2)
        
        else: #If the input was a mutated_sequences dictionary, the program will discard the flanks and return the dictionary with all the ΔHZEI scores
            if re.search(r"^Flank.*", name):
                continue
            hexplorer_results[name] = round(delta_score, 2)
            
    if csv_output is True:
        fieldnames = ["Variant", "ΔHZEI"]
        output_dir = Path(file_path) #The output path is defined
        output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
        file_path = output_dir / f"{file_name}.csv" #file name is concatenated
        with open(file_path, mode='w', newline='', encoding='utf-8-sig') as file:
            file.write("sep=,\n")
            writer = csv.writer(file)
            writer.writerow(fieldnames)
            for variant, hzei in hexplorer_results.items():
                    writer.writerow([variant, hzei])

    return hexplorer_results

def mes_canonical_analysis(self, genomic_coordinates = None, chromosome = None, strand = None, species = "human", coordinates = None, start_cdna = None, end_cdna = None, file_output=None, file_name=None, file_path="./results/mes_canonical_analysis"):
    """Evaluate canonical splice-site strength using MaxEntScan.
    
    The MaxEntScan score is a form of meassure for the splicing site strength

    Analyzes the acceptor and donor regions associated with the selected exons.
    For each analyzed splice site, calculates the percentage variation in the
    MaxEntScan score for every possible single-nucleotide polymorphism at the
    canonical positions.

    Exons whose splice sites fall within annotated 5' or 3' UTR regions are
    reported as being in UTR range rather than receiving the standard analysis.
    
    Optionally writes the results to a CSV file inside a "results" directory
    created locally.

    Parameters:
        genomic_coordinates (dict | None): Genomic coordinates of the
            exons, containing 'start' and 'end' coordinate lists. 
            If None, uses self.genomic_coordinates.
        chromosome (str | None): Chromosome containing the analyzed region.
            If None, uses self.chromosome.
        strand (int | None): Transcript strand (1 or -1). If None, uses
            self.strand or is attempted to be inferred from the exon coordinates.
        species (str): Species used for genomic sequence retrieval. Defaults
            to 'human', with the method using self.species when applicable.
        coordinates (str | None): Originally inputed coordinates, such as
            'ex.18', 'ex.18_20', or 'whole'. If None, uses self.coordinates.
            Only accepts exonic coordinates or 'whole', not cDNA coordinates.
        start_cdna (int | None): Starting cDNA coordinate, used to determine
            whether the selected interval includes a 5' UTR. If None, uses 
            self.start_cdna. If still None, the program assumes there is no 5'UTR
        end_cdna (int | str | None): Ending cDNA coordinate, used to determine
            whether the selected interval includes a 3' UTR. If None, uses 
            self.end_cdna. If still None, the program assumes there is no 3'UTR
        file_output (str | None): Optional csv output format
            May be introduced as a string (file_output='csv')
        file_name (str | None): Set the name to the csv. If not set, default can be:
            'MES_canonical_analysis' if self.transcript is empty
            "{self.transcript}_'MES_canonical_analysis'" if self.transcript not empty
            "{self.transcript}_{self.coordinates}_'MES_canonical_analysis'" if 
            neither self.transcript nor self.coordinates are empty
        file_path(str): Route where the output file will be created.
            Default to "./results/mes_canonical_analysis"
            This are directories, the file name will be concatenated after

    Returns:
        dict [str, list]: A dictionary whose keys are the exons acceptors/donnors
            and whose values are lists with the description of the splicing-strength
            variation for every possible variant. 
            {exon_acceptor/donnor: [Variant_nt1>nt2 = % of strength variation, variant_nt1>nt3...]}

    Raises:
        ValueError: If required genomic metadata is unavailable, the
            coordinate specification is unsupported, or the output format
            is invalid.

    Notes:
        The method retrieves genomic sequence from Ensembl API REST and uses
        MaxEntScan scoring functions for acceptor and donor sites. 
        
        Percentage variations are relative to the reference score.
        
        The method returns the results dictionary even when CSV output
        is requested.
    """
    csv_output = False #The program checks if the user wants a csv output
    if file_output is not None:
        
        if type(file_output) is str and file_output in ["csv"]:
            csv_output = True
            
        else:
            raise ValueError ("The only supported file output for mes_canonical_analysis() method is 'csv', and must be introduced as a string. Example file_output='csv'")
    
    if file_name is None:
        file_name = "MES_canonical_analysis"
        
        if self.transcript is not None:
            file_name = f"{self.transcript}_MES_canonical_analysis"
            
            if self.coordinates is not None:
                file_name = f"{self.transcript}_{self.coordinates}_MES_canonical_analysis"
        
    elif file_name is not None and re.search(r'[\\/:*?"<>|]', file_name):
        raise ValueError("The 'file_name' parameter contains invalid characters for file paths (e.g. /, \\, :, *, ?, \", <, >, |).")
        
    if genomic_coordinates is None: #This method will need to access Ensembl API REST to retrieve a sequence, so it needs the genomic coordinates
        if self.genomic_coordinates is not None:
            genomic_coordinates = self.genomic_coordinates
        else:
            raise ValueError ("Mes_canonical_analysis needs genomic coordinates to work. This genomic coordinates must be a dictionary with a key being 'start' containing every exon starting coordinate and another key called 'end' containing every exon ending coordinate")
    
    if chromosome is None: #This method will need to access Ensembl API REST to retrieve a sequence, so it needs the chromosome
        if self.chromosome is not None:
            chromosome = self.chromosome
        else:
            raise ValueError ("Mes_canonical_analysisi needs the transcript chromosome to work.")
        
    if strand is None: #This method will need to access Ensembl API REST to retrieve a sequence, so it needs the strand
        if self.strand is not None:
            strand = self.strand
        elif len(genomic_coordinates["start"]) > 1: #If more than one exon was inputed, the program may be able to deduce the strand
            if genomic_coordinates["start"][1] - genomic_coordinates["start"][0] < 0:
                strand = -1
            else:
                strand = 1
        else:
            raise ValueError ("Mes_canonical_analysis needs the transcript chromosome to work.")
    
    if species == "human": #As human is the default, it will take self.species in case the user changed it. Though if user didn't change it, it will still be human.
        species = self.species
    
    if coordinates is None: #MES_canonical_analysis needs to know what exons are to be analyzed
        if self.coordinates is not None:
            coordinates = self.coordinates
        else:
            raise ValueError ("""MES_canonical_analysis() needs the exonic coordinates to work. As it needs to know what exons to analize
The only coordinate inputs supported for MES_canonical_analysis() are whole exons (-c ex.18 or -c ex.18_20) or the whole transcript (not writing -c at all)""")    

    #Only exonic coordinates or whole transcript is allowed as coordinate input, as MaxEntScan is a whole-exon analysis (as it checks its acceptor and donnor)
    if match := re.search(r"^ex\.(\d+)(?:_(\d+))?$", coordinates, re.IGNORECASE): 
        current_exon = int(match.group(1))
    elif coordinates.lower() == "whole":
        current_exon = 1
    else:
        raise ValueError ("The only coordinate inputs supported for MES_canonical_analysis() are whole exons (-c ex.18 or -c ex.18_20) or the whole transcript (not writing -c at all)")
    
    if start_cdna is None: #start_cdna coordinate enables the program to know the pressence or abscence of 5'UTR
        if self.start_cdna is not None:
            start_cdna = self.start_cdna
        else:
            print("[ALERT] You have not inputed a start_cdna coordinate for MES_canonical_analysis, thus the program will asume there is NO 5'UTR region in your interval")
            start_cdna = 1 ##No 5'UTR nomenclature will not trigger the serch for the exon in which it beggins, which is the only use for start_cdna in this method. Thus 1 is selected as its generic nomenclature
    
    if end_cdna is None: #start_cdna coordinate enables the program to know the pressence or abscence of 3'UTR 
        if self.end_cdna is not None:
            end_cdna = self.end_cdna
        else:
            print("[ALERT] You have not inputed an end_cdna coordinate for MES_canonical_analysis, thus the program will asume there is NO 3'UTR region in your interval")
            end_cdna = 1 #No 3'UTR nomenclature will not trigger the serch for the exon in which it beggins, which is the only use for end_cdna in this method. Thus 1 is selected as its generic nomenclature
    
    #In this method, flank must be set to 0 and include_introns must be false
    flank = 0
    include_introns = False
    genomic_max_index = len(genomic_coordinates["start"]) - 1 #This variable will help the program to traverse through all the exons
    
    if start_cdna < 0: #There must be a 5'UTR
        accumulated_five_length = 0
        for i in range(0, genomic_max_index + 1): 
            accumulated_five_length += genomic_coordinates["end"][i] - genomic_coordinates["start"][i] + 1 #first exon length is accumulated
            if accumulated_five_length < -start_cdna: #If accumulation is smaller than the size of 5'UTR, there must be more exons inside 5'UTR
                continue
            else:
                last_five_exon = i + current_exon - 1 #This will be the last exon presenting 5'UTR
    else:
        last_five_exon = 0 #This will always be smaller than any current_exon
            
    if type(end_cdna) is str: #If end_cdna is a string, its nomenclature must be *(number) and thus, there is 3'UTR
        if match := re.search(r"^\*(\d+)", end_cdna):
            three_utr_length = int(match.group(1))
            accumulated_three_length = 0
            for i in range(0, genomic_max_index + 1):
                accumulated_three_length += genomic_coordinates["end"][genomic_max_index-i] - genomic_coordinates["start"][genomic_max_index-i] + 1
                if accumulated_three_length < three_utr_length:
                    continue
                else:
                    first_three_exon = genomic_max_index - i + current_exon + 1 #This will be the first exon presenting 3'UTR
    else:
        first_three_exon = current_exon + genomic_max_index + 2 #This will always be bigger than any current_exon
    
    #This variables will contain the starting and ending points of acceptors and donnors, distributed like [A,D,A,D...]
    new_start_list = []
    new_end_list= []

    if strand == 1:
        for start_coord, end_coord in zip(genomic_coordinates["start"], genomic_coordinates["end"]):
            new_start_list.append(start_coord - 20)
            new_end_list.append(start_coord + 2)
            new_start_list.append(end_coord - 2)
            new_end_list.append(end_coord + 6)
    elif strand == -1:
        for start_coord, end_coord in zip(genomic_coordinates["start"], genomic_coordinates["end"]):
            new_start_list.append(end_coord - 2)
            new_end_list.append(end_coord + 20)
            new_start_list.append(start_coord - 6)
            new_end_list.append(start_coord + 2)
        
    new_genomic_coordinates = {"start": new_start_list, "end": new_end_list} #The lists are merged into a dictionary suitable for get_sequence()
    #The program will obtained a concatenated sequence of acceptors and donnors of all the exons present in the interval
    concatenated_sequence = get_sequence(new_genomic_coordinates, chromosome, strand, flank, include_introns, species)
    sequence_position = 0
    mes_result = {} #This variable will store the result
    nts = ["A", "C", "G", "T"]
    current_exon -= 1 #This is performed, so that in next steps, when adding one per cycle, it returns the expected exon

    for i in range (0,len(new_genomic_coordinates["start"])): #Every instance of the loop will be an exon
        acceptor_list = []
        donor_list = []

        if i%2 == 0: #As acceptors go first, if the index is even, it must be an acceptor
            is_acceptor = True
            is_donor = False
            current_exon += 1 #Every time the program reaches an acceptor, it has advanced one exon
            mes_result[f"exon{current_exon}"] = {} #As the program reaches a new exon, it generates its key in the result dictionary
        if i%2 != 0:
            is_acceptor = False
            is_donor = True

        if (current_exon <= last_five_exon or current_exon > first_three_exon) and is_acceptor is True:
        #Which translates into, if the program is on any exon with 5'UTR or after the first exon presenting 3'UTR and is checking an acceptor
            acceptor_list.append("Acceptor is in UTR range")
            mes_result[f"exon{current_exon}"]["Acceptor strength variation"] = acceptor_list
            sequence_position += 23
        
        elif (current_exon < last_five_exon or current_exon >= first_three_exon) and is_donor is True:
        #Which translates into, if the program is on the last exon with 5'UTR or on any exon presenting 3'UTR and is checking a donor
            donor_list.append("Donor is in UTR range")
            mes_result[f"exon{current_exon}"]["Donor strength variation"] = donor_list
            sequence_position += 9
            
        elif current_exon > last_five_exon and current_exon <= first_three_exon and is_acceptor is True:
        #Which translates into, if program is after 5'UTR region and before or at the first exon presenting 3'UTR and is checking an acceptor
            instance_seq = concatenated_sequence[sequence_position : (sequence_position + 23)]
            mes_ref = maxent.score3(instance_seq) #The reference MES score is calculated
            for ind in range (18,20): #canonical positions will be mutated
                ref_base = instance_seq[ind]
                for alt_base in nts:
                    if alt_base == ref_base: continue
                    mutated_seq = snp_mutator(instance_seq, ind+1, alt_base)
                    mes_mut = maxent.score3(mutated_seq) #Mutated acceptor MES score is calculated
                    if mes_ref != 0: #The percentage of variation is calculated
                        percentage_variation = f"{round(((mes_mut - mes_ref)/mes_ref)*100, 2)}%"
                    else:
                        percentage_variation = f"MES score from reference was 0, an error may have occurred"
                    acceptor_list.append(f"Pos({ind-20}){ref_base}>{alt_base} = {percentage_variation}")
            mes_result[f"exon{current_exon}"]["Acceptor strength variation"] = acceptor_list #The analysis of the acceptor is added to the results
            sequence_position += 23
        
        elif current_exon >= last_five_exon and current_exon < first_three_exon and is_donor is True:
        #Which translates into, if program is after or at the last exon presenting 5'UTR region and before 3'UTR and is checking a donor 
            instance_seq = concatenated_sequence[sequence_position : (sequence_position + 9)]
            mes_ref = maxent.score5(instance_seq)
            for ind in range (3,5):
                ref_base = instance_seq[ind]
                for alt_base in nts:
                    if alt_base == ref_base: continue
                    mutated_seq = snp_mutator(instance_seq, ind+1, alt_base)
                    mes_mut = maxent.score5(mutated_seq)
                    if mes_ref != 0:
                        percentage_variation = f"{round(((mes_mut - mes_ref)/mes_ref)*100, 2)}%"
                    else:
                        percentage_variation = f"MES score from reference was 0, an error may have occurred"
                    donor_list.append(f"Pos({ind-2}){ref_base}>{alt_base} = {percentage_variation}")
            mes_result[f"exon{current_exon}"]["Donor strength variation"] = donor_list
            sequence_position += 9    

    if csv_output is True:
        fieldnames = ["Exon", "Site", "Variant", "Relative_strength_variation"]
        output_dir = Path(file_path) #The output path is defined
        output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
        file_concatenated_path = output_dir / f"{file_name}.csv" #file name is concatenated
        
        with open(file_concatenated_path, mode='w', newline='', encoding='utf-8-sig') as file:
            file.write("sep=,\n")
            writer = csv.writer(file)
            writer.writerow(fieldnames)
            for ex, ex_info in mes_result.items(): #The data for the different columns is extracted from the result dictionary
                for site, site_info in ex_info.items():
                    if site == "Acceptor strength variation":
                        site = "Acceptor"
                    elif site == "Donor strength variation":
                        site = "Donor"
                    for variant_info in site_info:
                        match = re.search(r"^([^=]+) = ([^=]+)$", variant_info)
                        variant = match.group(1)
                        strength_variation = match.group(2)
                        writer.writerow([ex, site, variant, strength_variation])
                    
    return mes_result