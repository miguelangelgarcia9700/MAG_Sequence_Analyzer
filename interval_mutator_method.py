import re
import csv
from auxiliar_functions import snp_mutator, exon_intron_boundaries
from pathlib import Path

#This module contains one Sequence-class method, extracted to reduce Sequence_class.py module size and make it more readable

def interval_mutator(self, transcript=None, start_cdna=None, end_cdna=None, genomic_coordinates=None, include_flanking=True, strand=None, intron_HGVS_limit=15, file_output=None, file_name=None, file_path="./results/interval_mutator"): 
    """Generate all possible single-nucleotide polymorphisms in a sequence.

    For every nucleotide in the selected sequence, generates the three
    alternative nucleotide substitutions and stores each mutated sequence
    in a dictionary. When sufficient transcript and genomic metadata are
    available, variants are labelled using cDNA/HGVS-style notation.
    Otherwise, generic position labels are used.

    Flanking regions can be excluded from the generated variants. Intronic
    positions are described using HGVS-style intronic notation up to the
    specified limit; positions farther inside an intron receive generic
    intronic labels.
    
    Optionally writes the results to a CSV file, a fasta file or both 
    inside a "results" directory created locally.
    (check parameter file output)

    Parameters:
        transcript (str | None): Transcript ensembl identifier.
            If None, uses self.transcript.
        start_cdna (int | None): Starting cDNA coordinate of the interval.
            If None, uses self.start_cdna. If negative, indicates 5'UTR pressence
        end_cdna (int | str | None): Ending cDNA coordinate of the interval.
            If None, uses self.end_cdna. A string such as '*100' indicates 3' UTR.
        genomic_coordinates (dict | None): Genomic coordinates of the
            exons, containing 'start' and 'end' coordinate lists. 
            If None, uses self.genomic_coordinates.
        include_flanking (bool): If True, includes available flanking
            nucleotides in the mutation analysis. Defaults to True
        strand (int | None): Transcript strand (1 or -1)
            If None, uses self.strand.
        intron_HGVS_limit (int): Maximum intronic distance from an exon
            boundary described using HGVS-style notation. Defaults to 15.
        file_output (str | list[str] | tuple[str, ...] | None): Optional
            output format: 'csv', 'fasta', or a list/tuple containing both.
            example: file_output='csv' or 'fasta' or ['csv', 'fasta']
        file_name (str | None): Set the name to the csv. If not set, default can be:
            'Interval_mutator' if self.transcript is empty
            "{self.transcript}_'Interval_mutator'" if self.transcript not empty
            "{self.transcript}_{self.coordinates}_'Interval_mutator'" if neither
            self.transcript nor self.coordinates are empty
        file_path(str): Route where the output file/s will be created.
            Default to "./results/interval_mutator"
            This are directories, the file name will be concatenated after

    Returns:
        dict[str, str]: A dictionary mapping variant labels to their
            corresponding mutated DNA sequences. {Variant: Mutated_sequence}

    Raises:
        ValueError: If the output format is unsupported or the intronic
            notation limit is less than 1.

    Notes:
        The method returns the results dictionary even when CSV and/or fasta
        output is requested.
        """
    #___INPUT CHECKERS___
    #Method parameters, if inputed, will have priority over their atribute equivalents.
    #Input is checked:
    csv_output = False
    fasta_output = False
    if file_output is not None:
        
        if type(file_output) is str:
            if file_output.lower() == "csv":
                csv_output = True
            elif file_output.lower() == "fasta":
                fasta_output = True
            else:
                raise ValueError ("The only supported file outputs for interval_mutator() method are 'csv' and 'fasta'. Example file_output='csv' or file_output=('csv', 'fasta') or file_output=['csv', 'fasta']")
        
        elif type(file_output) is list or type(file_output) is tuple:
            for output_type in file_output:
                if output_type.lower() == "csv":
                    csv_output = True
                elif output_type.lower() == "fasta":
                    fasta_output = True
                else:
                    raise ValueError ("The only supported file outputs for interval_mutator() method are 'csv' and 'fasta'. Example file_output='csv' or file_output=('csv', 'fasta') or file_output=['csv', 'fasta']") 
        
        else:
            raise ValueError ("The only supported file outputs for interval_mutator() method are 'csv' and 'fasta'. Example file_output='csv' or file_output=('csv', 'fasta') or file_output=['csv', 'fasta']")            
    
    if file_name is None:
        file_name = "Interval_mutator"
        
        if self.transcript is not None:
            file_name = f"{self.transcript}_Interval_mutator"
            
            if self.coordinates is not None:
                file_name = f"{self.transcript}_{self.coordinates}_Interval_mutator"
        
    elif file_name is not None and re.search(r'[\\/:*?"<>|]', file_name):
        raise ValueError("The 'file_name' parameter contains invalid characters for file paths (e.g. /, \\, :, *, ?, \", <, >, |).")
    
    if transcript is None:
        transcript = self.transcript
        
    if strand is None:
        strand = self.strand
        
    if start_cdna is None:
        start_cdna = self.start_cdna
    
    five_prime_utr = False
    if start_cdna is not None: #Either because start_cdna parameter was not None, or because self.start_cdna was not none.
        if start_cdna < 0: #The presence of a 5'UTR is checked
            five_prime_utr = True
            
    if end_cdna is None:
        end_cdna = self.end_cdna
        
    three_prime_utr = False
    if type(end_cdna) is str: #If end_cdna is a string, thus it is not None, while it also indicates de presence of 3'UTR
        match = re.search(r"^\*(\d+)",end_cdna) #This would mean that 3'UTR is present
        end_cdna = int(match.group(1))
        three_prime_utr = True
        three_prime_counter = 0
        
    if genomic_coordinates is None:
        genomic_coordinates = self.genomic_coordinates
    
    flank = self.flank #It would make no sense to include a flank different from the one inputed in SearchInfo, used to retrieve the sequence.
    if flank == 0: #The user may want to exclude the complete flanking though, thats why a include_flanking parameter exists
        include_flanking = False #When flank=0, include_flanking is automatically set to false
        
    if intron_HGVS_limit < 1:
        raise ValueError("intron_HGVS_limit must be 1 or higher")
    
    #___SEQUENCE MODIFIERS___
    #Flanks and/or UTR regions, previously included in extracting the sequence, but undesired in interval_mutator() have to be removed    
    if include_flanking is False: 
    #Thus, the program will cut off any flank that may have been included while acquiring the sequence.
        sequence = self.sequence[flank:(len(self.sequence) - flank)]
    else:
        sequence = self.sequence    
    
    #___EXTRA INFORMATION GETTER___
    # #If HGVS is possible to be written, the program keeps acquiring information. If all minimal input is present, HGVS is possible
    hgvsnom = True
    if transcript is None or start_cdna is None or end_cdna is None or genomic_coordinates is None:
        hgvsnom = False #The program will not be able to call the variants in an HGVS format and will name them in a generic way
        
    if hgvsnom is True:
        #===INTRONS INCLUDED?=== 
        #The program checks if introns where included while retrieving the sequence by comparing its length to the max possible genomic length 
        if len(genomic_coordinates["start"]) > 1 and (max(max(genomic_coordinates["start"]), max(genomic_coordinates["end"])) - min(min(genomic_coordinates["start"]), min(genomic_coordinates["end"])) + 1 + (flank * 2)) == len(self.sequence): 
        #If the sequence length is the biggest possible one (from the smallest GENOMIC coordinate to the biggest), it wolud include introns
            no_introns = False
        else: #Whereas if its shorter, then the introns must be absent
            no_introns = True
            
        #Exonic and intronic boundaries are only needed if there are introns, and strand is only needed to calculate those boundaries    
        if not no_introns: #If there are introns  
            #===STRAND?===
            #If the strand was not inputed the program may be able to calculate it
            if strand is None and len(genomic_coordinates["start"]) > 1:
                if genomic_coordinates["start"][1] - genomic_coordinates["start"][0] < 0:
                    strand = -1
                else:
                    strand = 1
                    
            #===EXONIC AND INTRONIC BOUNDARIES?===
            #The program uses the previously defined exon_intron_boundaries() auxiliar function to calculate it exonic and intronic lengths
            exon_intron_difference = exon_intron_boundaries(include_flanking, flank, genomic_coordinates, strand)
            
        else: #If there are no introns
            exon_intron_difference = [len(sequence) - (flank if include_flanking else 0)] #If there are no introns, exon_intron_difference just needs to exist, have 1 element and represent the whole cdna length.
    
    #___INDEPENDENT COMPONENTS NEEDED___
    mutated_sequences = {} #This dictionary will be returned, filled with the mutated sequences
    nts = ["A", "C", "G", "T"] 
    exon_intron_index = 0 #This will position the program in exon_intron_difference right index depending on sequence position (0 or even = exon and odd = intron)
    cdna_counter = -1
    intron_known_counter = 0
    intron_unknown_counter = 0
    intron_counter = 1 #If there were multiple introns, this counter will help name its unknown area (Between the intronic HGVS regions)
    
    #SEQUENCE MUTATOR
    #___POSITION LABELER___
    for i, ref_base in enumerate(sequence):
        #=====GENERIC LABELER===
        #Should any of this be absent, either not registered in the class atributes nor inputed in the method parameters, HGVS will not be able to be written.
        #Thus, only this conditional will exeute. Otherwise, the rest will.
        if hgvsnom is False:
            label = f"Position{i+1}" #A generic label is created for every variant
        
        #=====UPSTREAM FLANK=====
        elif include_flanking and i < flank: #If flanks have been inputed, the first labels to be defined are the upstream ones
            label = f"Flank{i - flank}"
        
        #=====EXON=====
        #When the upstream flank (if present) is finished, it comes the turn for the first included exon to be labeled.
        #Up from here, the exons and the introns will take turns to be labeled. 
        #The program will know if its currently on an exon or an intron thanks to exon_intron_index, as it happened with the get_sequence() function before
        elif exon_intron_index%2 == 0: #If exon_intron_index is even, the program is in an exon
            cdna_counter += 1
            cdna_pos = start_cdna + cdna_counter #The program starts on the first cDNA coordinate and will add 1 per iteration
            
            if five_prime_utr is True and cdna_pos == 0: #This would mark the transition between 5'UTR and codificant transcript, and c.0 does not exist
                cdna_pos += 1 #Current cycle is corrected
                cdna_counter += 1 #Next cycles are corrected (If not, next cycle would be c.1 again and would overwrite the result)
            if three_prime_utr is True and i + 1 >= len(sequence) - flank - end_cdna + 1: #This would mark the transition between the codificant transcript and 3'UTR
                three_prime_counter += 1
                cdna_pos = f"*{three_prime_counter}"

            label = f"{transcript}:c.{cdna_pos}"
            
            #Region transitioner
            if i + 1 == exon_intron_difference[exon_intron_index]: #If the length limit of the exon is reached, exon_intron_difference advances 1 index
                exon_intron_index += 1
                
        #=====INTRON=====
        elif exon_intron_index < len(exon_intron_difference) and exon_intron_index%2 != 0: #If exon_intron_index is odd, the program is in an intron
            
            #===DOWNSTREAM===
            #The program is in the downstream region after the exon that has just been labeled and before the intron_HGVS_limit currently set
            if i < exon_intron_difference[exon_intron_index - 1] + intron_HGVS_limit:
                intron_known_counter += 1
                intron_unknown_counter = 0
                label = f"{transcript}:c.{cdna_pos}+{intron_known_counter}"
            
            #===UNKNOWN===
            #The program is between both HGVS intronic regions
            elif exon_intron_difference[exon_intron_index - 1] + intron_HGVS_limit <= i < exon_intron_difference[exon_intron_index] - intron_HGVS_limit:
                intron_unknown_counter += 1
                label = f"Unknown intronic {intron_counter}:{intron_unknown_counter}"
                
            #===UPSTREAM===
            #The program is in the upstream region before the exon about to be labeled and after the intron_HGVS_limit currently set
            elif i >= exon_intron_difference[exon_intron_index] - intron_HGVS_limit:
                intron_known_counter -= 1
                label = f"{transcript}:c.{cdna_pos+1}-{intron_known_counter+1}"
            
            #Region transitioner
            if i + 1 == exon_intron_difference[exon_intron_index]: #If the length limit of the intron is reached, exon_intron_difference advances 1 index
                exon_intron_index += 1
                intron_known_counter = 0 #Intron known counter resets for next intron
                intron_counter += 1 
        
        #=====DOWNSTREAM FLANK=====
        #When the entire inputed interval has been labeled, its the turn of the downstream flank (if present).       
        elif include_flanking and i >= len(sequence) - flank: 
            downstream_flank_position = i - len(sequence) + flank + 1
            label = f"Flank{downstream_flank_position}"  
        
        #___BASE MUTATOR___
        #Now that the labels have been created, te program will actually mutate the sequence.
        for alt_base in nts:
            if alt_base == ref_base: #If the alternative base and the reference one are the same, the loop skips the itaration
                continue 
            else: #If not, the program executes the mutation
                mutated_seq = snp_mutator(sequence, i+1, alt_base)
                
            #___INTRODUCTION IN DICTIONARY___
            #The program now corresponds the current label with its mutation as key. And its respective mutated sequence as value
            mutated_sequences[f"{label}{ref_base}>{alt_base}"] = mutated_seq
    
    if csv_output is True: #If the user wants a csv output
        fieldnames = ["Variant", "Sequence"]
        output_dir = Path(file_path) #The output path is defined
        output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
        file_concatenated_path = output_dir / f"{file_name}.csv" #file name is concatenated
        
        with open(file_concatenated_path, mode='w', newline='', encoding='utf-8-sig') as file:
            writer = csv.writer(file)
            writer.writerow(fieldnames)
            for variant, instance_sequence in mutated_sequences.items():
                    writer.writerow([variant, instance_sequence])
                    
    if fasta_output is True: #If the user wants a fasta output
        output_dir = Path(file_path) #The output path is defined
        output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
        file_concatenated_path = output_dir / f"{file_name}.fasta" #file name is concatenated
        
        with open(file_concatenated_path, 'w', encoding='utf-8') as file:
            for variant, instance_sequence in mutated_sequences.items():
                file.write(f">{variant}\n")
                file.write(f"{instance_sequence}\n")

    #The program returns the dictionary with the mutated sequences and the wild-type sequence, which might not be the same
    #as the inputed one if the flanks have been excluded
    return mutated_sequences