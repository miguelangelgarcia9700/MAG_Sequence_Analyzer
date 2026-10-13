import re
from Bio.Restriction import RestrictionBatch, Analysis #This will be installed with biopython 
from Bio.Seq import Seq #This will be installed with biopython
from Bio import Entrez
import csv
from maxentpy import maxent
from auxiliar_functions import ensembl_request, get_sequence, snp_mutator, exon_intron_boundaries

class Sequence: 
    """This class includes several methods to analyze the sequence previously acquired by the SearchInfo class, or manually inputed (if used independently).
    This includes:
    """
    #In this program, Sequence class stores the sequence retrieved from SearchInfo, alongside the rest of its available information
    def __init__(self, sequence=None, transcript=None, start_cdna=None, end_cdna=None, chromosome = None, flank=0, coordinates = None, genomic_coordinates=None, strand = None, refseq_id=None, species="human"):
        if sequence is None:
            raise ValueError("Sequence is required to create a Sequence object")
        self.sequence = sequence
        self.transcript = transcript
        self.start_cdna = start_cdna
        self.end_cdna = end_cdna
        self.chromosome = chromosome
        self.flank = flank
        self.coordinates = coordinates
        self.genomic_coordinates = genomic_coordinates
        self.strand = strand
        self.refseq_id = refseq_id
        self.species = species
        self.mutated_sequences = None
    
    def get_mutated_sequences(self): #THIS IS AN AUXILIAR METHOD MEANT TO MAKE OTHER METHODS AHEAD MORE EFFICIENT (through lazy catch), THUS IT HAS NO REAL USE ON ITS OWN.
        """In long intervals, interval_mutator() method may take a long time to run. Thus, using it only once per run can be crucial
So, instead of writing directly self.mutated_sequences, self.modified_sequence = self.interval_mutator() on the atributes, this method will run the interval_mutator() method only when it's called"""
        
        if self.mutated_sequences is None: #Only call self.interval_mutator when self.mutated_sequences=None, which means it hasn't been already used
            self.mutated_sequences = self.interval_mutator()
            
        return self.mutated_sequences, self.modified_sequence
    
    def translate_dna_to_protein(self, dna_sequence=None, frame=0, table=1, stop_symbol="*", to_stop=False, gap="-"):
        """This method uses biopython's translate() function to translate a given DNA sequence to its protein form. 
    But, unlike biopython, by setting frame parameter to 0 (or not setting it at all), the method allows to automatically translate the sequence in all of its 3 possible frames

    Parameters:
        dna_sequence (str): The DNA sequence to translate.
        frame (int): The reading frame (0, 1, 2, or 3). Default is 0, returning all 3 possible frame translations. 1 Will start in the first nucleotide, 2 in the second and 3 in the third
        table (int): The translation table number. Default is 1. (More information on translation tables can be found in byopython documentation)
        stop_symbol (str): The character representing a stop codon. Default is "*".
        to_stop (bool): Whether to stop translation at the first stop codon. Default is False.
        gap (str): The character representing a gap in the protein sequence. Default is "-".
        traslate() function from biopython, may also include as parameter cds=bool, which evaluates certain characteristics in the sequence to determine if its a coding sequence. But as the aim of this method is to merely translate, it hasn't been included.
    Returns:
        dict: A dictionary containing the translated protein sequences for each frame.
    
Any value included in this method parameters will have priority over its class-atribute equivalents
        """
        
        #The program checks the input
        if dna_sequence is None:
            if self.sequence is not None:
                dna_sequence = self.sequence
            else:
                raise ValueError ("translate_dna_to_protein() method needs a dna_sequence, either in its parameters or in its object atributes (self.sequence)")
            
        if type(dna_sequence) is not str:
            raise ValueError("Invalid dna_sequence value. DNA sequence must be a string of nucleotides (A, T, C, G). Example: 'ATGCCTGATTGC'")
        
        if type(frame) is not int or frame not in [0, 1, 2, 3]:
            raise ValueError("Invalid frame value. Frame must be 0 (or not specified) for all frames, or 1, 2, or 3 for its respective frames.")
        
        if type(table) is not int or table < 1:
            raise ValueError("Invalid table value. Table must be a positive integer representing the translation table. Example, 1 for the standard genetic code.")
        
        if type(stop_symbol) is not str or len(stop_symbol) != 1:
            raise ValueError("Invalid stop_symbol value. Stop symbol must be a single character string. Default is '*'")
        
        if type (to_stop) is not bool:
            raise ValueError("Invalid to_stop value. to_stop must be a boolean value (True or False). Default is False")
        
        if type(gap) is not str or len(gap) != 1:
            raise ValueError("Invalid gap value. Gap must be a single character string. Default is '-'")
        
        dna_seq = Seq(dna_sequence) #Biopython needs the sequence to be a Seq-class object
        protein_translation = {} #Translation will be stored in this variable, which is a dictionary in case the user asked for the 3 possible translation frames
        
        if frame == 0: #If frame is not specified or set to 0, the program will translate in all possible frames
            for i in range(0,3):
                #The program uses the translate method from a Seq-class object from Biopython
                translation = dna_seq[i:].translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, gap=gap)
                protein_translation[f"frame_{i+1}"] = str(translation) #The translation is stored in the resulting dictionary
                
        elif frame != 0:
            translation = dna_seq[frame-1:].translate(table=table, stop_symbol=stop_symbol, to_stop=to_stop, gap=gap)
            protein_translation[f"frame_{frame}"] = str(translation  )

        return protein_translation
    
    def nucleotide_proportion(self, nucleotide, sequence=None, return_as_percentage=False):
        """This method will calculate the proportion of a given nucleotide in a given sequence, either as a percentage or in a scale to 1 (0,2 being 20% for example)
        """
        #The input is checked out
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its nucleotide_proportion, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence
                
        if nucleotide.upper() not in ["A", "C", "G", "T"]:
            raise ValueError("Invalid nucleotide. Must be A, C, G, or T.")
        
        if type(return_as_percentage) is not bool:
            raise ValueError("Invalid return_as_percentage value. Must be True or False.")
        
        #The program gets the number of the given nucleotide in the sequence and calculates the proportion
        nt_count = sequence.upper().count(nucleotide.upper())
        proportion = (nt_count / len(self.sequence))
        
        #The program returns the output in the desired form, always rounded by two decimals
        if return_as_percentage is not False:
            return f"{round(proportion * 100, 2)}%"
        else:
            return round(proportion, 2)
    
    def cg_content(self, sequence, return_as_percentage=False): 
        """This method calculates the CG content of the sequence, either as a percentage or in a scale to 1 (0,2 being 20% for example)
        """
        #The program checks the input
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence
                
        if type(return_as_percentage) is not bool:
            raise ValueError("Invalid return_as_percentage value. Must be True or False.")
        
        #The program gets the number of guanines and cytosines and calculate their proportion in the sequence
        g = sequence.lower().count("g")
        c = sequence.lower().count("c")
        proportion = ((g + c) / len(self.sequence))
        
        #The program returns the output in the desired form, always rounded by two decimals
        if return_as_percentage is not False:
            return f"{round(proportion * 100, 2)}%"
        else:
            return round(proportion, 2)
    
    def cpg_islands(self, sequence=None):
        """This method searches for CpG islands inside the inputed sequence, returning its position intervals locally (With the first nucleotide being position 1)
        """
    #According to 2013, Brenner's Encyclopedia of Genetics, a CpG island must have a length >= 200, C+G content >=50% and an observed/expected CpG ratio >= 0.6
    #CpG ratio might be diluted in very long sequences, thus, the program will take windows of 300 nucleotides (Arbitrary size > 200, user may change it at will)
    #Which will slide 50 nucleotides between them (Arbitrary size, user may change it at will)

        #The program checks the input
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence
                
        slide = 50
        window_start = 0
        window_end = 300
        startlist = []
        endlist = []
        
        while window_end <= len(self.sequence): #For each subsequence of 300nt, the program will check if it meets the CpG island criteria
            subseq = self.sequence[window_start:window_end]
            observed_cpg = subseq.lower().count("cg")
            expected_cpg =(subseq.lower().count("g")*subseq.lower().count("c"))/len(subseq)
            
            if expected_cpg == 0: 
                ratio = 0 #In case there wold be no C or G in the subsequence (which is extremely unlikely), the program avoids the ZeroDivisionError
            else: 
                ratio = observed_cpg / expected_cpg #The ratio is calculated
                
            #As the window is already >200 nt, if C+G content is >= 0.5 and the ratio is >= 0.6, the window meets the criteria to be a CpG island     
            if ((subseq.lower().count("g")+subseq.lower().count("c"))/len(subseq))>=0.5 and ratio >= 0.6:
                startlist.append(window_start + 1)
                endlist.append(window_end)
                
            window_start += 50 #The window slides
            window_end += 50
            
        if not startlist: #If no windows have been found, startlist will be empty
            return "No CpG island was found inside your sequence"
        
        else: #In case there are multiple consecutive windows that meet the criteria, the program will merge them into a single CpG island
            cleanstart = [startlist[0]]
            cleanend = []
            
            for i in range(1, len(startlist)):
                
                if startlist[i] - startlist[i-1] > slide:
                    cleanend.append(endlist[i-1])
                    cleanstart.append(startlist[i])
                         
            cleanend.append(endlist[-1])
            
            result = "CpG island found in position interval/s " #The result is given a more readable friendly format
            for i,cs in enumerate(cleanstart):
                result = result + str(cs) + "-" + str(cleanend[i]) + " "
            
            return result.strip()
    
    def inputed_enzymes_search(self, enzymes, file_output=None, file_name=None): 
        """This method BioPython to search inside the sequence for the targets of the restriction enzymes inputed by the user
The user may choose the output to be a .csv by setting the parameter file_output='csv' that will be created in the directory where this program is located, which can become rather helpfull the more enzymes are inputed (for a clearer disposition)
Even while outputing a csv, this method will still output the dictionary containing the result (which can still be printed in the terminal as well).
"""
        #The program checks the input
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence     
        
        csv_output = False
        if file_output is not None:
            
            if type(file_output) is str and file_output in ["csv"]:
                csv_output = True
                
            else:
                raise ValueError ("The only supported file output for inputed_enzymes_search() method is 'csv', and must be introduced as a string. Example file_output='csv'")
            
        
        if file_name is None:
            file_name = "Inputed_enzymes_search"
        
        if not enzymes:
            return "No enzymes were inputed. If you want a general search of the main restriction enzymes, try 'screen_enzymes()' method instead"
        
        else:
            
            try:
                enzyme_list = enzymes
                batch = RestrictionBatch(enzyme_list) #This line creates a RestrictionBatch object with the list of enzymes
                analysis = Analysis(batch, Seq(self.sequence)) #This line merges the batch of enzymes with the sequence, so that it may be analyzed
                enzyme_results = analysis.full()
                
                if csv_output is True: #If the user wants a csv output
                    fieldnames = ["Enzyme", "Positions_relative_to_the_sequence(first_nt=1)"] #The column names are set
                    
                    with open(f'{file_name}.csv', mode='w', newline='', encoding='utf-8-sig') as file:
                        writer = csv.writer(file)
                        writer.writerow(fieldnames) #The column names are written
                        
                        for enzyme, positions in enzyme_results.items():
                            position_string = ""
                            
                            #Positions for a single enzymes are merged into one cell
                            for i, p in enumerate(positions): 
                                if i == 0:
                                    position_string = str(p)
                                else:
                                    position_string = f"{position_string}, {p}"
                                    
                            writer.writerow([enzyme, position_string])
                
                return enzyme_results #The result consists on a dictionary with the enzymes as keys and their cutting positions as values
            
            except ValueError: return "An error happened while searching for your enzymes. Try checking if all your enzymes are separated by spaces and well written"
    
    def screen_enzymes(self, sequence=None): 
        """This method uses BioPython to search for 8 default restriction enzyme targets commonly used in molecular biology inside the inputed sequence
As the output is just the analysis of 8 restriction enzymes, the program only allows terminal output, unlike inputed_enzymes_search(), which also has csv output implemented
"""
        #The program checks the input
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence
        
        batch = RestrictionBatch(["EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "PstI", "KpnI", "NotI"])
        analysis = Analysis(batch, Seq(self.sequence))
        results = analysis.full()
        
        return results
    
    def search_seq(self, sequence=None, fragment=None): 
        """This method uses regex to search inside the main sequence, for the fragment inputed by the user(-s), returning its positions
"""
        #The program checks the input
        if sequence is None:
            if self.sequence is None:
                raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
            else:
                sequence = self.sequence
                
        if fragment is None:
            return "No fragment to search inside the main sequence was introduced"
        
        else:
            counter = 1
            matchdict = {}
            
            for match in re.finditer(fragment, self.sequence, re.IGNORECASE): #Iterating re.finditer() allows to create a dictionary with all matches
                matchdict[f"match{counter}"] = f"Positions {match.start()+1} to {match.end()}" #Result is given a readable format
                counter += 1
                
            if matchdict:
                return matchdict
            
            else:
                return "No matches were found inside your sequence"
            
    def interval_mutator(self, transcript=None, start_cdna=None, end_cdna=None, genomic_coordinates=None, include_flanking=True, strand=None, intron_HGVS_limit=15, file_output=None, file_name=None): 
    #intron_HGVS_limit was set to 15 because, usually, normal functioning introns wont be smaller than 30 nucleotides.
    #If file_output is set to None, the method will only output the dictionary. 
    #If set to "csv" or "fasta", the program will output said dictionary alongside the selected file inside the local directory
        """This method creates a dictionary with all possible SNPs within your sequence.
        
Transcript, start_cdna, end_cdna and genomic coordinates are essential for HGVS nomenclature. Thus, if any of those are absent, both in this method parameters and in its class atributes, variants will be named in a GENERIC way

Flanking inclusion may be turned off by writing 'include_flanking = False' inside the parenthesis, excluding flanks completely from results. In the other hand, if no flanks had been set at argparse, include_flanking will always be false

Introns will be written in HGVS nomenclature up to 15nt upstream and downstream, nucleotides between those limits will be named as 'Unknown intronic'. 
The amount of intronic nucleotides written with HGVS nomenclature can be changed by writing intron_HGVS_limit=number

The user may choose the output to be a .csv and or a .fasta by setting the parameter file_output='csv' or 'fasta' or ['csv', 'fasta']. The file/s will be created in the local directory
Setting the output to a .csv, a .fasta or both is recomended, as longer intervals will create huge outputs, which are difficult to read in the terminal
Even while outputing a csv, this method will still output the dictionary containing the result (which can be stored as a variable for further use or printed in the terminal as well).
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
            with open(f'{file_name}.csv', mode='w', newline='', encoding='utf-8') as file:
                writer = csv.writer(file)
                writer.writerow(fieldnames)
                for variant, instance_sequence in mutated_sequences.items():
                        writer.writerow([variant, instance_sequence])
                        
        if fasta_output is True: #If the user wants a fasta output
            with open(f"{file_name}.fasta", "w", encoding="utf-8") as file:
                for variant, instance_sequence in mutated_sequences.items():
                    file.write(f">{variant}\n")
                    file.write(f"{instance_sequence}\n")

        #The program returns the dictionary with the mutated sequences and the wild-type sequence, which might not be the same
        #as the inputed one if the flanks have been excluded
        return mutated_sequences
    
    def hexplorer_splicing_analysis(self, hex_score_matrix, reference_sequence=None, mutated_sequences=None, file_output=None, file_name=None):
        #If file_output is set to None, the method will only output the dictionary. 
        #If set to "csv" or "fasta", the program will output said dictionary alongside the selected file inside the local file
        """This method will allow to calculate the ΔHZEI(Delta score) from a HEXplorer analysis of splicing performed on every possible SNP of the inputed sequence (interval_mutator() dependant)        
This ΔHZEI can either be positive or negative. The more positive, the more likely is for the mutated sequence to have generated new enhancer sites around the mutation
The more negative it is, the more likely is for the mutated sequence to have generated new silencer sites around the mutation

If you want to study the whole inputed coordinates, INCLUDE A FLANKING OF AT LEAST 11(RECOMENDED). If not, the first and last 11 instances will not be reliable
Every flanking possition will be cut off the results, as its only use in here is to offer enough context for the interval to be correctly analyzed.

The mutated sequence/s must bue SNPs of the reference sequence, and thus have the same length, this method should not be used in in/del mutations.
If a mutated sequence is added as a single string instead of a dictionary (variant: sequence), the program will just return the ΔHZEI as a single value in the terminal. (.csv output will be disabled)

The score matrix for the hexamer scores is present in a dictionary form inside a different module called hex_matrix

The user may choose the output to be a .csv by setting the parameter file_output='csv'. The file will be created in the local directory
Setting the output to a .csv is recomended, as longer intervals will create huge outputs, which are difficult to read in the terminal
Even while outputing a .csv, this method will still output the dictionary containing the result (which can be stored as a variable for further use or printed in the terminal as well).
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
            with open(f'{file_name}.csv', mode='w', newline='', encoding='utf-8-sig') as file:
                file.write("sep=,\n")
                writer = csv.writer(file)
                writer.writerow(fieldnames)
                for variant, hzei in hexplorer_results.items():
                        writer.writerow([variant, hzei])
        
        return hexplorer_results
    
    def clinvar_screener(self, mutations=None, refseq_id=None, file_output=None, file_name=None, mail="your_email@example.com"):
        """The program will access the NCBI database for a clinvar screening of the variants. In this analysis, all the typically useful clinvar information of each variant will be retrieved

The user may choose the output to be a .csv by setting the parameter file_output='csv'. The file will be created in the local directory
Setting the output to a .csv is recomended, as longer intervals will create huge outputs, which are difficult to read in the terminal
Even while outputing a .csv, this method will still output the dictionary containing the result (which can be stored as a variable for further use or printed in the terminal as well).
"""
        
        csv_output = False
        if file_output is not None: #The program checks if the user wants a csv output
            
            if type(file_output) is str and file_output in ["csv"]:
                csv_output = True
                
            else:
                raise ValueError ("The only supported file output for clinvar_screener() method is 'csv', and must be introduced as a string. Example file_output='csv'")
        
        if file_name is None:
            file_name = "clinvar_screener"

        if mutations is None:
            mutations = self.get_mutated_sequences()
        
        #The program checks the mutation input
        mutation_list = []    
        if type(mutations) is str: #The program assumes tue user has inputed a single mutation in string form
            if match := re.search(r"(^ENST\d+):(c\.\d+[ACGT]>[ACGT]$)", mutations):
                ensembl_id = match.group(1)
                mutation_list.append(match.group(2))
            elif match := re.search(r"(^NM_.+):(c\.\d+[ACGT]>[ACGT]$)", mutations):
                refseq_id = match.group(1)
                mutation_list.append(match.group(2))
                ensembl_id = None
            else:
                raise ValueError ("Inputed mutation for clinvar_screener() must be written in HGVS nomenclatura, either with EnsemblID (ENST00000357654:c.5193G>A) or RefSeqID (NM_007294.4:c.5193G>A)")
        
        elif type(mutations) is list: #The program assumes tue user has inputed a list of mutations and will check every one of them
            for mutation in mutations: 
                if match := re.search(r"(^ENST\d+):(c\.\d+[ACGT]>[ACGT]$)", mutation):
                    ensembl_id = match.group(1)
                    mutation_list.append(match.group(2))
                elif match := re.search(r"(^NM_.+):(c\.\d+[ACGT]>[ACGT]$)", mutation):
                    refseq_id = match.group(1)
                    mutation_list.append(match.group(2))
                    ensembl_id = None
                else:
                    raise ValueError ("In an inputed mutation list for clinvar_screener(), every variant in it must be written in HGVS nomenclatura, either with EnsemblID (ENST00000357654:c.5193G>A) or RefSeqID (NM_007294.4:c.5193G>A)")
        
        elif type(mutations) is dict: #The program assumes tue user has inputed a dictionary(variants:sequences) of mutations and will check every one of the variants
            for mutation in list(mutations.keys()):
                if match := re.search(r"(^ENST\d+):(c\.\d+[ACGT]>[ACGT]$)", mutation):
                    ensembl_id = match.group(1)
                    mutation_list.append(match.group(2))
                elif match := re.search(r"(^NM_.+):(c\.\d+[ACGT]>[ACGT]$)", mutation):
                    refseq_id = match.group(1)
                    mutation_list.append(match.group(2))
                    ensembl_id = None
                else:
                    raise ValueError ("In an inputed mutation dictionary(variants:sequences) for clinvar_screener(), every variant in it must be written in HGVS nomenclatura, either with EnsemblID (ENST00000357654:c.5193G>A) or RefSeqID (NM_007294.4:c.5193G>A). Sequences will not be used")

        #If the user has not inputed a refseq_id, but a valid one exists in the atributes, the program will take it
        if refseq_id is None and self.refseq_id is not None and re.search(r"^NM_.+",self.refseq_id):
            refseq_id = self.refseq_id
            ensembl_id = None

        #If no refseq_id is available (this will happen if the inputed interval was in cDNA coordinates), the program will
        #access ensembl API REST to retrieve it
        if ensembl_id is not None:
            url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?expand=1;mane=1"
            headers = {"Content-Type": "application/json"}
            response = ensembl_request(url, headers)
            data= response.json()
            refseq_id = data["MANE"][0]["refseq_match"]
                 
        Entrez.email = mail #NCBI requires an identificative mail to access
        clinvar_result = []
        
        for mutation in mutation_list:
            hgvs = f"{refseq_id}:{mutation}"
            query = f'"{hgvs}"[All Fields]' #A query to search for the variant is generated
            search_handle = Entrez.esearch(db="clinvar", term=query) #In here, the programm accesses the NCBI, searching for the query in clinvar database
            #With Entrez.read, Biopython interprets the result from Entrez.esearch, returning a dictionary with its identifier called 'IdList'
            search_result = Entrez.read(search_handle) 
            search_handle.close()

            if search_result["IdList"]: #If the identifier is indeed present (which will probably be), the program will start a search by using it
                clinvar_id = search_result["IdList"][0] #The ID is stored in a variable
                #A summary from the variants info in clinvar is claimed by introducing its Id in esummary()
                summary_handle = Entrez.esummary(db="clinvar", id=clinvar_id)
                #The summary is interpreted by Biopython. Validate was set to false to avoid validation error (Asked to do so by the error itself) 
                summary = Entrez.read(summary_handle, validate=False) 
                summary_handle.close()
            else:
                print("Variant not found in clinvar")
                continue
            
            #The information is extracted from the retrieved summary
            document = summary["DocumentSummarySet"]["DocumentSummary"][0]
            if len(document["molecular_consequence_list"]) >= 1: #molecular_consequence_list may be empty. This will prevent an error if it is indeed empty
                molecular_consequence = document["molecular_consequence_list"][0]
            else:
                molecular_consequence = ""
            
            if len(document["germline_classification"]["trait_set"]) >= 1: #trait_set might be empty. This will prevent an error if it is indeed empty
                associated_illness = (document["germline_classification"]["trait_set"][0]["trait_name"]).split(",")[0]
            else:
                associated_illness = ""
            
            #The resulting dictionary is redacted
            variant_info = {
                "ID": document["title"],
                "Number_of_submissions": len(document["supporting_submissions"]["scv"]) - 1,
                "review_status": document["germline_classification"]["review_status"],
                "description": document["germline_classification"]["description"],
                "mutation_type": document["obj_type"],
                "molecular_consequence": molecular_consequence,
                "Associated_illness": associated_illness
            }

            clinvar_result.append(variant_info)
            
        if csv_output is True:
            fieldnames = [
                "ID", "Number_of_submissions", "review_status", 
                "description", "mutation_type", "molecular_consequence", "Associated_illness"]
            
            with open(f'{file_name}.csv', mode='w', newline='', encoding='utf-8-sig') as file:
                writer = csv.DictWriter(file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(clinvar_result)
        
        return clinvar_result
    
    def mes_canonical_analysis(self, genomic_coordinates = None, chromosome = None, strand = None, species = "human", coordinates = None, start_cdna = None, end_cdna = None, file_output=None, file_name=None):
        """This method performs a MaxEntScan analysis on the exons within inputed exonic coorddinates. 
The analysis returns the percentage of variation in the splicing strength of the Acceptor and Donnor regions of each exon, for any posible SNP variant.

The user may choose the output to be a .csv by setting the parameter file_output='csv'. The file will be created in the local directory
Setting the output to a .csv is recomended, as longer intervals will create huge outputs, which are difficult to read in the terminal
Even while outputing a .csv, this method will still output the dictionary containing the result (which can be stored as a variable for further use or printed in the terminal as well).
"""
        csv_output = False #The program checks if the user wants a csv output
        if file_output is not None:
            
            if type(file_output) is str and file_output in ["csv"]:
                csv_output = True
                
            else:
                raise ValueError ("The only supported file output for mes_canonical_analysis() method is 'csv', and must be introduced as a string. Example file_output='csv'")
        
        if file_name is None:
            file_name = "mes_canonical_analysis"
            
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
            with open(f'{file_name}.csv', mode='w', newline='', encoding='utf-8-sig') as file:
                file.write("sep=,\n") #This assures that decimals will be correctly formated in every language
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