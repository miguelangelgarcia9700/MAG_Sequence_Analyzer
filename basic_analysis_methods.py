import re
from Bio.Restriction import RestrictionBatch, Analysis #This will be installed with biopython 
from Bio.Seq import Seq #This will be installed with biopython
import csv
from pathlib import Path

def translate_dna_to_protein(self, dna_sequence=None, frame=0, table=1, stop_symbol="*", to_stop=False, gap="-"):
    """Translate a DNA sequence into one or more protein sequences.

    Uses Biopython's translation functionality to translate the supplied DNA
    sequence or, if none is supplied, the sequence stored in this object.

    If frame=0, translates the sequence in all three forward reading frames.
    Otherwise, translates only the requested frame. This method does not
    automatically calculate the reverse complement.

    Parameters:
        dna_sequence (str | None): DNA sequence to translate. If None, uses
            self.sequence.
        frame (int): Reading frame to use. Use 0 to translate all three
            forward frames, or 1, 2, or 3 to translate the respective frame.
        table (int): Biopython genetic code table identifier. Defaults to 1,
            which is the standard genetic code (genomic-human).
        stop_symbol (str): Character used to represent stop codons. (default='*')
        to_stop (bool): If True, stop translation at the first in-frame stop
            codon.
        gap (str): Character representing gaps during translation. (default='-')

    Returns:
        dict[str, str]: Protein translations keyed by frame, for example
        'frame_1', 'frame_2', and 'frame_3'. If a single frame is requested,
        the dictionary contains only that frame. {frame: translated_sequence}

    Raises:
        ValueError: If the DNA sequence is not a string, or if one of the
            translation parameters has an invalid value.
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

def nucleotide_proportion(self, nucleotide=None, sequence=None, return_as_percentage=False):
    """Calculate the proportion of a nucleotide in a DNA sequence.

    Counts the occurrences of the requested nucleotide and calculates its
    proportion relative to the sequence length.

    Parameters:
        nucleotide (str): Nucleotide to count. Must be A, C, G, or T.
        sequence (str | None): Sequence to analyze. If None, uses
            self.sequence.
        return_as_percentage (bool): If True, returns a percentage string;
            otherwise, returns a proportion between 0 and 1. Defaults to False.

    Returns:
        float | str: The proportion rounded to two decimal places, or a
            percentage string such as '20.0%'.

    Raises:
        ValueError: If the nucleotide is not supported or
            return_as_percentage is not a boolean.
    """
    #The input is checked out
    if sequence is None:
        if self.sequence is None:
            raise ValueError ("A sequence is required to calculate its nucleotide_proportion, either in this method parameters, or in the class atributes")
        else:
            sequence = self.sequence
            
    if nucleotide is None:
        raise ValueError ("""No nucleotide to be checked has been inputed in nucleotide_proportion() method.
Please, input the nucleotide whose proportion you want to check""")
            
    if nucleotide.upper() not in ["A", "C", "G", "T"]:
        raise ValueError("Invalid nucleotide. Must be A, C, G, or T.")
    
    if type(return_as_percentage) is not bool:
        raise ValueError("Invalid return_as_percentage value. Must be True or False.")
    
    #The program gets the number of the given nucleotide in the sequence and calculates the proportion
    nt_count = sequence.upper().count(nucleotide.upper())
    proportion = (nt_count / len(sequence))
    
    #The program returns the output in the desired form, always rounded by two decimals
    if return_as_percentage is not False:
        return f"{round(proportion * 100, 2)}%"
    else:
        return round(proportion, 2)

def cg_content(self, sequence=None, return_as_percentage=False): 
    """Calculate the GC content of a DNA sequence.

    Calculates the proportion of guanine (G) and cytosine (C) nucleotides
    relative to the sequence length.

    Parameters:
        sequence (str | None): Sequence to analyze. If None, uses
            self.sequence.
        return_as_percentage (bool): If True, returns a percentage string;
            otherwise, returns a proportion between 0 and 1. Defaults to False.

    Returns:
        float | str: The GC proportion rounded to two decimal places, or a
            percentage string such as '50.0%'.

    Raises:
        ValueError: If return_as_percentage is not a boolean or no sequence
            is available.
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
    proportion = ((g + c) / len(sequence))
    
    #The program returns the output in the desired form, always rounded by two decimals
    if return_as_percentage is not False:
        return f"{round(proportion * 100, 2)}%"
    else:
        return round(proportion, 2)

def cpg_islands(self, sequence=None):
    """Search for potential CpG islands in a DNA sequence.

    Scans the sequence using overlapping windows of 300 nucleotides, moving
    the window by 50 nucleotides each time. A window is considered a candidate
    CpG island when its GC content is at least 50% and its observed-to-expected
    CpG ratio is at least 0.6. Consecutive qualifying windows are merged.

    Parameters:
        sequence (str | None): Sequence to analyze. If None, uses
            self.sequence.

    Returns:
        str: A message containing the detected intervals in local, 1-based
            coordinates with the first nucleotide of the sequence being position 1. 
            Or a message indicating that no CpG island was found.

    Raises:
        ValueError: If no sequence is available.

    Notes:
        The window size, step size, and thresholds are currently defined
        inside the method.
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
    
    while window_end <= len(sequence): #For each subsequence of 300nt, the program will check if it meets the CpG island criteria
        subseq = sequence[window_start:window_end]
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

def inputed_enzymes_search(self, sequence=None, enzymes=None, file_output=None, file_name=None, file_path="./results/inputed_enzymes_search"): 
    """Search inside a DNA sequence for restriction sites of user-selected enzymes.

    Uses Biopython's RestrictionBatch and Analysis classes to identify the
    cutting positions of the specified restriction enzymes. 
    
    Optionally writes the results to a CSV file inside a "results" directory
    created locally.

    Parameters:
        enzymes (list[str] | tuple[str, ...]): Restriction enzyme names
            recognized by Biopython.
        file_output (str | None): Set to 'csv' (file_output='csv') to write 
            the results to a CSV file. If None, no external file is written, but
            the result can still be stored in a variable or printed in the teminal.
        file_name (str | None): Set the name to the csv. If not set, default can be:
            'Inputed_enzymes_search' if self.transcript is empty
            "{self.transcript}_'Inputed_enzymes_search'" if self.transcript not empty
            "{self.transcript}_{self.coordinates}_'Inputed_enzymes_search'" if neither
            self.transcript nor self.coordinates are empty
        file_path(str): Route where the output file/s will be created.
            Default to "./results/inputed_enzymes_search"
            This are directories, the file name will be concatenated after

    Returns:
        dict | str: A dictionary mapping enzyme names to their cutting
            positions, or a message if no enzymes were found, or an error
            message if Biopython raises a ValueError during the search.
        csv file: A csv file inside a 'results' directory created locally
            if file_output was set to 'csv'

    Raises:
        ValueError: If file_output has an unsupported value. Or if no
            restriction enzymes have been inputed in the parameters.

    Notes:
        The method returns the results dictionary even when CSV output 
        is requested.
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
        
        if self.transcript is not None:
            file_name = f"{self.transcript}_Inputed_enzymes_search"
            
            if self.coordinates is not None:
                file_name = f"{self.transcript}_{self.coordinates}_Inputed_enzymes_search"
        
    elif file_name is not None and re.search(r'[\\/:*?"<>|]', file_name):
        raise ValueError("The 'file_name' parameter contains invalid characters for file paths (e.g. /, \\, :, *, ?, \", <, >, |).")
    
    if enzymes is None:
        return "No enzymes were inputed. If you want a general search of the main restriction enzymes, try 'screen_enzymes()' method instead"
    
    else:
        
        try:
            enzyme_list = enzymes
            batch = RestrictionBatch(enzyme_list) #This line creates a RestrictionBatch object with the list of enzymes
            analysis = Analysis(batch, Seq(sequence)) #This line merges the batch of enzymes with the sequence, so that it may be analyzed
            enzyme_results = analysis.full()
            
            if csv_output is True: #If the user wants a csv output
                fieldnames = ["Enzyme", "Positions_relative_to_the_sequence(first_nt=1)"] #The column names are set
                output_dir = Path(file_path) #The output path is defined
                output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
                file_concatenated_path = output_dir / f"{file_name}.csv" #file name is concatenated
                
                with open(file_concatenated_path, mode='w', newline='', encoding='utf-8-sig') as file:
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
    """Search for restriction sites of a predefined set of enzymes.

    Uses Biopython to identify cutting positions for eight commonly used
    restriction enzymes: EcoRI, BamHI, HindIII, XhoI, XbaI, PstI, KpnI,
    and NotI.

    Parameters:
        sequence (str | None): Sequence to analyze. If None, uses
            self.sequence.

    Returns:
        dict: A dictionary mapping each restriction enzyme to its cutting
            positions, as returned by Biopython.

    Raises:
        ValueError: If no sequence is available.

    Notes:
        Results are returned in the terminal but are not exported to a file,
        unlike in inputed_enzymes_search
    """
    #The program checks the input
    if sequence is None:
        if self.sequence is None:
            raise ValueError ("A sequence is required to calculate its cg_content, either in this method parameters, or in the class atributes")
        else:
            sequence = self.sequence
    
    batch = RestrictionBatch(["EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "PstI", "KpnI", "NotI"])
    analysis = Analysis(batch, Seq(sequence))
    results = analysis.full()
    
    return results

def search_seq(self, sequence=None, fragment=None): 
    """Search for a user-supplied pattern within a DNA sequence.

    Uses Python regular expressions to locate matches without distinguishing
    between uppercase and lowercase letters.

    Parameters:
        sequence (str | None): Sequence to search into. If None, uses
            self.sequence.
        fragment (str | None): Regular expression representing the pattern
            to search for.

    Returns:
        dict[str, str] | str: A dictionary containing each match and its
            local, 1-based start and end positions, with the first nucleotide 
            of the sequence being position 1. {match_n: 1-based_interval}

    Raises:
        ValueError: If no sequence is available.
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
        
        for match in re.finditer(fragment, sequence, re.IGNORECASE): #Iterating re.finditer() allows to create a dictionary with all matches
            matchdict[f"match_{counter}"] = f"Positions {match.start()+1} to {match.end()}" #Result is given a readable format
            counter += 1
            
        if matchdict:
            return matchdict
        
        else:
            return "No matches were found inside your sequence"