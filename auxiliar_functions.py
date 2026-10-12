import argparse
import re
import requests
import time


def inputing():
    """Inputing function manages the data input through the library argparse
Information about any of the inputs is present inside its respective parser and can be retrieved by inputing in the terminal:
'python (module_name).py --help"""

    parser = argparse.ArgumentParser(description="Map cDNA coordinates to genomic sequences and analyze them.")
    parser.add_argument("-t", "--transcript", required=True, help="Ensembl Transcript ID ('ENST00000357654' or 'ENST00000357654.9' for example)")
    parser.add_argument("-c", "--coordinates", default="whole", help="Coordinates can be introduced in cDNA form (c.5172 or c.5172_5180) or in exonic form (ex.18 or ex.18_20) or not write -c at all, which will take the whole transcript by default")
    parser.add_argument("-f", "--flank", type=int, default=0, help="Flanking region in base pairs, before and after the main sequence (default: 0)")
    parser.add_argument("-s", "--fragment", default=None, help="Sequence fragment to be searched inside main sequence") 
    parser.add_argument("-e", "--enzymes", default=None, nargs="+", help="Restriction enzymes to be searched inside main sequence separated by spaces")
    parser.add_argument("-i", "--include_introns", action="store_true", help="If there are introns inside the inputed coordinates and you want to include them, write -i, else, if you only want cDNA, don't. Example: BRCA1:c.5152(end of ex 17)_5153(beggining ex 18) may include the intron inbetween")
    parser.add_argument("-m", "--mail", default="your_email@example.com", help="NCBI API REST to access Clinvar information needs a personal email of identification. Though, leaving it as default('your_email@example.com') by not writting -m at all, for some reason, works")
    args = parser.parse_args()
    
    #Before checking the info, input is stored in variables
    transcript = args.transcript.strip()
    coordinates = args.coordinates.strip()
    flank = args.flank
    fragment = args.fragment()
    enzymes = args.enzymes
    include_introns = args.include_introns
    mail = args.mail
    
    #Now the input passes a superficial check
    if match := re.match(r"^(ENST[\d]+)(\.[\d])?$", transcript):
        transcript = match.group(1)
    else:
        raise ValueError("Invalid Ensembl Transcript ID. Valid example = 'ENST00000357654' or 'ENST00000357654.9'")
    
    if fragment and not re.match(r"^[ACGT]+$", fragment, re.IGNORECASE):
        raise ValueError("Invalid sequence fragment. Only A, C, G and T are allowed")
    
    if enzymes: #This specific check is not very thorough, but can filter some incorrect inputs
        for enzyme in enzymes:
            if not re.match(r"^[A-Za-z0-9]+$", enzyme):
                raise ValueError(f"Invalid restriction enzyme name: {enzyme}. Only alphanumeric characters are allowed")
            
    if flank < 0:
        raise ValueError("Flanking region must be 0 or higher")
    
    #Keeping this order of return is crucial, so that variables are named correctly in the main() function
    return transcript, coordinates, flank, fragment, enzymes, include_introns, mail

def coordinate_extractor(coordinates):
    """This function checks if the coordinates are in a valid format through regex and separates start and end positions."""
    
    if match:= re.search(r"^c\.(\d+)(?:_(\d+))?$", coordinates, re.IGNORECASE): #Have coordinates been inputed in a cDNA format?
        
        coordtype = "cDNA"
        start = int(match.group(1))
        
        if match.group(2): #User may have inputed just one coordinate instead of an interval
            end = int(match.group(2)) 
        else: 
            end = start
        
        if end < start: 
            raise ValueError("Invalid cDNA coordinates. The end position can not be lower than the start position")
        
        if start < 1:
            raise ValueError("Invalid start coordinate. The start position can not be lower than 1")
        
        return start, end, coordtype
    
    elif match:= re.search(r"^ex\.(\d+)(?:_(\d+))?$", coordinates, re.IGNORECASE): #Have coordinates been inputed in a exonic format?
        
        coordtype = "exon"
        start = int(match.group(1))
        
        if match.group(2): #User may have inputed just one exon instead of an interval
            end = int(match.group(2))
        else: 
            end = start
        
        if end < start: 
            raise ValueError("Invalid exonic coordinates. The end position can not be lower than the start position")
        if start < 1:
            raise ValueError("Invalid start coordinate. The start position can not be lower than 1")
        
        return start, end, coordtype
    
    elif coordinates == "whole": #coordinates is set as "whole" by default if no coordinate has been inputed
        coordtype = "whole transcript"
        
        return None, None, coordtype
    
    #Should the input not be in valid HGVS format, the program will stop and ask the user to check the input
    else: raise ValueError("Invalid coordinate format. Valid cdna coordinates examples = c.5172 or c.5172_5180, Valid exonic examples = ex.18 or ex.18_20 and for Whole transcript don't write -c at all")

def ensembl_request(url, headers, max_retries=3):
    """This is an auxiliar function that, in case that ensembl is unresponsive, will retry another 3 times before it gives up and returns an error.
As an auxiliar function, this is NOT meant to be used independently"""

    for attempt in range(max_retries + 1):
        try:
            response = requests.get(url, headers=headers, timeout=30) #The program attempts to connect to the API REST

            if response.ok:
                return response

            #If Ensembl returns this responses, Ensembl servers are failing, and thus the program retries
            if response.status_code in [429, 500, 502, 503, 504]: 
                wait_time = 2 ** attempt
                print(
                    f"[WARN] Ensembl returned HTTP {response.status_code}. "
                    f"Retrying in {wait_time}s..."
                )
                time.sleep(wait_time)
                continue

            return response

        except requests.exceptions.RequestException as e:
            if attempt == max_retries: #If The program has retried the max number of retries permited, it returns the connection error
                raise RuntimeError(
                    f"Ensembl request failed after {max_retries + 1} attempts: {e}"
                ) from e

            wait_time = 2 ** attempt
            print(
                f"[WARN] Connection error with Ensembl: {e}. "
                f"Retrying in {wait_time}s..."
            )
            time.sleep(wait_time)

    raise RuntimeError("Ensembl request failed.")

def get_sequence(genomic_coords, chromosome, strand, flank = 0, include_introns = False, species = "human"): 
    """This method uses the Ensembl REST API, as well as some of the previously acquired information the get_info() method from a SearchInfo object, to extract the sequence
This function is used as an auxiliar function inside get_sequence() method in SearchInfo class and in mes_canonical_analysis() method in Sequence class.
Even while being an auxiliar function, with the correct input, it can be used independently"""
    
    #The program checks that the input is correct
    if type(genomic_coords) is not dict or set(genomic_coords.keys()) != {"start", "end"}:
        raise ValueError ("genomic_coords must be a dictionary with a key being 'start' containing every exon starting coordinate and another key called 'end' containing every exon ending coordinate")
    
    if type(flank) is not int:
        raise ValueError ("flank must be a positive integer or 0. Not writing a flank will set it to 0 by default")
    
    elif flank < 0:
        raise ValueError ("flank must be a positive integer or 0. Not writing a flank will set it to 0 by default")
    
    if type(include_introns) is not bool:
        raise ValueError ("include_introns is a boolean parameter, It can only be either True or False")
    
    #The program will take the bordering coordinates and apply the flanking, so that it can retrieve the whole sequence in the next steps
    startcoords = list(genomic_coords["start"])
    endcoords = list(genomic_coords["end"])  
    startcoords[startcoords.index(min(startcoords))] -= flank #Border coordinates are modified by flanking regions
    endcoords[endcoords.index(max(endcoords))] += flank  
    
    seq_with_introns = "" #A variable to contain the whole sequence is created
    url = f"https://rest.ensembl.org/sequence/region/{species}/{chromosome}:{min(startcoords)}..{max(endcoords)}:{strand}?"
    headers = {"Content-Type": "application/json"}
    response = ensembl_request(url, headers)
    seq_with_introns = seq_with_introns + response.json()["seq"]
    seq_no_introns = "" #In case the user didn't choose to include introns, this variable will be the container to the resulting sequence
    
    flanking_modified_coordinates = {"start": startcoords, "end": endcoords}
    #The next step will use an auxiliary function to calculate the length of the exons and the introns, out from the genomic coordinates
    #This will be useful in the next step, as the program will now know which parts of the sequence to cut if introns are to be excluded
    #exon_intron_boundaries() function is defined ahead in the code
    exon_intron_difference = exon_intron_boundaries(include_flanking = False, flank = flank, genomic_coordinates = flanking_modified_coordinates, strand = strand)
    exon_intron_index = 0 #This variable will be used as index for exon_intron_difference
    
    if include_introns is True: #No cutting is necessary
        
        return seq_with_introns
    
    else:
        
        for i,nt in enumerate(seq_with_introns):
            
            if exon_intron_index%2 == 0: #If the exon intron index is even, the loop is currently on an exon
                seq_no_introns = seq_no_introns + nt #Thus, the nucleotide is added to the resulting sequence
                
                if i + 1 == exon_intron_difference[exon_intron_index]: #If the exon has been completely traversed, program moves to next index
                    exon_intron_index += 1
                    
            elif exon_intron_index%2 != 0: #If the exon intron index is odd, the loop is currently on an intron (no nt will be added)
                
                if exon_intron_index < len(exon_intron_difference): #If the index is actually on exon_intron_difference (If not, the sequence has already been completely traversed)
                    
                    if i + 1 == exon_intron_difference[exon_intron_index]: #If the intron has been completely traversed, program moves to next index
                        exon_intron_index += 1
                        
        return seq_no_introns

def snp_mutator(sequence, position, new_base):
    """Mutates a single base at a designed position within a DNA sequence. (Starting from 1)
This function takes a sequence, a position inside it and the new desired nucleotide, and will return the mutated version of the original sequence
In this program, this function is used as an auxiliary in interval_mutator() method in Sequence-class objects. But it can very well be used independently"""
    
    if not (1 <= position <= len(sequence)): #If the given position is not inside the sequence, the program launches an error
        raise ValueError("Mutation position out of sequence bounds.")
      
    sequence = sequence.upper() #The program assures uniformity in its input
    new_base = new_base.upper()
    
    if new_base not in ["A", "C", "G", "T"]: #The program checks if the desired new_base is valid
        raise ValueError("Invalid replacement base. Must be A, C, G, or T.") 
       
    if sequence[position - 1] == new_base: #The program checks if the new base is identical to the one already present
        raise ValueError("The new base is identical to the reference base.")
    
    return sequence[:position - 1] + new_base + sequence[position:] #The mutated sequence is returned


def exon_intron_boundaries(include_flanking, flank, genomic_coordinates, strand):
    """This function calculates the length of exons and introns inside the sequence. 
This is being done in an accumulative way. For example, f two exons separated by an intron, all measured 50nt, this returns [50,100,150]
This function is NOT meant to be used independently (though it technically could), but instead is currently being used in interval_mutator() method in Sequence class
    """
    exon_intron_difference = [] #This variable will contain the resulting list
    accumulative = 0 #This variable will be used to generate a list with exon and intron endings locally inside the sequence        
    
    if include_flanking:
        accumulative += flank #5'flank will bias the genomic coordinates, thus, we add it to the accumulative variable to unbias it
        
    for i,(startcoord, endcoord) in enumerate(zip(genomic_coordinates["start"], genomic_coordinates["end"])):
        accumulative += endcoord-startcoord+1 #The program measures the exon
        exon_intron_difference.append(accumulative) #In this program, exons will always go first, as inputed coordinates always starts on an exon (due to them being cDNA, exonic or whole transcript)
        
        if i < len(genomic_coordinates["start"]) - 1: #If there are more exons ahead, thus, there must be an intron, else, iterations will end
            
            if strand == 1: #If the transcript is in the forward strand, the intron's length will be the start of the next exon minus the end of the current one
                accumulative += genomic_coordinates["start"][i+1] - endcoord - 1
            
            elif strand == -1: #If it's in the reverse strand, the intron's length will be the start of the current exon minus the end of the next one
                accumulative += startcoord - genomic_coordinates["end"][i+1] - 1
                
            exon_intron_difference.append(accumulative) #Program adds the intron next and continues iterating 
    
    return  exon_intron_difference