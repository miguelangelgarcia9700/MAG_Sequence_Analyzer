import sys
import re
from hex_matrix import hex_scores
from auxiliar_functions import inputing
from SearchInfo_class import SearchInfo
from Sequence_class import Sequence


def main():
    """Run the MAG_sequence_analyzer workflow.

    Retrieve user-provided arguments, obtain transcript information and its
    nucleotide sequence, initialize a Sequence object, and execute the
    available sequence analyses.

    Print selected results to the console and export analysis results to
    CSV and FASTA files where configured.

    Exit with an error message if invalid input raises a ValueError.
    
    Note:
        By default, every analysis will be executed. If you only want to
        execute several of them, just comment the ones you don't want
        
        By default, mes_canonical_analysis() method will be omited if
        the coordinates were inputed in a cDNA format (-c c.5172_5180)
        because mes_canonical_analysis() does not support that sort of
        input, only exonic (-c ex.18_21) or whole sequence (omitting -c
        or writting -c 'whole')
    """
    
    print("----/Executing your order, please wait/----")
    
    try: #This "try", alongside its "except" will make the exceptions more readable    
        #Variables are extracted from the inputed arguments through argparse   
        transcript, coordinates, flank, fragment, enzymes, include_introns, mail = inputing()
        
        #All the available information is inputed into a SearchInfo-class object, so that its sequence may be retrieved
        #The ONLY ESSENTIAL input for a SearchInfo-class object is the transcript.
        myinfo = SearchInfo(
            transcript = transcript, 
            coordinates = coordinates, 
            flank = flank, 
            fragment = fragment, 
            enzymes = enzymes, 
            include_introns = include_introns
            )
        
        #The sequence and any other useful information is inputed into a Sequence-class object, as this class methods are the 
        #ones to actually perform the sequence analysis.
        #The ONLY ESSENTIAL input to generate a Sequence-class object is the sequence, though some of the methods inside it
        #may have additional requierements. That's why its recomended to build Sequence-class objects alongside Search-info class ones,
        #as Search-Info class ones can provide Sequence-class objects with most of the information any method could need
        seq = Sequence(
            sequence = myinfo.sequence, 
            transcript = transcript, #The program could have taken the transcript from myinfo.transcript, but I thought it was safer to take it from the input
            start_cdna = myinfo.start_cdna, 
            end_cdna = myinfo.end_cdna, 
            flank = flank, #The program could have taken the flank from myinfo.flank, but I thought it was safer to take it from the input
            genomic_coordinates = myinfo.genomic_coordinates, 
            refseq_id = myinfo.refseq_id,
            chromosome = myinfo.chromosome,
            strand = myinfo.strand,
            coordinates = coordinates #The program could have taken the coordinates from myinfo.coordinates, but I thought it was safer to take it from the input
            )
        
        #----/OUTPUT/----
        
        print("=== SEQUENCE ===") #This kind of string prints are merely to make the output more readable
        print(myinfo.sequence)

        print("\n=== NUCLEOTIDE PROPORTION ===")
        print(f"The prooportion of Adenine inside the sequence is {seq.nucleotide_proportion('A', return_as_percentage=True)}")
        print(f"The prooportion of Cytosine inside the sequence is {seq.nucleotide_proportion('C', return_as_percentage=True)}")
        print(f"The prooportion of Guanine inside the sequence is {seq.nucleotide_proportion('G', return_as_percentage=True)}")
        print(f"The prooportion of Thymine inside the sequence is {seq.nucleotide_proportion('T', return_as_percentage=True)}")

        print("\n=== PROTEIN TRANSLATION ===")
        for frame, aminoacid_sequence in seq.translate_dna_to_protein().items():
            print(f"{frame}: {aminoacid_sequence}\n")

        print("\n=== CG CONTENT ===")
        print(seq.cg_content())
        
        print("\n=== CPG ISLANDS ===")
        print(seq.cpg_islands())
        
        print("\n=== INPUTED ENZYMES SEARCH ===")
        print(seq.inputed_enzymes_search(enzymes=enzymes, file_output="csv"))
        
        print("\n=== DEFAULT ENZYMES SEARCH ===")
        print(seq.screen_enzymes())
        
        print("\n=== INPUTED FRAGMENT SEARCH INSIDE MAIN SEQUENCE ===")
        print(seq.search_seq(fragment=fragment))
        
        #print("\n=== ALL POSSIBLE SNPs INSIDE MAIN SEQUENCE ===")
        seq.interval_mutator(file_output=["csv", "fasta"]) 

        #Uncoment the titles if you want to print the results in the terminal
        #print("\n=== HEXPLORER SPLICING ANALYSIS ===")
        seq.hexplorer_splicing_analysis(hex_score_matrix = hex_scores, file_output="csv")

        #print("\n=== CLINVAR VARIANTS SCREENING ===")
        seq.clinvar_screener(file_output="csv", mail=mail)

        #mes_canonical_analysis will not be executed if inputed coordinates are in cDNA form, as that would lead to an error
        if not re.search(r"^c\.(\d+)(?:_(\d+))?$", coordinates, re.IGNORECASE):
            #print("\n=== MaxEntScan(MES) SPLICING ANALYSIS ===")
            seq.mes_canonical_analysis(file_output="csv")

    except ValueError as e:
        sys.exit(f"[ERROR] {e}")
        
    except RuntimeError as e:
        sys.exit(f"[ERROR] Operation failed: {e}")

if __name__ == "__main__":
    main()
    