import sys
from hex_matrix import hex_scores
from auxiliar_functions import inputing
from SearchInfo_class import SearchInfo
from Sequence_class import Sequence


def main():
    """Main function integrates the rest of functions and classes, and will print every sort of genomic result by default as demonstration"""
    
    print("----/Executing your order, please wait/----")
    
    try: #This "try", alongside its "except" will make the exceptions more readable       
        transcript, coordinates, flank, fragment, enzymes, include_introns, mail = inputing() #Variables are extracted from the inputed arguments through argparse
        
        #All the available information is inputed into a SearchInfo-class object, so that its sequence may be retrieved
        #The ONLY ESSENTIAL input for a SearchInfo-class object is the transcript. As, if the program lacks any other information,
        #it will retrieve the sequence of the whole transcript (without introns) as default.
        #The rest of the information are modifiers that will alter the final sequence output (including flanks, introns...)
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

        print("\n=== PROTEIN TRANSLATION ===")
        print(seq.translate_dna_to_protein())

        print("\n=== CG CONTENT ===")
        print(seq.cg_content())
        
        print("\n=== CPG ISLANDS ===")
        print(seq.cpg_islands())
        
        print("\n=== INPUTED ENZYMES SEARCH ===")
        print(seq.inputed_enzymes_search(enzymes=enzymes, file_output="csv", file_name="Inputed_enzymes_search"))
        
        print("\n=== DEFAULT ENZYMES SEARCH ===")
        print(seq.screen_enzymes())
        
        print("\n=== INPUTED FRAGMENT SEARCH INSIDE MAIN SEQUENCE ===")
        print(seq.search_seq(fragment))
        
        # print("\n=== ALL POSSIBLE SNPs INSIDE MAIN SEQUENCE ===")
        seq.interval_mutator(file_output=["csv", "fasta"], file_name="Interval_mutator")
        # print(seq.interval_mutator())
        
        # print("\n=== HEXPLORER SPLICING ANALYSIS ===")
        seq.hexplorer_splicing_analysis(hex_score_matrix = hex_scores, file_output="csv", file_name="HEXplorer_splicing_analysis")
        # print(seq.hexplorer_splicing_analysis(hex_score_matrix = hex_scores))
        
        # print("\n=== CLINVAR VARIANTS SCREENING ===")
        seq.clinvar_screener(file_output="csv", file_name="Clinvar_variants_screening", mail=mail)
        # print(seq.clinvar_screener())
        
        # print("\n=== MaxEntScan(MES) SPLICING ANALYSIS ===")
        seq.mes_canonical_analysis(file_output="csv", file_name="MaxEntScan(MES)_splicing_analysis")
        # print(seq.mes_canonical_analysis(chromosome = myinfo.chromosome, strand = myinfo.strand, coordinates = myinfo.coordinates))

    except ValueError as e:
        sys.exit(f"[ERROR] {e}")

   



if __name__ == "__main__":
    main()
    