import re
from Bio import Entrez
import csv
from auxiliar_functions import ensembl_request
from pathlib import Path     
        
def clinvar_screener(self, mutations=None, refseq_id=None, file_output=None, file_name=None, file_path="./results/clinvar_screener", mail="your_email@example.com"):
    """Retrieve ClinVar information for variants specified in HGVS notation.

    Queries the NCBI ClinVar database for each supplied variant and retrieves
    available information about its classification, clinical significance,
    molecular consequence, associated condition, and supporting submissions.
    
    This screener can only work on one transcript at a time, so please, do not
    input multiple ensemble_IDs or refseq_IDs in the mutations parameter.
    (If executed from main, normally the program would use the transcript 
    inputed by -t)
    
    Optionally writes the results to a CSV file inside a "results" directory
    created locally.

    Parameters:
        mutations (str | list[str] | dict[str, str] | None): A single variant,
            a list of variants, or a dictionary whose keys are variant labels
            (like the one obtained from interval_mutator).
            Variants must use the supported transcript-prefixed HGVS format,
            such as 'ENST00000357654:c.5193G>A' or 'NM_007294.4:c.5193G>A'. 
            If None, the method uses self.mutated_sequences which is obtained 
            by default through interval_mutator method.
        refseq_id (str | None): RefSeq transcript identifier used for ClinVar
            queries. If None, the program uses self.refseq_id. If self.refseq_id
            is empty too, the program will access Ensembl API REST to get it.
        file_output (str | None): Optional csv output format
            May be introduced as a string (file_output='csv')
        file_name (str | None): Set the name to the csv. If not set, default can be:
            'Clinvar_screener' if self.transcript is empty
            "{self.transcript}_'Clinvar_screener'" if self.transcript not empty
            "{self.transcript}_{self.coordinates}_'Clinvar_screener'" if neither 
            self.transcript nor self.coordinates are empty
        file_path(str): Route where the output file will be created.
            Default to "./results/clinvar_screener"
            This are directories, the file name will be concatenated after
        mail (str): Email address supplied to NCBI through Biopython Entrez.
            (NCBI requires a mail to access its database. Though, the generic
            one set as default ("your_email@example.com") works for some reason)

    Returns:
        list[dict]: A list containing a dictionary per variant. Each dictionary
        contains its respective variant's ClinVar information. 
        Returns an empty list if no variants are found.

    Raises:
        ValueError: If the mutation input or output format is invalid.

    Notes:
        This method makes requests to NCBI. Results depend on the records
        available in ClinVar at the time of the query. Variants not found
        in ClinVar are reported to the terminal and skipped.
        
        The method returns the results dictionary even when CSV output
        is requested.
    """
    csv_output = False
    if file_output is not None: #The program checks if the user wants a csv output
        
        if type(file_output) is str and file_output in ["csv"]:
            csv_output = True
            
        else:
            raise ValueError ("The only supported file output for clinvar_screener() method is 'csv', and must be introduced as a string. Example file_output='csv'")
    
    if file_name is None:
        file_name = "Clinvar_screener"
        
        if self.transcript is not None:
            file_name = f"{self.transcript}_Clinvar_screener"
            
            if self.coordinates is not None:
                file_name = f"{self.transcript}_{self.coordinates}_Clinvar_screener"
        
    elif file_name is not None and re.search(r'[\\/:*?"<>|]', file_name):
        raise ValueError("The 'file_name' parameter contains invalid characters for file paths (e.g. /, \\, :, *, ?, \", <, >, |).")

    if mutations is None:
        mutations = []
        mutations_raw = self.get_mutated_sequences()
        for mutation in list(mutations_raw.keys()): #This will exclude flanking and any other sort of non-accepted variant
            if re.search(r"^ENST\d+:c\.\d+[ACGT]>[ACGT]$", mutation): 
                mutations.append(mutation)
            
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
    
    else:
        raise ValueError ("Invalid format for 'mutations' parameter in clinvar_screener, it should be either a string, a list or a dictionary")

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
            print(f"Variant {hgvs} not found in clinvar")
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
        
        output_dir = Path(file_path) #The output path is defined
        output_dir.mkdir(parents=True, exist_ok=True) #If the directories don't exist yet, the program creates them
        file_concatenated_path = output_dir / f"{file_name}.csv" #file name is concatenated
        
        with open(file_concatenated_path, mode='w', newline='', encoding='utf-8-sig') as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(clinvar_result)
    
    return clinvar_result