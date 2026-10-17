# MAG Sequence Analyzer
### Made by: Miguel Ángel García de Castro
### Contact mail: miguelangelgarcia9700@gmail.com
### Video tutorial of usage: here I'll put the link to the video

# MAG_Sequence_Analyzer

**MAG_Sequence_Analyzer** is a Python-based tool for transcript sequence analysis, variant exploration, splice-site evaluation, and clinical variant screening using Ensembl and NCBI ClinVar resources.

The program is designed to help researchers explore the potential molecular consequences of sequence variants, with a particular focus on nucleotide substitutions, splice-site alterations, and their possible clinical relevance.

The program counts with the following features to further automatize the analysis of the desired sequence:
    -Automatic information retrieval about the inputed transcript (or transcript fragment) from Ensembl database
    -Automatic DNA sequence retrieval of the inputed transcript (or transcript fragment) from Ensembl database

Between the different analysis the user may choose from, it can be found:
    -Translation of DNA to protein in any of the three possible frames (or the three at once, storing them in a dictionary).
    -Calculus of nucleotide proportion inside the inputed sequence.
    -Calculus of CG content inside the inputed sequence.
    -Location of CpG islands (if present).
    -Search of targets for the restriction enzymes inputed inside the inputed sequence (csv output available).
    -Search of targets for eight frequently used restriction enzymes (EcoRI, BamHI, HindIII, XhoI, XbaI, PstI, KpnI and NotI).
    -Search for a sequence fragment inside the main inputed sequence.
    -Automatic generation of every possible single nucleotide plymorphism of the inputed sequence (csv and fasta output available).
    -HEXplorer splicing analysis of the inputed sequence (csv output available).
    -Clinvar information retrieval of every possible single nucleotide plymorphism of the inputed sequence (csv output available).
    -MaxEntScan splicing analysis of every exon present in the inputed sequence.

It is intended for users with a background in biology or genetics. Thus, the program is user-friendly enough, so that even those who may have limited experience in bioinformatics may use it without any trouble.
The program combines several sequence-analysis methods in a single workflow and generates files that can be used for further analysis.

> **Important:** MAG_Sequence_Analyzer is a research and exploratory analysis tool. Its predictions and database results should not be interpreted as definitive evidence of pathogenicity or as a substitute for validated clinical genetic testing.

---

## Table of Contents

- [Overview](#overview)
- [Main Features](#main-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Running the Program](#running-the-program)
- [Input Parameters](#input-parameters)
- [Analysis Workflow](#analysis-workflow)
- [Output Files](#output-files)
- [Biological Interpretation](#biological-interpretation)
- [External Resources and Internet Connection](#external-resources-and-internet-connection)
- [Limitations and Considerations](#limitations-and-considerations)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements and References](#acknowledgements-and-references)

---

## Overview

MAG_Sequence_Analyzer integrates sequence retrieval, basic sequence characterization, variant generation, splice-related analyses, and clinical database screening into a single workflow.

The user provides a transcript identifier and selects the genomic or transcript region to investigate. The program retrieves the corresponding sequence information, performs the selected analyses, and either shows the results in the terminal, or saves them in files for subsequent inspection (only supported results may be saved in files, being csv or fasta).

The main objectives the program can achieve:
*(not necessary to perform all of them, only the ones the user may want)*

1. Retrieve transcript and sequence information from Ensembl.
2. Characterize the selected nucleotide sequence with several basic analysis.
3. Generate every possible single-nucleotide polymorphism across the selected region.
4. Evaluate potential effects on splicing using sequence-based prediction methods.
5. Search ClinVar for information associated with generated or user-provided variants.
6. Export results in formats suitable for further processing.

The analyses are organized into separate Python modules, making the project easier to maintain and extend.

---

## Main Features

### 1. Transcript and sequence retrieval

The program uses Ensembl to retrieve transcript information and the corresponding genomic sequence.

Depending on the selected input region, the analysis can be performed on a complete transcript or a specified region, indicated either in cDNA coordinates (c.5153_5260) or exonic coordinates (ex.3_9). The program also supports the use of flanking sequence and options related to the inclusion of intronic regions (as those are excluded by default).
*Further information in [Input Parameters](#input-parameters)*

### 2. Basic sequence analysis

The program provides several methods for characterizing a nucleotide sequence:

- **DNA translation:** translates nucleotide sequences into amino acid sequences in supported reading frames.
- **Nucleotide composition:** calculates the proportion of a specified nucleotide.
- **GC content:** calculates the fraction of guanine and cytosine nucleotides.
- **CpG island detection:** identifies candidate regions according to the thresholds implemented in the program.
- **Restriction enzyme analysis:** searches for recognition sites of user-selected enzymes and/or a predefined set of enzymes.
- **Sequence fragment search:** searches for a specified nucleotide fragment within the analyzed sequence.

The results can help characterize the sequence before investigating potential variant effects.

### 3. Systematic single-nucleotide polymorphism generation

The interval mutator generates alternative nucleotide substitutions across the analyzed region.

For each nucleotide position, the method considers the three alternative nucleotides that differ from the reference base. It also associates generated substitutions with positional information, including applicable cDNA and intronic annotations, naming each mutated sequence by its HGVS form (if possible).

The resulting collection of variants can be used as input for downstream analysis, like the HEXplorer, MaxEntScan and Clinvar ones.

This approach is useful for systematic exploration of a selected region. However, it does not mean that every possible variant is biologically relevant or clinically significant (this could be checked with the Clinvar screener).

### 4. Splicing analysis with HEXplorer

The program includes a HEXplorer-based analysis for evaluating changes in HEXplorer sequence scores (hzei) associated with nucleotide substitutions. 

The method compares the score of the reference sequence with the score of each corresponding mutated sequence. The difference between the scores is used to describe the predicted change associated with the substitution. A significant change towards the negative could mean the creation of new splicing silencing sites, whereas a significant change towards the possitive could indicate the pressence of newly created splicing enhancing sites.

These scores are computational predictions. They should not be interpreted as direct measurements of splicing efficiency or as definitive evidence that a variant alters splicing.

A single mutation only affects the hzei scores of 11 nucleotides around it. Thus, it should be taken in mind, that a variation in results between the wild type sequence and the mutated one, is only due to the affected nucleotide and other 5 nucleotides upsteram and downstream, and only there is where an enhancing or silencing region could be newly created, and not the whole sequence.
More information in the original article describing HEXplorer from Steffen Erkelenz et al. DOI: 10.1093/nar/gku736.

### 5. Canonical splice-site analysis with MaxEntScan

The program also implements a canonical splice-site analysis using MaxEntScan.

The method evaluates donor and acceptor splice-site sequences and examines selected nucleotide substitutions within the canonical region. The results include reference and alternative scores and the corresponding variation beween them, with a negative result being the corresponding percentual loss in strength of the analyzed splicing site and viceversa.

This analysis focuses on the canonical splice-site positions and should not be assumed to capture every possible mechanism of abnormal splicing.

### 6. ClinVar screening

The ClinVar screening module searches NCBI ClinVar for information associated with supported variant identifiers.

The module can process variants supplied directly or variants obtained from the program's own single-nucleotide polymorphism generation workflow (strongly recomended). It retrieves available database information, including fields such as the variant description, review status, molecular consequence, and associated conditions when present in the returned record.

**ClinVar screening operates on one transcript at a time.** To analyze several transcripts, create separate analyses for each transcript. (The program can automatically perform an analysis among any given interval within the given transcript, but cannot simultaneously analyze two different variants present in two different transcripts)

A ClinVar search result indicates that a database record was found and retrieved; it does not, by itself, establish that a variant is pathogenic or clinically relevant to a particular patient. This is performed through connection with NCBI Clinvar database

---

## Requirements

MAG_Sequence_Analyzer requires Python and the project's Python modules together in the same directory. Thus, it is recommended to download the whole MAG_Sequence_Analyzer folder from github and perform the analysis directly from said folder.

The code uses external libraries and online resources for sequence handling, splice prediction, and database queries. The exact dependency versions should be checked against the imports and tested environment of the project.

### Python environment

A dedicated virtual environment is recommended. This helps prevent conflicts with packages used by other projects, though it may not be strictly necessary.

The examples below use `venv`, which is included with standard Python installations.

### External Python libraries

The modules shared for this project use:

- **Biopython** (`Bio`) for biological sequence handling, translation, and restriction enzyme analysis.
- **maxentpy** for MaxEntScan-related calculations.
- **requests** for API calling (Ensembl)
- Python's standard library for file handling, regular expressions, CSV output, and other utility functions.

Additional dependencies may be required by the sequence-retrieval and supporting modules.

---

## Installation

### 1. Obtain the project

Download or clone the repository and open the project directory. 
*(Download the whole MAG_Sequence_Analyzer folder to assure all modules stay together)*

The project directory should contain `main.py` and the supporting Python modules described in the [Project Structure](#project-structure) section.

### 2. Create a virtual environment 
***(OPTIONAL but recomended for compatibility with future versions of python and libraries)***

From the project root, run:

```bash
python -m venv .venv
```

Activate the environment.

**Windows — Command Prompt:**

```bat
.venv\Scripts\activate.bat
```

**Windows — PowerShell:**

```powershell
.\.venv\Scripts\Activate.ps1
```

**macOS and Linux:**

```bash
source .venv/bin/activate
```

### 3. Install the dependencies

If the repository includes a `requirements.txt` file, install its dependencies with:

```bash
python -m pip install -r requirements.txt
```

### 4. Verify the project files

Make sure that all local modules are present in the expected directory. In particular, `main.py` imports other modules from the project, including `auxiliar_functions.py` and `hex_matrix.py`.

These are project files, not standard Python libraries. Installing third-party packages will not replace missing local modules.

---

## Project Structure

The project is organized into a main execution script, classes representing the transcript and sequence information, and modules containing specific analysis methods.

The following is the expected organization based on the modules used by the program:

```text
MAG_sequence_analyzer/
├── main.py
├── SearchInfo_class.py
├── Sequence_class.py
├── auxiliar_functions.py
├── basic_analysis_methods.py
├── interval_mutator_method.py
├── splicing_methods.py
├── clinvar_screening_method.py
├── hex_matrix.py
├── README.md
├── .gitignore
├── requirements.txt
└── LICENSE.txt

```

A `results/` directory is used for generated output files. Depending on the configured output paths, additional files or directories may be created during execution.

Python may also generate a `__pycache__/` directory containing cached bytecode files. These files are not required to run the project.

### Main files

| File | Purpose |
|---|---|
| `main.py` | Coordinates the workflow and invokes the analyses. |
| `SearchInfo_class.py` | Retrieves and organizes transcript, coordinate, and all the sequence information. |
| `Sequence_class.py` | Represents the analyzed sequence and exposes the analysis methods. |
| `auxiliar_functions.py` | Provides supporting functions for sequence retrieval, coordinate handling, mutation generation, and related tasks. |
| `basic_analysis_methods.py` | Implements basic sequence characterization, translation, restriction enzyme analysis, and fragment searching. (containing the funtioning code of some Sequence-class methods). |
| `interval_mutator_method.py` | Generates alternative single-nucleotide substitutions across the selected region. (containing the funtioning code of the Sequence-class method "interval_mutator()"). |
| `splicing_methods.py` | Implements HEXplorer-based and canonical splice-site analyses. (containing the funtioning code of the Sequence-class methods "hexplorer_splicing_analysis()" and "mes_canonical_analysis()")|
| `clinvar_screening_method.py` | Queries ClinVar for information associated with supported variants. (containing the funtioning code of the Sequence-class method "clinvar_screener()") |
| `hex_matrix.py` | Provides the score matrix used by the HEXplorer analysis in the "hexplorer_splicing_analysis()" method. A flanking of at least 11 is recommended|

The `main.py` script is the recommended entry point for running either the complete workflow, or specific analysis by commenting the undesired ones. 
Individual methods may also be called from outside of main() if their required objects and parameters are provided. (not recomended)

---

## Running the Program

### 1. Open the project

Open the project root in an IDE such as Visual Studio Code, or navigate to the directory using a terminal.

Ensure that the correct Python virtual environment is active.

### 2. Run the main script

From the project root, you will have to execute the main module, followed by the desired input.
Example with every kind of input (only the transcript is essential, the rest is optional):

```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11 -s GTGACCCCAGT -e NcoI NdeI SpeI SphI -i -m genericemail@gmail.com
```
*More information about input in [Input Parameters](#input-parameters)*


### 3. Wait for the analyses to finish

The program performs its workflow sequentially. Several analyses may require internet access to retrieve sequence information or query external databases.

Depending on the region length and number of variants generated, the analysis may take longer for larger regions.

### 4. Inspect the results

After execution, inspect both, the results printedd in the terminal, and the files generated in the configured output directory, normally under `./results/`.

The program supports CSV output for several analyses and FASTA output for the interval mutator, while other smaller analysis can only be printed in the terminal. The precise set of files depends on the selected parameters and the methods successfully completed.

---

## Input Parameters

The main script gathers the information required to construct the transcript and sequence objects and to configure the downstream analyses.

ONLY THE TRANSCRIPT IS ESSENTIAL, the rest of the input is optional, depending on the analysis you may want to perform

| Parameter | Purpose |
|---|---|
| Transcript | Identifies the transcript to analyze. | 
    -To input it just write -t followed by the ensembl identifier of the transcript while executing main:
```bash
python main.py -t ENST00000357654
```
| Parameter | Purpose |
|---|---|
| Coordinates | Specifies the region of interest, using a supported coordinate format. |
    -To input it just write -c followed by the cDNA (c.5172_5180) or exonic (ex.18_23) coordinate interval you want to analyze while executing main. Not writting -c at all or following it by 'whole' will trigger the analysis of the whole transcript:
```bash
python main.py -t ENST00000357654 -c ex.16_20
```
| Parameter | Purpose |
|---|---|
| Flank | Defines the amount of nucleotides flanking the main sequence |
    -To input it just write -f followed by the number of nucleotides you want to flank the interval with, while executing main:
```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11
```
| Parameter | Purpose |
|---|---|
| Fragment | Specifies a nucleotide sequence fragment to search for in the main sequence. |
    -To input it just write -s followed by the DNA sequence you want to search for inside the main sequence, while executing main:
```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11 -s GTGACCCCAGT
```
| Parameter | Purpose |
|---|---|
| Enzymes | Specifies restriction enzymes whose targets are to be searched inside the main sequence. |
    -To input it just write -e followed by the desired restriction enzymes, while executing main:
```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11 -s GTGACCCCAGT -e NcoI NdeI SpeI SphI
```
| Parameter | Purpose |
|---|---|
| Include introns | Controls whether intronic sequences inside the interval, are included in the sequence retrieval. |
    -To input it just write -i, while executing main:
```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11 -s GTGACCCCAGT -e NcoI NdeI SpeI SphI -i
```
| Parameter | Purpose |
|---|---|
| Email | Provides the email address, used for authentication in NCBI Entrez queries (for clinvar analysis method). |
    -To input it just write -m followed by the e-mail account you want to provide.
```bash
python main.py -t ENST00000357654 -c ex.16_20 -f 11 -s GTGACCCCAGT -e NcoI NdeI SpeI SphI -i -m genericemail@gmail.com
```

### Transcript identifier

Provide a transcript identifier compatible with the Ensembl lookup used by the program. Ensembl transcript identifiers commonly begin with `ENST`.

The identifier should refer to the intended transcript(ENST), not merely to the gene(ENSG). A single gene can have several transcripts with different exon structures and sequences. Thus, you should make sure your transcript identifier is preceded by ENST

### Coordinates and regions

The program supports several region-selection modes, including the complete transcript and exon-based selection. It also contains logic for cDNA coordinate inputs.

The precise coordinate syntax should be checked against the parser implemented in the current version.

Supported coordinate input::

- **Whole transcript:** analyze the full transcript region supported by the retrieval method.
- **Exon selection:** analyze one or more specified exons.
- **cDNA coordinates:** identify a region using positions in the transcript's coding DNA reference system.

For reproducible analyses, recording the transcript identifier and the exact coordinate expression used is recommended.

### Flanking sequence

Flanking sequence can provide context around a region of interest. It may be relevant to splice-site analysis and to the interpretation of variants near the selected interval.

A flank value of zero requests no additional flanking bases. The precise behavior of flanks depends on the analysis method and on whether the requested region includes intronic sequence.

### Restriction enzymes

The program accepts restriction enzyme information for its user-selected enzyme target search within the inputed transcript interval. It also performs a separate search for a predefined set of enzymes.

Use enzyme names recognized by the restriction enzyme library used by the program.

### Email address

The email address is passed to NCBI Entrez for ClinVar queries.

### Other inputs

The user may need to modify some method parameters in main, for example, to study nucleotide proportion, the user should indicate the nucleotide to analyze in the nucleotide_proportion() method parameter.
The next example would return the proportion of Adenine in a sequence stored in "seq" object, which is a Sequence-class object:
```python
seq.nucleotide_proportion('A', return_as_percentage=True)
```
Other example would be prompting for a file output in those methods who support it:
```python
seq.interval_mutator(file_output=["csv", "fasta"]) 
```
Or chosing frame in DNA translation to protein, in which not choosing a frame will translate in all three frames, while choosing one will only translate that specific frame:
```python
seq.translate_dna_to_protein(frame=1) 
```

Much more of less-impactful parameters may be modified in the functions and parameters present in main. For more information about them, you can either check those manually or use the command help(*the function/method you want help with*).

---

## Analysis Workflow

The main script coordinates the analyses in the following order.

### Step 1 — Input collection

The program gathers the transcript, region, flank, fragment, enzyme, intron-inclusion, and email settings.

### Step 2 — Transcript information retrieval

A `SearchInfo` object retrieves the information required for the selected transcript and region. This includes genomic and transcript-related coordinates, strand information, and sequence retrieval.

Then, the `SearchInfo` object retrieves the actual sequence by using the information retrieved in the previous step

### Step 3 — Sequence object creation

A `Sequence` object is created from the retrieved sequence, and can be complemented with the previously acquired information too. It serves as the central object for the subsequent analyses.

### Step 4 — Basic sequence characterization

The program reports sequence information and performs the configured basic analyses, including nucleotide proportions, translation to protein, GC content, CpG island detection, restriction enzyme searches, and fragment searching.

### Step 5 — Variant generation

The interval mutator generates alternative single-nucleotide polymorphisms for the selected region. Depending on the output configuration, results may be exported to CSV and/or FASTA files.

### Step 6 — HEXplorer analysis

The generated variants are evaluated using the HEXplorer-based score calculation, comparing each mutated sequence with the reference one. Depending on the output configuration, results may be exported to a CSV file.

### Step 7 — ClinVar screening

The ClinVar module searches for matching records for supported variants and retrieves available information from NCBI. Depending on the output configuration, results may be exported to a CSV file.

### Step 8 — Canonical splice-site analysis

The main script also invokes the MaxEntScan-based analysis when the coordinate input meets the conditions configured in the workflow.
 This analysis is skipped for cDNA-coordinate inputs, as this sort of input is not supported. This is due to MaxEntScan being an exonic type of analysis. On the other hand, exonic input and whole transcript input is analyzed without issue

The workflow therefore combines several complementary approaches. Not every analysis measures the same biological property, and the outputs should be interpreted separately before drawing an overall conclusion.

---

## Output Files

The program may write several supported analysis results to CSV and/or FASTA files, using output paths configured in the individual methods. 

By default, a `./results/` folder will be created in the local directory, and the output files will be stored in their respective folders inside the result diectory.

Default file names are generated from the analysis name and, where applicable, the transcript and analyzed coordinates. Supported methods allow the file name or output path to be customized.

### CSV files

CSV files are intended for bigger results that can be deifficult to read in the terminal and can be represented in a tabular way.
Those can be opened in spreadsheet applications like excell or imported into R, Python, and other data-analysis environments.

Depending on the method and the configured output, CSV files may contain:

- Variant identifiers and mutated sequences, obtained by the interval mutator.
- Splicing analysis tendency to switch towards a more silencing or enhancig behavior, obtained by HEXplorer analysis.
- Score differences or percentage variations in splicing, obtained by MaxEntScan analysis.
- ClinVar descriptions and associated database fields.
- Results from restriction enzyme searches.

The exact columns depend on the analysis method. Consult the corresponding module when interpreting individual columns.

### FASTA files

FASTA is a text format used to store biological sequences. It typically contains a header line beginning with `>` followed by the nucleotide or amino acid sequence.

The interval mutator supports FASTA output for the generated variant sequences. This format is useful for sequence inspection and for subsequent processing with compatible bioinformatics tools.

### Output management

The default output paths are relative to the working directory from which the script is run. For predictable results, run `main.py` from the project root unless the paths have been adapted.

---

## Biological Interpretation

The analyses implemented in MAG_Sequence_Analyzer provide complementary information, but each has a different interpretation, allowing the user to perform a good variety of analysis over the same sequence, depending on what the user could actually be looking for.

### Sequence composition and translation

Nucleotide proportions and GC content describe sequence composition. Translation predicts the amino acid sequence under the selected reading frame and translation table.

These analyses do not, by themselves, establish whether a variant affects gene expression, protein function, or phenotype.

### CpG islands

CpG island detection identifies regions that meet the thresholds defined in the implementation. The output should be interpreted as a computational classification according to those thresholds, rather than as experimental confirmation of a regulatory region.

### HEXplorer scores

HEXplorer scores are used to assess changes associated with the sequence context of potential splice-regulatory elements, with negative scores suggesting a possible switch towards becoming a silencing region and viceversa.

A score change suggests a difference according to the model; it does not prove that the variant causes exon skipping, intron retention, or another splicing abnormality.
More information in the original article describing HEXplorer from Steffen Erkelenz et al. DOI: 10.1093/nar/gku736.

### MaxEntScan scores

MaxEntScan evaluates sequence patterns associated with canonical splice sites. A difference between reference and alternative scores can help prioritize variants for further investigation, specially with negative scores, as those mean a decrease in strength of the splice site (represented in a percentual way)

The magnitude and direction of a score change should be interpreted in the context of the relevant splice site, the transcript, and any supporting experimental or computational evidence.

### ClinVar records

ClinVar aggregates information submitted by different sources about human genetic variants and their relationships to health conditions.

When reviewing a record, consider at least:

- The exact variant identifier and reference sequence.
- The transcript and coordinate system used.
- The clinical significance and review status reported by the database.
- The condition associated with the record.
- The supporting evidence and possible differences between submissions.

---

## External Resources and Internet Connection

MAG_Sequence_Analyzer relies on external resources for some operations.

### Ensembl

Ensembl is used to retrieve transcript metadata, coordinate mappings, and sequence information. This is done through Ensembl API REST and thus require internet connection. Thus, program features using this are also bound to Ensembl servers be working correctly, because if they are not, it could lead to connection errors, and the user may have to retry the analysis

Ensembl records and transcript annotations may change over time. For reproducibility, you may want to record the transcript identifier and the date of analysis.

### NCBI Entrez and ClinVar

The ClinVar screening module queries NCBI through Entrez.

Internet connectivity is therefore required for database-dependent operations. Query results can also vary if the database is updated.

For the email parameter, the user may follow NCBI's current guidance for automated access and be careful about making too many frequent requests.

### Local resources

The HEXplorer analysis depends on the score matrix supplied by the project. Ensure that all local modules and supporting data files are present and accessible inside the local directory before running the workflow.

---

## Limitations and Considerations

Before using the results, keep the following limitations in mind.

1. **One transcript per analysis:** the workflow is organized around a single transcript sequence that can be utterly analyzed in a completely automatic way, but it cannot annalyze multiple trancript in one go. Thus, the user would have to execute an analysis per desired transcript
2. **Predictions are not experimental results:** HEXplorer and MaxEntScan provide computational scores, not direct measurements of splicing.
3. **ClinVar is a database, not a diagnostic decision:** database matches and reported clinical significance require independent interpretation in context.
4. **Variant generation is systematic:** generating a substitution does not imply that it occurs in a sample, is evolutionarily plausible, or has a biological effect.
5. **External data can change:** Ensembl annotations and ClinVar records may be updated over time.
6. **Sequence length and computational cost:** analyzing larger regions generates more variants and can increase execution time and output size.
7. **Input validation:** provide valid transcript identifiers, supported coordinate formats, and compatible parameters. Invalid or incomplete inputs may prevent some analyses from running.
8. **Clinical use:** the program is not a substitute for expert clinical interpretation, validated diagnostic pipelines, or experimental confirmation.

---

## Troubleshooting

### Python reports that a module cannot be imported

Check that:

- The terminal is open in the project root.
- The expected local `.py` files are present in the same directory.
- The correct virtual environment is active (If it has been created).
- Required third-party libraries are installed.

### The program cannot retrieve transcript information

Check the transcript identifier and internet connection. Confirm that the identifier is recognized by Ensembl and that the requested region uses a supported coordinate format.

### A sequence or coordinate error occurs

Review the transcript identifier, coordinate expression, and flank value. Make sure that the requested positions correspond to the intended coordinate system and region.

### ClinVar returns no records

A search with no matching record does not necessarily indicate that the variant is benign or that it has never been studied.

Check the variant identifier, transcript reference, and coordinate system. Also confirm that the database query completed successfully and that the internet connection is working.

### The program does not generate an expected output file

Check the configured output directory, the output format selected for the method, and whether the analysis reached the relevant export step. Relative output paths depend on the working directory.

### An analysis takes longer than expected

The interval mutator generates multiple substitutions per nucleotide position, and downstream analyses can process these variants individually. Larger regions may therefore take longer and produce larger output files.

---

## Acknowledgements and References

MAG_Sequence_Analyzer builds on established biological databases, Python libraries, and sequence-analysis methods.

Users should consult the original documentation and publications for the tools relevant to their analyses:

- **Ensembl:** https://www.ensembl.org/
- **NCBI ClinVar:** https://www.ncbi.nlm.nih.gov/clinvar/
- **NCBI Entrez Programming Utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **Biopython:** https://biopython.org/
- **maxentpy:** https://github.com/kepbod/maxentpy
- **HEXplorer(web):** https://rna.hhu.de/HEXplorer/
- **MaxEntScan(web):** https://genebe.net/tools/maxentscan

The splice-site analysis uses MaxEntScan-related methods. Please refer
to the original publication:

Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs
with applications to RNA splicing signals. Journal of Computational
Biology. 2004;11:377–394.
https://pubmed.ncbi.nlm.nih.gov/15285897/


The HEXplorer analysis is based on the published HEXplorer method:

Erkelenz S, et al. Genomic HEXploring allows landscaping of novel
potential splicing regulatory elements. Nucleic Acids Research.
2014;42(16):10681–10691.
https://doi.org/10.1093/nar/gku736

---

## License

The original code developed for MAG_Sequence_Analyzer is distributed
under the MIT License. See the `LICENSE` file for the full license text.

Third-party software, data, and resources used by the project may be
subject to their own licenses or terms of use. These materials are not
necessarily covered by the project's MIT License.