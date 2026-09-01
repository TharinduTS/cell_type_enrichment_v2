# 8) Protein type data integration

Introduction
I found data on protein types from HPA:

"The human secretome comprises all proteins that are potentially secreted from the cell. These proteins were identified based on UniProt annotations of subcellular locations, along with computational predictions of signal peptides and transmembrane regions. The signal peptide is found in most secreted proteins, but also appears in some classes of membrane proteins. Therefore, the presence of transmembrane regions can be used to distinguish membrane proteins from secreted proteins. A whole-proteome scan of all Ensembl transcripts was performed using majority decision methods for signal peptide prediction (MDSEC) and membrane region prediction (MDM). Proteins with a predicted signal peptide by MDSEC and no predicted transmembrane region by MDM were considered secreted. Using this combined approach, the analysis predicts that 2520 genes (12% of all human protein-coding genes) encode at least one predicted secreted transcript or have at least one isoform with the subcellular location "Secreted" in UniProt."

HPA: https://www.proteinatlas.org/humanproteome/tissue/secretome#classification_of_the_secretome

And I downloaded the dataset from

```http
https://www.proteinatlas.org/search/sa_location%3ASecreted+-+unknown+location%2CSecreted+in+brain%2CSecreted+in+female+reproductive+system%2CSecreted+in+male+reproductive+system%2CSecreted+in+other+tissues%2CSecreted+to+blood%2CSecreted+to+digestive+system%2CSecreted+to+extracellular+matrix%2CIntracellular+and+membrane%2CImmunoglobulin+genes
```

It had protein type data in the column 'Protein type' as the first comma seperated value that starts with 'Predicted'. I extracted those protein type data with other needed columns

## 8-I Extract essential data

### CLI help

```py
USAGE
-----
    python extract_regex.py [options]

OPTIONS
-------

Input / Output
--------------
--input <path>           (Required) Path to input TSV/CSV file.
--output <path>          (Required) Path to output file.
--file-delim <char>      File delimiter. Default: \t
--output-delim <char>    Output file delimiter. Default: \t
--encoding <encoding>    File encoding. Default: utf‑8.

Column Selection
----------------
--id-cols <col1 col2...> (Required) Columns to keep in output (e.g., Gene Ensembl).
--target-col <column>    (Required) Column containing the delimited values to parse.

Token Parsing
-------------
--value-sep <char>       Separator inside the target column. Default: ,
--pattern <regex>        (Required) Regex pattern to match tokens.
                         Examples:
                           "Predicted|Secreted"
                           "^(Predicted|Secreted)"
                           ".*(Predicted|Secreted).*"
--nth <int>              Select the n‑th matching token (1‑based). Default: 1.
--all-matches			 Gets all matches in addition to n-th
--ignore-case            Enable case‑insensitive regex matching.
--match-anywhere         Use re.search() instead of re.match().
                         Default behavior uses re.match() (anchored at start).

Performance
-----------
--chunksize <int>        Process file in chunks (e.g., 200000) for large datasets.

BEHAVIOR SUMMARY
----------------
• The target column is split into tokens using the delimiter (default comma).
• Tokens are trimmed.
• Regex is applied to each token:
    - Default: re.match() → matches only at the start of the token.
    - With --match-anywhere: re.search() → matches anywhere in token.
• The script extracts the n‑th match (1‑based index).
• Output contains:
    - All selected --id-cols
    - A new column named:
        <target-col>__regex_nth<n>

EXAMPLES
--------

1) First token starting with “Predicted” or “Secreted”
------------------------------------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out.tsv \
        --file-delim $'\t' \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern 'Predicted|Secreted' \
        --nth 1 \
        --ignore-case

2) Second matching token
------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out2.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --pattern 'Predicted|Secreted' \
        --nth 2 \
        --ignore-case

3) Match anywhere, not just start
---------------------------------
    python extract_regex.py \
        --input genes.tsv \
        --output out_anywhere.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --value-sep "," \
        --pattern '(Predicted|Secreted)' \
        --match-anywhere \
        --ignore-case

4) Large files: chunked processing
----------------------------------
    python extract_regex.py \
        --input bigfile.tsv \
        --output out_big.tsv \
        --id-cols Gene Ensembl \
        --target-col "Protein class" \
        --pattern 'Predicted|Secreted' \
        --chunksize 250000

NOTES
-----
• Regex patterns must be quoted properly in shell environments.
• For TSV files, prefer --file-delim $'\t' to ensure correct parsing.
• Missing or empty target values produce empty output tokens.
```

### Run

```bash
python extract_regex.py \
  --input secretome_dataset.tsv \
  --output genes_with_protein_type.tsv \
  --file-delim $'\t' \
  --id-cols Gene Ensembl \
  --target-col "Protein class" \
  --value-sep "," \
  --pattern 'Predicted' \
  --ignore-case \
  --all-matches \
  --join-matches-sep ": "
```


