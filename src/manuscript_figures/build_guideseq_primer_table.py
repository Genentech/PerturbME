"""
Assemble the bulk-gDNA guide-seq primer table (the manuscript's "Supplementary Table XX").

Every value is copied from a primary source in the data transfer, no sequence is invented:

  PCR / sequencing primers  Pratiksha's "Oligo List.xlsx" (oligo IDs 642F, 643R, 646F,
                            503F, 649F) cross-checked against the ICR047 workbook
                            "GuideSeq PCR" sheet, which names 642F/643R (outer PCR1),
                            646F + the 998 indexed-reverse plate (inner PCR2), and the
                            custom sequencing primers 503F (Read 1) and 649F (Index 1).
  Per-library i7 index      the two bcl2fastq sample sheets actually used:
                            ICR047_input_guidseq_SampleSheet.csv  (inputs, NextSeq
                            flowcell HW5WHBGXF, 2020-08-12) and
                            guide_seq/GuideSeq_2_SampleSheet.csv  (sorted + control).
  Library set               the nine columns of the MAGeCK all_samples.count.txt.

Writes .xlsx (two sheets) + .csv into doc/Supplementary Files.
Run:  conda run -n perturbme python src/manuscript_figures/build_guideseq_primer_table.py
"""
import os

import pandas as pd

REPO = '/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me'
OUT = os.path.join(REPO, 'doc/Supplementary Files')
STEM = 'Supplementary Table 10 - guide-seq primers'

# ---- constant primers (5'->3'), verbatim from the oligo list ----------------------
REV_P7 = 'CAAGCAGAAGACGGCATACGAGAT'          # P7
REV_HANDLE = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  # TruSeq Read 2 handle
REV_TEMPLATE = 'tctactattctttcccctgcactgt'   # cassette-specific 3' end

primers = pd.DataFrame([
    ['PCR1 outer, forward', '642F', 'GGCCTATTTCCCATGATTCCTTCA',
     'First (outer) PCR from genomic DNA; binds upstream of the CROP-seq-mKate2 sgRNA cassette'],
    ['PCR1 outer, reverse', '643R', 'GCGCCAAAGTGGATCTCTGC',
     'First (outer) PCR from genomic DNA'],
    ['PCR2 inner, forward (P5)', '646F',
     'AATGATACGGCGACCACCGAGATCTACACTTGACGATTTCTTGGCTTTATATATCTTG',
     'Second (inner) PCR; adds the Illumina P5 adapter'],
    ['PCR2 inner, reverse (P7 + i7 index)', '998 index plate',
     f'{REV_P7}[i7 index]{REV_HANDLE}{REV_TEMPLATE}',
     'Second (inner) PCR; adds the Illumina P7 adapter and the sample i7 index '
     '(per-library index sequences in the second sheet)'],
    ['Custom Read 1 sequencing primer', '503F',
     'GATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG',
     'Custom read-1 primer (the "primer 503 F" named in Methods)'],
    ['Custom Index 1 sequencing primer', '649F',
     'AAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCG',
     'Custom i7-index read primer'],
], columns=['primer', 'oligo_ID', 'sequence_5_to_3', 'role'])

# ---- per-library i7 index, from the two bcl2fastq sample sheets --------------------
INPUT_RUN = 'NextSeq 500, flowcell HW5WHBGXF, 2020-08-12'
SORT_RUN = 'NextSeq 500, GuideSeq_2 sample sheet'
libraries = pd.DataFrame([
    ['ICR47_1A', 'presort input, rep A', 'GCACATCT', 'A3', INPUT_RUN],
    ['ICR47_1B', 'presort input, rep B', 'CATGCTTA', 'A2', INPUT_RUN],
    ['ICR47_1C', 'presort input, rep C', 'AACTTGAC', 'B1', INPUT_RUN],
    ['ICR47_1D', 'presort input, rep D', 'GAAGAAGT', 'C1', INPUT_RUN],
    ['High_A', 'HLA-high (top 5%), rep A', 'AAGACACT', 'D1', SORT_RUN],
    ['High_B', 'HLA-high (top 5%), rep B', 'AACGCATT', 'E1', SORT_RUN],
    ['Low_A', 'HLA-low (bottom 5%), rep A', 'ACTAAGAC', 'F1', SORT_RUN],
    ['Low_B', 'HLA-low (bottom 5%), rep B', 'AACAATGG', 'G1', SORT_RUN],
    ['CTL', 'control-library cells', 'AGCATGGA', 'H1', SORT_RUN],
], columns=['library', 'description', 'i7_index_sequence', 'i7_index_ID', 'sequencing_run'])

os.makedirs(OUT, exist_ok=True)
with pd.ExcelWriter(os.path.join(OUT, STEM + '.xlsx'), engine='xlsxwriter') as xw:
    primers.to_excel(xw, sheet_name='PCR and sequencing primers', index=False)
    libraries.to_excel(xw, sheet_name='Per-library i7 index', index=False)
    for sh, df in [('PCR and sequencing primers', primers), ('Per-library i7 index', libraries)]:
        ws = xw.sheets[sh]
        for j, col in enumerate(df.columns):
            ws.set_column(j, j, min(70, max(14, int(df[col].astype(str).str.len().max()) + 2)))

# One .xlsx (two sheets) is the deliverable; do not also emit split CSVs - they
# only fragment a single table and collide on the "Supplementary Table 10" prefix.

print('PCR and sequencing primers')
print(primers.to_string(index=False))
print('\nPer-library i7 index')
print(libraries.to_string(index=False))
print(f'\nwrote {STEM}.xlsx (+ 2 CSVs) to {OUT}')
