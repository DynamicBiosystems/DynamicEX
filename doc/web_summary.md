### Sequencing

---

- **Number of Reads**

Total number of read pairs that were assigned to this library in demultiplexing.

- **Valid Barcodes**

Fraction of reads with barcodes that match the whitelist after barcode correction.

- **Sequencing Saturation**

The fraction of reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that had a non-unique (cell-barcode, UMI, gene).

- **GC Content**

GC Content.

- **Q20 Bases in RNA Read**

Fraction of RNA read bases with Q-score >= 20, excluding very low quality/no-call (Q <= 2) bases from the denominator.

- **Q30 Bases in RNA Read**

Fraction of RNA read bases with Q-score >= 30, excluding very low quality/no-call (Q <= 2) bases from the denominator.

### Cells

---

- **Estimated Number of Cells**

The number of barcodes associated with at least one cell.

- **Fraction Reads in Cells**

The fraction of valid-barcode, confidently-mapped-to-exonic reads with cell-associated barcodes.

- **Mean Reads per Cell**

The total number of sequenced reads divided by the number of barcodes associated with cell-containing partitions.

- **Median UMI Counts per Cell**

The median number of UMI counts per %s cell-associated barcode.

- **Median Genes per Cell**

The median number of genes detected per cell-associated barcode. Detection is defined as the presence of at least 1 UMI count.

- **Total Genes Detected**

The number of genes with at least one UMI count in any cell.

### Mapping

---

- **Reads Mapped to Genome**

Fraction of reads that mapped to the genome.

- **Reads Mapped Confidently to Genome**

Fraction of reads that mapped uniquely to the genome.

- **Reads Mapped Confidently to Intergenic Regions**

Fraction of reads that mapped uniquely to an intergenic region of the genome.

- **Reads Mapped Confidently to Intronic Regions**

Fraction of reads that mapped uniquely to an intronic region of the genome.

- **Reads Mapped Confidently to Exonic Regions**

Fraction of reads that mapped uniquely to an exonic region of the genome.