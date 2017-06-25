### Beta 7 (6/22/16):
Breaking Changes:
1. Renamed "missing" to "missingGenos". "missing:fieldName" may be used for excluding fields during search.

### Beta 6 (6/18/16):
Breaking Changes:
1. 7th column, "missing" added. This contains any samples whose genotypes are missing ('.|.', '././', or 'N')
  - In multiallelic cases, any missing samples will be reported for each allele
2. When annotating snp files, Seqant not longer drops samples with <.95 confidence. This made denovo (or case -control) queries innacurate

Improvements:
1. seq-vcf uses a much more restrictive approach to calling indels
  - Single-variant sites can only have 1 base padding in case of insertions, and in case of deletions the alt must be 1 in length
  - In case of multiallelics where the reference is longer than 1bp (as happens when one of the variants is a deletion), the insertion must begin after the 1st base (padding)
  - In cases of multiallelics where the reference is longer than 1bp and the alt is longer than 1bp (but shorter than reference), the deletion must similarly occur beginning with the 2nd bse
  - If the reference is longer than 1bp, if the site is an insertion it will be skipped
  - If the alt is longer than 1bp and smaller than the reference (deletion, but with extra padding), it will be skipped

Bug fixes:
1. By the above, in highly repetitive regions, seq-vcf will no longer call the allele to the right of the common sequence shared between the ref and the allele
  - By VCF spec it is not clear that this is strictly a bug; in fact VCF 4.2 documentation opens with a multiallelic site (AC / ACT) that is called in this way, making it impossible to conclusively call alleles with repetitive sequences
  - However, looking at the gnomAD team's interpretation of multiallelic sites, the current approach seems more consistent.
  - I also believe that it is appropriate to indicate to people that they should never place an allele with more than 1 base padding, since this decreases the value of having a pos column, and is in general a ridiculous way to represent variants.

Performance:
1/2 CPU use while parsing vcf files containing many samples
Somewhat faster snp parsing of the same

Notes: seq-vcf now supports custom FILTER'ing (defaults to PASS and "." as being acceptable), based on either exclusion criteria (excludeFilter) and inclusion criteria (keepFilter), ala vcftools/bcftools. This will be exposed in the web interface at some point. Similarly seq-vcf now supports keeping id and info columns. See seq-vcf commit history for more detail (we'll be writing documetnation for seq-vcf in the future, time permitting).

### Beta 5:
Breaking Changes:
#### VCF handling
seq-vcf now allows only PASS or "." FILTER column values

### Beta 4:
Breaking Changes:
#### Search mapping 
refSeq.name2, refSeq.nearest.name2, and refSeq.spDisplayID are no longer indexed using edge n-grams. This means that by default, it is searched exactly (except capitalization doesn't matter).

To search refSeq.name2 (or any field) by prefix, simply type * (star) after the term (ex: dlg\*)

This was done because many users reported getting unexpected results, and to standardize quotes as being used for phrases, rather than exact matches.

### Beta 3:
Breaking Changes:
#### Clinvar
Normalized clinvar header names. This uses the new fieldMap property, which can be applied to any header column values found in any sparse track (.bed-like file), and which replaces the required_fields_map property, which worked only for chromStart, chromEnd, and chrom.

Sparse tracks still require chromStart, chromEnd, and chrom fields, either to be present, or to be provided via a fieldMap transformation.
```yaml
fieldMap:
    '#AlleleID': alleleID
    AlternateAllele: alternateAllele
    Chromosome: chrom
    ClinicalSignificance: clinicalSignificance
    Origin: origin
    OtherIDs: otherIDs
    PhenotypeIDS: phenotypeIDs
    NumberSubmitters: numberSubmitters
    PhenotypeList: phenotypeList
    ReferenceAllele: referenceAllele
    ReviewStatus: reviewStatus
    Start: chromStart
    Stop: chromEnd
    Type: type
```
Update clinvar track features 
```yaml 
features:
  - alleleID: number
  - phenotypeList
  - clinicalSignificance
  - type
  - origin
  - numberSubmitters: number
  - reviewStatus
  - referenceAllele
  - alternateAllele
```

#### RefSeq
Replaced geneSymbol with name2. name2 is provided by UCSC in the refGene database, and is the equivalent of the kgXref geneSymbol value, except is sometimes more complete. Since we are LEFT joining on refGene (aka providing the RefSeq transcripts), it makes sense to use the RefSeq/refGene value whenever possible

#### RefSeq.nearest
Replaced geneSymbol with name2. See RefSeq note above.

#### RefSeq.clinvar
Updated refSeq.clinvar overlap records: 
```yaml
join:
    features:
    - alleleID
    - phenotypeList
    - clinicalSignificance
    - type
    - origin
    - numberSubmitters
    - reviewStatus
    - chromStart
    - chromEnd
```
