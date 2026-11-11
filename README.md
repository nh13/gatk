# GATK: high-performance germline calling

This is a fork of [GATK](https://github.com/broadinstitute/gatk). Its default branch, `high_performance_germline_calling`, is GATK plus a set of performance changes to HaplotypeCaller and, from `hpgc-v5`, to joint calling. The first wave of changes is open as pull requests to GATK, and pull requests for the rest will follow once the first wave is reviewed. This branch carries all of them until they are merged upstream.

For GATK itself (documentation, tools, support), see the [upstream README](https://github.com/broadinstitute/gatk/blob/master/README.md).

> [!TIP]
> **HaplotypeCaller on this branch needs 3.2–3.4x less CPU than GATK 4.7.0.0 with `--max-effective-depth 100`, and 3.0–3.1x less without the depth cap.** 
> 
> It runs each whole-genome shard in 1 vCPU and 4 GB of memory, where GATK 4.7.0.0 as WARP runs it takes 2 vCPU and 10 GB.

## Results

Whole-genome HaplotypeCaller in DRAGEN-GATK mode, with WARP's arguments and WARP's 50-way sharding (see [How we benchmarked](#how-we-benchmarked)):

- GATK 4.7.0.0 runs as WARP runs it, with 2 vCPU, 10 GB, a 9 GB heap and 4 PairHMM threads
- This branch runs with 1 vCPU, 4 GB and a 3 GB heap, with and without the depth cap.

|  | HG00123 | HG00123 | HG00123 | NA18489 | NA18489 | NA18489 |
|---|---|---|---|---|---|---|
| GATK version | 4.7.0.0 | 4.7.0.0-hpgc-v4 | 4.7.0.0-hpgc-v4 | 4.7.0.0 | 4.7.0.0-hpgc-v4 | 4.7.0.0-hpgc-v4 |
| Mean depth | 35x | 35x | 35x | 65x | 65x | 65x |
| Depth cap | none | none | 100 | none | none | 100 |
| CPU hours | 18.1 | 6.1 | 5.6 | 33.2 | 10.7 | 9.8 |
| Allocated vCPU hours | 22.9 | 6.2 | 5.7 | 43.3 | 10.7 | 9.9 |
| Peak RSS | 6.7 GB | 3.3 GB | 1.7 GB | 7.1 GB | 3.3 GB | 1.8 GB |
| Shard wall time: min / median / max | 7.6 / 11.3 / 59.3 min | 4.8 / 6.9 / 17.4 min | 4.5 / 6.8 / 9.6 min | 16.2 / 22.1 / 89.7 min | 8.9 / 12.0 / 25.9 min | 7.3 / 11.8 / 15.9 min |
| Slowest ÷ median shard | 5.2 | 2.5 | **1.4** | 4.1 | 2.2 | **1.3** |
| Speedup vs 4.7.0.0 (CPU hours) | – | 3.0x | **3.2x** | – | 3.1x | **3.4x** |

HaplotypeCaller on this branch runs a shard in 1 vCPU and 4 GB or less, so its jobs pack onto the common 4 GB-per-vCPU machine shapes, such as AWS m-family instances and GCP n2-standard. With the depth cap a 1.5 GB heap is also enough, so a shard fits in 2 GB and jobs pack onto 2 GB-per-vCPU shapes such as AWS c-family instances. 

GATK 4.7.0.0 as WARP runs it reserves 2 vCPU and 10 GB per shard, so memory limits how many fit. A 16-vCPU, 64 GB machine runs 6 GATK 4.7.0.0 shards, with 4 vCPUs left idle, and 15–16 shards of this branch, depending on how much memory the scheduler holds back.

Without the depth cap the genotypes are ~identical to GATK's; with it they differ at a small number of sites, almost all outside the confidently callable genome (see [Concordance](#concordance)).

## How we benchmarked

**Samples.** Two 1000 Genomes high-coverage CRAMs (NYGC, GRCh38), read from the AWS Open Data copy:

- HG00123, about 35x: [`HG00123.final.cram`](https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3240197/HG00123.final.cram)
- NA18489, about 65x: [`NA18489.final.cram`](https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data/ERR3239336/NA18489.final.cram)

**Reference.** [`GRCh38_full_analysis_set_plus_decoy_hla.fa`](https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa) with its `.fai` and a sequence dictionary. The DRAGstr model for each sample comes from `ComposeSTRTableFile` and `CalibrateDragstrModel`.

**Shards.** WARP's calling regions, [`wgs_calling_regions.hg38.interval_list`](https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list), split 50 ways the way WARP splits them: `IntervalListTools` with `SCATTER_COUNT=50`, `SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW`, `UNIQUE=true`, `SORT=true` and `BREAK_BANDS_AT_MULTIPLES_OF=1000000`. Each shard is about 58 Mb.

**Arguments.** WARP's `HaplotypeCaller_GATK4_VCF` command in DRAGEN maximum-quality mode:

```
-contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation
--dragen-mode --dragstr-params-path <sample model>
-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF
```

**Resources per shard.**

| | vCPUs | Memory | Java heap | PairHMM threads |
|---|---|---|---|---|
| GATK 4.7.0.0 (WARP's task) | 2 | 10 GB | `-Xmx9000m -Xms9000m` | 4 (default) |
| This branch | 1 | 4 GB | `-Xmx3g -Xms3g` | 1 (kernels run on the calling thread) |

**Hardware.** AWS `r8a.4xlarge` instances: AMD EPYC 9R45 (Zen 5, AVX-512), 16 cores without SMT, 128 GB, running Amazon Linux 2023 and Amazon Corretto 17 JVM. Each shard is pinned to its own cores with `taskset`, so it sees exactly its allocation, as in a VM of that size. Shards are packed onto an instance up to its 16 cores. On this instance family a vCPU is a whole core. Under WARP on Google Cloud a vCPU is a hyperthread, so GATK 4.7.0.0 gets more CPU here than it would under WARP.

**Measurements.** For each shard, `/usr/bin/time` gives wall time, CPU time (user + system) and peak memory (maximum resident set size). Allocated vCPU hours are the sum over shards of wall time × vCPUs requested.

**Per-PR measurements.** Each PR is measured on top of the one before it, in table order, on four shards:

- two median 23 Mb chr1 shards;
- a 5.85 Mb segmental-duplication shard;
- a 5.85 Mb shard spanning the chr1 centromere and 1q21, the most expensive region of chr1.

Each shard is run on both samples with two repeats. The baseline and the candidate run side by side on single dedicated cores, with a 4 GB heap and one PairHMM thread, so the difference measures the code rather than threading. The gVCF records of each pair are compared for identity.

## What's Changed on this Branch?

Each batch is a tag on this branch (`4.7.0.0-hpgc-v1` to `4.7.0.0-hpgc-v5`), and a change's batch is the first tag that contains it. CPU change is measured on x86 as described [above](#how-we-benchmarked), against the row before, as HG00123 / NA18489; "within noise" means under about 1.5% or not consistent between the two samples.

| PR | Batch | Change | Changes gVCF? | Median shards | Segdup shard | chr1 centromere shard |
|---|---|---|---|---|---|---|
| [#9433](https://github.com/broadinstitute/gatk/pull/9433) | `hpgc-v1` | Adds `--max-effective-depth` for assembly-region walkers: reads are downsampled so no position in a region is covered by more than N reads, which bounds the cost of collapsed repeats and centromeres. Off unless set. | Only where depth exceeds the cap, when set | no change | −5.4% / −9.5% | −76.4% / −74.3% |
| [#9434](https://github.com/broadinstitute/gatk/pull/9434) | `hpgc-v1` | Stops installing Intel GKL's deflater and inflater, so BGZF compression uses htsjdk's libdeflate-backed defaults. | No (the same records; compressed bytes may differ) | within noise | within noise | within noise |
| [#9436](https://github.com/broadinstitute/gatk/pull/9436) | `hpgc-v1` | Caches `Haplotype.hashCode`. The PairHMM looks each haplotype up in a map once per read-haplotype pair, rehashing its bases every time. | No | within noise | within noise | within noise |
| [#9437](https://github.com/broadinstitute/gatk/pull/9437) | `hpgc-v1` | Replaces boxed `Integer` caches in the read adapter (soft start and end, adaptor boundary, cigar length) with primitives. The locus iterator reads them once per pileup element per locus. | No | within noise | within noise | within noise |
| [#9438](https://github.com/broadinstitute/gatk/pull/9438) | `hpgc-v1` | Presizes the kmer set used to find non-unique kmers for each read and kmer size, so it never rehashes. | No | within noise | within noise | within noise |
| [#9439](https://github.com/broadinstitute/gatk/pull/9439) | `hpgc-v1` | Keeps the locus iterator's per-sample read states in an `ArrayList`, dropping finished ones by in-place compaction, instead of a `LinkedList` with iterator removal. | No | within noise | within noise | within noise |
| [#9440](https://github.com/broadinstitute/gatk/pull/9440) | `hpgc-v1` | Holds the per-sample read-state managers in an array in sample order instead of a map keyed by sample name. | No | within noise | within noise | within noise |
| [#9441](https://github.com/broadinstitute/gatk/pull/9441) | `hpgc-v2` | Replaces Intel GKL's PairHMM, PD-PairHMM and Smith-Waterman with [fgkl](https://github.com/fulcrumgenomics/fgkl): Rust kernels (AVX-512, AVX2, NEON and scalar), picked at runtime and run on the calling thread. | PL and GQ in the last floating-point digits at some sites; see [Concordance](#concordance) | −7.3% / −21.3% | −28.2% / −41.5% | −50.7% / −52.1% |
| [`perf/read-major-ref-confidence`](https://github.com/tfenne/gatk/tree/perf/read-major-ref-confidence) | `hpgc-v3` | Reference confidence computed read by read: each read's cigar is walked once per region, instead of building a second pileup at every position. | No | −14.2% / −13.7% | −8.8% / −8.4% | −2.6% / −2.4% |
| [`perf/frd-primitive-arrays`](https://github.com/tfenne/gatk/tree/perf/frd-primitive-arrays) | `hpgc-v3` | DRAGEN foreign-read detection over primitive arrays: each read's support is computed once, and the three strand models share each pass over the reads. | No | within noise | within noise | −21.0% / −25.6% |
| [`perf/assembly-graph`](https://github.com/tfenne/gatk/tree/perf/assembly-graph) | `hpgc-v3` | Assembly graph: plain jgrapht adjacency with single-lookup queries, edges that store their own endpoints, reference-path pruning without streams, and no duplicate-edge search where none can exist. | No | −10.0% / −8.0% | −8.3% / −7.1% | −7.1% / −6.2% |
| [`perf/isactive-skip-no-alt-loci`](https://github.com/tfenne/gatk/tree/perf/isactive-skip-no-alt-loci) | `hpgc-v3` | Skips the activity-profile likelihoods at loci where no counted base is alt evidence, where the result is known to be exactly 0. | No | within noise | within noise | within noise |
| [`perf/mann-whitney-exact-dp`](https://github.com/tfenne/gatk/tree/perf/mann-whitney-exact-dp) | `hpgc-v3` | The exact Mann-Whitney U test counts rank sums with a small dynamic program instead of enumerating permutations. | No | −5.1% / −0.7% | −4.9% / −2.5% | −1.5% / −1.8% |
| [`perf/packed-nonunique-kmers`](https://github.com/tfenne/gatk/tree/perf/packed-nonunique-kmers) | `hpgc-v3` | Finds non-unique kmers by comparing two-bit-packed kmers in a primitive set. | No | −2.1% / −2.0% | −1.6% / −0.9% | within noise |
| [`perf/avoid-cigar-copies`](https://github.com/tfenne/gatk/tree/perf/avoid-cigar-copies) | `hpgc-v3` | Reads cigar elements in place instead of copying the cigar (locus iterator, read clipping). | No | −1.1% / −2.5% | −0.9% / −1.6% | within noise |
| [`perf/libs-per-locus-churn`](https://github.com/tfenne/gatk/tree/perf/libs-per-locus-churn) | `hpgc-v3` | Trims per-locus work in the locus iterator: a cached adaptor boundary, no partitioning cycle at loci where no read starts, and a presized pileup. | No | within noise | within noise | within noise |
| [`perf/identity-read-removal`](https://github.com/tfenne/gatk/tree/perf/identity-read-removal) | `hpgc-v3` | Removes reads from an assembly region by identity instead of `equals`. | No | −2.0% / −1.3% | −3.7% / −1.2% | within noise |
| [`perf/inactive-region-trims`](https://github.com/tfenne/gatk/tree/perf/inactive-region-trims) | `hpgc-v3` | Small trims in inactive regions: hard-clipped reads are tracked only for pileup detection, a single-sample path for overlapping read pairs, and the median depth from an `int[]`. | No | within noise | within noise | within noise |
| [`perf/lazy-called-haplotypes-message`](https://github.com/tfenne/gatk/tree/perf/lazy-called-haplotypes-message) | `hpgc-v3` | Builds the message for a validation check on the called haplotypes only when the check fails, instead of rendering every variant and haplotype of each region as a string. | No | within noise | within noise | within noise |
| [`perf/fgkl-0.2.0`](https://github.com/tfenne/gatk/tree/perf/fgkl-0.2.0) | `hpgc-v4` | Updates fgkl to 0.2.0: Smith-Waterman fill 1.3–2.2x faster, and a PairHMM that shares the computation of haplotype suffixes as well as prefixes, cutting its kernel time by more than half. | PL and GQ at some sites; see [Concordance](#concordance) | −3.8% / −7.5% | −6.6% / −13.9% | −30.8% / −31.9% |
| [`nh/hc-perf`](https://github.com/nh13/gatk/tree/nh/hc-perf) | — | Assembly graph vertices point at their own edges, so adjacency queries skip the vertex-map lookup; edge lists are indexed in read threading; cycle checks use Kahn's algorithm; edges keep their multiplicities in a primitive array. Reference confidence sweeps each read a cigar element at a time into flat per-position arrays, and skips the genotype allele check on its sites. IUPAC base intersection reads a primitive mask table. Not yet measured as above; on HG002 30x with an Apple M-series Mac, instructions retired fell 15.5% on a median shard (chr7:72–118 Mb) and 14.8% on an expensive one (chr1:92–156 Mb). | No | not yet measured | not yet measured | not yet measured |

### Joint calling

Changes to reblocking, GenomicsDB and GnarlyGenotyper, measured on 100 reblocked samples, or more where the row says so, as described [below](#how-we-measured-joint-calling), against the same run without the change.

| PR | Batch | Change | Changes output? | GenomicsDBImport | GnarlyGenotyper |
|---|---|---|---|---|---|
| [#9446](https://github.com/broadinstitute/gatk/pull/9446) | `hpgc-v5` | ReblockGVCF's header declares only the fields its records can carry. It no longer declares the final annotations it never computes (QD, FS, SOR, the allele-specific finals and others), the annotations it removes, DRAGEN's DRAGstr fields, or MIN_DP under `--floor-blocks`. GenomicsDB stores and processes every declared field for every record, and a reblocked gVCF now declares 24–28 fields instead of 44–48. | Header lines only; records are unchanged, and GnarlyGenotyper and GenotypeGVCFs make identical calls | −40% | −19% |
| [#9447](https://github.com/broadinstitute/gatk/pull/9447) | `hpgc-v5` | Adds `--genomicsdb-compression <codec>[:<level>]` to GenomicsDBImport, to compress the workspace's tiles with gzip, zstd or lz4 instead of GenomicsDB's default, gzip at level 6. lz4 doubles the workspace's size and zstd:1 adds about 10%. zstd needs the system's `libzstd` wherever the workspace is written or read. | No; off unless set, and GnarlyGenotyper's calls are identical for every codec | lz4 −25%, zstd:1 −20% | lz4 −13%, zstd:1 −6% |
| [`perf/gnarly-genotype-path`](https://github.com/tfenne/gatk/tree/perf/gnarly-genotype-path) | `hpgc-v5` | GnarlyGenotyper counts each sample's called alleles in an array indexed by allele instead of a map keyed by `Allele`, whose hash code is recomputed from its bases on every lookup. It also no longer builds and validates a throwaway record for each allele-specific annotation it finalises. | No; records are identical | — | −2% at 1,000 samples; −5% with allele-specific annotations |
| [`perf/genomicsdb-skip-non-variant-intervals`](https://github.com/tfenne/gatk/tree/perf/genomicsdb-skip-non-variant-intervals) | `hpgc-v5` | GnarlyGenotyper and GenotypeGVCFs ask GenomicsDB to skip the intervals that cannot produce a variant site: those where every sample is in a reference block, which are nearly all of a cohort's intervals, and those whose only alternate allele is a spanning deletion. The skip is off with `--keep-all-sites`, `--include-non-variant-sites` or `--force-output-intervals`, or with `--genomicsdb-skip-non-variant-intervals false`. It needs GenomicsDB export options that 1.5.5 lacks, so the branch builds against the fork's GenomicsDB, published as `com.tfenne:genomicsdb` 1.6.0. | No; records are identical | — | −59% at 1,000 samples, and a further −3% from the spanning deletions. GenotypeGVCFs: −44% at 100 samples, and a further −15% at 1,000 |
| [GenomicsDB 1.6.0](https://github.com/tfenne/GenomicsDB/tree/v1.6.0) | `hpgc-v5` | Builds against the fork's GenomicsDB, published as `com.tfenne:genomicsdb:1.6.0`, instead of GenomicsDB 1.5.5. Its import resolves each INFO and FORMAT field once per input file instead of by name in every record, casts its readers once instead of for every field, and fetches only the fields each record has. Its query keeps TileDB's files open between tile reads, reads tiles into reused buffers, and does less work per cell and per record when it combines samples. Its Linux libraries link their dependencies statically and need only glibc 2.28 or later and zlib. | No; workspaces and records are identical | −53% at 1,000 samples | −48% at 1,000 samples |
| [htsjdk 5.0.2](https://github.com/samtools/htsjdk/releases/tag/5.0.2) | `hpgc-v5` | Updates htsjdk from 5.0.0 to 5.0.2. 5.0.1 decodes BCF genotype data straight from each record's bytes instead of reading every byte through a synchronized stream ([#1856](https://github.com/samtools/htsjdk/pull/1856)), and writes VCF genotype fields with much less per-sample overhead, formatting numbers without `java.util.Formatter` ([#1857](https://github.com/samtools/htsjdk/pull/1857), [#1858](https://github.com/samtools/htsjdk/pull/1858)). 5.0.2 fixes a race in htsjdk's cached CRAM reference that occasionally failed a multi-threaded `CalibrateDragstrModel` on CRAM input with a reference MD5 mismatch ([#1868](https://github.com/samtools/htsjdk/pull/1868)). | No; records are identical | — | −31% at 3,000 samples |

#### How we measured joint calling

- **Samples.** 100 of the [1000 Genomes high-coverage](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) per-sample gVCFs that NYGC made with GATK 3.5 HaplotypeCaller, from the AnVIL workspace [`anvil-datastorage/1000G-high-coverage-2019`](https://anvil.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019), reblocked with WARP's arguments: `ReblockGVCF -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40`. For #9446, the same gVCFs reblocked with this branch and with `hpgc-v4`. Rows measured at 1,000 samples use 1,000 of the same gVCFs. "With allele-specific annotations" means those gVCFs with the allele-specific raw annotations this branch's HaplotypeCaller writes added to every variant record, built from the record's own values, since NYGC's gVCFs lack them. The htsjdk row uses 3,000 of the same 1000 Genomes samples called from their public CRAMs with this branch's HaplotypeCaller and WARP's arguments, reblocked as above, over chr20:4–8 Mb.
- **Region.** WARP's calling regions in chr20:1–16 Mb (15.9 Mb).
- **Import.** All samples in one batch through GenomicsDB's native reader (`--bypass-feature-reader`) with `--genomicsdb-shared-posixfs-optimizations`; the #9446 runs also used `--genomicsdb-compression lz4`. GenomicsDB 1.5.5 as released, except in the last four rows, whose runs used builds of the fork's GenomicsDB.
- **Genotyping.** GnarlyGenotyper with WARP's arguments (`-stand-call-conf 10 --max-alternate-alleles 5`) and `--genomicsdb-use-bcf-codec`.
- **Hardware.** Single-threaded wall time on an Apple M-series Mac, not the AWS instances used for HaplotypeCaller. The runs compared in each row were made in the same benchmark pass. The GenomicsDB and htsjdk rows compound one such comparison for each change they contain.

## Concordance

Each shard's gVCF was genotyped with single-sample `GenotypeGVCFs` from GATK 4.7.0.0, for every run alike. Variant genotypes were then compared at QUAL ≥ 30, which is the default calling threshold of both `GenotypeGVCFs` and `HaplotypeCaller` in VCF mode, as WARP runs it. A genotype counts as identical when its position, alleles and GT match.

**Without the depth cap**, 8 of 4.95 M variant sites differ in HG00123 and 13 of 5.99 M in NA18489. 11/21 are the same call with QUAL slightly above 30 on one side and below it on the other, from small floating-point differences between fgkl's PairHMM and Intel GKL's. Six are low-confidence genotypes, or the same repeat alleles written differently. Four are called by only one side, three of them in the chr1 and chr11 centromeres and 1q21.

**With `--max-effective-depth 100`**, 99.87% (HG00123) and 99.77% (NA18489) of variant genotypes are identical. The cap downsamples reads only where depth exceeds 100, which happens in collapsed repeats, centromeres and artifact pileups, and the sites that differ lie there: 0% (HG00123) and 1% (NA18489) of them fall inside the [GIAB HG001 v4.2.1 benchmark regions](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/), which hold 77% and 78% of all variant sites.

The alleles at these sites are mostly known to [gnomAD v4.1](https://gnomad.broadinstitute.org/) genomes, but gnomAD filters most of them itself: while 55–78% of them are in gnomAD only 7–25% with FILTER PASS, against 97–98% present and 89–91% PASSing for a sample of identical calls. This branch's own calls are no less often in gnomAD than GATK 4.7.0.0's (66% against 55–65%).

| Sites that differ with the cap | HG00123 | NA18489 |
|---|---|---|
| Called by this branch only | 2,851 | 5,322 |
| Called by GATK 4.7.0.0 only | 1,329 | 2,117 |
| The same call, with QUAL below 30 on one side | 1,904 | 4,767 |
| Different genotype at the same alleles | 1,187 | 3,379 |
| Different alleles at the same position | 961 | 1,413 |
| Other | 63 | 120 |
| **Total** (variant sites called by GATK 4.7.0.0: 4.95 M / 5.99 M) | **8,295** | **17,118** |

## Running it

Build with `./gradlew localJar` (Java 17), or download the jar attached to the latest [release](https://github.com/tfenne/gatk/releases). Use WARP's HaplotypeCaller arguments, add `--max-effective-depth 100`, and request 1 vCPU and 4 GB per shard with `-Xmx3g`, the heap the whole-genome results above were measured with. Test on your own data before using a smaller heap.

For joint calling, reblock with this branch's ReblockGVCF and WARP's arguments (`-do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40`; GnarlyGenotyper skips any variant without the `QUALapprox` that `-do-qual-approx` adds), and add `--genomicsdb-compression lz4` to GenomicsDBImport, or `zstd:1` if the workspace is copied between hosts that all have `libzstd`. From `hpgc-v5` the branch builds against the fork's GenomicsDB (`com.tfenne:genomicsdb:1.6.0`) and htsjdk 5.0.2, both on Maven Central.

## How the branch is maintained

The branch is upstream `master` with each PR branch merged in, plus fork-only commits (this README and CI settings). It only moves forward. Upstream is merged in periodically, and a PR branch that changes in review is merged in again. A PR that is merged upstream arrives through the next upstream merge.
