package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Records, on input records whose genotype priors were written by the caller for that record, the provenance of
 * those priors in {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY} (see {@link GenotypePriorProvenance}), and
 * counts input records whose priors cannot be used.
 *
 * <p>Priors on a record are taken as written for it only if the record has a single sample, the tool reads records as
 * the caller wrote them (not, e.g., from GenomicsDB, which merges samples), and the priors fit the record (see
 * {@link GenotypePriorProvenance#whyPriorsDoNotFitRecord}). A record that already has provenance keeps it. The priors
 * of other records are not used, and each kind is reported once and counted. Thread-safe.</p>
 */
public final class GenotypePriorProvenanceRecorder {

    /** The name of the argument of the tools that preserve genotype priors when merging GVCFs. */
    public static final String PRESERVE_GENOTYPE_PRIORS_LONG_NAME = "preserve-genotype-priors";

    /** The documentation of the argument of the tools that preserve genotype priors when merging GVCFs. */
    public static final String PRESERVE_GENOTYPE_PRIORS_DOC = "Record each sample's genotype priors (PG or PRI) with the " +
            "alleles they were written for, for GenotypeGVCFs --use-input-genotype-priors";

    private final boolean trustSingleSampleRecords;
    private final String remedy;
    private final Logger logger;
    private final AtomicLong recordsWithoutProvenance = new AtomicLong();
    private final AtomicLong recordsWithPriorsNotFittingThem = new AtomicLong();

    /**
     * @param trustSingleSampleRecords whether the tool reads records as the caller wrote them, so that the priors of a
     *                                 single-sample record may have been written for it
     * @param remedy what to do to make unusable priors usable, for warnings
     * @param logger where to log warnings
     */
    public GenotypePriorProvenanceRecorder(final boolean trustSingleSampleRecords, final String remedy, final Logger logger) {
        this.trustSingleSampleRecords = trustSingleSampleRecords;
        this.remedy = Utils.nonNull(remedy, "remedy");
        this.logger = Utils.nonNull(logger, "logger");
    }

    /**
     * @param vc an input record
     * @param count whether to count the record if its priors cannot be used; false when the same record was already
     *              given at an earlier position
     * @return the record with the provenance of its priors added, or {@code vc} itself if nothing was added
     */
    public VariantContext apply(final VariantContext vc, final boolean count) {
        Utils.nonNull(vc, "vc");
        if (vc.getGenotypes().stream().noneMatch(GenotypePriorProvenance::hasPriorsWithoutProvenance)) {
            return vc;
        }
        if (!trustSingleSampleRecords || vc.getNSamples() != 1) {
            if (count && recordsWithoutProvenance.getAndIncrement() == 0) {
                logger.warn(String.format("Input record at %s:%d has genotype priors (%s or %s) without %s, so they may not " +
                                "match its alleles and will not be used. %s Further such records are only counted.",
                        vc.getContig(), vc.getStart(), GATKVCFConstants.GENOTYPE_PRIOR_KEY, GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY,
                        GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, remedy));
            }
            return vc;
        }
        final Genotype genotype = vc.getGenotype(0);
        final double[] priors = GenotypePriorProvenance.readCallerPriors(vc, genotype);
        if (priors == null) {
            return vc;
        }
        final String misfit = GenotypePriorProvenance.whyPriorsDoNotFitRecord(vc, genotype, priors);
        if (misfit != null) {
            if (count && recordsWithPriorsNotFittingThem.getAndIncrement() == 0) {
                logger.warn(String.format("Sample %s at %s:%d has genotype priors that cannot have been written for its record (%s), " +
                                "as when a sample is taken out of a combined GVCF; they will not be used. Use GVCFs as written by " +
                                "the caller. Further such records are only counted.",
                        genotype.getSampleName(), vc.getContig(), vc.getStart(), misfit));
            }
            return vc;
        }
        final GenotypePriorProvenance provenance = new GenotypePriorProvenance(vc.getStart(), genotype.getPloidy(), vc.getAlleles(), priors);
        final Genotype withProvenance = new GenotypeBuilder(genotype)
                .attribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, provenance.encode())
                .make();
        return new VariantContextBuilder(vc).genotypes(GenotypesContext.create(withProvenance)).make();
    }

    /** @return the number of counted records whose priors had no provenance and were not used */
    public long getRecordsWithoutProvenance() { return recordsWithoutProvenance.get(); }

    /** @return the number of counted single-sample records whose priors did not fit them and were not used */
    public long getRecordsWithPriorsNotFittingThem() { return recordsWithPriorsNotFittingThem.get(); }

    /** Logs a warning with the number of records whose priors could not be used, if any. */
    public void logSummary() {
        if (getRecordsWithoutProvenance() > 0) {
            logger.warn(String.format("%d input records had genotype priors without %s, which were not used. %s",
                    getRecordsWithoutProvenance(), GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, remedy));
        }
        if (getRecordsWithPriorsNotFittingThem() > 0) {
            logger.warn(String.format("%d single-sample input records had genotype priors that cannot have been written for them, " +
                    "which were not used", getRecordsWithPriorsNotFittingThem()));
        }
    }
}
