package org.broadinstitute.hellbender.tools.walkers.genotyper;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.escape.Escaper;
import com.google.common.net.PercentEscaper;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.VariantContextGetters;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

/**
 * The genotype priors a caller applied to one sample's record, kept together with the allele basis they were
 * computed for: the record's start, the sample's ploidy and the record's alleles, reference first.
 *
 * <p>HaplotypeCaller, when it assigns genotypes from posteriors (e.g. in DRAGEN mode), writes each sample's
 * Phred-scaled genotype priors as the Number=G FORMAT field PG; Illumina DRAGEN writes them as PRI. Tools that merge
 * GVCFs remap PL to the merged alleles but copy these fields unchanged, so after a merge a vector of the right length
 * can refer to other alleles or another genotype order. A {@code GenotypePriorProvenance} keeps the priors with the
 * alleles they belong to, so a consumer can map them onto a later allele list exactly, or tell that they do not apply.
 * {@link GenotypePriorProvenanceRecorder} decides which input records' priors can be recorded.</p>
 *
 * <p>It is stored as the FORMAT field {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY}:
 * {@code 1|start|ploidy|allele/allele/...|prior/prior/...}. Characters reserved by VCF or by this encoding are
 * percent-encoded in alleles.</p>
 */
public final class GenotypePriorProvenance {

    /** The version of the serialized form written by this class, and the only one it reads. */
    public static final String VERSION = "1";

    private static final char FIELD_SEPARATOR = '|';
    private static final char VALUE_SEPARATOR = '/';
    private static final int NUMBER_OF_FIELDS = 5;

    /** Two prior vectors given for the same genotype (as PG and PRI) must agree to within this many Phred units. */
    private static final double MAX_PRIOR_DIFFERENCE_BETWEEN_KEYS = 1e-6;

    /**
     * The largest difference, in Phred units, between a caller's GP and its PL plus the priors (both scaled so the
     * best genotype is 0) for the priors to be taken as written for the record. HaplotypeCaller's GP is PL plus PG up
     * to rounding (at most 0.23 apart in one DRAGEN-mode genome); merged and re-extracted records whose PL was remapped
     * but whose PG and GP were not are usually much further apart.
     */
    @VisibleForTesting
    static final double MAX_POSTERIOR_DISCREPANCY = 1.0;

    /** Appended to errors about malformed priors, which are only read when an option asks for them. */
    private static final String PRIOR_OPTIONS_NOTE = " (Genotype priors are read only with --use-input-genotype-priors or --"
            + GenotypePriorProvenanceRecorder.PRESERVE_GENOTYPE_PRIORS_LONG_NAME + ".)";

    /** Besides alphanumerics, the characters left as they are in encoded alleles. */
    private static final String SAFE_ALLELE_CHARACTERS = "*<>_.-";
    private static final Escaper ALLELE_ESCAPER = new PercentEscaper(SAFE_ALLELE_CHARACTERS, false);

    /** Why a provenance's priors cannot be used for a site. */
    enum Inapplicability {
        /** The priors were written for a record starting at another position, e.g. a deletion spanning the site. */
        DIFFERENT_START,
        /** The site genotypes the sample at a different ploidy than the priors were written for. */
        DIFFERENT_PLOIDY,
        /** The site's reference allele does not extend the record's reference allele. */
        DIFFERENT_REFERENCE,
        /** An allele at the site was not among the record's alleles, so its priors are unknown. */
        NEW_ALLELE
    }

    private final int start;
    private final int ploidy;
    private final List<Allele> alleles;
    private final double[] phredPriors;

    /**
     * @param start       the start of the record the priors were written for
     * @param ploidy      the sample's ploidy at that record
     * @param alleles     the record's alleles, reference first
     * @param phredPriors the Phred-scaled genotype priors, in VCF genotype order for {@code alleles} and {@code ploidy}
     * @throws IllegalArgumentException if the priors do not fit the alleles and ploidy, or are not finite and non-negative
     */
    public GenotypePriorProvenance(final int start, final int ploidy, final List<Allele> alleles, final double[] phredPriors) {
        Utils.validateArg(start > 0, () -> "start must be positive but was " + start);
        Utils.validateArg(ploidy > 0, () -> "ploidy must be positive but was " + ploidy);
        Utils.nonEmpty(alleles, "alleles");
        Utils.validateArg(alleles.get(0).isReference(), "the first allele must be the reference allele");
        Utils.validateArg(alleles.stream().skip(1).noneMatch(Allele::isReference), "only the first allele may be the reference allele");
        Utils.validateArg(alleles.stream().allMatch(a -> a.isSymbolic() || a.isCalled()), "alleles may not be no-calls");
        Utils.validateArg(alleles.stream().map(Allele::getDisplayString).distinct().count() == alleles.size(), () -> "alleles must be distinct: " + alleles);
        Utils.nonNull(phredPriors, "phredPriors");
        final int expectedLength = GenotypeIndexCalculator.genotypeCount(ploidy, alleles.size());
        Utils.validateArg(phredPriors.length == expectedLength, () -> String.format(
                "expected %d genotype priors for %d alleles at ploidy %d but found %d", expectedLength, alleles.size(), ploidy, phredPriors.length));
        Utils.validateArg(Arrays.stream(phredPriors).allMatch(GenotypePriorProvenance::isValidPrior),
                () -> "genotype priors must be finite and non-negative Phred-scaled values: " + Arrays.toString(phredPriors));
        this.start = start;
        this.ploidy = ploidy;
        this.alleles = Collections.unmodifiableList(new ArrayList<>(alleles));
        this.phredPriors = phredPriors.clone();
    }

    int getStart() { return start; }

    int getPloidy() { return ploidy; }

    List<Allele> getAlleles() { return alleles; }

    /** @return a copy of the Phred-scaled priors */
    double[] getPhredPriors() { return phredPriors.clone(); }

    /**
     * Reads the priors that a caller wrote on a genotype (PG, or Illumina's PRI). Whether they were written for the
     * genotype's own record is for the caller to decide (see {@link #whyPriorsDoNotFitRecord}).
     *
     * @param vc the record the genotype belongs to
     * @param genotype the genotype
     * @return the Phred-scaled priors, or {@code null} if the genotype carries none or has no called ploidy
     * @throws UserException.BadInput if a prior is partly missing, not a number, negative or not finite, or PG and PRI disagree
     */
    static double[] readCallerPriors(final VariantContext vc, final Genotype genotype) {
        Utils.nonNull(vc, "vc");
        Utils.nonNull(genotype, "genotype");
        final double[] pg = readPhredValues(vc, genotype, GATKVCFConstants.GENOTYPE_PRIOR_KEY);
        final double[] pri = readPhredValues(vc, genotype, GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY);
        if (pg != null && pri != null && !agree(pg, pri)) {
            throw new UserException.BadInput(String.format("Sample %s at %s:%d has %s and %s genotype priors that disagree: %s and %s",
                    genotype.getSampleName(), vc.getContig(), vc.getStart(), GATKVCFConstants.GENOTYPE_PRIOR_KEY,
                    GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY, Arrays.toString(pg), Arrays.toString(pri)) + PRIOR_OPTIONS_NOTE);
        }
        final double[] priors = pg != null ? pg : pri;
        if (priors == null || genotype.getPloidy() == 0) {
            return null;
        }
        if (!Arrays.stream(priors).allMatch(GenotypePriorProvenance::isValidPrior)) {
            throw new UserException.BadInput(String.format("Sample %s at %s:%d has genotype priors that are not finite and " +
                    "non-negative Phred-scaled values: %s", genotype.getSampleName(), vc.getContig(), vc.getStart(), Arrays.toString(priors)) + PRIOR_OPTIONS_NOTE);
        }
        return priors;
    }

    /**
     * Checks that a caller's priors can have been written for this record: there is one per genotype of the record,
     * and, if the genotype has GP and PL, GP equals PL plus the priors (both scaled so the best genotype is 0) to within
     * {@link #MAX_POSTERIOR_DISCREPANCY}. A record merged with other samples' records and then taken out again (e.g. by
     * SelectVariants) has PL remapped to the merged alleles but PG and GP unchanged, so usually fails this check.
     *
     * @param vc the record
     * @param genotype the record's genotype carrying the priors
     * @param phredPriors the priors, from {@link #readCallerPriors}
     * @return why the priors do not fit the record, or {@code null} if they do
     * @throws UserException.BadInput if GP is malformed
     */
    static String whyPriorsDoNotFitRecord(final VariantContext vc, final Genotype genotype, final double[] phredPriors) {
        Utils.nonNull(vc, "vc");
        Utils.nonNull(genotype, "genotype");
        Utils.nonNull(phredPriors, "phredPriors");
        Utils.validateArg(genotype.getPloidy() > 0, "the genotype must have a called ploidy");
        final int genotypeCount = GenotypeIndexCalculator.genotypeCount(genotype.getPloidy(), vc.getNAlleles());
        if (phredPriors.length != genotypeCount) {
            return String.format("%d genotype priors for %d alleles at ploidy %d, which have %d genotypes",
                    phredPriors.length, vc.getNAlleles(), genotype.getPloidy(), genotypeCount);
        }
        final double[] posteriors = readPhredValues(vc, genotype, VCFConstants.GENOTYPE_POSTERIORS_KEY);
        if (posteriors == null || !genotype.hasPL()) {
            return null;  // nothing to check the priors against
        }
        final int[] likelihoods = genotype.getPL();
        if (posteriors.length != genotypeCount || likelihoods.length != genotypeCount) {
            return String.format("%d GP and %d PL values for %d genotypes", posteriors.length, likelihoods.length, genotypeCount);
        }
        final double[] expected = new double[genotypeCount];
        for (int i = 0; i < genotypeCount; i++) {
            expected[i] = likelihoods[i] + phredPriors[i];
        }
        final double best = Arrays.stream(expected).min().getAsDouble();
        for (int i = 0; i < genotypeCount; i++) {
            if (!(Math.abs(expected[i] - best - posteriors[i]) <= MAX_POSTERIOR_DISCREPANCY)) {
                return String.format("its GP %s is not its PL %s plus the priors %s", Arrays.toString(posteriors),
                        Arrays.toString(likelihoods), Arrays.toString(phredPriors));
            }
        }
        return null;
    }

    /**
     * Reads the provenance stored on a genotype in {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY}.
     *
     * @param genotype the genotype
     * @param contig the contig of the record holding the genotype, for error messages
     * @param position the position of the record holding the genotype, for error messages
     * @return the provenance, or {@code null} if the genotype carries none
     * @throws UserException.BadInput if the stored value cannot be decoded
     */
    static GenotypePriorProvenance fromGenotype(final Genotype genotype, final String contig, final int position) {
        Utils.nonNull(genotype, "genotype");
        final Object value = genotype.getExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY);
        if (value == null || VCFConstants.MISSING_VALUE_v4.equals(value.toString())) {
            return null;
        }
        try {
            return decode(value.toString());
        } catch (final IllegalArgumentException e) {
            throw new UserException.BadInput(String.format("Sample %s at %s:%d has a malformed %s value '%s': %s",
                    genotype.getSampleName(), contig, position, GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY, value, e.getMessage()) + PRIOR_OPTIONS_NOTE);
        }
    }

    /**
     * @param genotype a genotype
     * @return true if the genotype carries priors (PG or PRI) but no {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY}
     */
    static boolean hasPriorsWithoutProvenance(final Genotype genotype) {
        return !genotype.hasExtendedAttribute(GATKVCFConstants.GENOTYPE_PRIOR_PROVENANCE_KEY)
                && (hasValues(genotype, GATKVCFConstants.GENOTYPE_PRIOR_KEY) || hasValues(genotype, GATKVCFConstants.DRAGEN_GENOTYPE_PRIOR_KEY));
    }

    /**
     * Maps these priors onto the genotypes of a site.
     *
     * <p>The site's alleles must be alleles of the record, in any order and any subset, possibly with the reference
     * allele extended by bases appended to every allele (as when merging with a longer reference allele). An allele the
     * record did not have makes the priors inapplicable: the record's {@code <NON_REF>} prior covers all unseen alleles
     * together and does not determine the prior of any one of them.</p>
     *
     * @param siteStart the start of the site
     * @param siteAlleles the site's alleles, reference first
     * @param sitePloidy the sample's ploidy at the site
     * @return the log10 priors in VCF genotype order for {@code siteAlleles}, or the reason they do not apply
     */
    Mapping mapTo(final int siteStart, final List<Allele> siteAlleles, final int sitePloidy) {
        Utils.nonEmpty(siteAlleles, "siteAlleles");
        if (siteStart != start) {
            return Mapping.inapplicable(Inapplicability.DIFFERENT_START);
        }
        if (sitePloidy != ploidy) {
            return Mapping.inapplicable(Inapplicability.DIFFERENT_PLOIDY);
        }
        final byte[] siteReference = siteAlleles.get(0).getBases();
        final byte[] recordReference = alleles.get(0).getBases();
        if (!siteAlleles.get(0).isReference() || siteReference.length < recordReference.length
                || !Arrays.equals(siteReference, 0, recordReference.length, recordReference, 0, recordReference.length)) {
            return Mapping.inapplicable(Inapplicability.DIFFERENT_REFERENCE);
        }
        final byte[] referenceExtension = Arrays.copyOfRange(siteReference, recordReference.length, siteReference.length);

        final int[] siteToRecordAllele = new int[siteAlleles.size()];
        siteToRecordAllele[0] = 0;
        for (int i = 1; i < siteAlleles.size(); i++) {
            final int recordIndex = recordAlleleIndex(siteAlleles.get(i), referenceExtension);
            if (recordIndex < 0) {
                return Mapping.inapplicable(Inapplicability.NEW_ALLELE);
            }
            siteToRecordAllele[i] = recordIndex;
        }
        if (Arrays.stream(siteToRecordAllele).distinct().count() != siteToRecordAllele.length) {
            return Mapping.inapplicable(Inapplicability.NEW_ALLELE);  // two site alleles map to one record allele
        }

        final int[] siteToRecordGenotype = GenotypeIndexCalculator.newToOldGenotypeMap(ploidy, siteToRecordAllele);
        final double[] log10Priors = new double[siteToRecordGenotype.length];
        for (int i = 0; i < log10Priors.length; i++) {
            log10Priors[i] = -0.1 * phredPriors[siteToRecordGenotype[i]];
        }
        return Mapping.applicable(log10Priors);
    }

    /**
     * The index among this record's alternate alleles of the allele that a site allele represents, or -1 if none.
     * Called alleles match after the record allele is extended by the bases that extend the reference.
     */
    private int recordAlleleIndex(final Allele siteAllele, final byte[] referenceExtension) {
        if (siteAllele.isReference()) {
            return -1;
        }
        final boolean notExtended = siteAllele.isSymbolic() || siteAllele.equals(Allele.SPAN_DEL);
        int match = -1;
        for (int i = 1; i < alleles.size(); i++) {
            final Allele recordAllele = alleles.get(i);
            final boolean matches;
            if (notExtended) {
                // symbolic alleles and the spanning deletion are not extended, so they must match exactly
                matches = recordAllele.equals(siteAllele);
            } else {
                matches = recordAllele.isCalled() && !recordAllele.isSymbolic() && !recordAllele.equals(Allele.SPAN_DEL)
                        && Allele.extend(recordAllele, referenceExtension).equals(siteAllele);
            }
            if (matches) {
                if (match >= 0) {
                    return -1;  // ambiguous
                }
                match = i;
            }
        }
        return match;
    }

    /** @return the value to store in {@value GATKVCFConstants#GENOTYPE_PRIOR_PROVENANCE_KEY} */
    public String encode() {
        final String encodedAlleles = alleles.stream()
                .map(a -> percentEncode(a.getDisplayString()))
                .collect(Collectors.joining(String.valueOf(VALUE_SEPARATOR)));
        final String encodedPriors = Arrays.stream(phredPriors)
                .mapToObj(Double::toString)
                .collect(Collectors.joining(String.valueOf(VALUE_SEPARATOR)));
        return String.join(String.valueOf(FIELD_SEPARATOR), VERSION, Integer.toString(start), Integer.toString(ploidy), encodedAlleles, encodedPriors);
    }

    /**
     * Decodes a value written by {@link #encode()}.
     *
     * @throws IllegalArgumentException if the value is malformed or of another version
     */
    @VisibleForTesting
    static GenotypePriorProvenance decode(final String value) {
        Utils.nonNull(value, "value");
        final String[] fields = value.split("\\" + FIELD_SEPARATOR, -1);
        Utils.validateArg(fields.length == NUMBER_OF_FIELDS, () -> "expected " + NUMBER_OF_FIELDS + " '" + FIELD_SEPARATOR + "'-separated fields but found " + fields.length);
        Utils.validateArg(VERSION.equals(fields[0]), () -> "unsupported version '" + fields[0] + "'; expected " + VERSION);
        final int start = parseInt(fields[1], "start");
        final int ploidy = parseInt(fields[2], "ploidy");
        final String[] alleleStrings = fields[3].split(String.valueOf(VALUE_SEPARATOR), -1);
        final List<Allele> alleles = new ArrayList<>(alleleStrings.length);
        for (int i = 0; i < alleleStrings.length; i++) {
            final String allele = percentDecode(alleleStrings[i]);
            Utils.validateArg(!allele.isEmpty(), "alleles may not be empty");
            alleles.add(Allele.create(allele, i == 0));
        }
        final String[] priorStrings = fields[4].split(String.valueOf(VALUE_SEPARATOR), -1);
        final double[] priors = new double[priorStrings.length];
        for (int i = 0; i < priorStrings.length; i++) {
            try {
                priors[i] = Double.parseDouble(priorStrings[i]);
            } catch (final NumberFormatException e) {
                throw new IllegalArgumentException("genotype prior '" + priorStrings[i] + "' is not a number", e);
            }
        }
        return new GenotypePriorProvenance(start, ploidy, alleles, priors);
    }

    /**
     * Reads a Number=G Phred-scaled field from a genotype.
     *
     * @return the values, or {@code null} if the field is absent or entirely missing
     * @throws UserException.BadInput if only some values are missing or a value is not a number
     */
    private static double[] readPhredValues(final VariantContext vc, final Genotype genotype, final String key) {
        final List<Object> values = VariantContextGetters.attributeToList(genotype.getExtendedAttribute(key));
        if (values.isEmpty() || values.stream().allMatch(GenotypePriorProvenance::isMissing)) {
            return null;
        }
        final double[] result = new double[values.size()];
        for (int i = 0; i < result.length; i++) {
            final Object value = values.get(i);
            if (isMissing(value)) {
                throw new UserException.BadInput(String.format("Sample %s at %s:%d has missing values among its %s values: %s",
                        genotype.getSampleName(), vc.getContig(), vc.getStart(), key, values) + PRIOR_OPTIONS_NOTE);
            }
            try {
                result[i] = value instanceof Number ? ((Number) value).doubleValue() : Double.parseDouble(value.toString().trim());
            } catch (final NumberFormatException e) {
                throw new UserException.BadInput(String.format("Sample %s at %s:%d has a %s value that is not a number: '%s'",
                        genotype.getSampleName(), vc.getContig(), vc.getStart(), key, value) + PRIOR_OPTIONS_NOTE);
            }
        }
        return result;
    }

    private static boolean hasValues(final Genotype genotype, final String key) {
        return VariantContextGetters.attributeToList(genotype.getExtendedAttribute(key)).stream().anyMatch(v -> !isMissing(v));
    }

    private static boolean isMissing(final Object value) {
        return value == null || VCFConstants.MISSING_VALUE_v4.equals(value.toString().trim());
    }

    private static boolean isValidPrior(final double phredPrior) {
        return Double.isFinite(phredPrior) && phredPrior >= 0;
    }

    private static boolean agree(final double[] a, final double[] b) {
        if (a.length != b.length) {
            return false;
        }
        for (int i = 0; i < a.length; i++) {
            if (!(Math.abs(a[i] - b[i]) <= MAX_PRIOR_DIFFERENCE_BETWEEN_KEYS)) {
                return false;
            }
        }
        return true;
    }

    private static int parseInt(final String value, final String name) {
        try {
            return Integer.parseInt(value);
        } catch (final NumberFormatException e) {
            throw new IllegalArgumentException(name + " '" + value + "' is not an integer", e);
        }
    }

    /** Percent-encodes, as UTF-8, every character that is not safe in a VCF FORMAT value or this encoding. */
    @VisibleForTesting
    static String percentEncode(final String value) {
        return ALLELE_ESCAPER.escape(value);
    }

    private static boolean isSafe(final char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || SAFE_ALLELE_CHARACTERS.indexOf(c) >= 0;
    }

    @VisibleForTesting
    static String percentDecode(final String value) {
        final byte[] decoded = new byte[value.length()];
        int length = 0;
        for (int i = 0; i < value.length(); i++) {
            final char c = value.charAt(i);
            if (c == '%') {
                Utils.validateArg(i + 2 < value.length(), () -> "truncated percent-encoding in '" + value + "'");
                final int high = hexDigit(value.charAt(i + 1));
                final int low = hexDigit(value.charAt(i + 2));
                Utils.validateArg(high >= 0 && low >= 0, () -> "invalid percent-encoding in '" + value + "'");
                decoded[length++] = (byte) ((high << 4) | low);
                i += 2;
            } else {
                Utils.validateArg(isSafe(c), () -> "unexpected character '" + c + "' in '" + value + "'");
                decoded[length++] = (byte) c;
            }
        }
        return new String(decoded, 0, length, StandardCharsets.UTF_8);
    }

    /** @return the value of an ASCII hexadecimal digit, or -1 if {@code c} is not one */
    private static int hexDigit(final char c) {
        return c < 128 ? Character.digit(c, 16) : -1;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (!(o instanceof GenotypePriorProvenance)) return false;
        final GenotypePriorProvenance that = (GenotypePriorProvenance) o;
        return start == that.start && ploidy == that.ploidy && alleles.equals(that.alleles) && Arrays.equals(phredPriors, that.phredPriors);
    }

    @Override
    public int hashCode() {
        return 31 * Objects.hash(start, ploidy, alleles) + Arrays.hashCode(phredPriors);
    }

    @Override
    public String toString() {
        return encode();
    }

    /** The result of {@link #mapTo}: either log10 priors for the site's genotypes, or why there are none. */
    static final class Mapping {
        private final double[] log10Priors;
        private final Inapplicability inapplicability;

        private Mapping(final double[] log10Priors, final Inapplicability inapplicability) {
            this.log10Priors = log10Priors;
            this.inapplicability = inapplicability;
        }

        static Mapping applicable(final double[] log10Priors) {
            return new Mapping(log10Priors, null);
        }

        static Mapping inapplicable(final Inapplicability inapplicability) {
            return new Mapping(null, inapplicability);
        }

        boolean isApplicable() { return log10Priors != null; }

        /** @return the log10 priors in VCF genotype order for the site's alleles; only valid if {@link #isApplicable()} */
        double[] getLog10Priors() {
            Utils.validate(isApplicable(), "the priors do not apply to this site");
            return log10Priors.clone();
        }

        /** @return why the priors do not apply; only valid if not {@link #isApplicable()} */
        Inapplicability getInapplicability() {
            Utils.validate(!isApplicable(), "the priors apply to this site");
            return inapplicability;
        }
    }
}
