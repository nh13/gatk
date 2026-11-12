package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.ReadsPathDataSource;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.spark.AssemblyRegionArgumentCollection;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.broadinstitute.hellbender.testutils.IntegrationTestSpec;
import org.broadinstitute.hellbender.testutils.SamAssertionUtils;
import org.broadinstitute.hellbender.testutils.VariantContextTestUtils;
import org.broadinstitute.hellbender.tools.HaplotypeCallerSpark;
import org.broadinstitute.hellbender.tools.walkers.CombineGVCFs;
import org.broadinstitute.hellbender.tools.walkers.GenotypeGVCFs;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ReblockGVCF;
import org.broadinstitute.hellbender.tools.walkers.variantutils.ValidateVariants;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@Test(groups = {"variantcalling"})
public class HaplotypeCallerIntegrationTest extends CommandLineProgramTest {

    // If true, update the expected outputs in tests that assert an exact match vs. prior output,
    // instead of actually running the tests. Can be used with "./gradlew test -Dtest.single=HaplotypeCallerIntegrationTest"
    // to update all of the exact-match tests at once. After you do this, you should look at the
    // diffs in the new expected outputs in git to confirm that they are consistent with expectations.
    public static final boolean UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS = false;

    public static final String TEST_FILES_DIR = toolsTestDir + "haplotypecaller/";

    /*
     * Make sure that someone didn't leave the UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS toggle turned on
     */
    @Test
    public void assertThatExpectedOutputUpdateToggleIsDisabled() {
        Assert.assertFalse(UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS, "The toggle to update expected outputs should not be left enabled");
    }

    @DataProvider(name="HaplotypeCallerTestInputs")
    public Object[][] getHaplotypCallerTestInputs() {
        return new Object[][] {
                {NA12878_20_21_WGS_bam, b37_reference_20_21},
                {NA12878_20_21_WGS_cram, b37_reference_20_21}
        };
    }
    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that minimap2 data are supported and consistent with past results
     */
    @Test
    public void testVCFModeMinimap2IsConsistentWithPastResults() throws Exception {
        final String inputFileName = NA12878_20_21_WGS_mmp2_bam;
        final String referenceFileName = b37_reference_20_21;
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testVCFMode.mmp2.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in JunctionTree mode we're consistent with past JunctionTree results (over non-complicated data)
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testLinkedDebruijnModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testLinkedDebruijnModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testLinkedDebruijnMode.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--"+ReadThreadingAssemblerArgumentCollection.LINKED_DE_BRUIJN_GRAPH_LONG_NAME,
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /**
     * Test that exposure and modification of Smith-Waterman parameters that were previously constants is consistent with past results (over non-complicated data).
     * Expected outputs were generated by changing previously hardcoded constants in the previous commit to the exposed values below.
     * See https://github.com/broadinstitute/gatk/pull/6885.
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testExposureOfSmithWatermanParametersIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testExposureOfSmithWatermanParametersIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testExposureOfSmithWatermanParameters.HC.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_MATCH_VALUE_LONG_NAME, "1",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_MISMATCH_PENALTY_LONG_NAME, "-2",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_GAP_OPEN_PENALTY_LONG_NAME, "-3",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_DANGLING_END_GAP_EXTEND_PENALTY_LONG_NAME, "-4",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MATCH_VALUE_LONG_NAME, "5",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MISMATCH_PENALTY_LONG_NAME, "-6",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_OPEN_PENALTY_LONG_NAME, "-7",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_EXTEND_PENALTY_LONG_NAME, "-8",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_MATCH_VALUE_LONG_NAME, "9",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_MISMATCH_PENALTY_LONG_NAME, "-10",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_OPEN_PENALTY_LONG_NAME, "-11",
                "--" + AssemblyBasedCallerArgumentCollection.SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_EXTEND_PENALTY_LONG_NAME, "-12",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     *
     * Test currently throws an exception due to lack of support for allele-specific annotations in VCF mode
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = UserException.class)
    public void testVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        //NOTE: AlleleSpecific support in the VCF mode is implemented but bogus for now.
        // This test should not be treated as a strict check of correctness.
        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk4.alleleSpecific.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testVCFModeIsConcordantWithGATK3_8Results(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3Results", ".vcf");
        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in VCF mode is < 99% (" +  concordance + ")");
    }

    /*
     * Test that the this version of DRAGEN-GATK has not changed relative to the last version with the recommended arguments enabled
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testDRAGENDefaultArgIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testDragenSimpleModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk4.DRAGEN.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--dragen-mode",
                // STRE arguments
                "--dragstr-params-path", TEST_FILES_DIR+"example.dragstr-params.txt",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
        };

        runCommandLine(args);
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that the this version of DRAGEN-GATK has not changed relative to the last version with the recommended arguments enabled
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testDRAGENGATKModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testDRAGENGATKModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk4.DRAGEN.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                // FRD arguments
                "--apply-frd", "--transform-dragen-mapping-quality", "--mapping-quality-threshold-for-genotyping", "1", "--disable-cap-base-qualities-to-map-quality", "--minimum-mapping-quality", "1",
                // BQD arguments
                "--apply-bqd",  "--soft-clip-low-quality-ends",
                // Dynamic read disqualification arguments"
                "--enable-dynamic-read-disqualification-for-genotyping", "--expected-mismatch-rate-for-read-disqualification", "0.03",
                // Genotyper arguments
                "--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES",  "--standard-min-confidence-threshold-for-calling", "3", "--use-posteriors-to-calculate-qual",
                // STRE arguments
                "--dragstr-params-path", TEST_FILES_DIR+"example.dragstr-params.txt",
                // misc arguments
                "--enable-legacy-graph-cycle-detection", "--padding-around-indels", "150",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "1",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
        };

        runCommandLine(args);
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }


    /*
     * Test that in VCF mode we're >= 99% concordant with GATK3.8 results
     *
     * Test currently throws an exception due to lack of support for allele-specific annotations in VCF mode
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = UserException.class)
    public void testVCFModeIsConcordantWithGATK3_8ResultsAlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConcordantWithGATK3.8ResultsAlleleSpecificAnnotations", ".vcf");

        //Created by running
        //java -jar gatk.3.8-4-g7b0250253f.jar  -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // --out expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS VCF mode is < 99% (" +  concordance + ")");
    }

    /*
     * Test that in GVCF mode we're consistent with past GATK4 results
     */
    private File runHaplotypeCallerWithMaxEffectiveDepth(final String inputFileName, final String referenceFileName, final String maxEffectiveDepth, final String outputPrefix) {
        Utils.resetRandomGenerator();
        final File output = createTempFile(outputPrefix, ".vcf");
        final List<String> args = new ArrayList<>(Arrays.asList(
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"));
        if (maxEffectiveDepth != null) {
            args.addAll(Arrays.asList("--" + AssemblyRegionArgumentCollection.MAX_EFFECTIVE_DEPTH_LONG_NAME, maxEffectiveDepth,
                    "--" + AssemblyRegionArgumentCollection.MAX_EFFECTIVE_DEPTH_WINDOW_LONG_NAME, "500"));
        }
        runCommandLine(args);
        return output;
    }

    /**
     * A depth cap far above any depth in the input must leave the output bit-identical to the run without it:
     * no window can discard a read, so the downsampler never touches the read stream.
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testMaxEffectiveDepthAboveAllDepthsIsANoOp(final String inputFileName, final String referenceFileName) throws Exception {
        final File uncapped = runHaplotypeCallerWithMaxEffectiveDepth(inputFileName, referenceFileName, null, "testMaxEffectiveDepth.uncapped");
        final File capped = runHaplotypeCallerWithMaxEffectiveDepth(inputFileName, referenceFileName, "100000", "testMaxEffectiveDepth.highCap");

        IntegrationTestSpec.assertEqualTextFiles(capped, uncapped);
    }

    /**
     * A depth cap below the input's depth must be honored: the run completes, emits calls, and the capped
     * read set changes the output relative to the uncapped run.
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testMaxEffectiveDepthBelowInputDepthChangesOutput(final String inputFileName, final String referenceFileName) throws Exception {
        final File uncapped = runHaplotypeCallerWithMaxEffectiveDepth(inputFileName, referenceFileName, null, "testMaxEffectiveDepth.uncapped");
        final File capped = runHaplotypeCallerWithMaxEffectiveDepth(inputFileName, referenceFileName, "5", "testMaxEffectiveDepth.lowCap");

        final List<String> cappedCalls = VariantContextTestUtils.getVariantContexts(capped).stream()
                .map(vc -> vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles() + ":" + vc.getGenotype(0).getGenotypeString())
                .collect(Collectors.toList());
        final List<String> uncappedCalls = VariantContextTestUtils.getVariantContexts(uncapped).stream()
                .map(vc -> vc.getContig() + ":" + vc.getStart() + ":" + vc.getAlleles() + ":" + vc.getGenotype(0).getGenotypeString())
                .collect(Collectors.toList());
        Assert.assertFalse(cappedCalls.isEmpty(), "the capped run emitted no variants");
        Assert.assertNotEquals(cappedCalls, uncappedCalls, "a cap of 5 on ~30x data should change the calls");
    }

    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk4.g.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }


    /*
     * Minimal test that the non-seq graph haplotype detection code is equivalent using either seq graphs or kmer graphs
     *
     *  NOTE: --disable-sequence-graph-simplification is currently an experimental feature that does not directly match with
     *        the regular HaplotypeCaller. Specifically the haplotype finding code does not perform correctly at complicated
     *        sites, which is illustrated by this test. Use this mode at your own risk.
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", enabled = false)
    public void testGVCFModeIsConsistentWithPastResultsUsingKmerGraphs(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testGVCFMode.gatk4.g.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-ERC", "GVCF",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--disable-sequence-graph-simplification",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in GVCF mode we're consistent with past GATK4 results using AS_ annotations
     *
     * Updated on 09/01/17 to account for changes to AS_RankSum annotations the annotations were checked against GATK3
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testGVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults_AlleleSpecificAnnotations", ".g.vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk4.alleleSpecific.g.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in GVCF mode we're >= 99% concordant with GATK3 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs", enabled=false) //disabled after reference confidence change in #5172
    public void testGVCFModeIsConcordantWithGATK3_8Results(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3Results", ".g.vcf");
        //Created by running:
        //java -jar  gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf -G StandardHC -G Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.3.8-4-g7b0250253f.g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in GVCF mode is < 99% (" +  concordance + ")");
    }

    /*
     * Test that the this version of DRAGEN-GATK has not changed relative to the last version with the recommended arguments enabled
     */
    @Test
    public void testFRDBQDDRAGENGATKOnDRAGENProcessedFile() throws Exception {
        Utils.resetRandomGenerator();
        final String inputFileName = largeFileTestDir + "DRAGENExampleBamSites.bam";
        final String intervals = TEST_FILES_DIR + "DRAGENTestSites.bed";

        final File output = createTempFile("testFRDBQDDRAGENGATKOnDRAGENProcessedFile", ".vcf");
        final File expected = new File(TEST_FILES_DIR + "expected.testVCFMode.gatk4.FRDBQD.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", b37Reference,
                "-L", intervals,
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                // FRD arguments
                "--apply-frd", "--transform-dragen-mapping-quality", "--mapping-quality-threshold-for-genotyping", "1", "--disable-cap-base-qualities-to-map-quality", "--minimum-mapping-quality", "1",
                // BQD arguments
                "--apply-bqd",  "--soft-clip-low-quality-ends",
                // Dynamic read disqualification arguments"
                "--enable-dynamic-read-disqualification-for-genotyping", "--expected-mismatch-rate-for-read-disqualification", "0.03",
                // misc arguments
                "--enable-legacy-graph-cycle-detection", "--padding-around-indels", "150",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "1",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
        };

        runCommandLine(args);
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }


    @Test(dataProvider="HaplotypeCallerTestInputs", enabled=false) //disabled after reference confidence change in #5172
    public void testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();
        final File output = createTempFile("testGVCFModeIsConcordantWithGATK3_8AlelleSpecificResults", ".g.vcf");

        //Created by running:
        // java -jar gatk.3.8-4-g7b0250253f.jar -T HaplotypeCaller \
        // -I ./src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam \
        // -R src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10000000-10100000 \
        // -ERC GVCF --out expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf -G StandardHC -G Standard -G AS_Standard \
        // --disableDithering --no_cmdline_in_header  -dt NONE --maxReadsInRegionPerSample 100000000 --minReadsPerAlignmentStart 100000 \
        // -pairHMM VECTOR_LOGLESS_CACHING
        final File gatk3Output = new File(TEST_FILES_DIR + "expected.testGVCFMode.gatk3.8-4-g7b0250253f.alleleSpecific.g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-G", "StandardAnnotation",
                "-G", "StandardHCAnnotation",
                "-G", "AS_StandardAnnotation",
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-pairHMM", "AVX_LOGLESS_CACHING",
        };

        runCommandLine(args);

        final double concordance = calculateConcordance(output, gatk3Output);
        Assert.assertTrue(concordance >= 0.99, "Concordance with GATK 3.8 in AS GVCF mode is < 99% (" +  concordance + ")");
    }

    @Test
    public void testGVCFModeGenotypePosteriors() throws Exception {
        Utils.resetRandomGenerator();

        final String inputFileName = NA12878_20_21_WGS_bam;
        final String referenceFileName =b37_reference_20_21;

        final File output = createTempFile("testGVCFModeIsConsistentWithPastResults", ".g.vcf");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + GenotypeCalculationArgumentCollection.SUPPORTING_CALLSET_LONG_NAME,
                    largeFileTestDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf",
                "--" + GenotypeCalculationArgumentCollection.NUM_REF_SAMPLES_LONG_NAME, "2500",
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
        };

        runCommandLine(args);

        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        for (final VariantContext vc : results.getRight()) {
            final Genotype g = vc.getGenotype(0);
            if (g.hasDP() && g.getDP() > 0 && g.hasGQ() && g.getGQ() > 0) {
                Assert.assertTrue(g.hasExtendedAttribute(GATKVCFConstants.PHRED_SCALED_POSTERIORS_KEY));
            }
            if (isGVCFReferenceBlock(vc) ) {
                Assert.assertTrue(!vc.hasAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
            }
            else if (!vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE)){      //there are some variants that don't have non-symbolic alts
                Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.GENOTYPE_PRIOR_KEY));
            }
        }
    }
    
    /*
     * Test that GQs are correct when the --floor-blocks argument is supplied
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testFloorGVCFBlocks(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testFloorGVCFBlocks", ".vcf");
        
        final List<String> requestedGqBands = Arrays.asList("10","20","30","40","50","60");
        
        final ArgumentsBuilder args = new ArgumentsBuilder().addInput(new File(inputFileName))
        .addReference(new File(referenceFileName))
        .addInterval(new SimpleInterval("20:10009880-10012631"))
        .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false)
        .add(AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2")
        .add("pairHMM", "AVX_LOGLESS_CACHING")
        .addFlag("floor-blocks")
        .add("ERC", "GVCF")
        .addOutput(output);
        requestedGqBands.forEach(x -> args.add("GQB",x));
        runCommandLine(args);
        
        final List<String> allGqBands = new ArrayList<>(requestedGqBands);
        allGqBands.add("99");
        allGqBands.add("0");
        
        //The interval here is big, so use a FeatureDataSource to limit memory usage
        try (final FeatureDataSource<VariantContext> actualVcs = new FeatureDataSource<>(output)) {
            actualVcs.forEach(vc -> {
                //sometimes there are calls with alt alleles that are genotyped hom ref and those don't get floored
                if (vc.hasAttribute("END") && vc.getGenotype(0).hasGQ()) {
                    Assert.assertTrue(allGqBands.contains(Integer.toString(vc.getGenotype(0).getGQ())));
                }
            });
        }
    }

    // force calling bam (File) vcf (File) and intervals (String)
    @DataProvider(name="ForceCallingInputs")
    public Object[][] getForceCallingInputs() {
        return new Object[][] {
                {NA12878_20_21_WGS_bam, new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles.vcf"), "20:10000000-10010000"},
                {NA12878_20_21_WGS_bam, new File(toolsTestDir, "mutect/gga_mode.vcf"), "20:9998500-10010000"},
                {NA12878_20_21_WGS_bam, new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles_ExtremeLengthDeletion.vcf"), "20:9998500-10010000"} // This is designed to test https://github.com/broadinstitute/gatk/issues/8675, which stemmed from an edge case in the force calling logic where a deletion allele that is longer than the assembly window padding spans into the assembly window. This tests that we do not see an exception in this case.
        };
    }

    @Test(dataProvider = "ForceCallingInputs")
    public void testForceCalling(final String bamPath, final File forceCallingVcf, final String intervalString) throws IOException {
        Utils.resetRandomGenerator();
        final File output = createTempFile("testGenotypeGivenAllelesMode", ".vcf");

        final String[] args = {
                "-I", bamPath,
                "-R", b37Reference,
                "-L", intervalString,
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--" + AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, forceCallingVcf.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-ERC", "GVCF"
        };

        runCommandLine(args);

        final Map<Integer, List<Allele>> altAllelesByPosition = VariantContextTestUtils.streamVcf(output)
                .collect(Collectors.toMap(vc -> vc.getStart(), vc -> vc.getAlternateAlleles()));
        for (final VariantContext vc : new FeatureDataSource<VariantContext>(forceCallingVcf)) {
            final List<Allele> altAllelesAtThisLocus = altAllelesByPosition.get(vc.getStart());
            vc.getAlternateAlleles().stream().filter(a-> a.length() > 0 && BaseUtils.isNucleotide(a.getBases()[0])).forEach(a -> Assert.assertTrue(altAllelesAtThisLocus.contains(a)));
        }
    }

    @Test
    public void testForceCallingNotProducingNoCalls() throws IOException {
        Utils.resetRandomGenerator();
        final String bamPath = NA12878_20_21_WGS_bam;
        final File forceCallingVcf = new File(TEST_FILES_DIR, "testGenotypeGivenAllelesMode_givenAlleles.vcf");
        final String intervalString = "20:10000000-10010000";
        final File expected = new File(TEST_FILES_DIR, "expected.testForceCallingNotProducingNoCalls.gatk4.vcf");


        final File output = createTempFile("testGenotypeGivenAllelesMode", ".vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        // Prior to https://github.com/broadinstitute/gatk/pull/7740 the output for this test would be littered with 0/0 reference variants with empty alt alleles. This test is intended to enshrine the new behavior.
        // For an example of what the output used to look at see the issue https://github.com/broadinstitute/gatk/issues/7741.
        final String[] args = {
                "-I", bamPath,
                "-R", b37Reference,
                "-L", intervalString,
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--" + AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, forceCallingVcf.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + GenotypeCalculationArgumentCollection.CALL_CONFIDENCE_LONG_NAME, "1000" // This is intended to be an absurdly high number that forces most alleles to fail to demonstrate a bug that caused calls to be made in the vicinity of forced alleles
        };

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }


    // regression test for https://github.com/broadinstitute/gatk/issues/6495, where a mistake in assembly region trimming
    // caused variants in one-base or similarly short intervals to cause reads to be trimmed to one base long, yielding
    // no calls.
    @Test
    public void testSingleBaseIntervals() throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(b37Reference)
                .addInput(NA12878_20_21_WGS_bam)
                .addInterval("20:10000117")
                .addInterval("20:10000439")
                .addInterval("20:10000694")
                .addInterval("20:10001019")
                .addOutput(output);

        runCommandLine(args);

        Assert.assertEquals(VariantContextTestUtils.getVariantContexts(output).size(), 4);
    }

    @Test
    public void testBamoutProducesReasonablySizedOutput() {
        final Path bamOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".bam").toPath();
        innerTestBamoutProducesReasonablySizedOutput(bamOutput);
    }

    @Test(groups={"bucket"})
    public void testBamoutOnGcs() {
        final Path bamOutput = BucketUtils.getPathOnGcs(BucketUtils.getTempFilePath(
            getGCPTestStaging() + "testBamoutProducesReasonablySizedOutput", ".bam"));
        innerTestBamoutProducesReasonablySizedOutput(bamOutput);
    }

    private void innerTestBamoutProducesReasonablySizedOutput(Path bamOutput) {
        Utils.resetRandomGenerator();

        // We will test that when running with -bamout over the testInterval, we produce
        // a bam with a number of reads that is within 10% of what GATK3.5 produces with
        // -bamout over the same interval. This is just to test that we produce a reasonably-sized
        // bam for the region, not to validate the haplotypes, etc. We don't want
        // this test to fail unless there is a likely problem with -bamout itself (eg., empty
        // or truncated bam).
        final String testInterval = "20:10000000-10010000";
        final int gatk3BamoutNumReads = 5000;

        final File vcfOutput = createTempFile("testBamoutProducesReasonablySizedOutput", ".vcf");

        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(NA12878_20_21_WGS_bam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addOutput(new File(vcfOutput.getAbsolutePath()));
        argBuilder.add("L", testInterval);
        argBuilder.add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamOutput.toUri().toString());
        argBuilder.add("pairHMM", "AVX_LOGLESS_CACHING");
        argBuilder.add(AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME,2);

        runCommandLine(argBuilder.getArgsArray());

        try ( final ReadsDataSource bamOutReadsSource = new ReadsPathDataSource(bamOutput) ) {
            int actualBamoutNumReads = 0;
            for ( final GATKRead read : bamOutReadsSource ) {
                ++actualBamoutNumReads;
            }

            final int readCountDifference = Math.abs(actualBamoutNumReads - gatk3BamoutNumReads);
            Assert.assertTrue(((double)readCountDifference / gatk3BamoutNumReads) < 0.10,
                    "-bamout produced a bam with over 10% fewer/more reads than expected");
        }
    }

    @Test
    public void testSitesOnlyMode() {
        Utils.resetRandomGenerator();
        File out = createTempFile("GTStrippedOutput", "vcf");
        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", out.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + StandardArgumentDefinitions.SITES_ONLY_LONG_NAME,
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2"
        };

        runCommandLine(args);

        // Assert that the genotype field has been stripped from the file
        Pair<VCFHeader, List<VariantContext>> results = VariantContextTestUtils.readEntireVCFIntoMemory(out.getAbsolutePath());

        Assert.assertFalse(results.getLeft().hasGenotypingData());
        for (VariantContext v: results.getRight()) {
            Assert.assertFalse(v.hasGenotypes());
        }
    }

    @Test
    public void testForceActiveOption() throws Exception {
        Utils.resetRandomGenerator();
        final File out = createTempFile("GTStrippedOutput", "vcf");
        final File assemblyRegionOut = createTempFile("assemblyregions", ".igv");
        final File expectedAssemblyRegionOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionWithForceActiveRegions_assemblyregions.igv");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:1-5000",
                "-O", out.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + AssemblyRegionArgumentCollection.FORCE_ACTIVE_REGIONS_LONG_NAME, "true",
                "--" + AssemblyRegionArgumentCollection.ASSEMBLY_REGION_OUT_LONG_NAME, assemblyRegionOut.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };
        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(assemblyRegionOut, expectedAssemblyRegionOut);
    }

    @DataProvider(name="outputFileVariations")
    public Object[][] getOutputFileVariations() {
        return new Object[][]{
                // bamout index, bamout md5, vcf index, vcf md5
                { true, true, true, true },
                { true, false, true, false },
                { false, true, false, true },
                { false, false, false, false },
        };
    }

    @Test(dataProvider = "outputFileVariations")
    public void testOutputFileArgumentVariations(
            final boolean createBamoutIndex,
            final boolean createBamoutMD5,
            final boolean createVCFOutIndex,
            final boolean createVCFOutMD5) throws IOException {
        Utils.resetRandomGenerator();

        // run on small interval to test index/md5 outputs
        final String testInterval = "20:10000000-10001000";

        final File vcfOutput = createTempFile("testOutputFileArgumentVariations", ".vcf");
        final File bamOutput = createTempFile("testOutputFileArgumentVariations", ".bam");

        ArgumentsBuilder argBuilder = new ArgumentsBuilder();

        argBuilder.addInput(new File(NA12878_20_21_WGS_bam));
        argBuilder.addReference(new File(b37_reference_20_21));
        argBuilder.addOutput(new File(vcfOutput.getAbsolutePath()));
        argBuilder.add("L", testInterval);
        argBuilder.add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_SHORT_NAME, bamOutput.getAbsolutePath());
        argBuilder.add("pairHMM", "AVX_LOGLESS_CACHING");
        argBuilder.add(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_INDEX_LONG_NAME, createBamoutIndex);
        argBuilder.add(StandardArgumentDefinitions.CREATE_OUTPUT_BAM_MD5_LONG_NAME, createBamoutMD5);
        argBuilder.add(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_INDEX_LONG_NAME, createVCFOutIndex);
        argBuilder.add(StandardArgumentDefinitions.CREATE_OUTPUT_VARIANT_MD5_LONG_NAME, createVCFOutMD5);

        runCommandLine(argBuilder.getArgsArray());

        Assert.assertTrue(vcfOutput.exists(), "No VCF output file was created");

        // validate vcfout companion files
        final File vcfOutFileIndex = new File(vcfOutput.getAbsolutePath() + FileExtensions.TRIBBLE_INDEX);
        final File vcfOutFileMD5 = new File(vcfOutput.getAbsolutePath() + ".md5");
        Assert.assertEquals(vcfOutFileIndex.exists(), createVCFOutIndex, "The index file argument was not honored");
        Assert.assertEquals(vcfOutFileMD5.exists(), createVCFOutMD5, "The md5 file argument was not honored");

        // validate bamout companion files
        if (createBamoutIndex) {
            Assert.assertNotNull(SamFiles.findIndex(bamOutput));
        } else {
            Assert.assertNull(SamFiles.findIndex(bamOutput));
        }

        final File expectedBamoutMD5File = new File(bamOutput.getAbsolutePath() + ".md5");
        Assert.assertEquals(expectedBamoutMD5File.exists(), createBamoutMD5);

        // Check the output BAN header contains all of the inout BAM header Program Records (@PG)
        SamAssertionUtils.assertOutBamContainsInBamProgramRecords(new File(NA12878_20_21_WGS_bam), bamOutput);
    }

    @Test
    public void testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScores() throws IOException{
        final File testBAM = new File(TEST_FILES_DIR + "pretendTobeTetraPloidTetraAllelicSite.bam");
        final File output = createTempFile("testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScoresResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testHaplotypeCallerRemoveAltAlleleBasedOnHaptypeScores.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", testBAM.getAbsolutePath(),
                "-R", b37_reference_20_21,
                "-L", "20:11363580-11363600",
                "-O", outputPath,
                "-ploidy", "4",
                "--max-genotype-count", "15",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };
        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testCustomPloidyRegions(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("output", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testPloidyRegions.vcf");
        final File ploidyRegions = new File(TEST_FILES_DIR, "testPloidyRegions.bed");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ploidy", "3",
                "--ploidy-regions", ploidyRegions.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // Should fail due to non-integer values in name column of bed
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = IllegalArgumentException.class)
    public void testNonIntegerCustomPloidyRegions(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("output", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testPloidyRegions.vcf");
        final File ploidyRegions = new File(TEST_FILES_DIR, "testNonIntegerPloidyRegions.bed");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ploidy", "3",
                "--ploidy-regions", ploidyRegions.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // Should fail due to negative integer values in name column of bed
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = IllegalArgumentException.class)
    public void testNonPositiveCustomPloidyRegions(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("output", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testPloidyRegions.vcf");
        final File ploidyRegions = new File(TEST_FILES_DIR, "testNegativePloidyRegions.bed");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ploidy", "3",
                "--ploidy-regions", ploidyRegions.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // Should fail due to using chr20 contig name in ploidy-regions file when BAM sequence dictionary uses 20
    @Test(dataProvider="HaplotypeCallerTestInputs", expectedExceptions = UserException.class)
    public void testMismatchPloidyRegionsAndSequenceDictionary(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("output", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.testPloidyRegions.vcf");
        final File ploidyRegions = new File(TEST_FILES_DIR, "testMismatchPloidyRegions.bed");

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", output.getAbsolutePath(),
                "-ploidy", "3",
                "--ploidy-regions", ploidyRegions.getAbsolutePath(),
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected);
    }

    // Tests the somewhat ad hoc implementation of copying hcArgs instances in the engine, used for managing multiple ploidies used based on user inputs
    @Test
    public void testArgCollectionCopiesAreSynchronized() {
        Utils.resetRandomGenerator();

        HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
        hcArgs.applyBQD = !hcArgs.applyBQD;
        hcArgs.dontGenotype = !hcArgs.dontGenotype;
        hcArgs.standardArgs.annotateAllSitesWithPLs = !hcArgs.standardArgs.annotateAllSitesWithPLs;
        hcArgs.standardArgs.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = !hcArgs.standardArgs.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED;

        // Create copy with different ploidy and check tweaked values (and one control non-changed value) remain in sync
        HaplotypeCallerArgumentCollection hcArgsCopy = hcArgs.copyWithNewPloidy(3);
        Assert.assertEquals(hcArgs.applyBQD, hcArgsCopy.applyBQD);
        Assert.assertEquals(hcArgs.dontGenotype, hcArgsCopy.dontGenotype);
        Assert.assertEquals(hcArgs.standardArgs.annotateAllSitesWithPLs, hcArgsCopy.standardArgs.annotateAllSitesWithPLs);
        Assert.assertEquals(hcArgs.standardArgs.CONTAMINATION_FRACTION, hcArgsCopy.standardArgs.CONTAMINATION_FRACTION);
        Assert.assertEquals(hcArgs.standardArgs.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED, hcArgsCopy.standardArgs.genotypeArgs.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED);

        // Modify other values of copy and check they are now distinct, i.e. a "deep copy" was created above
        hcArgsCopy.applyBQD = !hcArgsCopy.applyBQD;
        Assert.assertNotEquals(hcArgs.applyBQD, hcArgsCopy.applyBQD);
    }

    // test that ReadFilterLibrary.NON_ZERO_REFERENCE_LENGTH_ALIGNMENT removes reads that consume zero reference bases
    // e.g. read name HAVCYADXX150109:1:2102:20528:2129 with cigar 23S53I
    @Test
    public void testReadsThatConsumeZeroReferenceReads() throws Exception {
        final String CONSUMES_ZERO_REFERENCE_BASES = publicTestDir + "org/broadinstitute/hellbender/tools/mutect/na12878-chr20-consumes-zero-reference-bases.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", CONSUMES_ZERO_REFERENCE_BASES,
                "-R", b37_reference_20_21,
                "-O", outputVcf.getAbsolutePath()
        };
        runCommandLine(args);
    }

    // Test fix for https://github.com/broadinstitute/gatk/issues/3466
    // This specifically tests the case of a read that ends up as empty
    // after a call to ReadClipper.hardClipSoftClippedBases()
    @Test
    public void testCompletelyClippedReadNearStartOfContig() throws Exception {
        final File testCaseFilesDir = new File(TEST_FILES_DIR, "issue3466_gatk_cigar_error");
        final File output = createTempFile("testCompletelyClippedReadNearStartOfContig", ".vcf");
        final File expected = new File(testCaseFilesDir, "expected_output_gatk3.vcf");

        final String[] args = {
                "-I", new File(testCaseFilesDir, "culprit.bam").getAbsolutePath(),
                "-R", new File(testCaseFilesDir, "GRCh37_MTonly.fa").getAbsolutePath(),
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertEquals(calculateConcordance(output, expected), 1.0);
    }

    // Test fix for https://github.com/broadinstitute/gatk/issues/3845
    // This specifically tests the case of a read at the start of a contig (position == 1)
    // that becomes completely clipped after a call to ReadClipper.revertSoftClippedBases()
    @Test
    public void testCompletelyClippedReadNearStartOfContig_revertSoftClipped() throws Exception {
        final File testCaseFilesDir = new File(TEST_FILES_DIR, "issue3845_revertSoftClip_bug");
        final File output = createTempFile("testCompletelyClippedReadNearStartOfContig_revertSoftClipped", ".vcf");
        final File expected = new File(testCaseFilesDir, "expected_testCompletelyClippedReadNearStartOfContig_revertSoftClipped_gatk4.vcf");

        final String[] args = {
                "-I", new File(testCaseFilesDir, "issue3845_bug.bam").getAbsolutePath(),
                "-R", new File(publicTestDir, "Homo_sapiens_assembly38_chrM_only.fasta").getAbsolutePath(),
                "-L", "chrM",
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        Assert.assertEquals(calculateConcordance(output, expected), 1.0);
    }

    @Test
    public void testAssemblyRegionAndActivityProfileOutput() throws Exception {
        final File output = createTempFile("testAssemblyRegionAndActivityProfileOutput", ".vcf");
        final File assemblyRegionOut = createTempFile("testAssemblyRegionAndActivityProfileOutput_assemblyregions", ".igv");
        final File expectedAssemblyRegionOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionAndActivityProfileOutput_assemblyregions.igv");
        final File expectedActivityProfileOut = new File(TEST_FILES_DIR, "expected.testAssemblyRegionAndActivityProfileOutput_activityprofile.igv");

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10003000",
                "-O", output.getAbsolutePath(),
                "-pairHMM", "AVX_LOGLESS_CACHING",
                "--" + AssemblyRegionArgumentCollection.ASSEMBLY_REGION_OUT_LONG_NAME, assemblyRegionOut.getAbsolutePath()
        };

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(assemblyRegionOut, expectedAssemblyRegionOut);
    }

    /*
    * Test that the min_base_quality_score parameter works
    */
    @Test
    public void testMinBaseQualityScore() throws Exception {
        Utils.resetRandomGenerator();

        final File outputAtLowThreshold = createTempFile("output", ".vcf");
        final File outputAtHighThreshold = createTempFile("output", ".vcf");

        final String[] lowThresholdArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", outputAtLowThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "20"
        };

        runCommandLine(lowThresholdArgs);

        final String[] highThresholdArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", outputAtHighThreshold.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.MIN_BASE_QUALITY_SCORE_LONG_NAME, "30"
        };

        runCommandLine(highThresholdArgs);

        try (final FeatureDataSource<VariantContext> lowThresholdSource = new FeatureDataSource<VariantContext>(outputAtLowThreshold);
             final FeatureDataSource<VariantContext> highThresholdSource = new FeatureDataSource<VariantContext>(outputAtHighThreshold)) {
            final List<VariantContext> variantsWithLowThreshold =
                    StreamSupport.stream(lowThresholdSource.spliterator(), false).collect(Collectors.toList());

            final List<VariantContext> variantsWithHighThreshold =
                    StreamSupport.stream(highThresholdSource.spliterator(), false).collect(Collectors.toList());

            final Set<Integer> lowStarts = variantsWithLowThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            final Set<Integer> highStarts = variantsWithHighThreshold.stream().map(VariantContext::getStart).collect(Collectors.toSet());
            lowStarts.removeAll(highStarts);
            final List<Integer> diff = lowStarts.stream().sorted().collect(Collectors.toList());
            Assert.assertEquals(diff, Arrays.asList(10002458, 10008758, 10009842));
        }

    }

    @DataProvider
    public Object[][] getContaminationCorrectionTestData() {
        // This bam is a snippet of GATKBaseTest.NA12878_20_21_WGS_bam artificially contaminated
        // at 15% with reads from another sample. We use such a high contamination fraction to ensure
        // that calls will be noticeably different vs. running without -contamination at all.
        final String contaminatedBam15Percent = largeFileTestDir + "contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam";

        final SimpleInterval traversalInterval = new SimpleInterval("20", 10100000, 10150000);

        // File specifying the per-sample contamination for contaminatedBam15Percent
        final String contaminationFile = TEST_FILES_DIR + "contamination_file_for_CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT";

        // Expected calls from GATK 3.x on the contaminated bam running with -contamination
        // Created in GATK 3.8-1-1-gdde23f56a6 using the command:
        // java -jar target/GenomeAnalysisTK.jar -T HaplotypeCaller -I ../hellbender/src/test/resources/large/contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam -R ../hellbender/src/test/resources/large/human_g1k_v37.20.21.fasta -L 20:10100000-10150000 -contamination 0.15 -o ../hellbender/src/test/resources/org/broadinstitute/hellbender/tools/haplotypecaller/expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.15PCT.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.vcf
        final String expectedGATK3ContaminationCorrectedCallsVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.15PCT.20.10100000-10150000.gatk3.8-1-1-gdde23f56a6.vcf";

        // Expected calls from GATK ~4.0.8.1 with indel reference confidence fix
        final String expectedGATK4ContaminationCorrectedCallsGVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.15PCT.20.10100000-10150000.postIndelRefConfUpdate.g.vcf";

        // Expected calls from GATK 4 on the uncontaminated bam (VCF mode)
        final String expectedGATK4UncontaminatedCallsVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.calls.20.10100000-10150000.vcf";

        // Expected calls from GATK 4 on the uncontaminated bam (GVCF mode)
        final String expectedGATK4UncontaminatedCallsGVCF = TEST_FILES_DIR + "expected.CEUTrio.HiSeq.WGS.b37.NA12878.calls.20.10100000-10150000.g.vcf";

        // bam,
        // known contamination fraction,
        // per-sample contamination file,
        // interval,
        // reference,
        // GVCF mode,
        // GATK 4 calls on the uncontaminated bam
        // GATK 3 calls on the contaminated bam with contamination correction
        return new Object[][] {
                { contaminatedBam15Percent,
                  0.15,
                  contaminationFile,
                  traversalInterval,
                  b37_reference_20_21,
                  false, // VCF mode
                  expectedGATK4UncontaminatedCallsVCF,
                  expectedGATK3ContaminationCorrectedCallsVCF
                },
                { contaminatedBam15Percent,
                  0.15,
                  contaminationFile,
                  traversalInterval,
                  b37_reference_20_21,
                  true, // GVCF mode
                  expectedGATK4UncontaminatedCallsGVCF,
                  expectedGATK4ContaminationCorrectedCallsGVCF
                }
        };
    }

    @Test
    public void testMnpsAreRepresentedAsSingleEvents() {
        final String bam = largeFileTestDir + "contaminated_bams/CEUTrio.HiSeq.WGS.b37.NA12878.CONTAMINATED.WITH.HCC1143.NORMALS.AT.15PERCENT.bam";
        final File outputVcf = createTempFile("output", ".vcf");
        final String[] args = {
                "-I", bam,
                "-R", b37_reference_20_21,
                "-L", "20:10100000-10150000",
                "-O", outputVcf.getAbsolutePath(),
                "--" + HaplotypeCallerArgumentCollection.MAX_MNP_DISTANCE_LONG_NAME, "1"
        };
        Utils.resetRandomGenerator();
        runCommandLine(args);

        final Map<Integer, Allele> altAllelesByPosition = VariantContextTestUtils.streamVcf(outputVcf)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc.getAlternateAllele(0)));

        final Map<Integer, Allele> expectedMnps = ImmutableMap.of(
                10102247, Allele.create("AC", false),
                10102530, Allele.create("TG", false),
                10103849, Allele.create("CA", false));

        expectedMnps.entrySet().forEach(entry -> {
            final int position = entry.getKey();
            Assert.assertEquals(altAllelesByPosition.get(position), entry.getValue());
        });
    }

    // test on an artificial bam with several contrived MNPs
    // this test is basically identical to a test in {@ link Mutect2IntegrationTest}
    @Test
    public void testMnps() throws Exception {
        Utils.resetRandomGenerator();
        final File bam = new File(toolsTestDir, "mnp.bam");

        for (final int maxMnpDistance : new int[] {0, 1, 2, 3, 5}) {
            final File outputVcf = createTempFile("unfiltered", ".vcf");

            final List<String> args = Arrays.asList("-I", bam.getAbsolutePath(),
                    "-R", b37_reference_20_21,
                    "-L", "20:10019000-10022000",
                    "-O", outputVcf.getAbsolutePath(),
                    "-" + HaplotypeCallerArgumentCollection.MAX_MNP_DISTANCE_SHORT_NAME, Integer.toString(maxMnpDistance));
            runCommandLine(args);

            checkMnpOutput(maxMnpDistance, outputVcf);
        }
    }

    // this is particular to our particular artificial MNP bam -- we extract a method in order to use it for HaplotypeCaller
    private static void checkMnpOutput(int maxMnpDistance, File outputVcf) {
        // note that for testing HaplotypeCaller GVCF mode we will always have the symbolic <NON REF> allele
        final Map<Integer, List<String>> alleles = VariantContextTestUtils.streamVcf(outputVcf)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).map(Allele::getBaseString).collect(Collectors.toList())));

        // phased, two bases apart
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10019970), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10019968), Arrays.asList("GAG"));
            Assert.assertTrue(!alleles.containsKey(10019970));
        }

        // adjacent and out of phase
        Assert.assertEquals(alleles.get(10020229), Arrays.asList("A"));
        Assert.assertEquals(alleles.get(10020230), Arrays.asList("G"));

        // 4-substitution MNP w/ spacings 2, 3, 4
        if (maxMnpDistance < 2) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020432), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 3) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAG"));
            Assert.assertEquals(alleles.get(10020435), Arrays.asList("G"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else if (maxMnpDistance < 4) {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTG"));
            Assert.assertEquals(alleles.get(10020439), Arrays.asList("G"));
        } else {
            Assert.assertEquals(alleles.get(10020430), Arrays.asList("GAGTTGTCTG"));
        }

        // two out of phase DNPs that overlap and have a base in common
        if (maxMnpDistance > 0) {
            Assert.assertEquals(alleles.get(10020680), Arrays.asList("TA"));
            Assert.assertEquals(alleles.get(10020681), Arrays.asList("AT"));
        }
    }

    @Test(dataProvider = "getContaminationCorrectionTestData")
    public void testContaminationCorrection( final String contaminatedBam,
                                   final double contaminationFraction,
                                   final String contaminationFile,
                                   final SimpleInterval interval,
                                   final String reference,
                                   final boolean gvcfMode,
                                   final String gatk4UncontaminatedCallsVCF,
                                   final String expectedContaminationCorrectedCallsVCF ) throws Exception {
        final File uncorrectedOutput = createTempFile("testContaminationCorrectionUncorrectedOutput", gvcfMode ? ".g.vcf" : ".vcf");
        final File correctedOutput = createTempFile("testContaminationCorrectionCorrectedOutput", gvcfMode ? ".g.vcf" : ".vcf");
        final File correctedOutputUsingContaminationFile = createTempFile("testContaminationCorrectionCorrectedOutputUsingContaminationFile", gvcfMode ? ".g.vcf" : ".vcf");

        // Generate raw uncorrected calls on the contaminated bam, for comparison purposes
        // Note that there are a huge number of MNPs in this bam, and that in {@code gatk4UncontaminatedCallsVCF} and
        // {@code expectedContaminationCorrectedCallsVCF} these are represented as independent consecutive SNPs
        // Thus if we ever turn on MNPs by default, this will fail
        final String[] noContaminationCorrectionArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", uncorrectedOutput.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, (gvcfMode ? ReferenceConfidenceMode.GVCF.toString() : ReferenceConfidenceMode.NONE.toString()),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        Utils.resetRandomGenerator();
        runCommandLine(noContaminationCorrectionArgs);

        // Run with contamination correction using -contamination directly
        final String[] contaminationCorrectionArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", correctedOutput.getAbsolutePath(),
                "-contamination", Double.toString(contaminationFraction),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, (gvcfMode ? ReferenceConfidenceMode.GVCF.toString() : ReferenceConfidenceMode.NONE.toString()),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        Utils.resetRandomGenerator();
        runCommandLine(contaminationCorrectionArgs);

        // Run with contamination correction using a contamination file
        final String[] contaminationCorrectionFromFileArgs = {
                "-I", contaminatedBam,
                "-R", reference,
                "-L", interval.toString(),
                "-O", correctedOutputUsingContaminationFile.getAbsolutePath(),
                "-contamination-file", contaminationFile,
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, (gvcfMode ? ReferenceConfidenceMode.GVCF.toString() : ReferenceConfidenceMode.NONE.toString()),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        Utils.resetRandomGenerator();
        runCommandLine(contaminationCorrectionFromFileArgs);

        // Calculate concordance vs. the calls on the uncontaminated bam. Count only actual calls
        // (ignoring GVCF reference blocks, if present) for the purposes of the concordance calculation.
        final double uncorrectedCallsVSUncontaminatedCallsConcordance = calculateConcordance(uncorrectedOutput, new File(gatk4UncontaminatedCallsVCF), gvcfMode);
        final double correctedCallsVSUncontaminatedCallsConcordance = calculateConcordance(correctedOutput, new File(gatk4UncontaminatedCallsVCF), gvcfMode);
        final double correctedCallsFromFileVSUncontaminatedCallsConcordance = calculateConcordance(correctedOutputUsingContaminationFile, new File(gatk4UncontaminatedCallsVCF), gvcfMode);

        // Calculate concordance vs. the contamination-corrected calls from GATK3. Here we count
        // all records in the output, including reference blocks.
        final double correctedCallsVSExpectedCorrectedCallsConcordance = calculateConcordance(correctedOutput, new File(expectedContaminationCorrectedCallsVCF));
        final double correctedCallsFromFileVSExpectedCorrectedCallsConcordance = calculateConcordance(correctedOutputUsingContaminationFile, new File(expectedContaminationCorrectedCallsVCF));

        // Sanity checks: the concordance when running with -contamination should be the same as the
        // concordance when running with -contamination-file
        assertEqualsDoubleSmart(correctedCallsVSUncontaminatedCallsConcordance, correctedCallsFromFileVSUncontaminatedCallsConcordance, 0.001,
                "concordance running with -contamination and -contamination-file should be identical, but wasn't");
        assertEqualsDoubleSmart(correctedCallsVSExpectedCorrectedCallsConcordance, correctedCallsFromFileVSExpectedCorrectedCallsConcordance, 0.001,
                "concordance running with -contamination and -contamination-file should be identical, but wasn't");

        // With -contamination correction on, concordance vs. the calls on the uncontaminated bam
        // should improve by at least 1% (actual improvement tends to be much larger, but we just
        // want to make sure that the calls don't actually get worse here)
        Assert.assertTrue(correctedCallsVSUncontaminatedCallsConcordance > uncorrectedCallsVSUncontaminatedCallsConcordance + 0.01,
                "running with -contamination should have improved concordance vs. calls on the uncontaminated bam by at least 1%, but it didn't");

        // With -contamination-file correction on, concordance vs. the calls on the uncontaminated bam
        // should improve by at least 1% (actual improvement tends to be much larger, but we just
        // want to make sure that the calls don't actually get worse here)
        Assert.assertTrue(correctedCallsFromFileVSUncontaminatedCallsConcordance > uncorrectedCallsVSUncontaminatedCallsConcordance + 0.01,
                "running with -contamination-file should have improved concordance vs. calls on the uncontaminated bam by at least 1%, but it didn't");

        // -contamination corrected output in GATK4 should have >= 99% concordance
        // vs. -contamination corrected output in GATK3
        Assert.assertTrue(correctedCallsVSExpectedCorrectedCallsConcordance >= 0.99,
                "output with -contamination correction should have >= 99% concordance with contamination-corrected calls from GATK3, but it didn't");

        // -contamination-file corrected output in GATK4 should have >= 99% concordance
        // vs. -contamination corrected output in GATK3
        Assert.assertTrue(correctedCallsFromFileVSExpectedCorrectedCallsConcordance >= 0.99,
                "output with -contamination-file correction should have >= 99% concordance with contamination-corrected calls from GATK3, but it didn't");
    }

    // With 100% contamination in VCF mode, make sure that we get no calls
    @Test
    public void test100PercentContaminationNoCallsInVCFMode() throws Exception {
        final File output = createTempFile("test100PercentContaminationNoCallsInVCFMode", ".vcf");

        final String[] contaminationArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-contamination", "1.0"
        };
        runCommandLine(contaminationArgs);

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        // Check that we get no calls, but a non-empty VCF header.
        Assert.assertTrue(result.getLeft().getMetaDataInInputOrder().size() > 0, "There should be a non-empty header present");
        Assert.assertEquals(result.getRight().size(), 0, "There should be no calls with 100% contamination in VCF mode");
    }

    // With 100% contamination in GVCF mode, make sure that we get no calls (only reference blocks)
    @Test
    public void test100PercentContaminationNoCallsInGVCFMode() throws Exception {
        final File output = createTempFile("test100PercentContaminationNoCallsInGVCFMode", ".g.vcf");

        final String[] contaminationArgs = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000000-10010000",
                "-O", output.getAbsolutePath(),
                "-contamination", "1.0",
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, ReferenceConfidenceMode.GVCF.toString(),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        runCommandLine(contaminationArgs);

        final Pair<VCFHeader, List<VariantContext>> result = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath());

        // Check that we get a non-empty VCF header.
        Assert.assertTrue(result.getLeft().getMetaDataInInputOrder().size() > 0, "There should be a non-empty header present");

        // Check that we get only GVCF reference blocks, and no actual calls
        final List<VariantContext> vcfRecords = result.getRight();
        Assert.assertTrue(! vcfRecords.isEmpty(), "VCF should be non-empty in GVCF mode");
        for ( final VariantContext vc : vcfRecords ) {
            Assert.assertTrue(isGVCFReferenceBlock(vc), "Expected only GVCF reference blocks (no actual calls)");
        }
    }

    // check that a single, errorful, read that induces a cycle does not cause an assembly region to lose a real variant
    @Test
    public void testPrunedCycle() throws Exception {
        final File output = createTempFile("output", ".vcf");

        final String[] args = {
                "-I", TEST_FILES_DIR + "pruned_cycle.bam",
                "-R", b37Reference,
                "-L", "1:169510380",
                "--" + IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME, "100",
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
                "-O", output.getAbsolutePath()
        };
        runCommandLine(args);

        final Optional<VariantContext> het = VariantContextTestUtils.streamVcf(output)
                .filter(vc -> vc.getStart() == 169510380)
                .findFirst();

        Assert.assertTrue(het.isPresent());
        Assert.assertTrue(het.get().getGenotype(0).getAD()[1] > 100);
    }

    @DataProvider
    public Object[][] getMaxAlternateAllelesData() {
        return new Object[][] {
                // bam, reference, interval string, max alternate alleles, GVCF mode toggle
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10008000-10010000", 1, false },
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10002000-10011000", 1, true },
                { NA12878_20_21_WGS_bam, b37_reference_20_21, "20:10002000-10011000", 2, true }
        };
    }

    /*
     * Test for the --max-alternate-alleles argument
     */
    @Test(dataProvider = "getMaxAlternateAllelesData")
    public void testMaxAlternateAlleles(final String bam, final String reference, final String intervalString,
                                        final int maxAlternateAlleles, final boolean gvcfMode) {
        final File outputNoMaxAlternateAlleles = createTempFile("testMaxAlternateAllelesNoMaxAlternateAlleles", (gvcfMode ? ".g.vcf" : ".vcf"));
        final File outputWithMaxAlternateAlleles = createTempFile("testMaxAlternateAllelesWithMaxAlternateAlleles", (gvcfMode ? ".g.vcf" : ".vcf"));

        // Run both with and without --max-alternate-alleles over our interval, so that we can
        // prove that the argument is working as intended.
        final String[] argsNoMaxAlternateAlleles = {
                "-I", bam,
                "-R", reference,
                "-L", intervalString,
                "-O", outputNoMaxAlternateAlleles.getAbsolutePath(),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, (gvcfMode ? ReferenceConfidenceMode.GVCF.toString() : ReferenceConfidenceMode.NONE.toString()),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        runCommandLine(argsNoMaxAlternateAlleles);

        final String[] argsWithMaxAlternateAlleles = {
                "-I", bam,
                "-R", reference,
                "-L", intervalString,
                "-O", outputWithMaxAlternateAlleles.getAbsolutePath(),
                "--max-alternate-alleles", Integer.toString(maxAlternateAlleles),
                "--" + AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, (gvcfMode ? ReferenceConfidenceMode.GVCF.toString() : ReferenceConfidenceMode.NONE.toString()),
                "--" + AssemblyBasedCallerArgumentCollection.ALLELE_EXTENSION_LONG_NAME, "2",
        };
        runCommandLine(argsWithMaxAlternateAlleles);

        final List<VariantContext> callsNoMaxAlternateAlleles = VariantContextTestUtils.readEntireVCFIntoMemory(outputNoMaxAlternateAlleles.getAbsolutePath()).getRight();
        final List<VariantContext> callsWithMaxAlternateAlleles = VariantContextTestUtils.readEntireVCFIntoMemory(outputWithMaxAlternateAlleles.getAbsolutePath()).getRight();

        // First, find all calls in the VCF produced WITHOUT --max-alternate-alleles that have
        // more than maxAlternateAlleles alt alleles, excluding NON_REF. For each call, calculate
        // and store the expected list of subsetted alleles:
        final Map<SimpleInterval, List<Allele>> expectedSubsettedAllelesByLocus = new HashMap<>();
        for ( final VariantContext vc : callsNoMaxAlternateAlleles ) {
            if ( getNumAltAllelesExcludingNonRef(vc) > maxAlternateAlleles ) {
                final List<Allele> mostLikelyAlleles = AlleleSubsettingUtils.calculateMostLikelyAlleles(vc, HomoSapiensConstants.DEFAULT_PLOIDY, maxAlternateAlleles, false);
                expectedSubsettedAllelesByLocus.put(new SimpleInterval(vc), mostLikelyAlleles);
            }
        }

        // Then assert that we saw at least one call with more than maxAlternateAlleles alt alleles
        // when running without --max-alternate-alleles (otherwise, the tests below won't be meaningful):
        Assert.assertTrue(! expectedSubsettedAllelesByLocus.isEmpty(),
                "Without --max-alternate-alleles, there should be at least one call in the output with more than " + maxAlternateAlleles +
                        " alt alleles in order for this test to be meaningful");

        // Now assert that in the VCF produced WITH --max-alternate-alleles, there are no calls with
        // more than maxAlternateAlleles alt alleles, excluding NON_REF. Also check each call that would
        // have had more than maxAlternateAlleles alleles against the expected list of subsetted alleles,
        // to ensure that we selected the most likely alleles:
        for ( final VariantContext vc : callsWithMaxAlternateAlleles ) {

            // No call should have more than the configured number of alt alleles (excluding NON_REF)
            Assert.assertTrue(getNumAltAllelesExcludingNonRef(vc) <= maxAlternateAlleles,
                    "Number of alt alleles exceeds --max-alternate-alleles " + maxAlternateAlleles + " for VariantContext: " + vc);

            // If there's an entry for this locus in our table of expected alleles post-subsetting, assert
            // that we selected the right alleles during subsetting.
            List<Allele> alleleSubsettingExpectedResult = expectedSubsettedAllelesByLocus.get(new SimpleInterval(vc));
            if ( alleleSubsettingExpectedResult != null ) {

                // CollectionUtils.isEqualCollection() will compare the lists of Alleles without
                // regard to ordering
                Assert.assertTrue(CollectionUtils.isEqualCollection(vc.getAlleles(), alleleSubsettingExpectedResult),
                        "For call " + vc + " expected alleles after subsetting were: " + alleleSubsettingExpectedResult +
                                 " but instead found alleles: " + vc.getAlleles());
            }

            // For completeness sake, also check the genotypes to ensure that no genotypes reference
            // an allele not present in the VC:
            for ( final Genotype genotype : vc.getGenotypes() ) {
                if ( genotype.isAvailable() ) {
                    for ( final Allele genotypeAllele : genotype.getAlleles() ) {
                        if ( genotypeAllele.isCalled() ) {
                            Assert.assertTrue(vc.hasAllele(genotypeAllele),
                                    "Allele " + genotypeAllele + " was present in genotype " + genotype +
                                            " but not in the VariantContext itself");
                        }
                    }
                }
            }
        }
    }

    @Test
    public void testContaminatedHomVarDeletions() {
        final String bam = toolsTestDir + "haplotypecaller/deletion_sample.snippet.bam";
        final String intervals = "chr3:46373452";

        final File outputContaminatedHomVarDeletions = createTempFile("testContaminatedHomVarDeletions", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder().addInput(new File(bam))
                .addReference(hg38Reference)
                .addInterval(new SimpleInterval(intervals))
                .addOutput(outputContaminatedHomVarDeletions)
                .add(IntervalArgumentCollection.INTERVAL_PADDING_LONG_NAME, 50);
        runCommandLine(args);

        List<VariantContext> vcs = VariantContextTestUtils.getVariantContexts(outputContaminatedHomVarDeletions);

        //check known homozygous deletion for correct genotype
        for (final VariantContext vc : vcs) {
            final Genotype gt = vc.getGenotype(0);
            if (gt.hasAD()) {
                final int[] ads = gt.getAD();
                final double alleleBalance = ads[1] / (ads[0] + ads[1]);
                if (alleleBalance > 0.9) {
                    Assert.assertTrue(vc.getGenotype(0).isHomVar());
                }
            }
        }
    }

    // the order of adjacent deletion and insertion events is arbitary, hence the locus of an insertion can equally well
    // be assigned to the beginning or end of an adjacent deletion.  However, assigning to the beginning of the deletion
    // creates a single event that looks like a MNP or complicated event eg ACG -> TTGT.  There is nothing wrong with this -- indeed
    // it is far preferable to calling two events -- however GenomicsDB as of March 2020 does not support MNPs.  Once it does,
    // this test will not be necessary.
    @Test
    public void testAdjacentIndels() {
        final File bam = new File(TEST_FILES_DIR, "issue_6473_adjacent_indels.bam");
        final String interval = "17:7578000-7578500";

        final File output = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(bam)
                .addReference(b37Reference)
                .addInterval(interval)
                .addOutput(output);
        runCommandLine(args);

        Assert.assertTrue(VariantContextTestUtils.streamVcf(output)
                .allMatch(vc -> vc.isBiallelic() && (vc.getReference().length() == 1 || vc.getAlternateAllele(0).length() == 1)));
    }

    /*
    Prior to IUPAC ReadTransformer fix, this test yields
    java.lang.IllegalArgumentException: Unexpected base in allele bases 'ATTCATTTCACAAGGGTAAAGCTTTCTTTGGATTCAGCAGGTTGGAAAATCTGTTTTTCACCTTTCTGTGAATGGACGTTTGGGAGCTCATTGAGGCCAGTGRCAATAAAGGAGATATCTCAGGGTGAAAAATAAAAGACAGGAATGTGAGAATTGGCTTTGTGATGTGAGCATTCATTTCACAAAGTTAAACCTTTCTTTTCATTCAGCAGTTAGAAATCACTGGTTTTGTAGAATCTG'
    where there's an R in that big long string representing the assembled haplotype.  The R is in the hg38 reference and
    gets propagated into cram reads.
     */
    @Test
    public void testAdjacentIUPACBasesinReads() {
        final File bam = new File(TEST_FILES_DIR, "cramWithR.cram");
        final String interval = "chr10:39239400-39239523";

        final File output = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(bam)
                .addReference(hg38Reference)
                .addInterval(interval)
                .addOutput(output)
                //We can check the warning in the test log provided we're explicit about the logging level
                .add(StandardArgumentDefinitions.VERBOSITY_NAME, Log.LogLevel.WARNING.name());
        runCommandLine(args);

        final List<VariantContext> outputVCs = VariantContextTestUtils.readEntireVCFIntoMemory(output.getAbsolutePath()).getRight();
        Assert.assertEquals(outputVCs.size(), 1);
        Assert.assertEquals(outputVCs.get(0).getStart(), 39239403);
        Assert.assertEquals(outputVCs.get(0).getAlternateAllele(0).getBaseString(), "A");
    }

    // this test has a reference with 8 repeats of a 28-mer and an alt with 7 repeats.  This deletion left-aligns to the
    // beginning of the padded assembly region, and an exception occurs if we carelessly drop the leading deletion from
    // the alt haplotype's cigar.  This is a regression test for https://github.com/broadinstitute/gatk/issues/6533.
    // In addition to testing that no error is thrown, we also check that the 28-base deletion is called
    @Test
    public void testLeadingDeletionInAltHaplotype() {
        final File bam = new File(TEST_FILES_DIR, "alt-haplotype-leading-deletion.bam");

        final File output = createTempFile("output", ".vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(bam)
                .addReference(b37Reference)
                .addInterval("10:128360000-128362000")
                .addOutput(output);
        runCommandLine(args);

        final boolean has28BaseDeletion = VariantContextTestUtils.streamVcf(output)
                .anyMatch(vc -> vc.getStart() == 128361367 && vc.getAlternateAlleles().stream().anyMatch(a -> a.length() == vc.getReference().length() - 28));

        Assert.assertTrue(has28BaseDeletion);
    }

    @Test
    public void testGvcfFlowConsistentWithPreviousResults() throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testGvcfBeforeRebase.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-likelihood-optimized-comp", true)
                .add("ERC", "GVCF")
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testGvcfUsingFlowModeStandardConsistentWithPreviousResults() throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testGvcfBeforeRebaseUsingFlowModeStandard.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("flow-mode", "STANDARD")
                .add("ERC", "GVCF")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testGvcfUsingFlowModeAdvancedConsistentWithPreviousResults() throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testGvcfBeforeRebaseUsingFlowModeAdvanced.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("flow-mode", "ADVANCED")
                .add("ERC", "GVCF")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testGvcfKeepLoneAlleles() throws Exception {
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testGvcfKeepLoneAlleles.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();


        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-likelihood-optimized-comp", true)
                .add("ERC", "GVCF")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( !UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testVcfFliowModeConsistentWithPreviousResults() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testVcfBeforeRebase.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testVcfWithAssemblyComplexityAnnotation() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("testVcfWithAssemblyComplexityAnnotation", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "testVcfWithAssemblyComplexityAnnotation.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-filter-lone-alleles", true)
                .add("A", "AssemblyComplexity")
                .addFlag("assembly-complexity-reference-mode")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testGvcfWithAssemblyComplexityAnnotationRevamp() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("testVcfWithAssemblyComplexityAnnotation", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "testGvcfWithAssemblyComplexityAnnotationRevamp.expected.flowbased.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("mbq", 0).
                addFlag("flow-filter-alleles")
                .add("flow-filter-alleles-sor-threshold","3")
                .add("flow-assembly-collapse-hmer-size ","12")
                .add("flow-matrix-mods","10,12,11,12")
                .add("bam-writer-type","CALLED_HAPLOTYPES_NO_READS")
                .addFlag("apply-frd")
                .add("ERC","GVCF")
                .add("minimum-mapping-quality","1")
                .add("mapping-quality-threshold-for-genotyping","1")
                .addFlag("override-fragment-softclip-check")
                .add("flow-likelihood-parallel-threads","2")
                .addFlag("flow-likelihood-optimized-comp")
                .addFlag("adaptive-pruning")
                .add("pruning-lod-threshold","3.0")
                .addFlag("enable-dynamic-read-disqualification-for-genotyping")
                .add("dynamic-read-disqualification-threshold","10")
                .add("G","StandardAnnotation")
                .add("G","StandardHCAnnotation")
                .add("G","AS_StandardAnnotation")
                .add("GQB","10").add("GQB","20").add("GQB","30").add("GQB","40").add("GQB","50").add("GQB","60").add("GQB","70").add("GQB","80").add("GQB","90")
                .add("likelihood-calculation-engine","FlowBased")
                .add("A", "AssemblyComplexity")
                .addFlag("assembly-complexity-reference-mode")
                .add("verbosity","INFO")
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    /**
     * Helper method for testMaxAlternateAlleles
     *
     * @param vc VariantContext to check
     * @return number of alt alleles in vc, excluding NON_REF (if present)
     */
    private int getNumAltAllelesExcludingNonRef( final VariantContext vc ) {
        final List<Allele> altAlleles = vc.getAlternateAlleles();
        int numAltAllelesExcludingNonRef = 0;

        for ( final Allele altAllele : altAlleles ) {
            if ( ! altAllele.equals(Allele.NON_REF_ALLELE) ) {
                ++numAltAllelesExcludingNonRef;
            }
        }

        return numAltAllelesExcludingNonRef;
    }

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     */
    public static double calculateConcordance( final File actual, final File expected ) {
        return calculateConcordance(actual, expected, false);
    }

    /*
     * Calculate rough concordance between two vcfs, comparing only the positions, alleles, and the first genotype.
     * Can optionally ignore GVCF blocks in the concordance calculation.
     */
    public static double calculateConcordance( final File actual, final File expected, final boolean ignoreGVCFBlocks ) {
        final Set<String> actualVCFKeys = new HashSet<>();
        final Set<String> expectedVCFKeys = new HashSet<>();
        int concordant = 0;
        int discordant = 0;

        // unzip to avoid https://github.com/broadinstitute/gatk/issues/4224
        File actualUnzipped = IOUtils.gunzipToTempIfNeeded(actual);
        try ( final FeatureDataSource<VariantContext> actualSource = new FeatureDataSource<>(actualUnzipped);
              final FeatureDataSource<VariantContext> expectedSource = new FeatureDataSource<>(expected) ) {

            for ( final VariantContext vc : actualSource ) {
                if ( ! ignoreGVCFBlocks || ! isGVCFReferenceBlock(vc) ) {
                    actualVCFKeys.add(keyForVariant(vc));
                }
            }

            for ( final VariantContext vc : expectedSource ) {
                if ( ! ignoreGVCFBlocks || ! isGVCFReferenceBlock(vc) ) {
                    expectedVCFKeys.add(keyForVariant(vc));
                }
            }

            for ( final String vcKey : actualVCFKeys ) {
                if ( ! expectedVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
                else {
                    ++concordant;
                }
            }

            for ( final String vcKey : expectedVCFKeys ) {
                if ( ! actualVCFKeys.contains(vcKey) ) {
                    ++discordant;
                }
            }
        }

        return (double)concordant / (double)(concordant + discordant);
    }

    private static String keyForVariant( final VariantContext variant ) {
        Genotype genotype = variant.getGenotype(0);
        if (genotype.isPhased()) { // unphase it for comparisons, since we rely on comparing the genotype string below
            genotype = new GenotypeBuilder(genotype).phased(false).make();
        }
        return String.format("%s:%d-%d %s %s", variant.getContig(), variant.getStart(), variant.getEnd(),
                variant.getAlleles(), genotype.getGenotypeString(false));
    }

    private static boolean isGVCFReferenceBlock( final VariantContext vc ) {
        return vc.hasAttribute(VCFConstants.END_KEY) &&
               vc.getAlternateAlleles().size() == 1 &&
               vc.getAlternateAllele(0).equals(Allele.NON_REF_ALLELE);
    }

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testPileupCallingDefaultsConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.pileupCallerDefaults.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",

                // NOTE: These arguments are intended to force the assembly to fail at all sites in order to test what gets recovered by the pileup code
                "--" + ReadThreadingAssemblerArgumentCollection.KMER_SIZE_LONG_NAME, "1",
                "--" + ReadThreadingAssemblerArgumentCollection.DONT_INCREASE_KMER_SIZE_LONG_NAME,

                // Pileup Caller args
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME,
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ENABLE_INDELS,
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    /*
     * Test that in VCF mode we're consistent with past GATK4 results
     */
    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testPileupCallingDRAGENModeConsistentWithPastResultsWithIndels(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.pileupCallerDRAGEN.WithIndels.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",

                // NOTE: These arguments are intended to force the assembly to fail at all sites in order to test what gets recovered by the pileup code
                "--" + ReadThreadingAssemblerArgumentCollection.KMER_SIZE_LONG_NAME, "1",
                "--" + ReadThreadingAssemblerArgumentCollection.DONT_INCREASE_KMER_SIZE_LONG_NAME,

                // Pileup Caller args
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME,
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH, "0",
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME, "0.01", // an excessive strict read badness filter designed to amplify the effect (the default here is 40%)
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ENABLE_INDELS,
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testPileupCallingDRAGENModeConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.pileupCallerDRAGEN.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "-pairHMM", "AVX_LOGLESS_CACHING",

                // NOTE: These arguments are intended to force the assembly to fail at all sites in order to test what gets recovered by the pileup code
                "--" + ReadThreadingAssemblerArgumentCollection.KMER_SIZE_LONG_NAME, "1",
                "--" + ReadThreadingAssemblerArgumentCollection.DONT_INCREASE_KMER_SIZE_LONG_NAME,

                // Pileup Caller args
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME,
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH, "0",
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME, "0.01", // an excessive strict read badness filter designed to amplify the effect (the default here is 40%)
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false"
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(dataProvider="HaplotypeCallerTestInputs")
    public void testPileupCallingDRAGEN378ModeConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.pileupCallerDRAGEN.378.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "--dragen-378-concordance-mode",
                "--dragstr-params-path", TEST_FILES_DIR+"example.dragstr-params.txt",

                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @Test(dataProvider="HaplotypeCallerTestInputs", groups={"bucket"})
    public void testPileupCallingDRAGEN378OptimizedModeConsistentWithPastResults(final String inputFileName, final String referenceFileName) throws Exception {
        Utils.resetRandomGenerator();

        final File output = createTempFile("testVCFModeIsConsistentWithPastResults", ".vcf");
        final File expected = new File(TEST_FILES_DIR, "expected.pileupCallerDRAGEN.378.gatk4.vcf");

        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final String[] args = {
                "-I", inputFileName,
                "-R", referenceFileName,
                "-L", "20:10000000-10100000",
                "-O", outputPath,
                "--dragen-378-concordance-mode",
                "--dragstr-params-path", TEST_FILES_DIR+"example.dragstr-params.txt",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",

                // This optimization in PDHMM mode should have essentially no bearing on the output, this test deliberately
                // uses the same expected output as the above test to enforce that over this short test range the optimization
                // doesn't impact the results in any meaningful way.
                "--" + PileupDetectionArgumentCollection.PDHMM_READ_OVERLAP_OPTIMIZATION
        };

        runCommandLine(args);

        // Test for an exact match against past results
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected);
        }
    }

    @DataProvider(name="PairHMMResultsModes")
    public Object[][] PairHMMResultsModes() {
        return new Object[][] {
                // The native PairHMM is checked against the Java results within a tolerance, see testPairHMMResultsFile
                {PairHMM.Implementation.AVX_LOGLESS_CACHING, new File(TEST_FILES_DIR, "expected.Java.hmmresults.txt")},
                {PairHMM.Implementation.LOGLESS_CACHING, new File(TEST_FILES_DIR, "expected.Java.hmmresults.txt")},
                {PairHMM.Implementation.ORIGINAL, new File(TEST_FILES_DIR, "expected.Original.hmmresults.txt")},
                {PairHMM.Implementation.EXACT, new File(TEST_FILES_DIR, "expected.Exact.hmmresults.txt")},

        };
    }

    @Test(dataProvider = "PairHMMResultsModes")
    public void testPairHMMResultsFile(final PairHMM.Implementation implementation, final File expected) throws Exception {
        Utils.resetRandomGenerator();

        final File vcfOutput = createTempFile("hmmResultFileTest", ".vcf");
        final File hmmOutput = createTempFile("hmmResult", ".txt");

        final boolean isNative = implementation == PairHMM.Implementation.AVX_LOGLESS_CACHING;
        final String hmmOutputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS && !isNative ? expected.getAbsolutePath() : hmmOutput.getAbsolutePath();

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10000117",
                "-ip", "100",
                "-O", vcfOutput.getAbsolutePath(),
                "--pair-hmm-results-file", hmmOutputPath,
                "-pairHMM", implementation.name()
        };

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            // The native PairHMM runs in single precision, so its likelihoods are compared to the Java PairHMM's
            // within a tolerance. The Java PairHMMs must reproduce their recorded output exactly.
            if (isNative) {
                assertPairHMMResultsMatch(hmmOutput, expected, NATIVE_PAIR_HMM_TOLERANCE);
            } else {
                IntegrationTestSpec.assertEqualTextFiles(hmmOutput, expected);
            }
        }
    }

    /** Largest allowed difference in log10 likelihood between the native and Java PairHMMs, as in VectorPairHMMUnitTest. */
    private static final double NATIVE_PAIR_HMM_TOLERANCE = 1e-5;

    /**
     * Asserts that two PairHMM results files record the same read/haplotype pairs, in any order, with log10
     * likelihoods within {@code tolerance} of each other.
     */
    private static void assertPairHMMResultsMatch(final File actual, final File expected, final double tolerance) throws IOException {
        final List<String[]> actualResults = readPairHMMResults(actual);
        final List<String[]> expectedResults = readPairHMMResults(expected);
        Assert.assertEquals(actualResults.size(), expectedResults.size(), "number of recorded read/haplotype pairs");
        for (int i = 0; i < actualResults.size(); i++) {
            final String[] actualResult = actualResults.get(i);
            final String[] expectedResult = expectedResults.get(i);
            final String pair = String.join(" ", Arrays.copyOf(expectedResult, 6));
            Assert.assertEquals(String.join(" ", Arrays.copyOf(actualResult, 6)), pair, "read/haplotype pair");
            Assert.assertEquals(Double.parseDouble(actualResult[6]), Double.parseDouble(expectedResult[6]), tolerance, "likelihood of " + pair);
        }
    }

    /**
     * Reads a PairHMM results file into its space-separated tokens per pair, sorted so that two files recording the
     * same pairs can be compared line by line.
     */
    private static List<String[]> readPairHMMResults(final File file) throws IOException {
        return Files.readAllLines(file.toPath()).stream()
                .filter(line -> !line.isEmpty() && !line.startsWith("#"))
                .map(line -> line.trim().split(" "))
                .sorted(Comparator.comparing(tokens -> String.join(" ", tokens)))
                .collect(Collectors.toList());
    }

    @Test()
    // Since the PDHMM isn't currently slot-in compared to the pairHMM modes this has to be computed seperately...
    public void testPDPairHMMResultsFile() throws Exception {
        Utils.resetRandomGenerator();

        final File vcfOutput = createTempFile("hmmResultFileTest", ".vcf");
        final File hmmOutput = createTempFile("hmmResult", ".txt");
        final File expected = new File(largeFileTestDir, "expected.PDHMM.hmmresults.txt");

        final String hmmOutputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : hmmOutput.getAbsolutePath();

        final String[] args = {
                "-I", NA12878_20_21_WGS_bam,
                "-R", b37_reference_20_21,
                "-L", "20:10009875", // site selected as one of the more complex in our standard test inputs (which means lots of PD events in the haplotypes)
                "-ip", "300",
                "-O", vcfOutput.getAbsolutePath(),
                "--dragen-mode",
                "--dragstr-params-path", TEST_FILES_DIR+"example.dragstr-params.txt",

                // Pileup Caller args
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME,
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ABSOLUTE_ALT_DEPTH, "0",
                "--" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, "false",

                // Pileup caller and PDHMM arguments
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_BAD_READ_RATIO_LONG_NAME, "0.40",
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ENABLE_INDELS,
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_ACTIVE_REGION_PHRED_THRESHOLD_LONG_NAME, "3.0",
                "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_FILTER_ASSEMBLY_HAPS_THRESHOLD, "0.4",
                "--" + PileupDetectionArgumentCollection.GENERATE_PARTIALLY_DETERMINED_HAPLOTYPES_LONG_NAME,

                "--pdhmm-results-file", hmmOutputPath
        };

        runCommandLine(args);

        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(hmmOutput, expected);
        }
    }

    @Test
    public void testVcfFlowLikelihoodOptimizedCompConsistentWithPreviousResults() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testVcfBeforeRebase.expected.flowbased.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(output)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-likelihood-optimized-comp", true)
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
    }

    @Test
    public void testVcfFlowLikelihoodParralelThreadsConsistentWithPreviousResults() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "testVcfBeforeRebase.expected.flowbased.vcf");

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(output)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBased")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-likelihood-optimized-comp", true)
                .add("flow-likelihood-parallel-threads", 4)
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);

        IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
    }

    @Test
    public void testVcfFlowBasedHMMConsistentWithPreviousResults() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "test_flowBasedHMM.expected.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBasedHMM")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    @Test
    public void testVcfFlowBasedHMMStepwiseConsistentWithPreviousResults() throws Exception {
        Utils.resetRandomGenerator();
        final File input = new File(largeFileTestDir, "input_jukebox_for_test.bam");
        final File output = createTempFile("output", ".vcf");

        final File expected = new File(TEST_FILES_DIR, "test_flowBasedHMM_Stepwise.expected.vcf");
        final String outputPath = UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ? expected.getAbsolutePath() : output.getAbsolutePath();

        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addReference(hg38Reference)
                .addInterval("chr9:81149486-81177047")
                .addOutput(outputPath)
                .addInput(input)
                .add("smith-waterman", "FASTEST_AVAILABLE")
                .add("likelihood-calculation-engine", "FlowBasedHMM")
                .add("mbq", "0")
                .add("kmer-size", 10)
                .add("flow-filter-alleles", true)
                .add("flow-filter-alleles-sor-threshold", 40)
                .add("flow-assembly-collapse-hmer-size", 12)
                .add("flow-matrix-mods", "10,12,11,12")
                .addFlag("use-flow-aligner-for-stepwise-hc-filtering")
                .add("flow-filter-lone-alleles", true)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);

        runCommandLine(args);
        if ( ! UPDATE_EXACT_MATCH_EXPECTED_OUTPUTS ) {
            IntegrationTestSpec.assertEqualTextFiles(output, expected, "#");
        }
    }

    // 20 kb of chromosome 21 where many of the NA12878 test BAM's reads have MAPQ 0
    private static final String MQ0_RESCUE_INTERVAL = "21:10000000-10020000";

    /** @return the arguments of a DRAGEN-mode run on the rescue interval (or the given ones), the mode the rescue requires */
    private ArgumentsBuilder mq0RescueArgs(final File output, final ReferenceConfidenceMode mode, final String... intervals) {
        return mq0RescueArgsWithoutDragen(output, mode, intervals).addFlag(HaplotypeCallerArgumentCollection.DRAGEN_3412_GATK_MODE_LONG_NAME);
    }

    /** @return the value --dragen-mode gives an argument */
    private static String dragenModeValue(final String name) {
        final String[] pairs = new HaplotypeCallerArgumentCollection().getDragenVersion3412NameValuePairs();
        for (int i = 0; i < pairs.length; i += 2) {
            if (pairs[i].equals(name)) {
                return pairs[i + 1];
            }
        }
        throw new IllegalArgumentException(name + " is not set by --dragen-mode");
    }

    private ArgumentsBuilder mq0RescueArgsWithoutDragen(final File output, final ReferenceConfidenceMode mode, final String... intervals) {
        final ArgumentsBuilder args = new ArgumentsBuilder()
                .addInput(NA12878_20_21_WGS_bam)
                .addReference(b37_reference_20_21)
                .addOutput(output)
                .add(AssemblyBasedCallerArgumentCollection.EMIT_REF_CONFIDENCE_LONG_NAME, mode)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false);
        for (final String interval : intervals.length == 0 ? new String[]{MQ0_RESCUE_INTERVAL} : intervals) {
            args.addInterval(interval);
        }
        return args;
    }

    // outputs of the runs several tests share, by mode and whether rescue is on
    private final Map<String, File> mq0RescueRuns = new HashMap<>();

    private synchronized File mq0RescueRun(final ReferenceConfidenceMode mode, final boolean rescue) {
        return mq0RescueRuns.computeIfAbsent(mode + "." + rescue, key -> {
            final File output = createTempFile("mq0Rescue." + key, mode == ReferenceConfidenceMode.NONE ? ".vcf" : ".g.vcf");
            final ArgumentsBuilder args = mq0RescueArgs(output, mode);
            if (rescue) {
                args.addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME);
            }
            runCommandLine(args);
            return output;
        });
    }

    private static final String GVCF_BLOCK_ALT = "<NON_REF>";

    private static boolean isMq0RescueRecord(final VariantContext vc) {
        return vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_KEY);
    }

    private static List<String> dataLines(final File vcf) throws IOException {
        return Files.readAllLines(vcf.toPath()).stream().filter(line -> !line.startsWith("#")).collect(Collectors.toList());
    }

    /** @return the number of reads with MAPQ 0 (or above 0) that HaplotypeCaller's other read filters keep, overlapping a position */
    private static int countBamReads(final String contig, final int pos, final boolean mq0) {
        try (final ReadsDataSource reads = new ReadsPathDataSource(IOUtils.getPath(NA12878_20_21_WGS_bam))) {
            int n = 0;
            for (final Iterator<GATKRead> it = reads.query(new SimpleInterval(contig, pos, pos)); it.hasNext(); ) {
                final GATKRead r = it.next();
                if (!r.isUnmapped() && !r.isSecondaryAlignment() && !r.isDuplicate() && !r.failsVendorQualityCheck()
                        && (mq0 == (r.getMappingQuality() == 0))) {
                    n++;
                }
            }
            return n;
        }
    }

    @DataProvider
    public Object[][] mq0RescueModes() {
        return new Object[][]{{ReferenceConfidenceMode.GVCF}, {ReferenceConfidenceMode.NONE}, {ReferenceConfidenceMode.BP_RESOLUTION}};
    }

    /**
     * With --mq0-rescue the output is the regular output with the rescue records merged in: non-reference calls whose
     * alleles MAPQ-0 reads carry, marked, away from the regular calls, overlapping no other record.
     */
    @Test(dataProvider = "mq0RescueModes")
    public void testMq0Rescue(final ReferenceConfidenceMode mode) throws IOException {
        assertMq0RescueOutput(mq0RescueRun(mode, false), mq0RescueRun(mode, true), mode);
    }

    /**
     * @return whether a record is a call, worked out here rather than with the code under test: a non-reference,
     * non-symbolic alternate allele, a called genotype carrying one, and QUAL at least the rescue's minimum
     */
    private static boolean isCallForTest(final VariantContext vc) {
        final double minCallQual = new Mq0RescueArgumentCollection().minCallQual;
        return vc.getPhredScaledQual() >= minCallQual && vc.getAlternateAlleles().stream().anyMatch(a -> !a.isSymbolic() && !a.equals(Allele.SPAN_DEL))
                && vc.getGenotypes().stream().anyMatch(g -> g.isCalled() && g.getAlleles().stream()
                        .anyMatch(a -> a.isCalled() && a.isNonReference() && !a.isSymbolic() && !a.equals(Allele.SPAN_DEL)));
    }

    private static int distanceForTest(final VariantContext a, final VariantContext b) {
        return a.getEnd() < b.getStart() ? b.getStart() - a.getEnd() : b.getEnd() < a.getStart() ? a.getStart() - b.getEnd() : 0;
    }

    /**
     * Checks a rescue run's output against the same run without rescue.
     *
     * @return the rescue records
     */
    private static List<VariantContext> assertMq0RescueOutput(final File plain, final File rescued, final ReferenceConfidenceMode mode) throws IOException {
        return assertMq0RescueOutput(plain, rescued, mode, 10);
    }

    /**
     * Checks a rescue run's output against the same run without rescue.
     *
     * @param minRescue the fewest rescue records expected
     * @return the rescue records
     */
    private static List<VariantContext> assertMq0RescueOutput(final File plain, final File rescued, final ReferenceConfidenceMode mode,
                                                              final int minRescue) throws IOException {
        final Pair<VCFHeader, List<VariantContext>> plainVcf = VariantContextTestUtils.readEntireVCFIntoMemory(plain.getAbsolutePath());
        final Pair<VCFHeader, List<VariantContext>> rescuedVcf = VariantContextTestUtils.readEntireVCFIntoMemory(rescued.getAbsolutePath());

        // the header is the regular one plus the rescue lines
        final VCFHeader header = rescuedVcf.getLeft();
        Assert.assertNotNull(header.getFilterHeaderLine(GATKVCFConstants.MQ0_RESCUE_FILTER_NAME));
        for (final String key : GATKVCFConstants.MQ0_RESCUE_INFO_KEYS) {
            Assert.assertTrue(header.hasInfoLine(key), key);
            Assert.assertFalse(plainVcf.getLeft().hasInfoLine(key), key);
        }
        Assert.assertEquals(header.getMetaDataInInputOrder().size(), plainVcf.getLeft().getMetaDataInInputOrder().size() + 4);

        final List<VariantContext> records = rescuedVcf.getRight();
        final List<VariantContext> rescue = records.stream().filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).collect(Collectors.toList());
        Assert.assertTrue(rescue.size() >= minRescue, "expected rescue records, got " + rescue.size());

        // Apart from the rescue records, the records are those of the regular output, unchanged (in VCF mode, QD
        // included, so the rescue draws nothing from the shared random generator), except for the reference records
        // and the regular non-calls under rescue records, which these replace. GVCF blocks are left out of the
        // comparison, since rescue records split them.
        final List<String> plainLines = dataLines(plain);
        final List<VariantContext> plainRecords = plainVcf.getRight();
        Assert.assertEquals(plainLines.size(), plainRecords.size());
        final List<String> expected = new ArrayList<>();
        for (int i = 0; i < plainLines.size(); i++) {
            final VariantContext vc = plainRecords.get(i);
            final boolean block = vc.getAlternateAlleles().equals(List.of(Allele.NON_REF_ALLELE));
            final boolean replaced = !isCallForTest(vc) && rescue.stream().anyMatch(r -> r.overlaps(vc));
            if (!replaced && !(mode == ReferenceConfidenceMode.GVCF && block)) {
                expected.add(plainLines.get(i));
            }
        }
        final List<String> actual = dataLines(rescued).stream()
                .filter(line -> !line.split("\t", 8)[6].equals(GATKVCFConstants.MQ0_RESCUE_FILTER_NAME))
                .filter(line -> mode != ReferenceConfidenceMode.GVCF || !line.split("\t", 6)[4].equals(GVCF_BLOCK_ALT))
                .collect(Collectors.toList());
        Assert.assertEquals(actual, expected);

        final List<VariantContext> regularCalls = records.stream()
                .filter(vc -> !isMq0RescueRecord(vc) && isCallForTest(vc)).collect(Collectors.toList());
        Assert.assertFalse(regularCalls.isEmpty());
        final SimpleInterval interval = new SimpleInterval(MQ0_RESCUE_INTERVAL);
        int lastStart = 0;
        for (final VariantContext vc : records) {
            Assert.assertTrue(vc.getStart() > lastStart, "out of order or at the same position: " + vc);
            lastStart = vc.getStart();
        }
        final int distance = new Mq0RescueArgumentCollection().distance;
        for (final VariantContext vc : rescue) {
            Assert.assertTrue(interval.contains(new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart())), vc.toString());
            Assert.assertEquals(vc.getFilters(), Set.of(GATKVCFConstants.MQ0_RESCUE_FILTER_NAME), vc.toString());
            // a call: a non-reference genotype with QUAL at least the minimum, and no '*' allele
            Assert.assertTrue(isCallForTest(vc), vc.toString());
            Assert.assertFalse(vc.hasAllele(Allele.SPAN_DEL), vc.toString());
            // no other record overlaps it
            Assert.assertTrue(records.stream().noneMatch(other -> other != vc && other.overlaps(vc)), vc + " overlaps another record");
            // MAPQ-0 reads cover it, but no more than the BAM has; and no more reads with MAPQ above 0 than the BAM has
            final int mq0Depth = vc.getAttributeAsInt(GATKVCFConstants.MQ0_RESCUE_MQ0_DEPTH_KEY, -1);
            Assert.assertTrue(mq0Depth > 0 && mq0Depth <= countBamReads(vc.getContig(), vc.getStart(), true), vc.toString());
            final int depth = vc.getAttributeAsInt(GATKVCFConstants.MQ0_RESCUE_DEPTH_KEY, -1);
            Assert.assertTrue(depth >= 0 && depth <= countBamReads(vc.getContig(), vc.getStart(), false), vc.toString());
            // no regular call within the rescue distance
            final int nearest = regularCalls.stream().mapToInt(fp -> distanceForTest(fp, vc)).min().getAsInt();
            Assert.assertTrue(nearest > distance, vc + " is " + nearest + " bp from a regular call");
            // a SNP's called alt base is carried by a MAPQ-0 read
            if (vc.isSNP()) {
                final Set<String> calledAlts = vc.getGenotype(0).getAlleles().stream().filter(Allele::isNonReference)
                        .map(Allele::getBaseString).collect(Collectors.toSet());
                Assert.assertTrue(calledAlts.stream().anyMatch(alt -> mq0BasesAt(vc.getContig(), vc.getStart()).contains(alt.charAt(0))), vc.toString());
            }
        }
        return rescue;
    }

    /** @return the bases MAPQ-0 reads (primary, not duplicates) carry at a position */
    private static Set<Character> mq0BasesAt(final String contig, final int pos) {
        final Set<Character> bases = new HashSet<>();
        try (final ReadsDataSource reads = new ReadsPathDataSource(IOUtils.getPath(NA12878_20_21_WGS_bam))) {
            for (final Iterator<GATKRead> it = reads.query(new SimpleInterval(contig, pos, pos)); it.hasNext(); ) {
                final GATKRead r = it.next();
                if (r.getMappingQuality() == 0 && !r.isSecondaryAlignment() && !r.isDuplicate()) {
                    final int offset = ReadUtils.getReadIndexForReferenceCoordinate(r, pos).getLeft();
                    if (offset >= 0 && offset < r.getLength()) {
                        bases.add((char) r.getBase(offset));
                    }
                }
            }
        }
        return bases;
    }

    /**
     * The rescue records of GVCF and VCF output differ only as Mq0RescueArgumentCollection documents. A GVCF rescue
     * record with QUAL at least the calling threshold is in the VCF, with the same genotype but for phasing (the GVCF
     * record may list more alleles), and one below it is not,
     * since VCF mode does not make such calls. Each VCF rescue record the GVCF lacks is explained by a regular GVCF
     * record the VCF does not hold (a non-call record, or a call below the threshold) that overlaps it or, for a call,
     * lies within the rescue distance of it.
     */
    @Test
    public void testMq0RescueSitesInVcfAndGvcfModes() {
        // the calling threshold --dragen-mode sets
        final double threshold = Double.parseDouble(dragenModeValue(GenotypeCalculationArgumentCollection.CALL_CONFIDENCE_LONG_NAME));
        final int distance = new Mq0RescueArgumentCollection().distance;
        final List<VariantContext> gvcfRecords = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.GVCF, true)).collect(Collectors.toList());
        final List<VariantContext> vcfRecords = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.NONE, true)).collect(Collectors.toList());
        final Map<Integer, VariantContext> gvcf = gvcfRecords.stream().filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        final Map<Integer, VariantContext> vcf = vcfRecords.stream().filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        for (final VariantContext vc : gvcf.values()) {
            final VariantContext inVcf = vcf.get(vc.getStart());
            if (vc.getPhredScaledQual() >= threshold) {
                Assert.assertNotNull(inVcf, vc.toString());
                Assert.assertTrue(inVcf.getGenotype(0).sameGenotype(vc.getGenotype(0), true), vc + " " + inVcf);
            } else {
                Assert.assertNull(inVcf, vc.toString());
            }
        }
        final Set<Integer> vcfRegularStarts = vcfRecords.stream().filter(vc -> !isMq0RescueRecord(vc)).map(VariantContext::getStart).collect(Collectors.toSet());
        final List<VariantContext> gvcfOnlyRegular = gvcfRecords.stream()
                .filter(vc -> !isMq0RescueRecord(vc) && !vc.getAlternateAlleles().equals(List.of(Allele.NON_REF_ALLELE)) && !vcfRegularStarts.contains(vc.getStart()))
                .collect(Collectors.toList());
        final List<VariantContext> vcfOnly = vcf.values().stream().filter(vc -> !gvcf.containsKey(vc.getStart())).collect(Collectors.toList());
        Assert.assertFalse(vcf.isEmpty());
        for (final VariantContext vc : vcfOnly) {
            Assert.assertTrue(gvcfOnlyRegular.stream().anyMatch(regular -> regular.overlaps(vc) || (isCallForTest(regular) && distanceForTest(regular, vc) <= distance)),
                    vc + " is not explained by a GVCF-only regular record");
        }
    }

    /**
     * Records force-called from --alleles win over rescue records: no rescue record overlaps a given allele, whatever
     * the genotype of the record force-called there. A filtered --alleles record, which is not force-called without
     * --force-call-filtered-alleles, does not keep a rescue record out.
     */
    @Test
    public void testMq0RescueGivenAlleles() throws IOException {
        // given alleles at sites the rescue calls: at two of them an alt base other than the rescued one, which is
        // force-called 0/0, at the third the rescued alt itself, and at the fourth the rescued alt, filtered
        final List<VariantContext> rescuedSnps = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.GVCF, true))
                .filter(vc -> isMq0RescueRecord(vc) && vc.getReference().length() == 1
                        && vc.getAlternateAlleles().stream().filter(a -> !a.isSymbolic()).allMatch(a -> a.length() == 1))
                .limit(4).collect(Collectors.toList());
        Assert.assertEquals(rescuedSnps.size(), 4);
        final File given = createTempFile("mq0Rescue.given", ".vcf");
        final VCFHeader givenHeader = new VCFHeader();
        try (final ReferenceDataSource reference = ReferenceDataSource.of(IOUtils.getPath(b37_reference_20_21))) {
            givenHeader.setSequenceDictionary(reference.getSequenceDictionary());
        }
        givenHeader.addMetaDataLine(new VCFFilterHeaderLine("given", "a filtered given allele"));
        try (final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(given.toPath(), givenHeader.getSequenceDictionary(), false, Options.INDEX_ON_THE_FLY)) {
            writer.writeHeader(givenHeader);
            for (int i = 0; i < rescuedSnps.size(); i++) {
                final VariantContext vc = rescuedSnps.get(i);
                final Allele rescuedAlt = vc.getAlternateAllele(0);
                final Allele alt = i >= 2 ? rescuedAlt : Allele.create(List.of("A", "C", "G", "T").stream()
                        .filter(b -> !b.equals(vc.getReference().getBaseString()) && vc.getAlleles().stream().noneMatch(a -> a.getBaseString().equals(b)))
                        .findFirst().get());
                final VariantContextBuilder builder = new VariantContextBuilder("given", vc.getContig(), vc.getStart(), vc.getStart(), List.of(vc.getReference(), alt));
                writer.add((i == 3 ? builder.filter("given") : builder.passFilters()).make());
            }
        }
        final File output = createTempFile("mq0Rescue.givenAlleles", ".g.vcf");
        runCommandLine(mq0RescueArgs(output, ReferenceConfidenceMode.GVCF).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(AssemblyBasedCallerArgumentCollection.FORCE_CALL_ALLELES_LONG_NAME, given));
        final List<VariantContext> records = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        Assert.assertTrue(records.stream().filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).count() >= 10);
        Assert.assertTrue(records.stream().anyMatch(r -> isMq0RescueRecord(r) && r.getStart() == rescuedSnps.get(3).getStart()), "filtered given allele");
        for (final VariantContext vc : rescuedSnps.subList(0, 3)) {
            final List<VariantContext> atSite = records.stream().filter(r -> r.overlaps(vc)).collect(Collectors.toList());
            Assert.assertTrue(atSite.stream().noneMatch(HaplotypeCallerIntegrationTest::isMq0RescueRecord), atSite.toString());
            Assert.assertTrue(atSite.stream().anyMatch(r -> r.getStart() == vc.getStart() && !r.getAlternateAlleles().equals(List.of(Allele.NON_REF_ALLELE))),
                    "no force-called record at " + vc.getStart() + ": " + atSite);
        }
    }

    /**
     * The rescue draws its random numbers from its own generators: with --max-reads-per-alignment-start low enough that
     * the read stream and the rescue both downsample at random, the regular records are the same as without rescue.
     */
    @Test
    public void testMq0RescueRandomIsolation() throws IOException {
        final Map<Boolean, File> outputs = new HashMap<>();
        for (final boolean rescue : new boolean[]{false, true}) {
            outputs.put(rescue, createTempFile("mq0Rescue.random." + rescue, ".vcf"));
            final ArgumentsBuilder args = mq0RescueArgs(outputs.get(rescue), ReferenceConfidenceMode.NONE)
                    .add(AssemblyRegionArgumentCollection.MAX_STARTS_LONG_NAME, 2);
            if (rescue) {
                args.addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME);
            }
            Utils.resetRandomGenerator();
            runCommandLine(args);
        }
        // the regular records are the same, but for the non-call records rescue records replace
        assertMq0RescueOutput(outputs.get(false), outputs.get(true), ReferenceConfidenceMode.NONE, 1);
    }

    /**
     * The rescue works with the options that change how regions are called: --pileup-detection, --max-mnp-distance,
     * and --dragen-378-concordance-mode in place of --dragen-mode.
     */
    @DataProvider
    public Object[][] mq0RescueOptions() {
        return new Object[][]{
                {List.of("--" + HaplotypeCallerArgumentCollection.DRAGEN_378_GATK_MODE_LONG_NAME), 1},
                {List.of("--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME), 1},
                {List.of("--" + AssemblyBasedCallerArgumentCollection.MAX_MNP_DISTANCE_LONG_NAME, "2"), 1},
                {List.of("--" + Mq0RescueArgumentCollection.MQ0_RESCUE_MODE_LONG_NAME, "ASSEMBLE", "--" + PileupDetectionArgumentCollection.PILEUP_DETECTION_LONG_NAME), 1}};
    }

    @Test(dataProvider = "mq0RescueOptions")
    public void testMq0RescueWithOptions(final List<String> options, final int minRescue) throws IOException {
        // --dragen-378-concordance-mode takes the place of --dragen-mode, in VCF mode: in GVCF mode its regular calls
        // fail on this interval without the rescue too
        final boolean dragen378 = options.contains("--" + HaplotypeCallerArgumentCollection.DRAGEN_378_GATK_MODE_LONG_NAME);
        final ReferenceConfidenceMode mode = dragen378 ? ReferenceConfidenceMode.NONE : ReferenceConfidenceMode.GVCF;
        final Map<Boolean, File> outputs = new HashMap<>();
        for (final boolean rescue : new boolean[]{false, true}) {
            outputs.put(rescue, createTempFile("mq0Rescue.options." + rescue, dragen378 ? ".vcf" : ".g.vcf"));
            final ArgumentsBuilder args = dragen378 ? mq0RescueArgsWithoutDragen(outputs.get(rescue), mode, "21:10000000-10005000") :
                    mq0RescueArgs(outputs.get(rescue), mode, "21:10000000-10005000");
            options.forEach(args::addRaw);
            if (rescue) {
                args.addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME);
            }
            runCommandLine(args);
        }
        assertMq0RescueOutput(outputs.get(false), outputs.get(true), mode, minRescue);
    }

    /**
     * A rescue record at the last position of the run, decided after the last region (here, the only one).
     */
    @Test(dataProvider = "mq0RescueModes")
    public void testMq0RescueAtLastPosition(final ReferenceConfidenceMode mode) {
        final File output = createTempFile("mq0Rescue.lastPosition", ".vcf");
        runCommandLine(mq0RescueArgs(output, mode, "21:10017000-10017628").addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME));
        final List<VariantContext> records = VariantContextTestUtils.streamVcf(output).collect(Collectors.toList());
        final VariantContext last = records.get(records.size() - 1);
        Assert.assertTrue(isMq0RescueRecord(last) && last.getStart() == 10017628, last.toString());
    }

    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = ".*cannot be used with --population-callset.*")
    public void testMq0RescueRejectsPopulationCallset() {
        runCommandLine(mq0RescueArgs(createTempFile("mq0Rescue", ".g.vcf"), ReferenceConfidenceMode.GVCF)
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(GenotypeCalculationArgumentCollection.SUPPORTING_CALLSET_LONG_NAME, largeFileTestDir + "1000G.phase3.broad.withGenotypes.chr20.10100000.vcf"));
    }

    /**
     * The rescue writes nothing to the side outputs: the bamout, the assembled-variants output, the haplotype
     * histogram, the assembly regions and the genotyper debug output are the same as without rescue. (The assembly
     * graph output names vertices by object identity, so it differs from run to run; INJECT-mode rescue calls do not
     * assemble.)
     */
    @Test
    public void testMq0RescueSideOutputs() throws IOException {
        final List<String> textOutputs = List.of(ReadThreadingAssemblerArgumentCollection.DEBUG_ASSEMBLY_VARIANTS_OUT_LONG_NAME,
                ReadThreadingAssemblerArgumentCollection.HAPLOTYPE_DEBUG_HISTOGRAM_OUTPUT_LONG_NAME,
                AssemblyRegionArgumentCollection.ASSEMBLY_REGION_OUT_LONG_NAME, HaplotypeCallerArgumentCollection.DEBUG_GENOTYPER_OUTPUT_LONG_NAME);
        final Map<Boolean, File> bams = new HashMap<>();
        final Map<Boolean, Map<String, File>> texts = new HashMap<>();
        for (final boolean rescue : new boolean[]{false, true}) {
            bams.put(rescue, createTempFile("mq0Rescue.bamout." + rescue, ".bam"));
            final ArgumentsBuilder args = mq0RescueArgs(createTempFile("mq0Rescue.side." + rescue, ".g.vcf"), ReferenceConfidenceMode.GVCF, "21:10000000-10005000")
                    .add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_LONG_NAME, bams.get(rescue));
            final Map<String, File> files = new LinkedHashMap<>();
            for (final String name : textOutputs) {
                files.put(name, createTempFile("mq0Rescue." + name + "." + rescue, name.equals(ReadThreadingAssemblerArgumentCollection.DEBUG_ASSEMBLY_VARIANTS_OUT_LONG_NAME) ? ".vcf" : ".txt"));
                args.add(name, files.get(name));
            }
            texts.put(rescue, files);
            if (rescue) {
                args.addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME);
            }
            runCommandLine(args);
        }
        for (final String name : textOutputs) {
            final List<String> withRescue = Files.readAllLines(texts.get(true).get(name).toPath());
            Assert.assertFalse(withRescue.isEmpty(), name);
            Assert.assertEquals(withRescue, Files.readAllLines(texts.get(false).get(name).toPath()), name);
        }
        final List<String> withRescue = samRecords(bams.get(true));
        Assert.assertFalse(withRescue.isEmpty());
        Assert.assertEquals(withRescue, samRecords(bams.get(false)));
    }

    private static List<String> samRecords(final File bam) {
        try (final ReadsDataSource reads = new ReadsPathDataSource(bam.toPath())) {
            return StreamSupport.stream(reads.spliterator(), false).map(r -> r.convertToSAMRecord(reads.getHeader()).getSAMString())
                    .collect(Collectors.toList());
        }
    }

    /**
     * With --sample-name on a multi-sample input, the rescue uses only that sample's reads, MAPQ-0 and regular: its depths
     * are those of a run that filters the other samples' reads out.
     */
    @Test
    public void testMq0RescueSampleName() {
        final String bam = largeFileTestDir + "multiSampleSubsetted.bam";
        final File sampleName = createTempFile("mq0Rescue.sampleName", ".vcf");
        runCommandLine(new ArgumentsBuilder().addInput(bam).addReference(hg38Reference).addInterval("chr1:1590000-1600000")
                .addFlag(HaplotypeCallerArgumentCollection.DRAGEN_3412_GATK_MODE_LONG_NAME)
                .addOutput(sampleName).add("sample-name", "NA12878").addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false));
        final File filtered = createTempFile("mq0Rescue.sampleFilter", ".vcf");
        runCommandLine(new ArgumentsBuilder().addInput(bam).addReference(hg38Reference).addInterval("chr1:1590000-1600000")
                .addFlag(HaplotypeCallerArgumentCollection.DRAGEN_3412_GATK_MODE_LONG_NAME)
                .addOutput(filtered).add("read-filter", "SampleReadFilter").add("sample", "NA12878")
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME).add(StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE, false));
        final Map<Integer, VariantContext> expected = VariantContextTestUtils.streamVcf(filtered)
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        final List<VariantContext> actual = VariantContextTestUtils.streamVcf(sampleName).filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .filter(vc -> expected.containsKey(vc.getStart())).collect(Collectors.toList());
        Assert.assertFalse(actual.isEmpty());
        for (final VariantContext vc : actual) {
            for (final String key : List.of(GATKVCFConstants.MQ0_RESCUE_DEPTH_KEY, GATKVCFConstants.MQ0_RESCUE_MQ0_DEPTH_KEY)) {
                Assert.assertEquals(vc.getAttributeAsInt(key, -1), expected.get(vc.getStart()).getAttributeAsInt(key, -2), vc.toString());
            }
        }
    }

    /**
     * Known rescue calls: heterozygous SNPs where MAPQ-0 reads (MQ0R_DP0) make up a good part of the depth next to the
     * reads with MAPQ above 0 (MQ0R_DP1).
     */
    @Test
    public void testMq0RescueCalls() {
        final Map<Integer, VariantContext> rescue = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.GVCF, true))
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        // position, called alt, minimum QUAL, MQ0R_DP1, MQ0R_DP0
        final Object[][] expected = {
                {10003018, "G", 25.0, 37, 29},
                {10005921, "T", 25.0, 66, 162},
                {10007597, "G", 25.0, 65, 41},
                {10011449, "T", 15.0, 33, 26}};
        for (final Object[] e : expected) {
            final VariantContext vc = rescue.get((Integer) e[0]);
            Assert.assertNotNull(vc, "no rescue record at " + e[0]);
            final Genotype g = vc.getGenotype(0);
            Assert.assertTrue(g.isHet(), vc.toString());
            Assert.assertTrue(g.getAlleles().stream().anyMatch(a -> a.getBaseString().equals(e[1])), vc.toString());
            Assert.assertTrue(vc.getPhredScaledQual() >= (Double) e[2], vc.toString());
            Assert.assertEquals(vc.getAttributeAsInt(GATKVCFConstants.MQ0_RESCUE_DEPTH_KEY, -1), e[3], vc.toString());
            Assert.assertEquals(vc.getAttributeAsInt(GATKVCFConstants.MQ0_RESCUE_MQ0_DEPTH_KEY, -1), e[4], vc.toString());
        }
    }

    /**
     * The rescue genotyping does not cap base qualities at MAPQ even if the cap that --dragen-mode turns off is turned
     * back on: without that, every MAPQ-0 base would score at the minimum base quality. The rescue records the two runs
     * share have the same QUAL (the capped regular calls differ, so the runs can rescue different sites).
     */
    @Test
    public void testMq0RescueKeepsBaseQualitiesUncapped() {
        final File capped = createTempFile("mq0Rescue.capped", ".g.vcf");
        runCommandLine(mq0RescueArgs(capped, ReferenceConfidenceMode.GVCF).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(LikelihoodEngineArgumentCollection.DISABLE_CAP_BASE_QUALITIES_TO_MAP_QUALITY_LONG_NAME, false));
        final Map<Integer, VariantContext> uncapped = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.GVCF, true))
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).collect(Collectors.toMap(VariantContext::getStart, vc -> vc));
        final List<VariantContext> shared = VariantContextTestUtils.streamVcf(capped).filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .filter(vc -> uncapped.containsKey(vc.getStart())).collect(Collectors.toList());
        Assert.assertTrue(shared.size() >= 10, shared.size() + " shared rescue records");
        shared.forEach(vc -> Assert.assertEquals(vc.getPhredScaledQual(), uncapped.get(vc.getStart()).getPhredScaledQual(), 0.01, vc.toString()));
    }

    /**
     * With --ploidy-regions, rescue calls use the region's ploidy; with several intervals on two contigs (rescue records on
     * the second only, as the test BAM's MAPQ-0 reads on 20 make no DRAGEN-mode rescue call), the records of both come out
     * in order; and the same run gives the same output twice.
     */
    @Test
    public void testMq0RescuePloidyContigsAndDeterminism() throws IOException {
        final File ploidyBed = createTempFile("mq0Rescue.ploidy", ".bed");
        Files.writeString(ploidyBed.toPath(), "21\t9999999\t10020000\t1\n");
        final File haploid = createTempFile("mq0Rescue.haploid", ".g.vcf");
        runCommandLine(mq0RescueArgs(haploid, ReferenceConfidenceMode.GVCF).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(HaplotypeCallerArgumentCollection.PLOIDY_REGIONS_NAME, ploidyBed));
        final List<VariantContext> haploidRescue = VariantContextTestUtils.streamVcf(haploid)
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).collect(Collectors.toList());
        Assert.assertFalse(haploidRescue.isEmpty());
        haploidRescue.forEach(vc -> Assert.assertEquals(vc.getGenotype(0).getPloidy(), 1, vc.toString()));

        final File twoContigs = createTempFile("mq0Rescue.twoContigs", ".g.vcf");
        runCommandLine(mq0RescueArgs(twoContigs, ReferenceConfidenceMode.GVCF, "20:10130000-10150000", "21:10000000-10005000")
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME));
        final List<VariantContext> records = VariantContextTestUtils.streamVcf(twoContigs).collect(Collectors.toList());
        final List<String> rescueContigs = records.stream().filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord)
                .map(VariantContext::getContig).distinct().collect(Collectors.toList());
        Assert.assertEquals(rescueContigs, List.of("21"));
        Assert.assertEquals(records.stream().map(VariantContext::getContig).distinct().collect(Collectors.toList()), List.of("20", "21"));
        for (int i = 1; i < records.size(); i++) {
            final VariantContext a = records.get(i - 1);
            final VariantContext b = records.get(i);
            Assert.assertTrue(Integer.parseInt(a.getContig()) < Integer.parseInt(b.getContig())
                    || (a.getContig().equals(b.getContig()) && a.getStart() < b.getStart()), b.toString());
        }

        final File again = createTempFile("mq0Rescue.again", ".g.vcf");
        runCommandLine(mq0RescueArgs(again, ReferenceConfidenceMode.GVCF).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME));
        IntegrationTestSpec.assertEqualTextFiles(again, mq0RescueRun(ReferenceConfidenceMode.GVCF, true));
    }

    /**
     * ASSEMBLE mode makes rescue records INJECT cannot (here, a SNP that INJECT does not insert), with the same
     * guarantees, in every reference-confidence mode (all in DRAGEN mode, which the rescue requires): non-reference calls whose alleles MAPQ-0 reads carry, overlapping
     * no regular call. Its assemblies add nothing to the side outputs: the rescue assembler writes nothing to the
     * haplotype histogram, and rescue calls nothing to the assembled-variants output or the bamout. (The assembly graph
     * output names vertices by object identity, so it differs from run to run and is not compared.) The same run gives
     * the same output, and --mq0-rescue-max-event-sets changes nothing in ASSEMBLE mode.
     */
    @Test
    public void testMq0RescueAssemble() throws IOException {
        final Map<Boolean, File> outputs = new HashMap<>();
        final Map<Boolean, List<File>> sideOutputs = new HashMap<>();
        for (final boolean assemble : new boolean[]{false, true}) {
            outputs.put(assemble, createTempFile("mq0Rescue.assemble." + assemble, ".g.vcf"));
            sideOutputs.put(assemble, List.of(createTempFile("mq0Rescue.assembledVariants." + assemble, ".vcf"),
                    createTempFile("mq0Rescue.histogram." + assemble, ".txt"), createTempFile("mq0Rescue.bamout." + assemble, ".bam")));
            final ArgumentsBuilder args = mq0RescueArgs(outputs.get(assemble), ReferenceConfidenceMode.GVCF)
                    .add(ReadThreadingAssemblerArgumentCollection.DEBUG_ASSEMBLY_VARIANTS_OUT_LONG_NAME, sideOutputs.get(assemble).get(0))
                    .add(ReadThreadingAssemblerArgumentCollection.HAPLOTYPE_DEBUG_HISTOGRAM_OUTPUT_LONG_NAME, sideOutputs.get(assemble).get(1))
                    .add(AssemblyBasedCallerArgumentCollection.BAM_OUTPUT_LONG_NAME, sideOutputs.get(assemble).get(2));
            if (assemble) {
                args.addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME).add(Mq0RescueArgumentCollection.MQ0_RESCUE_MODE_LONG_NAME, Mq0RescueArgumentCollection.Mode.ASSEMBLE);
            }
            runCommandLine(args);
        }
        final File assembled = outputs.get(true);
        final List<VariantContext> assembleRescue = assertMq0RescueOutput(mq0RescueRun(ReferenceConfidenceMode.GVCF, false), assembled, ReferenceConfidenceMode.GVCF);
        final Set<Integer> injectStarts = VariantContextTestUtils.streamVcf(mq0RescueRun(ReferenceConfidenceMode.GVCF, true))
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).map(VariantContext::getStart).collect(Collectors.toSet());
        final VariantContext assembleOnly = assembleRescue.stream().filter(vc -> vc.getStart() == 10007110).findFirst().orElse(null);
        Assert.assertNotNull(assembleOnly);
        Assert.assertFalse(injectStarts.contains(10007110));
        Assert.assertEquals(assembleOnly.getReference().getBaseString(), MQ0_RESCUE_ASSEMBLE_ONLY_REF, assembleOnly.toString());
        Assert.assertTrue(assembleOnly.getGenotype(0).isHet() && assembleOnly.getGenotype(0).getAlleles().stream()
                .anyMatch(a -> a.getBaseString().equals(MQ0_RESCUE_ASSEMBLE_ONLY_ALT)), assembleOnly.toString());
        Assert.assertTrue(assembleOnly.getPhredScaledQual() >= 15, assembleOnly.toString());

        Assert.assertEquals(dataLines(sideOutputs.get(true).get(0)), dataLines(sideOutputs.get(false).get(0)));
        IntegrationTestSpec.assertEqualTextFiles(sideOutputs.get(true).get(1), sideOutputs.get(false).get(1));
        final List<String> bamout = samRecords(sideOutputs.get(true).get(2));
        Assert.assertFalse(bamout.isEmpty());
        Assert.assertEquals(bamout, samRecords(sideOutputs.get(false).get(2)));

        final File again = createTempFile("mq0Rescue.assemble.again", ".g.vcf");
        runCommandLine(mq0RescueArgs(again, ReferenceConfidenceMode.GVCF).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                .add(Mq0RescueArgumentCollection.MQ0_RESCUE_MODE_LONG_NAME, Mq0RescueArgumentCollection.Mode.ASSEMBLE)
                .add(Mq0RescueArgumentCollection.MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME, 0));
        IntegrationTestSpec.assertEqualTextFiles(again, assembled);

        for (final ReferenceConfidenceMode mode : List.of(ReferenceConfidenceMode.NONE, ReferenceConfidenceMode.BP_RESOLUTION)) {
            final File output = createTempFile("mq0Rescue.assemble." + mode, ".vcf");
            runCommandLine(mq0RescueArgs(output, mode).addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME)
                    .add(Mq0RescueArgumentCollection.MQ0_RESCUE_MODE_LONG_NAME, Mq0RescueArgumentCollection.Mode.ASSEMBLE));
            assertMq0RescueOutput(mq0RescueRun(mode, false), output, mode);
        }
    }

    // the SNP at 10007110 that ASSEMBLE mode rescues and INJECT mode does not insert
    private static final String MQ0_RESCUE_ASSEMBLE_ONLY_REF = "A";
    private static final String MQ0_RESCUE_ASSEMBLE_ONLY_ALT = "G";

    @Test(expectedExceptions = CommandLineException.class, expectedExceptionsMessageRegExp = ".*mq0-rescue-mode.*must be INJECT or ASSEMBLE.*")
    public void testMq0RescueRejectsNullMode() {
        runCommandLine(mq0RescueArgs(createTempFile("mq0Rescue", ".g.vcf"), ReferenceConfidenceMode.GVCF)
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME).add(Mq0RescueArgumentCollection.MQ0_RESCUE_MODE_LONG_NAME, "null"));
    }

    /**
     * The rescue GVCF is valid as a GVCF and is accepted by the GVCF tools downstream.
     */
    @Test
    public void testMq0RescueGvcfIsValid() {
        final File gvcf = mq0RescueRun(ReferenceConfidenceMode.GVCF, true);
        runCommandLine(new ArgumentsBuilder().addReference(b37_reference_20_21).addVCF(gvcf).addInterval(MQ0_RESCUE_INTERVAL)
                .addFlag(ValidateVariants.GVCF_VALIDATE), ValidateVariants.class.getSimpleName());
        final File genotyped = createTempFile("mq0Rescue.genotyped", ".vcf");
        runCommandLine(new ArgumentsBuilder().addReference(b37_reference_20_21).addVCF(gvcf).addOutput(genotyped),
                GenotypeGVCFs.class.getSimpleName());
        final File combined = createTempFile("mq0Rescue.combined", ".g.vcf");
        runCommandLine(new ArgumentsBuilder().addReference(b37_reference_20_21).addVCF(gvcf).addVCF(mq0RescueRun(ReferenceConfidenceMode.GVCF, false))
                .addOutput(combined), CombineGVCFs.class.getSimpleName());
        final File reblocked = createTempFile("mq0Rescue.reblocked", ".g.vcf");
        runCommandLine(new ArgumentsBuilder().addReference(b37_reference_20_21).addVCF(gvcf).addOutput(reblocked),
                ReblockGVCF.class.getSimpleName());
        // The rescue calls come out of joint genotyping as calls. GenotypeGVCFs drops their FILTER and MQ0R flag but
        // keeps the depths MQ0R_DP0 and MQ0R_DP1, the markers left.
        final List<VariantContext> rescue = VariantContextTestUtils.readEntireVCFIntoMemory(gvcf.getAbsolutePath()).getRight().stream()
                .filter(HaplotypeCallerIntegrationTest::isMq0RescueRecord).collect(Collectors.toList());
        final Map<Integer, VariantContext> genotypedByStart = VariantContextTestUtils.streamVcf(genotyped)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc, (a, b) -> a));
        final List<VariantContext> genotypedRescue = rescue.stream().map(vc -> genotypedByStart.get(vc.getStart()))
                .filter(Objects::nonNull).collect(Collectors.toList());
        Assert.assertFalse(genotypedRescue.isEmpty());
        genotypedRescue.forEach(vc -> {
            Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_MQ0_DEPTH_KEY) && vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_DEPTH_KEY), vc.toString());
            Assert.assertFalse(vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_KEY), vc.toString());
            Assert.assertFalse(vc.isFiltered(), vc.toString());
        });
        // CombineGVCFs keeps the depths too, and drops the FILTER and the flag
        final Map<Integer, VariantContext> combinedByStart = VariantContextTestUtils.streamVcf(combined)
                .collect(Collectors.toMap(VariantContext::getStart, vc -> vc, (a, b) -> a));
        rescue.forEach(r -> {
            final VariantContext vc = combinedByStart.get(r.getStart());
            Assert.assertNotNull(vc, r.toString());
            Assert.assertTrue(vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_MQ0_DEPTH_KEY) && vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_DEPTH_KEY), vc.toString());
            Assert.assertFalse(vc.hasAttribute(GATKVCFConstants.MQ0_RESCUE_KEY), vc.toString());
            Assert.assertFalse(vc.isFiltered(), vc.toString());
        });
    }

    /** The rescue is validated only in DRAGEN mode, so it requires --dragen-mode (or --dragen-378-concordance-mode). */
    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = ".*mq0-rescue is validated only in DRAGEN mode and requires --dragen-mode.*")
    public void testMq0RescueRequiresDragenMode() {
        runCommandLine(mq0RescueArgsWithoutDragen(createTempFile("mq0Rescue", ".g.vcf"), ReferenceConfidenceMode.GVCF)
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME));
    }

    /** Without --mq0-rescue, HaplotypeCaller runs without --dragen-mode as before, warning about rescue-only arguments. */
    @Test
    public void testMq0RescueArgumentsWithoutRescueOrDragen() {
        runCommandLine(mq0RescueArgsWithoutDragen(createTempFile("mq0Rescue.off", ".vcf"), ReferenceConfidenceMode.NONE, "21:10000000-10001000")
                .add(Mq0RescueArgumentCollection.MQ0_RESCUE_DISTANCE_LONG_NAME, 20));
    }

    @Test(expectedExceptions = UserException.class, expectedExceptionsMessageRegExp = ".*not supported by RampedHaplotypeCaller.*")
    public void testMq0RescueRejectedByRampedHaplotypeCaller() {
        runCommandLine(mq0RescueArgs(createTempFile("mq0Rescue", ".g.vcf"), ReferenceConfidenceMode.GVCF)
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME), RampedHaplotypeCaller.class.getSimpleName());
    }

    @Test(expectedExceptions = CommandLineException.class)
    public void testMq0RescueNotAnArgumentOfHaplotypeCallerSpark() {
        runCommandLine(mq0RescueArgs(createTempFile("mq0Rescue", ".g.vcf"), ReferenceConfidenceMode.GVCF)
                .addFlag(Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME), HaplotypeCallerSpark.class.getSimpleName());
    }
}
