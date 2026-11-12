package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Collections;
import java.util.List;

/**
 * Unit tests for {@link Mq0RescueArgumentCollection}, whose rescue-only arguments HaplotypeCaller warns about when they are
 * given without --mq0-rescue.
 */
public class Mq0RescueArgumentCollectionUnitTest extends GATKBaseTest {

    private static CommandLineArgumentParser parse(final String... args) {
        final CommandLineArgumentParser parser = new CommandLineArgumentParser(new Mq0RescueArgumentCollection(),
                Collections.emptyList(), Collections.emptySet());
        Assert.assertTrue(parser.parseArguments(System.err, args));
        return parser;
    }

    @Test
    public void testGivenRescueOnlyArguments() {
        Assert.assertEquals(Mq0RescueArgumentCollection.givenRescueOnlyArguments(parse()), List.of());
        // given on the command line counts, even at the default value; --mq0-rescue itself is not a tuning argument
        Assert.assertEquals(Mq0RescueArgumentCollection.givenRescueOnlyArguments(parse("--" + Mq0RescueArgumentCollection.MQ0_RESCUE_LONG_NAME,
                        "--" + Mq0RescueArgumentCollection.MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME, "2",
                        "--" + Mq0RescueArgumentCollection.MQ0_RESCUE_DISTANCE_LONG_NAME, "10")),
                List.of(Mq0RescueArgumentCollection.MQ0_RESCUE_DISTANCE_LONG_NAME, Mq0RescueArgumentCollection.MQ0_RESCUE_MAX_EVENT_SETS_LONG_NAME));
        // every rescue-only argument is an argument of the collection
        for (final String name : Mq0RescueArgumentCollection.RESCUE_ONLY_ARGUMENTS) {
            Assert.assertNotNull(parse().getNamedArgumentDefinitionByAlias(name), name);
        }
    }
}
