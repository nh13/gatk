package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.MultiDeBruijnVertex;
import org.jgrapht.alg.CycleDetector;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Tests that a vertex's pointer to its edges never answers for the wrong graph: vertices shared by cloned graphs,
 * removed and re-added vertices, and lookups through a vertex equal to, but not the same as, the one a graph holds.
 */
public final class AssemblyGraphSpecificsUnitTest extends GATKBaseTest {

    @Test
    public void aVertexSharedByAClonedGraphAnswersForEachGraph() {
        final SeqGraph g = new SeqGraph(11);
        final SeqVertex a = new SeqVertex("A");
        final SeqVertex b = new SeqVertex("C");
        final SeqVertex c = new SeqVertex("G");
        g.addVertices(a, b, c);
        g.addEdge(a, b);
        final SeqGraph copy = g.clone();
        copy.addEdge(a, c);

        // Alternate between the graphs, so the vertices' pointers keep moving from one graph to the other.
        for (int i = 0; i < 3; i++) {
            Assert.assertEquals(g.outDegreeOf(a), 1);
            Assert.assertEquals(copy.outDegreeOf(a), 2);
            Assert.assertEquals(g.inDegreeOf(c), 0);
            Assert.assertEquals(copy.inDegreeOf(c), 1);
            Assert.assertTrue(g.containsEdge(a, b));
            Assert.assertFalse(g.containsEdge(a, c));
            Assert.assertTrue(copy.containsEdge(a, c));
        }

        copy.removeVertex(b);
        Assert.assertFalse(copy.containsVertex(b));
        Assert.assertTrue(g.containsVertex(b));
        Assert.assertEquals(g.inDegreeOf(b), 1);
        Assert.assertThrows(IllegalArgumentException.class, () -> copy.inDegreeOf(b));
    }

    @Test
    public void aRemovedVertexIsNotInTheGraphUntilAddedAgain() {
        final SeqGraph g = new SeqGraph(11);
        final SeqVertex a = new SeqVertex("A");
        final SeqVertex b = new SeqVertex("C");
        g.addVertices(a, b);
        g.addEdge(a, b);
        Assert.assertTrue(g.removeVertex(b));
        Assert.assertFalse(g.containsVertex(b));
        Assert.assertFalse(g.vertexSet().contains(b));
        Assert.assertThrows(IllegalArgumentException.class, () -> g.outgoingEdgesOf(b));
        Assert.assertEquals(g.outDegreeOf(a), 0);
        Assert.assertFalse(g.removeVertex(b));

        Assert.assertTrue(g.addVertex(b));
        Assert.assertEquals(g.inDegreeOf(b), 0);
        g.addEdge(b, a);
        Assert.assertEquals(g.inDegreeOf(a), 1);
        Assert.assertEquals(new ArrayList<>(g.vertexSet()), Arrays.asList(a, b));
    }

    @Test
    public void aVertexEqualToOneInTheGraphFindsThatVertexsEdges() {
        // Vertices that merge identical nodes are equal when their sequences are.
        final TestGraph g = new TestGraph(3);
        final MultiDeBruijnVertex held = new MultiDeBruijnVertex("ACG".getBytes(), true);
        final MultiDeBruijnVertex next = new MultiDeBruijnVertex("CGT".getBytes(), true);
        final MultiDeBruijnVertex equalToHeld = new MultiDeBruijnVertex("ACG".getBytes(), true);
        Assert.assertEquals(equalToHeld, held);
        g.addVertices(held, next);
        g.addEdge(held, next);

        Assert.assertTrue(g.containsVertex(equalToHeld));
        Assert.assertFalse(g.addVertex(equalToHeld));
        Assert.assertEquals(g.outDegreeOf(equalToHeld), 1);
        Assert.assertTrue(g.containsEdge(equalToHeld, next));
        Assert.assertSame(g.getEdgeSource(g.getEdge(equalToHeld, next)), held);

        // Removing the vertex through an equal object takes the held vertex out of the graph.
        Assert.assertTrue(g.removeVertex(equalToHeld));
        Assert.assertFalse(g.containsVertex(held));
        Assert.assertThrows(IllegalArgumentException.class, () -> g.outDegreeOf(held));
        Assert.assertEquals(g.inDegreeOf(next), 0);
    }

    @Test
    public void aVertexEqualToOneInAnotherGraphFindsItsOwnEdges() {
        final TestGraph g1 = new TestGraph(3);
        final TestGraph g2 = new TestGraph(3);
        final MultiDeBruijnVertex x1 = new MultiDeBruijnVertex("ACG".getBytes(), true);
        final MultiDeBruijnVertex x2 = new MultiDeBruijnVertex("ACG".getBytes(), true);
        final MultiDeBruijnVertex y = new MultiDeBruijnVertex("CGT".getBytes(), true);
        g1.addVertices(x1, y);
        g1.addEdge(x1, y);
        g2.addVertex(x2);

        Assert.assertEquals(g1.outDegreeOf(x2), 1);
        Assert.assertEquals(g2.outDegreeOf(x1), 0);
        Assert.assertEquals(g2.outDegreeOf(x2), 0);
        Assert.assertEquals(g1.outDegreeOf(x1), 1);
        g2.removeVertex(x1);
        Assert.assertFalse(g2.containsVertex(x2));
        Assert.assertTrue(g1.containsVertex(x1));
        Assert.assertEquals(g1.outDegreeOf(x2), 1);
    }

    /**
     * Random edits applied to a graph and to its clone, compared after every step with an adjacency model kept in
     * lists, which also checks that vertices and edges stay in insertion order.
     */
    @Test
    public void randomEditsOfAGraphAndItsCloneMatchAListModel() {
        final Random random = new Random(7);
        final List<SeqVertex> pool = new ArrayList<>();
        for (int i = 0; i < 12; i++) {
            pool.add(new SeqVertex("A"));
        }
        final SeqGraph original = new SeqGraph(11);
        final GraphModel originalModel = new GraphModel();
        for (int i = 0; i < 6; i++) {
            original.addVertex(pool.get(i));
            originalModel.addVertex(pool.get(i));
        }
        SeqGraph copy = null;
        GraphModel copyModel = null;
        for (int step = 0; step < 3000; step++) {
            if (step == 1000) {
                copy = original.clone();
                copyModel = originalModel.copy();
            }
            final boolean editCopy = copy != null && random.nextBoolean();
            final SeqGraph g = editCopy ? copy : original;
            final GraphModel model = editCopy ? copyModel : originalModel;
            final SeqVertex u = pool.get(random.nextInt(pool.size()));
            final SeqVertex v = pool.get(random.nextInt(pool.size()));
            switch (random.nextInt(4)) {
                case 0:
                    Assert.assertEquals(g.addVertex(u), model.addVertex(u));
                    break;
                case 1:
                    Assert.assertEquals(g.removeVertex(u), model.removeVertex(u));
                    break;
                case 2:
                    if (model.vertices.contains(u) && model.vertices.contains(v)) {
                        Assert.assertEquals(g.addEdge(u, v) != null, model.addEdge(u, v));
                    }
                    break;
                default:
                    if (model.vertices.contains(u) && model.vertices.contains(v)) {
                        Assert.assertEquals(g.removeEdge(u, v) != null, model.removeEdge(u, v));
                    }
                    break;
            }
            model.assertMatches(g, pool);
            if (copy != null) {
                (editCopy ? originalModel : copyModel).assertMatches(editCopy ? original : copy, pool);
            }
        }
    }

    @Test
    public void theEdgeListsAreReadOnlyViewsInIterationOrder() {
        final SeqGraph g = new SeqGraph(11);
        final SeqVertex a = new SeqVertex("A");
        final SeqVertex b = new SeqVertex("C");
        final SeqVertex c = new SeqVertex("G");
        g.addVertices(a, b, c);
        g.addEdge(a, c);
        g.addEdge(a, b);
        g.addEdge(b, c);
        Assert.assertEquals(g.outgoingEdgeList(a), new ArrayList<>(g.outgoingEdgesOf(a)));
        Assert.assertEquals(g.incomingEdgeList(c), new ArrayList<>(g.incomingEdgesOf(c)));
        Assert.assertThrows(UnsupportedOperationException.class, () -> g.outgoingEdgeList(a).remove(0));
        Assert.assertThrows(UnsupportedOperationException.class, () -> g.incomingEdgeList(c).add(new BaseEdge(false, 1)));
        // a view, so it follows later edits
        final List<BaseEdge> out = g.outgoingEdgeList(b);
        g.addEdge(b, a);
        Assert.assertEquals(out.size(), 2);
    }

    /** Vertices and edges as lists of vertices and of (source, target) pairs, in insertion order. */
    private static final class GraphModel {
        final List<SeqVertex> vertices = new ArrayList<>();
        final List<SeqVertex[]> edges = new ArrayList<>();

        GraphModel copy() {
            final GraphModel copy = new GraphModel();
            copy.vertices.addAll(vertices);
            copy.edges.addAll(edges);
            return copy;
        }

        boolean addVertex(final SeqVertex v) {
            if (vertices.contains(v)) {
                return false;
            }
            vertices.add(v);
            return true;
        }

        boolean removeVertex(final SeqVertex v) {
            if (!vertices.remove(v)) {
                return false;
            }
            edges.removeIf(e -> e[0] == v || e[1] == v);
            return true;
        }

        boolean addEdge(final SeqVertex u, final SeqVertex v) {
            if (edges.stream().anyMatch(e -> e[0] == u && e[1] == v)) {
                return false;
            }
            edges.add(new SeqVertex[]{u, v});
            return true;
        }

        boolean removeEdge(final SeqVertex u, final SeqVertex v) {
            return edges.removeIf(e -> e[0] == u && e[1] == v);
        }

        void assertMatches(final SeqGraph g, final List<SeqVertex> pool) {
            Assert.assertEquals(new ArrayList<>(g.vertexSet()), vertices);
            Assert.assertEquals(g.edgeSet().size(), edges.size());
            for (final SeqVertex v : pool) {
                Assert.assertEquals(g.containsVertex(v), vertices.contains(v));
                if (!vertices.contains(v)) {
                    Assert.assertThrows(IllegalArgumentException.class, () -> g.outDegreeOf(v));
                    continue;
                }
                final List<SeqVertex> targets = new ArrayList<>();
                final List<SeqVertex> sources = new ArrayList<>();
                for (final SeqVertex[] e : edges) {
                    if (e[0] == v) {
                        targets.add(e[1]);
                    }
                    if (e[1] == v) {
                        sources.add(e[0]);
                    }
                }
                final List<SeqVertex> graphTargets = new ArrayList<>();
                g.outgoingEdgesOf(v).forEach(e -> graphTargets.add(g.getEdgeTarget(e)));
                final List<SeqVertex> graphSources = new ArrayList<>();
                g.incomingEdgesOf(v).forEach(e -> graphSources.add(g.getEdgeSource(e)));
                Assert.assertEquals(graphTargets, targets);
                Assert.assertEquals(graphSources, sources);
                Assert.assertEquals(g.outDegreeOf(v), targets.size());
                Assert.assertEquals(g.inDegreeOf(v), sources.size());
                // a self-loop is one of the vertex's edges once, but adds to both its in- and out-degree
                Assert.assertEquals(g.degreeOf(v), targets.size() + sources.size());
                final List<SeqVertex[]> touching = new ArrayList<>();
                for (final SeqVertex[] e : edges) {
                    if (e[0] == v || e[1] == v) {
                        touching.add(e);
                    }
                }
                final Set<BaseEdge> edgesOfV = g.edgesOf(v);
                Assert.assertEquals(edgesOfV.size(), touching.size());
                for (final SeqVertex[] e : touching) {
                    Assert.assertTrue(edgesOfV.stream().anyMatch(edge -> g.getEdgeSource(edge) == e[0] && g.getEdgeTarget(edge) == e[1]));
                }
                for (final SeqVertex w : pool) {
                    final Set<BaseEdge> between = g.getAllEdges(v, w);
                    if (!vertices.contains(w)) {
                        Assert.assertNull(between);
                        continue;
                    }
                    final boolean joined = edges.stream().anyMatch(e -> e[0] == v && e[1] == w);
                    Assert.assertEquals(between.size(), joined ? 1 : 0);
                    for (final BaseEdge edge : between) {
                        Assert.assertSame(g.getEdgeSource(edge), v);
                        Assert.assertSame(g.getEdgeTarget(edge), w);
                    }
                }
            }
            for (final SeqVertex v : pool) {
                if (!vertices.contains(v)) {
                    for (final SeqVertex w : pool) {
                        Assert.assertNull(g.getAllEdges(v, w));
                    }
                }
            }
        }
    }

    /** The cycle check agrees with jgrapht's cycle detector on random graphs, self-loops and removals included. */
    @Test
    public void hasCyclesAgreesWithJgraphtsCycleDetector() {
        final Random random = new Random(11);
        int cyclic = 0;
        for (int trial = 0; trial < 3000; trial++) {
            final SeqGraph g = new SeqGraph(11);
            final int n = 1 + random.nextInt(12);
            final List<SeqVertex> vertices = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                vertices.add(new SeqVertex("A"));
                g.addVertex(vertices.get(i));
            }
            final int nEdges = random.nextInt(2 * n);
            for (int e = 0; e < nEdges; e++) {
                final int from = random.nextInt(n);
                // mostly forward edges, so that both acyclic and cyclic graphs are common
                final int to = random.nextInt(8) == 0 ? random.nextInt(n) : from + random.nextInt(n - from);
                g.addEdge(vertices.get(from), vertices.get(to));
            }
            if (random.nextInt(4) == 0) {
                g.removeVertex(vertices.get(random.nextInt(n)));
            }
            final boolean expected = new CycleDetector<>(g).detectCycles();
            Assert.assertEquals(g.hasCycles(), expected);
            cyclic += expected ? 1 : 0;
        }
        Assert.assertTrue(cyclic > 100 && cyclic < 2900, "cyclic graphs: " + cyclic);
    }
}
