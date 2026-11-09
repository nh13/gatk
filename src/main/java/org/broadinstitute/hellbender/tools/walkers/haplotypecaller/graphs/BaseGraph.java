package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import com.google.common.collect.Sets;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.jgrapht.EdgeFactory;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.specifics.Specifics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Common code for graphs used for local assembly.
 *
 * A graph never holds two edges between the same pair of vertices, and an edge object joins one pair of vertices for
 * life (see {@link BaseEdge}). The graph refuses parallel edges itself rather than through jgrapht, so jgrapht's
 * {@link #isAllowingMultipleEdges()} and {@link #getType()} report that parallel edges are allowed.
 *
 * Each vertex keeps a pointer to its edges in the graph that most recently added it or looked it up (see
 * {@link BaseVertex}), so even read-only queries such as {@link #outDegreeOf} write to the vertex. A graph, its clones
 * and any other graph sharing its vertices must therefore be used by one thread at a time. The pointer also keeps the
 * last graph that used a vertex reachable for as long as the vertex is, so vertices should not outlive their graphs.
 */
public abstract class BaseGraph<V extends BaseVertex, E extends BaseEdge> extends AbstractBaseGraph<V, E> {
    private static final long serialVersionUID = 1l;
    protected final int kmerSize;

    /** The specifics jgrapht stores this graph in, kept so adjacency queries can reach them directly. */
    private AssemblyGraphSpecifics<V, E> assemblySpecifics;

    /**
     * Construct a TestGraph with kmerSize
     * @param kmerSize
     */
    protected BaseGraph(final int kmerSize, final EdgeFactory<V,E> edgeFactory) {
        // A directed, unweighted graph with loops. jgrapht is told parallel edges are allowed only so that it never
        // searches for one: the addEdge methods below refuse them themselves, and addEdgeWhereNoneExists lets a caller
        // that has ruled one out skip the search.
        super(edgeFactory, true, true, true, false);
        Utils.validateArg(kmerSize > 0, () -> "kmerSize must be > 0 but got " + kmerSize);
        this.kmerSize = kmerSize;
    }

    /**
     * How big of a kmer did we use to create this graph?
     * @return
     */
    public final int getKmerSize() {
        return kmerSize;
    }

    /**
     * Stores the graph in {@link AssemblyGraphSpecifics}: plain directed adjacency, as in jgrapht's
     * {@code DirectedSpecifics}, rather than jgrapht's default fast-lookup variant, which keeps an extra map from every
     * (source, target) pair to its edges and allocates a pair per lookup. Assembly graphs have few edges per vertex, so
     * scanning a vertex's outgoing edges is cheaper. Vertices and edges are kept in insertion order, and a vertex points
     * at its own edges, so adjacency queries need no map lookup. Assembly graphs are always directed.
     *
     * jgrapht's constructor and clone call this before this class's field initializers would run, so
     * {@link #assemblySpecifics} has none.
     */
    @Override
    protected Specifics<V, E> createSpecifics(final boolean directed) {
        assemblySpecifics = new AssemblyGraphSpecifics<>(this);
        return assemblySpecifics;
    }

    @Override
    public boolean containsVertex(final V v) {
        return assemblySpecifics.containsVertex(v);
    }

    /**
     * The edges leaving {@code v}, in the order {@link #outgoingEdgesOf} iterates them, as a read-only list for hot
     * loops that index it rather than allocate an iterator. The list is a view of the graph's storage, so the graph's
     * edges must not be changed while it is in use.
     *
     * @throws IllegalArgumentException if {@code v} is not in the graph
     */
    protected final List<E> outgoingEdgeList(final V v) {
        return assemblySpecifics.edges(v).outgoingList();
    }

    /**
     * The edges entering {@code v}, in the order {@link #incomingEdgesOf} iterates them, as a read-only list for hot
     * loops that index it rather than allocate an iterator. The list is a view of the graph's storage, so the graph's
     * edges must not be changed while it is in use.
     *
     * @throws IllegalArgumentException if {@code v} is not in the graph
     */
    protected final List<E> incomingEdgeList(final V v) {
        return assemblySpecifics.edges(v).incomingList();
    }

    @Override
    public Set<E> outgoingEdgesOf(final V v) {
        return assemblySpecifics.outgoingEdgesOf(v);
    }

    @Override
    public Set<E> incomingEdgesOf(final V v) {
        return assemblySpecifics.incomingEdgesOf(v);
    }

    @Override
    public int outDegreeOf(final V v) {
        return assemblySpecifics.outDegreeOf(v);
    }

    @Override
    public int inDegreeOf(final V v) {
        return assemblySpecifics.inDegreeOf(v);
    }

    /**
     * Adds an edge from {@code source} to {@code target}, unless the graph already has an edge between them.
     *
     * An edge stores its endpoints (see {@link BaseEdge}), so an edge object that already joins other vertices, in
     * this graph or another, is rejected rather than silently re-pointed wherever it is held.
     *
     * @return true if the edge was added, false if the graph already had an edge from {@code source} to {@code target}
     * @throws NullPointerException if {@code e} is null
     * @throws IllegalArgumentException if {@code e} already joins vertices other than {@code source} and {@code target},
     *                                  in this graph or another
     */
    @Override
    public boolean addEdge(final V source, final V target, final E e) {
        Objects.requireNonNull(e);
        rejectEdgeJoiningOtherVertices(source, target, e);
        return !containsEdge(source, target) && super.addEdge(source, target, e);
    }

    /**
     * Adds an edge from {@code source} to {@code target} when the caller has already established that the graph has no
     * edge between them, skipping the search for one that {@link #addEdge(BaseVertex, BaseVertex, BaseEdge)} makes.
     *
     * @return true if the edge was added, false if the graph already contains {@code e}
     * @throws NullPointerException if {@code e} is null
     * @throws IllegalArgumentException if {@code e} already joins vertices other than {@code source} and {@code target},
     *                                  in this graph or another
     */
    protected final boolean addEdgeWhereNoneExists(final V source, final V target, final E e) {
        Objects.requireNonNull(e);
        rejectEdgeJoiningOtherVertices(source, target, e);
        // The precondition is checked only where assertions are enabled, as in tests.
        assert !containsEdge(source, target) : "an edge from " + source + " to " + target + " already exists";
        return super.addEdge(source, target, e);
    }

    /**
     * Adds a new edge from {@code source} to {@code target}, made by the graph's edge factory, unless the graph already
     * has an edge between them.
     *
     * @return the new edge, or null if the graph already had an edge from {@code source} to {@code target}
     */
    @Override
    public E addEdge(final V source, final V target) {
        return containsEdge(source, target) ? null : super.addEdge(source, target);
    }

    private static void rejectEdgeJoiningOtherVertices(final BaseVertex source, final BaseVertex target, final BaseEdge e) {
        if (e.joinsOtherVertices(source, target)) {
            throw new IllegalArgumentException("edge " + e + " already joins " + e.describeJoinedVertices()
                    + "; add its duplicate() to join " + source + " -> " + target);
        }
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference node (meaning that it appears on the reference path in the graph)
     */
    public final boolean isReferenceNode( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");

        if (hasRefEdge(incomingEdgeList(v)) || hasRefEdge(outgoingEdgeList(v))) {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref node, otherwise it's not
        return vertexSet().size() == 1;
    }

    /** @return whether any of {@code edges} is a reference edge */
    private boolean hasRefEdge(final List<E> edges) {
        for (int i = 0; i < edges.size(); i++) {
            if (edges.get(i).isRef()) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a source node (in degree == 0)
     */
    public final boolean isSource( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");
        return inDegreeOf(v) == 0;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a sink node (out degree == 0)
     */
    public final boolean isSink( final V v ) {
        Utils.nonNull(v, "Attempting to test a null vertex.");
        return outDegreeOf(v) == 0;
    }

    /**
     * Get the set of source vertices of this graph
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @return a non-null set
     */
    public final LinkedHashSet<V> getSources() {
        return vertexSet().stream().filter(v -> isSource(v)).collect(Collectors.toCollection(LinkedHashSet::new));
    }

    /**
     * Get the set of sink vertices of this graph
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @return a non-null set
     */
    public final LinkedHashSet<V> getSinks() {
        return vertexSet().stream().filter(v -> isSink(v)).collect(Collectors.toCollection(LinkedHashSet::new));
    }

    /**
     * Convert this kmer graph to a simple sequence graph.
     *
     * Each kmer suffix shows up as a distinct SeqVertex, attached in the same structure as in the kmer
     * graph.  Nodes that are sources are mapped to SeqVertex nodes that contain all of their sequence
     *
     * @return a newly allocated SequenceGraph
     */
    public SeqGraph toSequenceGraph() {
        final SeqGraph seqGraph = new SeqGraph(kmerSize);
        final Map<V, SeqVertex> vertexMap = new HashMap<>();

        // create all of the equivalent seq graph vertices
        for ( final V dv : vertexSet() ) {
            final SeqVertex sv = new SeqVertex(dv.getAdditionalSequence(isSource(dv)));
            sv.setAdditionalInfo(dv.getAdditionalInfo());
            vertexMap.put(dv, sv);
            seqGraph.addVertex(sv);
        }

        // walk through the nodes and connect them to their equivalent seq vertices. This graph has no parallel
        // edges and every sequence vertex was created above for exactly one vertex here, so no copied edge can
        // duplicate another: the copies are added without a duplicate-edge search.
        for( final E e : edgeSet() ) {
            final SeqVertex seqInV = vertexMap.get(getEdgeSource(e));
            final SeqVertex seqOutV = vertexMap.get(getEdgeTarget(e));
            seqGraph.addEdgeWhereNoneExists(seqInV, seqOutV, e.copy());
        }

        return seqGraph;
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @return  non-null byte array
     */
    public final byte[] getAdditionalSequence( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");
        return v.getAdditionalSequence(isSource(v));
    }

    /**
     * Pull out the additional sequence implied by traversing this node in the graph
     * @param v the vertex from which to pull out the additional base sequence
     * @param isSource if true, treat v as a source vertex regardless of in-degree
     * @return  non-null byte array
     */
    public static final byte[] getAdditionalSequence( final BaseVertex v, final boolean isSource) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");
        return v.getAdditionalSequence(isSource);
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference source
     */
    public final boolean isRefSource( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");

        // confirm that no incoming edges are reference edges
        if (hasRefEdge(incomingEdgeList(v))) {
            return false;
        }

        // confirm that there is an outgoing reference edge
        if (hasRefEdge(outgoingEdgeList(v))) {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref source, otherwise it's not
        return vertexSet().size() == 1;
    }

    /**
     * @param v the vertex to test
     * @return  true if this vertex is a reference sink
     */
    public final boolean isRefSink( final V v ) {
        Utils.nonNull(v, "Attempting to pull sequence from a null vertex.");

        // confirm that no outgoing edges are reference edges
        if (hasRefEdge(outgoingEdgeList(v))) {
            return false;
        }

        // confirm that there is an incoming reference edge
        if (hasRefEdge(incomingEdgeList(v))) {
            return true;
        }

        // edge case: if the graph only has one node then it's a ref sink, otherwise it's not
        return vertexSet().size() == 1;
    }

    /**
     * @return the reference source vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public V getReferenceSourceVertex( ) {
        for (final V v : vertexSet()) {
            if (isRefSource(v)) {
                return v;
            }
        }
        return null;
    }

    /**
     * @return the reference sink vertex pulled from the graph, can be null if it doesn't exist in the graph
     */
    public V getReferenceSinkVertex( ) {
        for (final V v : vertexSet()) {
            if (isRefSink(v)) {
                return v;
            }
        }
        return null;
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the next reference vertex if it exists, otherwise null
     */
    public final V getNextReferenceVertex( final V v ) {
        return getNextReferenceVertex(v, false, Optional.<E>empty());
    }

    /**
     * Traverse the graph and get the next reference vertex if it exists
     * @param v the current vertex, can be null
     * @param allowNonRefPaths if true, allow sub-paths that are non-reference if there is only a single outgoing edge
     * @param blacklistedEdge optional edge to ignore in the traversal down; useful to exclude the non-reference dangling paths
     * @return the next vertex (but not necessarily on the reference path if allowNonRefPaths is true) if it exists, otherwise null
     */
    public final V getNextReferenceVertex( final V v, final boolean allowNonRefPaths, final Optional<E> blacklistedEdge ) {
        if( v == null ) { return null; }

        final Set<E> outgoingEdges = outgoingEdgesOf(v);

        if (outgoingEdges.isEmpty()){
            return null;
        }

        for( final E edgeToTest : outgoingEdges ) {
            if( edgeToTest.isRef() ) {
                return getEdgeTarget(edgeToTest);
            }
        }

        if (!allowNonRefPaths){
            return null;
        }

        //singleton or empty set
        final Set<E> blacklistedEdgeSet = blacklistedEdge.isPresent() ? Collections.singleton(blacklistedEdge.get()) : Collections.emptySet();

        // if we got here, then we aren't on a reference path
        final List<E> edges = outgoingEdges.stream().filter(e -> !blacklistedEdgeSet.contains(e)).limit(2).collect(Collectors.toList());
        return edges.size() == 1 ? getEdgeTarget(edges.get(0)) : null;
    }

    /**
     * Traverse the graph and get the previous reference vertex if it exists
     * @param v the current vertex, can be null
     * @return  the previous reference vertex if it exists or null otherwise.
     */
    public final V getPrevReferenceVertex( final V v ) {
        if( v == null ) { return null; }
        return incomingEdgesOf(v).stream().map(e -> getEdgeSource(e)).filter(vrtx -> isReferenceNode(vrtx)).findFirst().orElse(null);
    }

    /**
     * Walk along the reference path in the graph and pull out the corresponding bases
     * @param fromVertex    starting vertex
     * @param toVertex      ending vertex
     * @param includeStart  should the starting vertex be included in the path
     * @param includeStop   should the ending vertex be included in the path
     * @return              byte[] array holding the reference bases, this can be null if there are no nodes between the starting and ending vertex (insertions for example)
     */
    public byte[] getReferenceBytes( final V fromVertex, final V toVertex, final boolean includeStart, final boolean includeStop ) {
        Utils.nonNull(fromVertex, "Starting vertex in requested path cannot be null.");
        Utils.nonNull(toVertex, "From vertex in requested path cannot be null.");

        byte[] bytes = null;
        V v = fromVertex;
        if( includeStart ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        v = getNextReferenceVertex(v); // advance along the reference path
        while( v != null && !v.equals(toVertex) ) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
            v = getNextReferenceVertex(v); // advance along the reference path
        }
        if( includeStop && v != null && v.equals(toVertex)) {
            bytes = ArrayUtils.addAll(bytes, getAdditionalSequence(v));
        }
        return bytes;
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    @SafeVarargs
    @SuppressWarnings("varargs")
    public final void addVertices(final V... vertices) {
        Utils.nonNull(vertices);
        addVertices(Arrays.asList(vertices));
    }

    /**
     * Convenience function to add multiple vertices to the graph at once
     * @param vertices one or more vertices to add
     */
    public final void addVertices(final Collection<V> vertices) {
        Utils.nonNull(vertices);
        vertices.forEach(v -> addVertex(v));
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    @SafeVarargs
    public final void addEdges(final V start, final V... remaining) {
        Utils.nonNull(start, "start vertex");
        if (remaining == null || remaining.length == 0){
            return;
        }
        V prev = start;
        for ( final V next : remaining ) {
            Utils.nonNull(next, "null vertex");
            addEdge(prev, next);
            prev = next;
        }
    }

    /**
     * Convenience function to add multiple edges to the graph
     * @param start the first vertex to connect
     * @param remaining all additional vertices to connect
     */
    @SafeVarargs
    public final void addEdges(final Supplier<E> template, final V start, final V... remaining) {
        Utils.nonNull(template, "template edge");
        Utils.nonNull(start, "start vertex");

        V prev = start;
        for ( final V next : remaining ) {
            Utils.nonNull(next, "null vertex");
            addEdge(prev, next, template.get());
            prev = next;
        }
    }

    /**
     * Get the set of vertices connected by outgoing edges of V
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @param v a non-null vertex
     * @return a set of vertices connected by outgoing edges from v
     */
    public final Set<V> outgoingVerticesOf(final V v) {
        Utils.nonNull(v);
        final Set<V> targets = new LinkedHashSet<>();
        for (final E e : outgoingEdgesOf(v)) {
            targets.add(getEdgeTarget(e));
        }
        return targets;
    }

    /**
     * Get the set of vertices connected to v by incoming edges
     * NOTE: We return a LinkedHashSet here in order to preserve the determinism in the output order of VertexSet(),
     *       which is deterministic in output due to the underlying sets all being LinkedHashSets.
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v
     */
    public final Set<V> incomingVerticesOf(final V v) {
        Utils.nonNull(v);
        final Set<V> sources = new LinkedHashSet<>();
        for (final E e : incomingEdgesOf(v)) {
            sources.add(getEdgeSource(e));
        }
        return sources;
    }

    /**
     * Get the set of vertices connected to v by incoming or outgoing edges
     * @param v a non-null vertex
     * @return a set of vertices {X} connected X -> v or v -> Y
     */
    public final Set<V> neighboringVerticesOf(final V v) {
        Utils.nonNull(v);
        return Sets.union(incomingVerticesOf(v), outgoingVerticesOf(v));
    }

    /**
     * Print out the graph in the dot language for visualization
     * @param destination File to write to
     */
    public final void printGraph(final File destination, final int pruneFactor) {
        try (PrintStream stream = new PrintStream(new FileOutputStream(destination))) {
            printGraph(stream, true, pruneFactor);
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(destination.getAbsolutePath(), e);
        }
    }

    public final void printGraph(final PrintStream graphWriter, final boolean writeHeader, final int pruneFactor) {
        if ( writeHeader ) {
            graphWriter.println("digraph assemblyGraphs {");
        }

        for( final E edge : edgeSet() ) {
            final String edgeString =  String.format("\t%s -> %s ", getEdgeSource(edge).toString(), getEdgeTarget(edge).toString());
            final String edgeLabelString;
            if (edge.getMultiplicity() > 0 && edge.getMultiplicity() < pruneFactor){
                edgeLabelString = String.format("[style=dotted,color=grey,label=\"%s\"];", edge.getDotLabel());
            } else {
                edgeLabelString = String.format("[label=\"%s\"];", edge.getDotLabel());
            }
            graphWriter.print(edgeString);
            graphWriter.print(edgeLabelString);
            if( edge.isRef() ) {
                graphWriter.println(edgeString + " [color=red];");
            }
        }

        for( final V v : vertexSet() ) {
            graphWriter.println(String.format("\t%s [label=\"%s\",shape=box]", v.toString(),
                    new String(getAdditionalSequence(v)) + " (" + v.hashCode() + ") " + v.getAdditionalInfo()) );
        }

        getExtraGraphFileLines().forEach(graphWriter::println);

        if ( writeHeader ) {
            graphWriter.println("}");
        }
    }

    // Extendable method intended to allow for adding extra material to the graph
    public List<String> getExtraGraphFileLines() {
        return Collections.emptyList();
    }

    /**
     * Remove edges that are connected before the reference source and after the reference sink
     *
     * Also removes all vertices that are orphaned by this process
     */
    public final void cleanNonRefPaths() {
        // Only non-reference edges are removed below, which cannot change which vertices are the reference source and sink.
        final V refSource = getReferenceSourceVertex();
        final V refSink = getReferenceSinkVertex();
        if( refSource == null || refSink == null ) {
            return;
        }

        // Remove non-ref edges connected before and after the reference path
        final Collection<E> edgesToCheck = new HashSet<>();
        edgesToCheck.addAll(incomingEdgesOf(refSource));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( incomingEdgesOf(getEdgeSource(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        edgesToCheck.addAll(outgoingEdgesOf(refSink));
        while( !edgesToCheck.isEmpty() ) {
            final E e = edgesToCheck.iterator().next();
            if( !e.isRef() ) {
                edgesToCheck.addAll( outgoingEdgesOf(getEdgeTarget(e)) );
                removeEdge(e);
            }
            edgesToCheck.remove(e);
        }

        removeSingletonOrphanVertices();
    }

    /**
     * Remove all vertices in the graph that have in and out degree of 0
     */
    public void removeSingletonOrphanVertices() {
        // Run through the graph and clean up singular orphaned nodes
        //Note: need to collect nodes to remove first because we can't directly modify the list we're iterating over
        final List<V> toRemove = vertexSet().stream().filter(v -> isSingletonOrphan(v)).collect(Collectors.toList());
        removeAllVertices(toRemove);
    }

    private boolean isSingletonOrphan(final V v) {
        Utils.nonNull(v);
        return inDegreeOf(v) == 0 && outDegreeOf(v) == 0 && !isRefSource(v);
    }

    /**
     * Remove all vertices on the graph that cannot be accessed by following any edge,
     * regardless of its direction, from the reference source vertex
     */
    public final void removeVerticesNotConnectedToRefRegardlessOfEdgeDirection() {
        final V refV = getReferenceSourceVertex();
        final Set<V> connected = refV == null ? Collections.emptySet() : verticesReachableFrom(refV, true, true);
        removeAllVertices(verticesNotIn(connected));
    }

    /**
     * Remove all vertices in the graph that aren't on a path from the reference source vertex to the reference sink vertex
     *
     * More aggressive reference pruning algorithm than removeVerticesNotConnectedToRefRegardlessOfEdgeDirection,
     * as it requires vertices to not only be connected by a series of directed edges but also prunes away
     * paths that do not also meet eventually with the reference sink vertex
     */
    public final void removePathsNotConnectedToRef() {
        final V refSource = getReferenceSourceVertex();
        final V refSink = getReferenceSinkVertex();
        if ( refSource == null || refSink == null ) {
            throw new IllegalStateException("Graph must have ref source and sink vertices");
        }

        // keep only the vertices reachable both forward from the ref source and backward from the ref sink
        final Set<V> onPathFromRefSource = verticesReachableFrom(refSource, false, true);
        onPathFromRefSource.retainAll(verticesReachableFrom(refSink, true, false));
        removeAllVertices(verticesNotIn(onPathFromRefSource));

        // simple sanity checks that this algorithm is working.
        if ( getSinks().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference sink, but found " + getSinks());
        }

        if ( getSources().size() > 1 ) {
            throw new IllegalStateException("Should have eliminated all but the reference source, but found " + getSources());
        }
    }

    /**
     * Semi-lenient comparison of two graphs, truing true if g1 and g2 have similar structure
     *
     * By similar this means that both graphs have the same number of vertices, where each vertex can find
     * a vertex in the other graph that's seqEqual to it.  A similar constraint applies to the edges,
     * where all edges in g1 must have a corresponding edge in g2 where both source and target vertices are
     * seqEqual
     *
     * @param g1 the first graph to compare
     * @param g2 the second graph to compare
     * @param <T> the type of the nodes in those graphs
     * @return true if g1 and g2 are equals
     */
    public static <T extends BaseVertex, E extends BaseEdge> boolean graphEquals(final BaseGraph<T,E> g1, final BaseGraph<T,E> g2) {
        Utils.nonNull(g1, "g1");
        Utils.nonNull(g2, "g2");
        final Set<T> vertices1 = g1.vertexSet();
        final Set<T> vertices2 = g2.vertexSet();
        final Set<E> edges1 = g1.edgeSet();
        final Set<E> edges2 = g2.edgeSet();

        if ( vertices1.size() != vertices2.size() || edges1.size() != edges2.size() ) {
            return false;
        }

        //for every vertex in g1 there is a vertex in g2 with an equal getSequenceString
        final boolean ok= vertices1.stream().map(v1 -> v1.getSequenceString()).allMatch(v1seqString -> vertices2.stream().anyMatch(v2 -> v1seqString.equals(v2.getSequenceString())));
        if (! ok){
            return false;
        }

        //for every edge in g1 there is an equal edge in g2
        final boolean okG1 = edges1.stream().allMatch(e1 -> edges2.stream().anyMatch(e2 -> g1.seqEquals(e1, e2, g2)));
        if (! okG1){
            return false;
        }
        //for every edge in g2 there is an equal edge in g1
        return edges2.stream().allMatch(e2 -> edges1.stream().anyMatch(e1 -> g2.seqEquals(e2, e1, g1)));
    }

    // For use when comparing edges across graphs!
    private boolean seqEquals( final E edge1, final E edge2, final BaseGraph<V,E> graph2 ) {
        return (getEdgeSource(edge1).seqEquals(graph2.getEdgeSource(edge2))) && (getEdgeTarget(edge1).seqEquals(graph2.getEdgeTarget(edge2)));
    }


    /**
     * Get the incoming edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single incoming edge to v, or null if none exists
     */
    public final E incomingEdgeOf(final V v) {
        Utils.nonNull(v);
        return getSingletonEdge(incomingEdgesOf(v));
    }

    /**
     * Get the outgoing edge of v.  Requires that there be only one such edge or throws an error
     * @param v our vertex
     * @return the single outgoing edge from v, or null if none exists
     */
    public final E outgoingEdgeOf(final V v) {
        Utils.nonNull(v);
        return getSingletonEdge(outgoingEdgesOf(v));
    }

    /**
     * Helper function that gets the a single edge from edges, null if edges is empty, or
     * throws an error is edges has more than 1 element
     * @param edges a set of edges
     * @return a edge
     */
    private E getSingletonEdge(final Collection<E> edges) {
        Utils.validateArg(edges.size() <= 1, () -> "Cannot get a single incoming edge for a vertex with multiple incoming edges " + edges);
        return edges.isEmpty() ? null : edges.iterator().next();
    }

    /**
     * Add edge between source -> target if none exists, or add e to an already existing one if present
     *
     * @param source source vertex
     * @param target vertex
     * @param e edge to add
     */
    public final void addOrUpdateEdge(final V source, final V target, final E e) {
        Utils.nonNull(source, "source");
        Utils.nonNull(target, "target");
        Utils.nonNull(e, "edge");

        final E prev = getEdge(source, target);
        if ( prev != null ) {
            prev.add(e);
        } else {
            addEdge(source, target, e);
        }
    }

    @Override
    public String toString() {
        return "BaseGraph{" +
                "kmerSize=" + kmerSize +
                '}';
    }

    /**
     * Get the set of vertices within distance edges of source, regardless of edge direction
     *
     * @param source the source vertex to consider
     * @param distance the distance
     * @return a set of vertices within distance of source
     */
    private Set<V> verticesWithinDistance(final V source, final int distance) {
        if ( distance == 0 ) {
            return Collections.singleton(source);
        }

        final Set<V> found = new HashSet<>();
        found.add(source);
        for ( final V v : neighboringVerticesOf(source) ) {
            found.addAll(verticesWithinDistance(v, distance - 1));
        }

        return found;
    }

    /**
     * Get a graph containing only the vertices within distance edges of target
     * @param target a vertex in graph
     * @param distance the max distance
     * @return a non-null graph
     */
    public final BaseGraph<V,E> subsetToNeighbors(final V target, final int distance) {
        Utils.nonNull(target, "Target cannot be null");
        Utils.validateArg(containsVertex(target), () -> "Graph doesn't contain vertex " + target);
        Utils.validateArg(distance >= 0, () -> "Distance must be >= 0 but got " + distance);

        final Set<V> toKeep = verticesWithinDistance(target, distance);
        final Collection<V> toRemove = new HashSet<>(vertexSet());
        toRemove.removeAll(toKeep);

        final BaseGraph<V,E> result = clone();
        result.removeAllVertices(toRemove);

        return result;
    }

    /**
     * Get a subgraph of graph that contains only vertices within a given number of edges of the ref source vertex
     * @return a non-null subgraph of this graph
     */
    public final BaseGraph<V,E> subsetToRefSource(final int refSourceNeighborhood) {
        Utils.validateArg(refSourceNeighborhood > 0, () -> "refSourceNeighborhood needs to be positive but was " + refSourceNeighborhood);
        return subsetToNeighbors(getReferenceSourceVertex(), refSourceNeighborhood);
    }

    /**
     * Checks whether the graph contains all the vertices in a collection.
     *
     * @param vertices the vertices to check. Must not be null and must not contain a null.
     *
     * @throws IllegalArgumentException if {@code vertices} is {@code null}.
     *
     * @return {@code true} if all the vertices in the input collection are present in this graph.
     * Also if the input collection is empty. Otherwise it returns {@code false}.
     */
    public final boolean containsAllVertices(final Collection<? extends V> vertices) {
        Utils.nonNull(vertices, "the input vertices collection cannot be null");
        Utils.containsNoNull(vertices, "null vertex");
        return vertices.stream().allMatch(v -> containsVertex(v));
    }

    /**
     * Checks for the presence of directed cycles in the graph.
     *
     * @return {@code true} if the graph has cycles, {@code false} otherwise.
     */
    public final boolean hasCycles() {
        return new CycleDetector<>(this).detectCycles();
    }

    @Override
    @SuppressWarnings("unchecked")
    public BaseGraph<V,E> clone()  {
        return (BaseGraph<V,E>) super.clone();
    }

    /**
     * Breadth-first search from one vertex along edges in the chosen directions.
     *
     * @param start the vertex to start from, which must be in this graph
     * @param followIncomingEdges whether to follow edges backward, from target to source
     * @param followOutgoingEdges whether to follow edges forward, from source to target
     * @return every vertex reachable from start by following edges in the given directions, start included
     */
    private Set<V> verticesReachableFrom(final V start, final boolean followIncomingEdges, final boolean followOutgoingEdges) {
        final Set<V> reached = new HashSet<>();
        final Deque<V> toVisit = new ArrayDeque<>();
        reached.add(start);
        toVisit.add(start);
        while ( ! toVisit.isEmpty() ) {
            final V v = toVisit.poll();
            if ( followIncomingEdges ) {
                for ( final E e : incomingEdgesOf(v) ) {
                    final V source = getEdgeSource(e);
                    if ( reached.add(source) ) {
                        toVisit.add(source);
                    }
                }
            }
            if ( followOutgoingEdges ) {
                for ( final E e : outgoingEdgesOf(v) ) {
                    final V target = getEdgeTarget(e);
                    if ( reached.add(target) ) {
                        toVisit.add(target);
                    }
                }
            }
        }
        return reached;
    }

    /** @return the vertices of this graph that are not in keep, in the graph's vertex order */
    private List<V> verticesNotIn(final Set<V> keep) {
        final List<V> others = new ArrayList<>();
        for ( final V v : vertexSet() ) {
            if ( ! keep.contains(v) ) {
                others.add(v);
            }
        }
        return others;
    }
}
