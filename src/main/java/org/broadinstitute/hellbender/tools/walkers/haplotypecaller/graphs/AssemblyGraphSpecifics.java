package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs;

import org.jgrapht.graph.AbstractBaseGraph;
import org.jgrapht.graph.specifics.Specifics;
import org.jgrapht.util.ArrayUnenforcedSet;

import java.io.Serializable;
import java.util.AbstractSet;
import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Directed adjacency storage for assembly graphs, equivalent to jgrapht's plain {@code DirectedSpecifics} except in
 * how a vertex's edges are found.
 *
 * As there, vertices are kept in a {@link LinkedHashMap} in insertion order, each with its incoming and outgoing edges
 * in {@link ArrayUnenforcedSet}s in insertion order, so iteration orders are unchanged. The difference is that a
 * vertex also points at its entry, through {@link BaseVertex#graphEdges}, so the adjacency queries that dominate
 * assembly find a vertex's edges without hashing it and probing the map.
 *
 * A vertex can be in several graphs at once (a cloned graph shares its vertices), but has one pointer, which holds the
 * entry of the graph that most recently added or looked the vertex up through its map. An entry records the graph
 * that owns it and whether its vertex is still in that graph, so the pointer is used only when it belongs to this graph
 * and is current; otherwise the vertex is looked up in the map as before, which also handles a vertex object that is
 * equal to, but not the same as, the one the graph holds.
 *
 * Unlike {@code DirectedSpecifics}, a query for a vertex that is not in the graph throws rather than adding the vertex:
 * {@link NullPointerException} for null and {@link IllegalArgumentException} otherwise, the same exceptions jgrapht's
 * membership assertions throw.
 *
 * Because lookups re-point a vertex at the entry they find, even read-only queries such as {@code outDegreeOf} write
 * to the vertex. Graphs that share vertices (a graph and its clones) must therefore be used by one thread at a time,
 * as assembly uses them. A vertex's pointer also keeps the graph that last touched it reachable for as long as the
 * vertex is, so vertices should not outlive the graphs built from them by much; assembly drops both together.
 */
final class AssemblyGraphSpecifics<V extends BaseVertex, E> implements Specifics<V, E>, Serializable {
    private static final long serialVersionUID = 1L;

    private final AbstractBaseGraph<V, E> graph;
    private final Map<V, VertexEdges<V, E>> vertexMap = new LinkedHashMap<>();
    private final Set<V> vertexSet = new VertexSet();

    /** One vertex's incoming and outgoing edges in one graph. */
    static final class VertexEdges<V, E> implements Serializable {
        private static final long serialVersionUID = 1L;

        private final AssemblyGraphSpecifics<?, ?> owner;
        private final V vertex;
        /** False once the vertex has been removed from {@link #owner}. */
        private boolean inGraph = true;
        /** Scratch count of incoming edges not yet removed, used only within {@link #hasCycles}. */
        private int pendingIncoming;
        final ArrayUnenforcedSet<E> incoming = new ArrayUnenforcedSet<>(1);
        final ArrayUnenforcedSet<E> outgoing = new ArrayUnenforcedSet<>(1);
        private transient Set<E> unmodifiableIncoming = null;
        private transient Set<E> unmodifiableOutgoing = null;
        private transient List<E> incomingList = null;
        private transient List<E> outgoingList = null;

        private VertexEdges(final AssemblyGraphSpecifics<?, ?> owner, final V vertex) {
            this.owner = owner;
            this.vertex = vertex;
        }

        /** The incoming edges as a read-only set, made on first use and then reused. */
        Set<E> unmodifiableIncoming() {
            if (unmodifiableIncoming == null) {
                unmodifiableIncoming = Collections.unmodifiableSet(incoming);
            }
            return unmodifiableIncoming;
        }

        /** The outgoing edges as a read-only set, made on first use and then reused. */
        Set<E> unmodifiableOutgoing() {
            if (unmodifiableOutgoing == null) {
                unmodifiableOutgoing = Collections.unmodifiableSet(outgoing);
            }
            return unmodifiableOutgoing;
        }

        /** The incoming edges as a read-only list, for loops that index them. */
        List<E> incomingList() {
            if (incomingList == null) {
                incomingList = Collections.unmodifiableList(incoming);
            }
            return incomingList;
        }

        /** The outgoing edges as a read-only list, for loops that index them. */
        List<E> outgoingList() {
            if (outgoingList == null) {
                outgoingList = Collections.unmodifiableList(outgoing);
            }
            return outgoingList;
        }
    }

    AssemblyGraphSpecifics(final AbstractBaseGraph<V, E> graph) {
        this.graph = graph;
    }

    /**
     * @return the entry of a vertex in this graph
     * @throws NullPointerException if {@code vertex} is null
     * @throws IllegalArgumentException if {@code vertex} is not in this graph
     */
    @SuppressWarnings("unchecked")
    VertexEdges<V, E> edges(final V vertex) {
        final VertexEdges<?, ?> pointed = vertex.graphEdges;
        if (pointed != null && pointed.owner == this && pointed.inGraph) {
            return (VertexEdges<V, E>) pointed;
        }
        final VertexEdges<V, E> found = vertexMap.get(vertex);
        if (found == null) {
            throw new IllegalArgumentException("no such vertex in graph: " + vertex);
        }
        if (found.vertex == vertex) {
            vertex.graphEdges = found;
        }
        return found;
    }

    /** Whether the graph holds a vertex equal to {@code o}. */
    boolean containsVertex(final Object o) {
        if (o instanceof BaseVertex) {
            final VertexEdges<?, ?> pointed = ((BaseVertex) o).graphEdges;
            if (pointed != null && pointed.owner == this && pointed.inGraph) {
                return true;
            }
        }
        return vertexMap.containsKey(o);
    }

    /**
     * Whether the graph has a directed cycle, a self-loop included. Kahn's algorithm: repeatedly remove vertices with
     * no remaining incoming edges; the graph is acyclic exactly when that removes every vertex.
     */
    boolean hasCycles() {
        final ArrayDeque<VertexEdges<V, E>> ready = new ArrayDeque<>();
        for (final VertexEdges<V, E> entry : vertexMap.values()) {
            entry.pendingIncoming = entry.incoming.size();
            if (entry.pendingIncoming == 0) {
                ready.add(entry);
            }
        }
        int removed = 0;
        while (!ready.isEmpty()) {
            final VertexEdges<V, E> entry = ready.poll();
            removed++;
            final ArrayUnenforcedSet<E> outgoing = entry.outgoing;
            for (int i = 0; i < outgoing.size(); i++) {
                final VertexEdges<V, E> target = edges(graph.getEdgeTarget(outgoing.get(i)));
                if (--target.pendingIncoming == 0) {
                    ready.add(target);
                }
            }
        }
        return removed != vertexMap.size();
    }

    @Override
    public void addVertex(final V vertex) {
        // The graph adds a vertex only when it holds no equal one, so this always makes a new entry.
        final VertexEdges<V, E> entry = new VertexEdges<>(this, vertex);
        vertexMap.put(vertex, entry);
        vertex.graphEdges = entry;
    }

    /**
     * Removes the vertex equal to {@code o}, marks its entry as no longer in this graph, and clears the vertex's pointer
     * if it points at that entry.
     *
     * @return whether the graph held such a vertex
     */
    private boolean removeVertex(final Object o) {
        final VertexEdges<V, E> removed = vertexMap.remove(o);
        if (removed == null) {
            return false;
        }
        removed.inGraph = false;
        if (removed.vertex.graphEdges == removed) {
            removed.vertex.graphEdges = null;
        }
        return true;
    }

    /**
     * The vertex set the graph reads and removes vertices through. It iterates in insertion order and cannot be
     * modified except by {@link #remove}, so every removal marks the vertex's entry as no longer in the graph.
     */
    private final class VertexSet extends AbstractSet<V> implements Serializable {
        private static final long serialVersionUID = 1L;

        @Override
        public Iterator<V> iterator() {
            final Iterator<V> keys = vertexMap.keySet().iterator();
            return new Iterator<V>() {
                @Override
                public boolean hasNext() {
                    return keys.hasNext();
                }

                @Override
                public V next() {
                    return keys.next();
                }
            };
        }

        @Override
        public int size() {
            return vertexMap.size();
        }

        @Override
        public boolean contains(final Object o) {
            return containsVertex(o);
        }

        @Override
        public boolean remove(final Object o) {
            return removeVertex(o);
        }
    }

    @Override
    public Set<V> getVertexSet() {
        return vertexSet;
    }

    @Override
    public Set<E> getAllEdges(final V sourceVertex, final V targetVertex) {
        Set<E> edges = null;
        if (containsVertex(sourceVertex) && containsVertex(targetVertex)) {
            edges = new ArrayUnenforcedSet<>();
            final ArrayUnenforcedSet<E> outgoing = edges(sourceVertex).outgoing;
            for (int i = 0; i < outgoing.size(); i++) {
                final E e = outgoing.get(i);
                if (graph.getEdgeTarget(e).equals(targetVertex)) {
                    edges.add(e);
                }
            }
        }
        return edges;
    }

    @Override
    public E getEdge(final V sourceVertex, final V targetVertex) {
        if (containsVertex(sourceVertex) && containsVertex(targetVertex)) {
            final ArrayUnenforcedSet<E> outgoing = edges(sourceVertex).outgoing;
            for (int i = 0; i < outgoing.size(); i++) {
                final E e = outgoing.get(i);
                if (graph.getEdgeTarget(e).equals(targetVertex)) {
                    return e;
                }
            }
        }
        return null;
    }

    @Override
    public void addEdgeToTouchingVertices(final E e) {
        edges(graph.getEdgeSource(e)).outgoing.add(e);
        edges(graph.getEdgeTarget(e)).incoming.add(e);
    }

    @Override
    public int degreeOf(final V vertex) {
        return inDegreeOf(vertex) + outDegreeOf(vertex);
    }

    @Override
    public Set<E> edgesOf(final V vertex) {
        final VertexEdges<V, E> entry = edges(vertex);
        if (entry.incoming.isEmpty() && entry.outgoing.isEmpty()) {
            // pruning removes many vertices left without edges, each after asking for its edges
            return Collections.emptySet();
        }
        final ArrayUnenforcedSet<E> inAndOut = new ArrayUnenforcedSet<>(entry.incoming);
        inAndOut.addAll(entry.outgoing);

        // we have two copies for each self-loop - remove one of them.
        if (graph.isAllowingLoops()) {
            final Set<E> loops = getAllEdges(vertex, vertex);
            for (int i = 0; i < inAndOut.size();) {
                final E e = inAndOut.get(i);
                if (loops.contains(e)) {
                    inAndOut.remove(i);
                    loops.remove(e); // so we remove it only once
                } else {
                    i++;
                }
            }
        }
        return Collections.unmodifiableSet(inAndOut);
    }

    @Override
    public int inDegreeOf(final V vertex) {
        return edges(vertex).incoming.size();
    }

    @Override
    public Set<E> incomingEdgesOf(final V vertex) {
        return edges(vertex).unmodifiableIncoming();
    }

    @Override
    public int outDegreeOf(final V vertex) {
        return edges(vertex).outgoing.size();
    }

    @Override
    public Set<E> outgoingEdgesOf(final V vertex) {
        return edges(vertex).unmodifiableOutgoing();
    }

    @Override
    public void removeEdgeFromTouchingVertices(final E e) {
        edges(graph.getEdgeSource(e)).outgoing.remove(e);
        edges(graph.getEdgeTarget(e)).incoming.remove(e);
    }
}
