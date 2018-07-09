/*
Name: Rebecca Hu
PID: A13999923
 */

import java.util.*;

/**
 * Implementation of Vertex class
 * @author Rebecca Hu
 * @since 6/4/2018
 */
public class Vertex {

    private static final double INFINITY = Double.POSITIVE_INFINITY;

    public String name; // the name of this vertex
    public int x; // the x coordinates of this vertex on map
    public int y; // the y coordinates of this vertex on map
    public ArrayList<Edge> sentEdges; //list of edges that are directed from this vertex
    public ArrayList<Edge> receivedEdges; //list of edges that are directed to this vertex
    public Vertex prev; // list of shortest path to another vertex
    public double cost; //initialize to infinity


    public Vertex(String name, int x, int y) {
        this.name = name;
        this.x = x;
        this.y = y;
        this.sentEdges = new ArrayList<>();
        this.receivedEdges = new ArrayList<>();
        this.prev = null;
        this.cost = INFINITY;
    }


    /**
     * adds an edge to the vertex's sentEdges ArrayList
     *
     * @param e, edge to be added
     */
    public void sendEdge(Edge e){
        this.sentEdges.add(e);
    }

    /**
     * adds an edge to the vertex's recievedEdges ArrayList
     *
     * @param e, edge to be received
     */
    public void receiveEdge (Edge e){
        this.receivedEdges.add(e);
    }

    @Override
    public int hashCode() {
        // we assume that each vertex has a unique name
        return name.hashCode();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null) {
            return false;
        }
        if (!(o instanceof Vertex)) {
            return false;
        }
        Vertex oVertex = (Vertex) o;

        return name.equals(oVertex.name) && x == oVertex.x && y == oVertex.y;
    }

    public String toString() {
        return name + " (" + x + ", " + y + ")";
    }

}