/*
Name: Rebecca Hu
PID: A13999923
 */

import java.lang.reflect.Array;
import java.util.*;

/**
 * Implementation of Graph class
 * @author Rebecca Hu
 * @since 6/4/2018
 */
public class Graph {

    private static final double INFINITY = Double.POSITIVE_INFINITY;
    private static final double HALF = 0.5;
    private static final int SQUARED = 2;

    private HashMap<Integer, Vertex> vertices;


    /**
     * Constructor for Graph
     */
    public Graph() {
        this.vertices = new HashMap<>();
    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if given vertex
     * already exist in the graph.
     *
     * @param v vertex to be added to the graph
     * @throws IllegalArgumentException if two vertices with the same name are added.
     */
    public void addVertex(Vertex v) throws IllegalArgumentException {
        int hashKey = v.hashCode();
        if (this.vertices.containsKey(hashKey)){
            throw new IllegalArgumentException();
        } else {
            this.vertices.put(hashKey, v);
        }
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    public Collection<Vertex> getVertices() {
        return this.vertices.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        int hashKey = name.hashCode();
        return this.vertices.get(hashKey);
    }

    /**
     * Adds a directed edge from vertex u to vertex v, Throws IllegalArgumentException if one of
     * the vertex does not exist
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight weight of the edge between vertex u and v
     * @throws IllegalArgumentException if one of the vertex does not exist
     */
    public void addEdge(String nameU, String nameV, Double weight) throws IllegalArgumentException {
        int keyU = nameU.hashCode();
        int keyV = nameV.hashCode();
        if (!this.vertices.containsKey(keyU) || !this.vertices.containsKey(keyV)){
            throw new IllegalArgumentException();
        } else {
            Edge newEdge = new Edge(this.getVertex(nameU), this.getVertex(nameV), weight);
            this.getVertex(nameU).sendEdge(newEdge);
            this.getVertex(nameV).receiveEdge(newEdge);
        }
    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight  weight of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double weight) {
        this.addEdge(nameU, nameV, weight);
        this.addEdge(nameV, nameU, weight);
    }

    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
        double xDiff = ux-vx;
        double yDiff = uy-vy;
        //distance formula
        return Math.pow( Math.pow(xDiff, SQUARED) + Math.pow(yDiff, SQUARED), HALF);
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    public void computeAllEuclideanDistances() {
        Iterator iterator = this.getVertices().iterator();

        while (iterator.hasNext()){    // for each vertex in the hash table

            Vertex v = (Vertex) iterator.next();

            for (int i = 0; i < v.sentEdges.size(); i++){ //each edge in vertex's sentEdges list
                Edge e = v.sentEdges.get(i);
                //set new distance
                e.distance = computeEuclideanDistance(
                        e.source.x, e.source.y, e.target.x, e.target.y);
            }
        }
    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {
        Iterator iterator = this.getVertices().iterator();
        while (iterator.hasNext()){
            Vertex currV = (Vertex) iterator.next();
            currV.prev = null;    //clears the last shortest path found
            currV.cost = INFINITY;
        }
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using DFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void DFS(String s, String t) {
        resetAllVertices();

        ArrayList<Vertex> explored = new ArrayList<>();
        Stack<Vertex> frontier = new Stack<>();
        frontier.push(this.getVertex(s));

        while (!frontier.isEmpty()) {
            Vertex curr = frontier.pop();

            if (!explored.contains(curr)) {     //add to explored if has not been explored yet
                explored.add(curr);
            } else {
                continue;       //skip to next vertex if has been explored already
            }

            //add all neighbors to frontier (if not already explored)
            for (int i = 0; i < curr.sentEdges.size(); i++) {
                Vertex neighbor = curr.sentEdges.get(i).target;
                if (!explored.contains(neighbor)) {
                    frontier.push(neighbor);     //if vertex not been explored, add to frontier
                    neighbor.prev = curr;
                }
            }

                if (curr == this.getVertex(t)) {       //if reached target vertex, terminate
                return;
            }
        }
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using BFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void BFS(String s, String t) {
        resetAllVertices();

        ArrayList<Vertex> explored = new ArrayList<>();
        LinkedList<Vertex> frontier = new LinkedList<>();
        frontier.add(this.getVertex(s));

        while (!frontier.isEmpty()) {
            Vertex curr = frontier.remove();

            if (!explored.contains(curr)) {     //add to explored if has not been explored yet
                explored.add(curr);
            } else {
                continue;       //skip to next vertex if has been explored already
            }

            for (int i = 0; i < curr.sentEdges.size(); i++) {
                Vertex neighbor = curr.sentEdges.get(i).target;
                //add all neighbors to frontier (if not already explored/in frontier)
                if (!explored.contains(neighbor) && !frontier.contains(neighbor)) {
                    //if vertex not been explored & not in already frontier, add to frontier
                    frontier.add(neighbor);
                    neighbor.prev = curr;
                }
            }

            if (curr == this.getVertex(t)) {        //if reached target vertex, terminate
                return;
            }
        }
    }

    /**
     * Helper class for Dijkstra and A*, used in priority queue
     */
    private class CostVertex implements Comparable<CostVertex> {
        double cost;
        Vertex vertex;

        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
        }

        public int compareTo(CostVertex o) {
            return Double.compare(cost, o.cost);
        }
    }

    /**
     * Find the shortest path from vertex with name s to vertex with name t, using Dijkstra
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void Dijkstra(String s, String t) {
        resetAllVertices();

        PriorityQueue<CostVertex> frontier = new PriorityQueue<>();

        //cost is initialized to infinity for all vertices

        //set cost for start to 0
        Vertex startV = this.getVertex(s);
        startV.cost = 0;
        CostVertex start = new CostVertex(startV.cost, startV);
        frontier.add(start);

        while (!frontier.isEmpty()){
            // Visit vertex with minimum distance from startV
            CostVertex current = frontier.poll();

            if (current.vertex.name.equals(t)){
                break;
            }

            //for each neighbor of current vertex
            for (int i = 0; i < current.vertex.sentEdges.size(); i++){
                Vertex neighbor = current.vertex.sentEdges.get(i).target;
                double edgeWeight = current.vertex.sentEdges.get(i).distance;
                double altPathDistance = current.cost + edgeWeight;

                // If shorter path from startV to adjV is found,
                // update adjV's distance and predecessor
                if (altPathDistance < neighbor.cost){
                    neighbor.cost = altPathDistance;
                    neighbor.prev = current.vertex;
                    CostVertex next = new CostVertex(neighbor.cost, neighbor);
                    frontier.add(next);
                }
            }
        }
    }

    /**
     * Helper method to calculate the h value in A*
     *
     * @param cur the current vertex being explored
     * @param goal the goal vertex to reach
     * @return the h value of cur and goal vertices
     */
    private double hValue(String cur, String goal) {
        Vertex start = this.getVertex(cur);
        Vertex end = this.getVertex(goal);
        return computeEuclideanDistance(start.x, start.y, end.x, end.y);
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    public void AStar(String s, String t) {
        resetAllVertices();

        PriorityQueue<CostVertex> frontier = new PriorityQueue<>();

        //cost is initialized to infinity for all vertices

        //set cost for start to 0
        Vertex startV = this.getVertex(s);
        startV.cost = 0;
        CostVertex start = new CostVertex(startV.cost, startV);
        frontier.add(start);

        while (!frontier.isEmpty()){
            // Visit vertex with minimum distance from startV
            CostVertex current = frontier.poll();

            if (current.vertex.name.equals(t)){
                break;
            }

            //for each neighbor of current vertex
            for (int i = 0; i < current.vertex.sentEdges.size(); i++){
                Vertex neighbor = current.vertex.sentEdges.get(i).target;
                double edgeWeight = current.vertex.sentEdges.get(i).distance;
                double altPathDistance = current.cost + edgeWeight + hValue(neighbor.name, t);

                // If shorter path from startV to adjV is found,
                // update adjV's distance and predecessor
                if (altPathDistance < neighbor.cost){
                    neighbor.cost = altPathDistance;
                    neighbor.prev = current.vertex;
                    CostVertex next = new CostVertex(neighbor.cost, neighbor);
                    frontier.add(next);
                }
            }
        }
    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) {
        ArrayList<Edge> path = new ArrayList<>(); // list of edges that connect vertices
        ArrayList<Edge> finalPath = new ArrayList<>();

        Vertex curr = this.getVertex(t);

        //follow the path of previous vertices starting from the end vertex
        while (curr.prev != null){
            for (int i = 0; i < curr.prev.sentEdges.size() ; i++){
                if (curr.prev.sentEdges.get(i).target == curr){
                    path.add(curr.prev.sentEdges.get(i));
                    curr = curr.prev;
                    break;
                }
            }
        }

        //reserve the order
        for (int i = path.size() - 1; i >=0; i--){
            finalPath.add(path.get(i));
        }

        return finalPath;
    }

}