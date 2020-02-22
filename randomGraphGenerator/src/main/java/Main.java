import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.Collectors;

public class Main {
    public List<Integer> getActualComponentList() {
        return actualComponentList;
    }

    public void setActualComponentList(List<Integer> actualComponentList) {
        this.actualComponentList = actualComponentList;
    }

    private List<Integer> actualComponentList = new ArrayList<>();

    public static void main(String[] args) {

        final int NUMBER_OF_NODES = 550;


        TreeMap<Integer, List<Integer>> graphOfVertices = generateRandomGraph(NUMBER_OF_NODES);
        TreeMap<Integer, List<Integer>> graphOfVerticesWithEdges = generateRandomEdges(graphOfVertices, 0.0007);
        System.out.println("for global average degree equal lower than 1: " + calcGlobalAverage(graphOfVerticesWithEdges));
        int[][] matrixWithOfProbLessThanZero = floydWarshalAlgorithm(graphOfVerticesWithEdges);
        int avgOfProbLessThanZero = findAverageOfGraph(matrixWithOfProbLessThanZero);
        System.out.println("avg of graph for equal lower than 1 : " + avgOfProbLessThanZero);
        System.out.println("avg distance of nodes is : " + averagePathNode(matrixWithOfProbLessThanZero));
        DFS dds = new DFS(NUMBER_OF_NODES, graphOfVerticesWithEdges);
        ccsOfNodes(graphOfVerticesWithEdges, matrixWithOfProbLessThanZero);
        TreeMap<Integer, Integer> kk = findNodesWithSameDegree(graphOfVerticesWithEdges, avgOfProbLessThanZero);
        dds.connectedComponents();

        XYChart ccOfDegreesChart = new XYChartBuilder().width(600).height(500).title("Degree distrubution of nodes' degrees").xAxisTitle("X - nodes").yAxisTitle("Y - amount").build();
        List<Integer> vv = kk.values().stream().collect(Collectors.toList());
        List<Integer> nodes = createNodesUpToMax(avgOfProbLessThanZero);
        ccOfDegreesChart.addSeries("Degree distribution of nodes' degrees with lower prob", nodes, vv);
        ccOfDegreesChart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        ccOfDegreesChart.getStyler().setChartTitleVisible(false);
        ccOfDegreesChart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        ccOfDegreesChart.getStyler().setMarkerSize(16);
        new SwingWrapper(ccOfDegreesChart).displayChart();


        System.out.println();
        TreeMap<Integer, List<Integer>> graphOfVerticesWithProbEqual = generateRandomGraph(NUMBER_OF_NODES);
        TreeMap<Integer, List<Integer>> graphOfVerticesWithEdgesWithProbEqual = generateRandomEdges(graphOfVerticesWithProbEqual, 0.001); //0.001842
        System.out.println("for global average degree equal to 1: " + calcGlobalAverage(graphOfVerticesWithEdgesWithProbEqual));
        int[][] matrixWithProbEqual = floydWarshalAlgorithm(graphOfVerticesWithEdgesWithProbEqual);
        int avgWithProbEqual = findAverageOfGraph(matrixWithProbEqual);
        System.out.println("avg of graph for degree equal to 1 : " + avgWithProbEqual);
        System.out.println("avg distance of nodes is : " + averagePathNode(matrixWithProbEqual));
        ccsOfNodes(graphOfVerticesWithEdgesWithProbEqual, matrixWithProbEqual);
        DFS dd = new DFS(NUMBER_OF_NODES, graphOfVerticesWithEdgesWithProbEqual);

        TreeMap<Integer, Integer> kkequal = findNodesWithSameDegree(graphOfVerticesWithEdgesWithProbEqual, avgWithProbEqual);
        XYChart distDegreesChartEqual = new XYChartBuilder().width(600).height(500).title("Degree distrubution of nodes' degrees").xAxisTitle("X - nodes").yAxisTitle("Y - amount").build();
        List<Integer> vvEq = kkequal.values().stream().collect(Collectors.toList());
        List<Integer> nodesEq = createNodesUpToMax(avgWithProbEqual);
        distDegreesChartEqual.addSeries("Degree distribution of nodes' degrees with equal prob", nodesEq, vvEq);
        distDegreesChartEqual.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        distDegreesChartEqual.getStyler().setChartTitleVisible(false);
        distDegreesChartEqual.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        distDegreesChartEqual.getStyler().setMarkerSize(16);
        new SwingWrapper(distDegreesChartEqual).displayChart();

        dd.connectedComponents();

        int matrixArray[][] = new int[graphOfVerticesWithEdgesWithProbEqual.size()][graphOfVerticesWithEdgesWithProbEqual.size()];
        for (int i = 0; i < matrixArray.length; i++) {
            for (int j = 0; j < matrixArray.length; j++) {
                matrixArray[i][j] = 0;
            }
        }
        for (int i = 0; i < graphOfVerticesWithEdgesWithProbEqual.keySet().size(); i++) {
            for (int j = 0; j < graphOfVerticesWithEdgesWithProbEqual.get(i).size(); j++) {
                if (j == graphOfVerticesWithEdgesWithProbEqual.get(i).get(j)) {
                    matrixArray[i][graphOfVerticesWithEdgesWithProbEqual.get(i).get(j)] = 1;
                    matrixArray[graphOfVerticesWithEdgesWithProbEqual.get(i).get(j)][i] = 1;
                }
            }
        }
//        System.out.println();
//        for (int i = 1; i < matrixArray.length; i++) {
//            System.out.print(i < 10 ? i + "=>  " : i + "=> ");
//            for (int j = 1; j < matrixArray.length; j++) {
//                System.out.print(matrixArray[i][j] + " ");
//
//            }
//            System.out.println();
//        }


        System.out.println();
        TreeMap<Integer, List<Integer>> graphOfVerticesWithProbHigher = generateRandomGraph(NUMBER_OF_NODES);
        TreeMap<Integer, List<Integer>> graphOfVerticesWithEdgesWithProbHigher = generateRandomEdges(graphOfVerticesWithProbHigher, 0.002); //0.0028
        System.out.println("for global average degree equal higher than 1: " + calcGlobalAverage(graphOfVerticesWithEdgesWithProbHigher));
        int[][] matrixWithProbHigher = floydWarshalAlgorithm(graphOfVerticesWithEdgesWithProbHigher);
        int avgWithProbHigher = findAverageOfGraph(matrixWithProbHigher);
        System.out.println("avg of graph for degree higher than 1 : " + avgWithProbHigher);
        System.out.println("avg distance of nodes is : " + averagePathNode(matrixWithProbHigher));
        ccsOfNodes(graphOfVerticesWithEdgesWithProbHigher, matrixWithProbHigher);
        DFS ddh = new DFS(NUMBER_OF_NODES, graphOfVerticesWithEdgesWithProbHigher);
        ddh.connectedComponents();

        TreeMap<Integer, Integer> kkHigher = findNodesWithSameDegree(graphOfVerticesWithEdgesWithProbHigher, avgWithProbHigher);
        XYChart distDegreesChartHigher = new XYChartBuilder().width(600).height(500).title("Degree distribution of nodes' degrees with higher prob").xAxisTitle("X - nodes").yAxisTitle("Y - amount").build();
        List<Integer> vvHigh = kkHigher.values().stream().collect(Collectors.toList());
        List<Integer> nodesHigh = createNodesUpToMax(avgWithProbHigher);
        distDegreesChartHigher.addSeries("Degree distribution of nodes' degrees with higher prob", nodesHigh, vvHigh);
        distDegreesChartHigher.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
        distDegreesChartHigher.getStyler().setChartTitleVisible(false);
        distDegreesChartHigher.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        distDegreesChartHigher.getStyler().setMarkerSize(16);
        new SwingWrapper(distDegreesChartHigher).displayChart();

    }


    private static TreeMap<Integer, List<Integer>> generateRandomGraph(int numberOfVertices) {
        TreeMap<Integer, List<Integer>> graphOfVertices = new TreeMap<>();
        for (int i = 0; i < numberOfVertices; i++) {
            graphOfVertices.put(i, new ArrayList<>());
        }

        return graphOfVertices;
    }

    private static TreeMap<Integer, List<Integer>> generateRandomEdges(TreeMap<Integer, List<Integer>> generatedGraph, double probability) {
        for (int i = 0; i < generatedGraph.keySet().size(); i++) {
            for (int j = 0; j < generatedGraph.keySet().size(); j++) {
                double rand = Math.random();
//
                if (i != j) {
                    if (!generatedGraph.get(i).contains(j) && !generatedGraph.get(j).contains(i)) {
                        if (rand < probability) {
                            List<Integer> verticesOfCurrentNode = generatedGraph.get(i);
                            List<Integer> verticesOfOppositeNode = generatedGraph.get(j);
                            if (verticesOfCurrentNode == null) {
                                verticesOfCurrentNode = new ArrayList<>();
                            }
                            if (verticesOfOppositeNode == null) {
                                verticesOfOppositeNode = new ArrayList<>();
                            }
                            verticesOfOppositeNode.add(i);
                            generatedGraph.put(j, verticesOfOppositeNode);
                            verticesOfCurrentNode.add(j);
                            generatedGraph.put(i, verticesOfCurrentNode);
                        }
                    }
                }
            }
        }
        return generatedGraph;
    }

    public static double calcGlobalAverage(TreeMap<Integer, List<Integer>> generatedGraph) {
        int totalAmountOfVertices = 0;
        for (int i = 0; i < generatedGraph.keySet().size(); i++) {
            totalAmountOfVertices += generatedGraph.get(i).size();
        }
        System.out.println("totalAmountOfVertices:" + totalAmountOfVertices / 2);
        return (double) totalAmountOfVertices / (double) 550;
    }

    private static int[][] floydWarshalAlgorithm(TreeMap<Integer, List<Integer>> graphMap) {
        int V = graphMap.keySet().size();
        int dist[][] = new int[V][V];
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                dist[i][j] = !graphMap.get(i).contains(j) ? 65336 : 1;
            }
        }
        for (int k = 1; k < V; k++) {
            // Pick all vertices as source one by one
            for (int i = 1; i < V; i++) {
                // Pick all vertices as destination for the
                // above picked source
                for (int j = 1; j < V; j++) {
                    // If vertex k is on the shortest path from
                    // i to j, then update the value of dist[i][j]
                    if (dist[i][k] + dist[k][j] < dist[i][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
        return dist;
    }

    private static int findAverageOfGraph(int[][] pathArray) {
        int maxExcentricity = 0;
        for (int[] arr : pathArray) {
            int maxExcentricityOfNode = 0;
            if (Arrays.stream(arr).filter(value -> value != 65336).max().isPresent()) {
                maxExcentricityOfNode = Arrays.stream(arr).filter(value -> value != 65336).max().getAsInt();
            }
            maxExcentricity = maxExcentricityOfNode > maxExcentricity ? maxExcentricityOfNode : maxExcentricity;
        }
        return maxExcentricity;
    }

    private static double averagePathNode(int[][] pathArray) {
        int totalNumberOfshortesDistance = 0;

        for (int i = 1; i < pathArray.length; i++) {
            for (int j = 1; j < pathArray.length; j++) {
                totalNumberOfshortesDistance += pathArray[i][j] == 65336 ? 0 : pathArray[i][j];
            }
        }
//        double number = (double) 1 / pathArray.length - 1;
        int V = pathArray.length;
        double number = (double) 2 / ((V) * (V - 1));
        return (number * totalNumberOfshortesDistance) / 2;
    }

    private static void ccsOfNodes(TreeMap<Integer, List<Integer>> neighboursOfAllNodes, int[][] matrixArray) {
        TreeMap<Integer, Double> nodesWithCoeficients = new TreeMap<>();
        for (int i = 0; i < matrixArray.length; i++) {
            nodesWithCoeficients.put(i, 0.0);
        }
        double totalCoefNumber = 0.0;
        for (int node : neighboursOfAllNodes.keySet()) {

            double dsdsn = 0;
            List<Integer> listOfNeighbours = neighboursOfAllNodes.get(node);
            for (int i = 0; i < neighboursOfAllNodes.get(node).size(); i++) {
                dsdsn += calculateNumberOfVerticesOfNeighbours(neighboursOfAllNodes.get(node).get(i), listOfNeighbours, neighboursOfAllNodes);
            }
            int numberOfMaxNeighoubrs = (neighboursOfAllNodes.get(node).size() * ((neighboursOfAllNodes.get(node)).size() - 1));
            double result = ((double) (dsdsn)) / (double) numberOfMaxNeighoubrs;
            if (Double.isNaN(result)) {
                nodesWithCoeficients.put(node, 0.0);


            } else {
                nodesWithCoeficients.put(node, result);
                totalCoefNumber += result;
            }


//            System.out.println("For node " + node + ", the CC is: " + nodesWithCoeficients.get(node));
        }
        System.out.println("Clustering effect : " + totalCoefNumber / neighboursOfAllNodes.size());
    }

    private static double calculateNumberOfVerticesOfNeighbours(int node, List<Integer> neighboursList, TreeMap<Integer, List<Integer>> neighboursOfAllNodes) {
        int totalNumberOfVertices = 0;
        for (int i = 0; i < neighboursOfAllNodes.get(node).size(); i++) {
            if (neighboursOfAllNodes.get(node).contains(i) && neighboursList.contains(i)) {
                totalNumberOfVertices++;
            }
        }
        return (double) (totalNumberOfVertices);
    }


    private static TreeMap<Integer, Integer> findNodesWithSameDegree(TreeMap<Integer, List<Integer>> neighboursOfAllNodes, int maxNode) {
        TreeMap<Integer, Integer> nodesWithSameDegree = new TreeMap<>();
        for (int i = 0; i <= maxNode; i++) {
            nodesWithSameDegree.put(i, 0);
        }

        for (int i = 0; i < neighboursOfAllNodes.keySet().size(); i++) {
            for (int j = 0; j <= maxNode; j++) {
                if (j == neighboursOfAllNodes.get(i).size()) {
//                        List<Integer> listOfNodes = nodesWithSameDegree.get(j);
                    int actual_size = nodesWithSameDegree.get(j);
                    actual_size++;
//                    listOfNodes.add(i);
                    nodesWithSameDegree.put(j, actual_size);
                }
            }

        }
        return nodesWithSameDegree;

    }

    private static List<Integer> createNodesUpToMax(int maxNode) {
        List<Integer> nodes = new ArrayList<>();
        for (int i = 0; i <= maxNode; i++) {
            nodes.add(i);
        }
        return nodes;
    }

//
//    private List<Integer> DFSUtil(int v, boolean[] visited, TreeMap<Integer, List<Integer>> nodesWithVertices) {
//        // Mark the current node as visited and print it
//        visited[v] = true;
//        List<Integer> componentList = new ArrayList<>();
//        componentList.add(v);
//        System.out.print(v + " ");
//        this.actualComponentList.add(v);
//        // Recur for all the vertices
//        // adjacent to this vertex
//        for (int x : nodesWithVertices.get(v)) {
//            if (!visited[x]) {
//                DFSUtil(x, visited, nodesWithVertices);
//                componentList.add(x);
//                this.actualComponentList.add(x);
//
//            }
//        }
//        return componentList;
//    }
//
//      static void connectedComponents(int V, TreeMap<Integer, List<Integer>> nodesWithVertices) {
//        // Mark all the vertices as not visited
//        boolean[] visited = new boolean[V];
//        for (int v = 0; v < V; ++v) {
//            if (!visited[v]) {
//
//                // print all reachable vertices
//                // from v
//                List<Integer> componentList = DFSUtil(v, visited, nodesWithVertices);
//                if(componentList.size()>1){
////                    System.out.print("Component  of size " + componentList.size() + " with nodes: "+componentList.toString());
////                    System.out.println();
//                }
//                System.out.println();
//                setActualComponentList(new ArrayList<>());
////                DFSUtil(v, visited, nodesWithVertices);
////                System.out.print("Component " + v + " of graph : ");
////                System.out.println();
//
//            }
//        }
//    }


}
