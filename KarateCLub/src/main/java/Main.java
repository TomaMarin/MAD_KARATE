import com.sun.prism.paint.Color;
import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategorySeries;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class Main {

    public static void main(String[] args) throws IOException {
        List<List<String>> records = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader("KarateClub.csv"))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] values = line.split(";");
                records.add(Arrays.asList(values));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        int maxArrayNumber = Integer.parseInt(records.get(records.size() - 1).get(0));


        //Matrix array init
        int[][] matrixArray = new int[maxArrayNumber + 1][maxArrayNumber + 1];
        for (int i = 1; i < maxArrayNumber; i++) {
            for (int j = 1; j < maxArrayNumber; j++) {
                matrixArray[i][j] = 0;
            }
        }

//        Matrix array value init
        for (int i = 0; i <= maxArrayNumber; i++) {
            for (int j = 0; j < records.size(); j++) {
                if (i == Integer.parseInt(records.get(j).get(0))) {
                    matrixArray[i][Integer.parseInt(records.get(j).get(1))] = 1;
                    matrixArray[Integer.parseInt(records.get(j).get(1))][i] = 1;
                }
            }
        }

        int[][] matrixArrayWithoutZero = new int[matrixArray.length - 1][matrixArray.length - 1];
        for (int i = 1; i < matrixArray.length; i++) {
            for (int j = 1; j < matrixArray.length; j++) {
                matrixArrayWithoutZero[i - 1][j - 1] = matrixArray[i][j];
            }
        }
        System.out.println(matrixArray.length);
        System.out.println();

        for (int i = 1; i < matrixArray.length; i++) {
            System.out.print(i < 10 ? i + "=>  " : i + "=> ");
            for (int j = 1; j < matrixArray.length; j++) {
                System.out.print(matrixArray[i][j] + " ");

            }
            System.out.println();
        }

        System.out.println();
        System.out.println("max degree of node: " + findMax(matrixArray));
        System.out.println("min degree size of node: " + findMin(matrixArray));
        System.out.println("global degree avg: " + calculateAvg(matrixArray));

        System.out.println("-----------------------------------------------------");
        int[][] pathMatrix = floydWarshalAlgorithm(matrixArray);
        for (int i = 1; i < pathMatrix.length; i++) {
            System.out.print(i < 10 ? i + "=>  " : i + "=> ");
            for (int j = 1; j < pathMatrix.length; j++) {
                System.out.print(pathMatrix[i][j] + " ");

            }
            System.out.println();
        }


        double aveg = averagePathNode(pathMatrix);
        HashMap<Integer, Integer> histogramData = calculateFrequencyOfDegrees(matrixArray);

        System.out.println("avg of distances is " + aveg);
        System.out.println("avg of graph is " + findAverageOfGraph(pathMatrix));
        closenessCentrality(pathMatrix);
        CategoryChart chart = new CategoryChart(600, 500);
        int[] yData = calculateFrequency(histogramData, findMax(matrixArray));
        int[] xData = new int[findMax(matrixArray) + 1];
        for (int i = 0; i < xData.length; i++) {
            xData[i] = i;
        }

        chart.addSeries("Frequency of matrix", (xData), yData);
        chart.getStyler().setDefaultSeriesRenderStyle(CategorySeries.CategorySeriesRenderStyle.Bar);
        chart.getStyler().setChartTitleVisible(false);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        chart.getStyler().setMarkerSize(16);

//        new SwingWrapper(chart).displayChart();

        double[] yDataofRelativeFreq = calculateRelativeFrequency(yData, maxArrayNumber);
        double[] xDataofRelativefreq = new double[xData.length];
        for (int i = 0; i < xDataofRelativefreq.length; i++) {
            xDataofRelativefreq[i] = xData[i];
        }

        CategoryChart relFrequencyChart = new CategoryChart(600, 500);
        relFrequencyChart.addSeries("Relative Frequency of matrix", (xDataofRelativefreq), yDataofRelativeFreq);
        relFrequencyChart.getStyler().setDefaultSeriesRenderStyle(CategorySeries.CategorySeriesRenderStyle.Bar);
        relFrequencyChart.getStyler().setChartTitleVisible(false);
        relFrequencyChart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        relFrequencyChart.getStyler().setMarkerSize(16);
//        new SwingWrapper(relFrequencyChart).displayChart();

        TreeMap<Integer, List<Integer>> neighboursOfAllNodes = new TreeMap<>();
        for (int i = 1; i < matrixArray.length; i++) {
            neighboursOfAllNodes.put(i, new ArrayList<>());
        }
        findAllNeighboursOfNode(1, matrixArray);
        for (int i = 1; i <= neighboursOfAllNodes.keySet().size(); i++) {

            neighboursOfAllNodes.put(i, findAllNeighboursOfNode(i, matrixArray));
        }

//        int node = 2;
        TreeMap<Integer, Double> nodesWithCoeficients = new TreeMap<>();
        for (int i = 1; i < 35; i++) {
            nodesWithCoeficients.put(i, 0.0);
        }

        double totalCoefNumber = 0.0;
        System.out.println();
        for (Integer node : nodesWithCoeficients.keySet()) {


            int dsdsn = 0;
            List<Integer> listOfNeighbours = neighboursOfAllNodes.get(node);
            for (int i = 0; i < neighboursOfAllNodes.get(node).size(); i++) {
                dsdsn += calculateNumberOfVerticesOfNeighbours(neighboursOfAllNodes.get(node).get(i), neighboursOfAllNodes.get(node), matrixArray);
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

        totalCoefNumber = (totalCoefNumber) / 34;
        System.out.println("Total average CC of nodes is: " + totalCoefNumber);
        TreeMap<Integer, List<Integer>> nodesWithSameDegree = findNodesWithSameDegree(histogramData, 17);
        double[] ccOfDegrees = calculateCCOfDegreesOfNodes(nodesWithSameDegree, nodesWithCoeficients);
        double[] degrees = new double[18];
        for (int i = 1; i < nodesWithSameDegree.keySet().size(); i++) {
            if (i == 0)
                continue;
            if (i != 1 && ccOfDegrees[i] == 0)
                continue;
            degrees[i] = (double) (i);
        }
        ccOfDegrees[0] = Double.NaN;
        degrees[0] = Double.NaN;

        double[][] similarityMatrix = createSimilarityMatrix(matrixArrayWithoutZero);
        System.out.println(similarityMatrix);
        writeValuesToFile(similarityMatrix, "similarityMatrix.csv");
        ClusterRow[] clusteringResults = findMostSimilarCluster(similarityMatrix);
        double k = 0.0;
        Graph graph = new SingleGraph("Tutorial 1");
        graph.addAttribute("ui.stylesheet", "url('file:E:/Documents/skola/ING/1.semester/MADI/git_repo/MAD_KARATE/KarateCLub/src/main/resources/ui.stylesheet.css')");

        for (int i = 0; i < clusteringResults.length; i++) {
            if (!clusteringResults[i].getIndexes().isEmpty()) {
                for (int j = 0; j < clusteringResults[i].getIndexes().size(); j++) {

                    graph.addNode(Integer.toString(clusteringResults[i].getIndexes().get(j))).setAttribute("ui.color", k);
                    graph.getNode(Integer.toString(clusteringResults[i].getIndexes().get(j))).addAttribute("ui.label", Integer.toString(clusteringResults[i].getIndexes().get(j) +1));
                }
                k+=0.5;
            }
        }

//        for (int i = 0; i < matrixArray.length; i++) {
//            graph.addNode(Integer.toString(i)).setAttribute("ui.color", 0);
//        }
        for (int i = 0; i < matrixArrayWithoutZero.length; i++) {
            for (int j = i; j < matrixArrayWithoutZero.length; j++) {
                if (matrixArrayWithoutZero[i][j] == 1) {

                    graph.addEdge((Integer.toString(i).concat(Integer.toString(j))), i, j);
                }
            }
        }
        graph.display();
    }

    private static void writeValuesToFile(double[][] vals, String fileName) throws IOException {
        FileWriter writer = new FileWriter(fileName, false);
        for (int i = 0; i < vals.length; i++) {
            for (int j = 0; j < vals.length; j++) {
                if (j == vals.length - 1) {
                    writer.write(vals[i][j] + "");
                } else {
                    writer.write(vals[i][j] + ",");
                }
            }
            writer.write("\r\n");   // write new line
        }
        writer.close();
    }

    private static int returnIndexOfClusterByElementIndex(int index, ClusterRow[] rows) {
        for (int i = 0; i < rows.length; i++) {
            if (rows[i].getIndexes().contains(index)) {
                return i;
            }
        }
        return Integer.MAX_VALUE;
    }

    private static boolean areMoreZeroedRowsThanNumber(int number, ClusterRow[] rows) {
        int numberOfNonZeroedRows = 0;
        for (ClusterRow r : rows) {
            if (!r.isMergedRow()) {
                numberOfNonZeroedRows++;
            }
        }
        return numberOfNonZeroedRows >= number;
    }

    private static ClusterRow[] findMostSimilarCluster(double[][] similarityMatrix) {
        for (int i = 0; i < similarityMatrix.length; i++) {
            for (int j = 0; j < similarityMatrix.length; j++) {
                if (i == j) {
                    similarityMatrix[i][j] = 0.0;
                }
            }
        }

        double mostSimilarClusterValue = 0.0;
        int indexI = 0;
        int indexJ = 0;
        double[][] clusterMatrix = similarityMatrix;
        List<Integer> ssss = new ArrayList<>();
        int counter = 0;
        final int N = 3;

        ClusterRow[] rows = new ClusterRow[similarityMatrix.length];
        for (int i = 0; i < rows.length; i++) {
            rows[i] = new ClusterRow(similarityMatrix[i], new ArrayList<>(Collections.singletonList(i)));
        }
        while (areMoreZeroedRowsThanNumber(N, rows)) {
            mostSimilarClusterValue = 0.0;
            indexI = 0;
            indexJ = 0;
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows.length; j++) {
                    if (rows[i].getProbabilities()[j] != 0.0 && i != j && !rows[i].isMergedRow()) {
                        if (mostSimilarClusterValue < rows[i].getProbabilities()[j]) {
                            mostSimilarClusterValue = rows[i].getProbabilities()[j];
                            indexI = i;
                            indexJ = j;
                        }
                    }
                }
            }
            if (returnIndexOfClusterByElementIndex(indexJ, rows) != Integer.MAX_VALUE) {
                indexJ = returnIndexOfClusterByElementIndex(indexJ, rows);
            }

            rows[indexI].setProbabilities(createRowFrom2Rows(indexI, indexJ, rows));
            rows[indexI].setIndexes(setIndexesFrom2Rows(indexI, indexJ, rows));
            double[] probsToZeroOut = rows[indexJ].getProbabilities();
            Arrays.fill(probsToZeroOut, 0.0);
            if (indexI != indexJ) {
                rows[indexJ].setProbabilities(probsToZeroOut);
                rows[indexJ].setIndexes(new ArrayList<>());
            }
        }
//        rows = Arrays.stream(rows).filter(clusterRow -> clusterRow.getIndexes().size() == 0).collect(Collectors.toList()).toArray(rows);
        return rows;
    }

    private static List<Integer> setIndexesFrom2Rows(int iIndexOfRowToMerge, int jIndexOfRowToMerge, ClusterRow[] rows) {
        List<Integer> indexesOfCluster = rows[iIndexOfRowToMerge].getIndexes();
        indexesOfCluster.sort(Integer::compareTo);
        if (!indexesOfCluster.contains(jIndexOfRowToMerge)) {
            indexesOfCluster.add(jIndexOfRowToMerge);
        }
        for (Integer index : rows[jIndexOfRowToMerge].getIndexes()) {
            if (!indexesOfCluster.contains(index)) {
                indexesOfCluster.add(index);
            }
        }
        return indexesOfCluster;
    }

    private static double[] createRowFrom2Rows(int iIndexOfRowToMerge, int jIndexOfRowToMerge, ClusterRow[] rows) {
        double[] mergedRow = new double[rows.length];
        for (int i = 0; i < rows.length; i++) {
            if (!rows[iIndexOfRowToMerge].getIndexes().contains(i)) {
                if (rows[iIndexOfRowToMerge].getProbabilities()[i] > rows[jIndexOfRowToMerge].getProbabilities()[i]) {
                    mergedRow[i] = rows[iIndexOfRowToMerge].getProbabilities()[i];
                } else {
                    mergedRow[i] = rows[jIndexOfRowToMerge].getProbabilities()[i];
                }
            }
        }
        mergedRow[jIndexOfRowToMerge] = 0.0;
        return mergedRow;
    }

    private static int calculateDegreeOfNode(int[] node) {
        int numberOfEdges = 0;
        for (int value : node) {
            if (value == 1) {
                numberOfEdges++;
            }
        }
        return numberOfEdges;
    }

    private static double[][] createSimilarityMatrix(int[][] adjacencyMatrix) {
        double[][] similarityMatrix = new double[adjacencyMatrix.length][adjacencyMatrix.length];
        for (int i = 0; i < similarityMatrix.length; i++) {
            for (int j = 0; j < similarityMatrix.length; j++) {
                int numberOfNeighbors = findNeighbors(adjacencyMatrix[i], adjacencyMatrix[j]);
                int degreeofI = calculateDegreeOfNode(adjacencyMatrix[i]);
                int degreeofJ = calculateDegreeOfNode(adjacencyMatrix[j]);

                similarityMatrix[i][j] = calculateCosineSimilarity(numberOfNeighbors, degreeofI, degreeofJ);
            }
        }
        return similarityMatrix;
    }

    private static double calculateCosineSimilarity(int numberOfNeighbors, double rowADegree, double rowBDegree) {
        return (double) numberOfNeighbors / Math.sqrt(rowADegree * rowBDegree);
    }

    private static int findNeighbors(int[] rowA, int[] rowB) {
        int numberOfNeighbors = 0;
        for (int i = 0; i < rowA.length; i++) {
            if (rowA[i] == 1 && rowB[i] == 1) {
                numberOfNeighbors++;
            }
        }
        return numberOfNeighbors;
    }

    private static HashMap<Integer, Integer[]> findKNearestNeighborsForDataSet(double[][] matrixArray, int k) {
        double[][] neigbors = new double[matrixArray.length][k];
        HashMap<Integer, Integer[]> nearestsNeighbors = new HashMap<>();
        for (int i = 0; i < matrixArray.length; i++) {
            nearestsNeighbors.put(i, findKNearestNeighbors(matrixArray[i], k, i));
        }
        return nearestsNeighbors;
    }


    private static Integer[] findKNearestNeighbors(double[] row, int k, int rowIndex) {
        Integer[] nearestNeighbors = new Integer[k];
        Integer[] example = new Integer[k];
        List<Integer> minIndexes = new LinkedList<>();
        for (int i = 0; i < k; i++) {
            double min = Double.MAX_VALUE;
            int minIndex = 0;


            for (int j = 0; j < row.length; j++) {
                if (row[j] < min && !minIndexes.contains(j) && j != rowIndex) {
                    min = row[j];
                    minIndex = j;
                }
            }
            minIndexes.add(minIndex);
        }
        nearestNeighbors = minIndexes.toArray(example);
        return nearestNeighbors;
    }

    private static double[] calculateCCOfDegreesOfNodes(TreeMap<Integer, List<Integer>> data, TreeMap<Integer, Double> nodesWithCoeficients) {
        double[] ccOfNodes = new double[data.keySet().size()];

        for (int i = 1; i < ccOfNodes.length; i++) {
            double ccValOfAllNdodes = 0.0;
            for (int j = 0; j < data.get(i).size(); j++) {
                ccValOfAllNdodes += nodesWithCoeficients.get(data.get(i).get(j));
            }
            if (Double.isNaN(ccValOfAllNdodes / data.get(i).size()))
                continue;
            ccOfNodes[i] = Double.isNaN(ccValOfAllNdodes / data.get(i).size()) ? 0.0 : ccValOfAllNdodes / data.get(i).size();
        }
        ccOfNodes[0] = Double.NaN;
        return ccOfNodes;
    }

    private static int calculateNumberOfVerticesOfNeighbours(int node, List<Integer> neighboursList, int[][] matrix) {
        int totalNumberOfVertices = 0;
        for (int i = 1; i < matrix[node].length; i++) {
            if (matrix[node][i] == 1 && neighboursList.contains(i)) {
                totalNumberOfVertices++;
            }
        }
        return (totalNumberOfVertices);
    }

    private static List<Integer> findAllNeighboursOfNode(int node, int[][] matrixArray) {
        List<Integer> neighboursOfNode = new ArrayList<>();

        for (int i = 1; i < matrixArray[node].length; i++) {
            if (matrixArray[node][i] == 1) {
                neighboursOfNode.add(i);
            }
        }
        return neighboursOfNode;
    }

    private static int findMax(int[][] matrixArray) {
        int max = 0;
        for (int i = 0; i < matrixArray.length; i++) {
            int actualValue = 0;
            for (int j = 0; j < matrixArray.length; j++) {
                actualValue += matrixArray[i][j];
            }
            if (actualValue > max) {
                max = actualValue;
            }
        }

        return max;

    }

    private static int findMin(int[][] matrixArray) {
        int min = findMax(matrixArray);
        for (int i = 1; i < matrixArray.length; i++) {
            int actualValue = 0;
            for (int j = 0; j < matrixArray.length; j++) {
                actualValue += matrixArray[i][j];
            }
            if (actualValue < min) {
                min = actualValue;
            }
        }

        return min;
    }

    private static double calculateAvg(int[][] matrixArray) {
        double result;
        int correctNodes = 0;
        int numberOfNodes = 0;
        for (int i = 0; i < matrixArray.length; i++) {

            for (int j = 0; j < matrixArray.length; j++) {
                if (matrixArray[i][j] == 1) {
                    numberOfNodes++;
                }
            }
            correctNodes++;
        }

        result = ((double) numberOfNodes) / ((double) correctNodes - 1);
        return result;
    }

    private static HashMap<Integer, Integer> calculateFrequencyOfDegrees(int[][] matrix) {
        HashMap<Integer, Integer> data = new HashMap<>();
        for (int i = 0; i < matrix.length; i++) {
            data.put(i, 0);
        }

        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix.length; j++) {
                if (matrix[i][j] == 1) {
                    data.replace(j, data.get(j) + 1);
                    data.replace(i, data.get(i) + 1);

                }
            }
        }
        return data;
    }

    private static int[] calculateFrequency(HashMap<Integer, Integer> data, int maxVal) {
        int[] freqs = new int[maxVal + 1];
        for (int i = 1; i < freqs.length; i++) {
            freqs[i] = 0;
            for (int j = 1; j < data.size(); j++) {
                if (i == data.get(j) / 2) {
                    freqs[i]++;
                }
            }
        }
        return freqs;
    }

    private static TreeMap<Integer, List<Integer>> findNodesWithSameDegree(HashMap<Integer, Integer> data, int maxNode) {
        TreeMap<Integer, List<Integer>> nodesWithSameDegree = new TreeMap<>();
        for (int i = 0; i <= maxNode; i++) {
            nodesWithSameDegree.put(i, new ArrayList<>());
        }

        for (int i = 0; i < data.keySet().size(); i++) {
            for (int j = 1; j <= maxNode; j++) {
                if (j == data.get(i) / 2) {
                    List<Integer> listOfNodes = nodesWithSameDegree.get(j);
                    listOfNodes.add(i);
                    nodesWithSameDegree.put(j, listOfNodes);
                }
            }

        }
        return nodesWithSameDegree;

    }


    private static double[] calculateRelativeFrequency(int[] data, int N) {
        double[] freqs = new double[data.length];
        for (int i = 1; i < freqs.length; i++) {
            freqs[i] = ((double) data[i]) / (double) N;
        }

        return freqs;
    }


    private static int[][] floydWarshalAlgorithm(int[][] adjacencyMatrix) {
        int V = adjacencyMatrix.length;
        int dist[][] = new int[V][V];
        for (int i = 1; i < V; i++) {
            for (int j = 1; j < V; j++) {
                dist[i][j] = adjacencyMatrix[i][j] == 0 ? 65336 : adjacencyMatrix[i][j];
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
            int maxExcentricityOfNode = Arrays.stream(arr).max().getAsInt();
            maxExcentricity = maxExcentricityOfNode > maxExcentricity ? maxExcentricityOfNode : maxExcentricity;
        }
        return maxExcentricity;
    }

    private static double averagePathNode(int[][] pathArray) {
        int totalNumberOfshortesDistance = 0;

        for (int i = 1; i < pathArray.length; i++) {
            for (int j = 1; j < pathArray.length; j++) {
                totalNumberOfshortesDistance += pathArray[i][j];
            }
        }
//        double number = (double) 1 / pathArray.length - 1;
        int V = pathArray.length - 1;
        double number = (double) 2 / ((V) * (V - 1));
        return (number * totalNumberOfshortesDistance) / 2;
    }

    private static void closenessCentrality(int[][] pathArray) {
        double centrality = 0;
        for (int i = 1; i < pathArray.length; i++) {
            int dOfI = 0;

            for (int j = 1; j < pathArray.length; j++) {
                if (i != j) {
                    dOfI += pathArray[i][j];
                }
            }
            System.out.println("for node " + i + " is closeness " + (double) (pathArray.length - 1) / (double) dOfI);
        }
//         centrality = Arrays.stream(pathArray[1]).sum();
//        return  (double)(pathArray.length)/centrality;
    }

}




