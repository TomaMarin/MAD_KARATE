import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategorySeries;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Main {

    public static void main(String[] args) {
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
                if (i == data.get(j)) {
                    freqs[i]++;
                }
            }
        }
        return freqs;
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
                if(i!=j){
                    dOfI += pathArray[i][j];
                }
            }
            System.out.println("for node "+ i + " is closeness "+ (double)(pathArray.length-1)/(double) dOfI);
        }
//         centrality = Arrays.stream(pathArray[1]).sum();
//        return  (double)(pathArray.length)/centrality;
    }

}




