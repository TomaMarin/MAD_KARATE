import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategorySeries;
import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.style.Styler;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class ListMain {
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
        HashMap<Integer, List<Integer>> listHashMap = new HashMap<>();

        for (int i = 0; i <= maxArrayNumber; i++) {
            listHashMap.put(i, new ArrayList<>());
        }


        for (int i = 0; i <= maxArrayNumber; i++) {
            for (int j = 0; j < records.size(); j++) {
                if (i == Integer.parseInt(records.get(j).get(0))) {
                    listHashMap.get(i).add(Integer.parseInt(records.get(j).get(1)));
                    listHashMap.get(Integer.parseInt(records.get(j).get(1))).add(i);
                }
            }
        }

        System.out.println();
        System.out.println("max size of node: " + findMax(listHashMap));
        System.out.println("min size of node: " + findMin(listHashMap));
        System.out.println("global avg: " + calculateAvg(listHashMap));
        System.out.println();

        for (int  i = 1 ;i < listHashMap.keySet().size();i++) {
            System.out.println("Vertice => " + i + " has neighbors : " + listHashMap.get(i));
        }
        HashMap<Integer, Integer> histogramData = calculateFrequencyOfDegrees(maxArrayNumber,listHashMap);

        CategoryChart chart = new CategoryChart(600, 500);
        int[] yData = calculateFrequency(histogramData, findMax(listHashMap));
        int[] xData = new int[findMax(listHashMap) + 1];
        for (int i = 0; i < xData.length; i++) {
            xData[i] = i;
        }
        chart.addSeries("Frequency of list ", (xData), yData);
        chart.getStyler().setDefaultSeriesRenderStyle(CategorySeries.CategorySeriesRenderStyle.Bar);
        chart.getStyler().setChartTitleVisible(false);
        chart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        chart.getStyler().setMarkerSize(16);

        new SwingWrapper(chart).displayChart();

        double[] yDataofRelativeFreq = calculateRelativeFrequency(yData, maxArrayNumber);
        double[] xDataofRelativefreq = new double[xData.length];
        for (int i = 0; i < xDataofRelativefreq.length; i++) {
            xDataofRelativefreq[i] = xData[i];
        }

        CategoryChart relFrequencyChart = new CategoryChart(600, 500);
        relFrequencyChart.addSeries("Relative Frequency of list", (xDataofRelativefreq), yDataofRelativeFreq);
        relFrequencyChart.getStyler().setDefaultSeriesRenderStyle(CategorySeries.CategorySeriesRenderStyle.Bar);
        relFrequencyChart.getStyler().setChartTitleVisible(false);
        relFrequencyChart.getStyler().setLegendPosition(Styler.LegendPosition.OutsideE);
        relFrequencyChart.getStyler().setMarkerSize(16);
        new SwingWrapper(relFrequencyChart).displayChart();

    }

    private static double calculateAvg(HashMap<Integer, List<Integer>> data) {
        int avg = 0;
        double result;
        int correctNodes = 0;
        for (List<Integer> integerList : data.values()) {
            avg += integerList.size();
            correctNodes++;
        }
        result = ((double)  avg) / ((double) correctNodes - 1);
        return result;
    }

    private static int findMax(HashMap<Integer, List<Integer>> data) {
        int max = 0;
        for (List<Integer> integerList : data.values()) {
            max = max > integerList.size() ? max : integerList.size();
        }
        return max;

    }

    private static int findMin(HashMap<Integer, List<Integer>> data) {
        int min = data.get(1).size();
            for (int i = 1; i < data.values().size(); i++) {
                min = min < data.get(i).size() ? min : data.get(i).size();
            }
        return min;

    }


    //
    private static HashMap<Integer, Integer> calculateFrequencyOfDegrees(int maxNumber, HashMap<Integer, List<Integer>> matrix) {
        HashMap<Integer, Integer> data = new HashMap<>();
        for (int i = 0; i < maxNumber+1; i++) {
            data.put(i, 0);
        }

//        for (int i = 1; i < maxNumber ; i++) {
//            for (int j = 1; j <maxNumber ; j++) {
//                if (matrix.get(i).get(j)){
//                    data.replace(j,data.get(j)+1);
//                    data.replace(i,data.get(i)+1);
//
//                }
//            }
//        }


        for (Integer key:
                matrix.keySet()) {
            for (int i = 0; i < matrix.get(key).size(); i++) {
                int  verticeValue = matrix.get(key).get(i);
                data.replace(key,data.get(key)+1);
                data.replace(verticeValue,data.get(verticeValue)+1);

            }
        }

        return data;
    }



    private static int[] calculateFrequency(HashMap<Integer, Integer> data, int maxVal){
        int [] freqs = new int[maxVal+1];
        for (int i = 1; i <freqs.length; i++) {
            freqs[i]= 0;
            for (int j = 1; j < data.size(); j++) {
                if (i == data.get(j)){
                    freqs[i]++;
                }
            }
        }
        return  freqs;
    }

    private static double[] calculateRelativeFrequency(int[] data, int N) {
        double[] freqs = new double[data.length];
        for (int i = 1; i < freqs.length; i++) {
            freqs[i] = ((double) data[i]) / (double) N;
        }

        return freqs;
    }

}
