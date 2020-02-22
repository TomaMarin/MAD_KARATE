import java.util.Arrays;
import java.util.List;
import java.util.Objects;

public class ClusterRow {
    private double[] probabilities;
    private List<Integer> indexes;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ClusterRow that = (ClusterRow) o;
        return Arrays.equals(probabilities, that.probabilities) &&
                Objects.equals(indexes, that.indexes);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(indexes);
        result = 31 * result + Arrays.hashCode(probabilities);
        return result;
    }

    public double[] getProbabilities() {
        return probabilities;
    }

    public void setProbabilities(double[] probabilities) {
        this.probabilities = probabilities;
    }

    public List<Integer> getIndexes() {
        return indexes;
    }

    public void setIndexes(List<Integer> indexes) {
        this.indexes = indexes;
    }

    public ClusterRow() {
    }

    public ClusterRow(double[] probabilities, List<Integer> indexes) {
        this.probabilities = probabilities;
        this.indexes = indexes;
    }

    public boolean isMergedRow() {
        for (double prob : getProbabilities()) {
            if (prob > 0.0000) {
                return false;
            }
        }
        return true;
    }
}
