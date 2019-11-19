import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class DFS {
    public List<Integer> getActualComponentList() {
        return actualComponentList;
    }

    public void setActualComponentList(List<Integer> actualComponentList) {
        this.actualComponentList = actualComponentList;
    }

    private List<Integer> actualComponentList = new ArrayList<>();


    private int V;
    private TreeMap<Integer, List<Integer>> nodesWithVertices;

    public DFS(int v, TreeMap<Integer, List<Integer>> nodesWithVertices) {
        this.V = v;
        this.nodesWithVertices = nodesWithVertices;
    }

    private List<Integer> DFSUtil(int v, boolean[] visited, TreeMap<Integer, List<Integer>> nodesWithVertices) {
        // Mark the current node as visited and print it
        visited[v] = true;
        List<Integer> componentList = new ArrayList<>();
        componentList.add(v);
//        System.out.print(v + " ");
        if(!this.actualComponentList.contains(v)) {
            this.actualComponentList.add(v);
        }
        // Recur for all the vertices
        // adjacent to this vertex
        for (int x : nodesWithVertices.get(v)) {
            if (!visited[x]) {
                DFSUtil(x, visited, nodesWithVertices);
                componentList.add(x);
                if(!this.actualComponentList.contains(x)) {
                    this.actualComponentList.add(x);
                }

            }
        }
        return componentList;
    }

    void connectedComponents() {
        // Mark all the vertices as not visited
        int V = this.V;
        boolean[] visited = new boolean[V];
        for (int v = 0; v < V; ++v) {
            if (!visited[v]) {

                // print all reachable vertices
                // from v
                List<Integer> componentList = DFSUtil(v, visited, this.nodesWithVertices);

                if (this.actualComponentList.size() > 1) {
                    System.out.print("Component  of size " + this.actualComponentList.size() + " with nodes: "+this.actualComponentList.toString());
                    System.out.println();
                }
//                System.out.println();
                setActualComponentList(new ArrayList<>());
//                DFSUtil(v, visited, nodesWithVertices);
//                System.out.print("Component " + v + " of graph : ");
//                System.out.println();

            }
        }
    }
}
