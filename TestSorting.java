import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

public class TestSorting {
	public List<int[]> test(int[][] buildings) {
		List<int[]> height = new ArrayList<>();
	    for(int[] b:buildings) {
	        height.add(new int[]{b[0], -b[2]});
	        height.add(new int[]{b[1], b[2]});
	    }
		Collections.sort(height, (a, b) -> {
	        if(a[0] != b[0]) 
	            return a[0] - b[0];
	        return a[1] - b[1];
	});
		for(int[] h:height) {
			for(int i:h) {
		    	System.out.print(i);
			}
			System.out.println();
		}

		return height;
		
	}
    public static void main(String[] args) {
    	TestSorting hi  = new TestSorting();
    	int[][] buildings = {{2,9,10},{3,7,15},{5,12,12}};
    	hi.test(buildings);
    }
	
}

