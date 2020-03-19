	import java.util.*;
public class House {
	// you can also use imports, for example:

	// you can write to stdout for debugging purposes, e.g.
	// System.out.println("this is a debug message");

	    public int[] solution(int[] store, int[] house) {
	        // write your code in Java SE 8
	        int[] rst = new int[house.length];
		        Arrays.sort(store);
		        for(int i =0; i< house.length; i++) {
		        	int idx = Arrays.binarySearch(store, house[i]);
		        	if(idx<0) {
		        		idx = -(idx+1);
		        		if(idx == 0) rst[i] = store[0];
		        		else if(idx == store.length) rst[i] = store[store.length-1];
		        		else if(store[idx]-house[i]>= house[i] - store[idx-1]) {
		        			rst[i] = store[idx-1];
		        		}else {
		        			rst[i] = store[idx];
		        		}
		        	}else {
		        		rst[i]= store[idx];
		        	}
		        }
		        return rst;
	    }

	    public int solution(int[] A) {
	        // write your code in Java SE 8
	        int n = A.length;
	        int sum = 0;
	        for(int i =0; i< A.length; i++){
	            sum+= A[i];
	        }
	        int half = sum/2;
	        int[][] dp = new int[n][half+1];
	        for(int j = 0; j<=half; j++){
	            if(A[0]<=j){
	                dp[0][j] = A[0];
	            }else  dp[0][j] = 0;
	        }
	        // dp[i][j] stores the maximum knapsack value using the items 0,...,i, and with total weight <= j
	        for(int i = 1; i<n; i++){
	            for (int j = 0; j<=half; j++){
	                if(A[i]>j){
	                    dp[i][j] = dp[i-1][j];
	                }else{
	                    dp[i][j] = Math.max(dp[i-1][j], A[i] + dp[i-1][j-A[i]]);
	                }
	            }
	        }
	        int val = dp[n-1][half];
	        System.out.println("val "+val);

	        if(val*2 == sum) return 0;
	        if(val*2 > sum) return val - (sum-val);
	        return sum - val -val;


	    }

}
