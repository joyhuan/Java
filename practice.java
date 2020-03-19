//import java.util.Map;
//import java.util.PriorityQueue;
//import java.util.ArrayDeque;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.Collections;
//import java.util.Comparator;
//import java.util.Deque;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.Iterator;
//import java.util.LinkedHashMap;
//import java.util.LinkedList;
//import java.util.Queue;
//import java.util.Random;
//import java.util.Set;
//import java.util.Stack;
//import java.util.Vector;
//import java.util.stream.Collectors;
import static org.junit.Assert.fail;

import java.time.LocalDate;
import java.util.*;
import java.util.stream.Collectors;

public class practice {

    
    private boolean canJumpHelper(int[] nums,  int end){         
        int len = nums.length; 
        if (len<=1 || end<=1) return true; 
        for (int i = end-1; i>=0; i--){
            if (nums[i]>=end-i){
                if(i==0){
                    return true; 
                }else{
                    if(canJumpHelper(nums,i-1)) 
                    	return true; 
                }
            }
        }
        return false; 
    }
    
    public boolean canJump(int[] nums) {
        return canJumpHelper(nums, nums.length-1);
    }
//    public List<List<Integer>> levelOrder(TreeNode root) {
//        List<List<Integer>> all = new ArrayList<>();
//    	if (root == null) return all;
////    	int level = 0;
//    	Queue<TreeNode> myQ = new LinkedList<TreeNode>();	
//    	myQ.add(root);
//        all.add(new ArrayList<>(root.val));
//    	Queue<TreeNode> newQ = new LinkedList<TreeNode>();	
//    	List<Integer> list = new ArrayList<>();
//    	while (!myQ.isEmpty()) {
//    		while (!myQ.isEmpty()) {
//        		TreeNode curr = myQ.poll(); 
//        		if (curr.left != null)
//        			list.add(curr.left.val);
//        			newQ.add(curr.left);
//        		if(curr.right != null) {
//        			list.add(curr.right.val);
//        			newQ.add(curr.right);
//        		}
//        	}
//        	all.add(new ArrayList<>(list));
//        	myQ = newQ; 
//            newQ = new LinkedList<TreeNode>();
//        	list = new ArrayList<>();
//    	}
//    	return all;
//    }
    public List<List<Integer>> levelOrder(TreeNode root) {
       	List<List<Integer>> all = new LinkedList<List<Integer>>();
    	if (root == null) return all;
//    	int level = 0;
    	Queue<TreeNode> myQ = new LinkedList<TreeNode>();	
    	myQ.offer(root);
    		while (!myQ.isEmpty()) {
    			int levelNum = myQ.size();
    	    	List<Integer> list = new LinkedList<Integer>();
    	    	for (int i =0; i<levelNum; i++) {
            		TreeNode curr = myQ.poll(); 
    	    		if (curr.left != null)
            			myQ.add(curr.left);
            		if(curr.right != null) {
            			myQ.add(curr.right);
            		}
            		list.add(curr.val);
    	    	}
        		all.add(list);
        	}
    	return all;
    }
//    private int n;
//    private int m; 
//    public int numIslands(char[][] grid) {
//        int count = 0;
//        n = grid.length;
//        if (n==0) return 0;
//        m = grid[0].length;
//        for (int i =0; i<n; i++){
//            for (int j =0; i<m; j++){
//                if (grid[i][j]=='1'){
//                    DFSMarking(grid, i, j);
//                    count += 1; 
//                }
//            }
//        }
//        return count;
//    }
//    
//    private void DFSMarking(char[][] grid, int i, int j){
//        if (i<0 || j< 0 || i>=n || j>=m || grid[i][j] != '1') return ;
//        grid[i][j] = '0';
//        DFSMarking(grid, i+1, j);
//        DFSMarking(grid, i-1, j);
//        DFSMarking(grid, i, j+1);
//        DFSMarking(grid, i, j-1);
//    }
    
//    public boolean canPartition(int[] nums) {
//        int length = nums.length;
//        double sum = 0; 
//        for (int i =0 ; i<length; i++) {
//        	sum += nums[i];
//        }
//        double half = sum/2; 
//        if (half - (int)half != 0) return false; 
////        return twoGoals(nums, 0, (int)half, (int)half);
//        return DP1(nums, (int)half);
//    }
//    
//    private boolean DP1(int[] nums, int goal) {
//    	int m = nums.length;
//    	int n =goal*2; 
//    	boolean[][] dp = new boolean[m][n];
//    	for (int i = 0; i< n; i++) {
//    		dp[i][0] = true;
//    	}
//    	for (int j = 1; j<m; j++) {
//    		dp[0][j] = false;
//    	}
//    	
//    	for (int i = 1; i< m; i++) {
//    		for (int j=1; j<n; j++) {
//    			dp[i][j] = dp[i-1][j] || dp[i-1][j-nums[i-1]];
//
//    			}
//    		}
//    	return dp[m][n];
//    }
    public boolean canPartition(int[] nums) {
        int length = nums.length;
     int sum = 0;
    
    for (int num : nums) {
        sum += num;
    }
    
    if ((sum & 1) == 1) {
        return false;
    }
    sum /= 2;
    
        return DP1(nums, sum);
    }
    
    private boolean DP1(int[] nums, int goal) {
    	int n = nums.length+1;
    	int m =goal+1; 
    	boolean[][] dp = new boolean[n][m];
           for (int i = 0; i < dp.length; i++) {
        Arrays.fill(dp[i], false);
    }
    	for (int i = 0; i< n; i++) {
    		dp[i][0] = true;
    	}
    	for (int j = 1; j<m; j++) {
    		dp[0][j] = false;
    	}
    	
    	for (int i = 1; i< n; i++) {
    		for (int j=1; j<m; j++) {
    			dp[i][j] = dp[i-1][j];
    			if (j>=nums[i-1]) {
    				dp[i][j] = dp[i][j]|| dp[i-1][j-nums[i-1]];
    			}

    			}
    		}
    	return dp[nums.length][goal];
    }

    
//    private boolean twoGoals(int[] nums, int start, int goal1, int goal2) {
////    	while (start < nums.length) {
//    	if (goal1 < 0 || goal2 <0 || start >= nums.length) return false; 
//    		return twoGoals(nums, 1, goal1-nums[start], goal2) || twoGoals(nums, 1, goal1, goal2-nums[start]);
//    }
//    Efficient Solution 
//    public List<List<Integer>> levelOrder(TreeNode root) {
//        List<List<Integer>> levelList = new ArrayList<>();
//        helper(root, levelList, 0);
//        return levelList;
//    }
//    
//    private void helper(TreeNode node, List<List<Integer>> levelList, int level) {
//        if (node == null) return;
//        if (levelList.size() <= level) {
//            levelList.add(new ArrayList<Integer>());
//        }
//        levelList.get(level).add(node.val);
//        helper(node.left, levelList, level+1);
//        helper(node.right, levelList, level+1);
//    }
    
// 用Queue来解决的valideBST （还没有搞对）   
//    public boolean isValidBST(TreeNode root) {
//        Queue<TreeNode> myQ = new LinkedList<TreeNode>();
//        myQ.offer(root);
//        while (!myQ.isEmpty()) {
//        	int compare = (myQ.peek()).val; 
//        	TreeNode left =myQ.peek().left; 
//        	if (left!=null) {
//        		if(left.val >= compare) return false; 
//        		myQ.offer(left);
//        	}
//        	TreeNode right =myQ.peek().right; 
//        	if (right!=null) {
//        		if(right.val <=compare) return false;
//        		myQ.offer(right);
//        	}
//        	myQ.poll();
//        }
//        return true; 
//    }
//    
//    
//    public boolean isValidBST2(TreeNode root) {
//        if (root == null) return true; 
//        Queue<TreeNode> myQ = new LinkedList<TreeNode>();
//        myQ.offer(root);
//        while (!myQ.isEmpty()) {
//        	int Lcompare = (myQ.peek()).val; 
//            int Rcompare = (myQ.peek()).val; 
//        	TreeNode left =myQ.peek().left; 
//        	if (left!=null) {
//        		if(left.val >= Lcompare || left.val <= Rcompare) return false; 
//                Rcompare = Math.max(Rcompare,left.val);
//        		myQ.offer(left);
//        	}
//        	TreeNode right =myQ.peek().right; 
//        	if (right!=null) {
//        		if(right.val <=Rcompare || right.val >= Lcompare) return false;
//                Lcompare = Math.min(Rcompare,right.val);
//        		myQ.offer(right);
//        	}
//        	myQ.poll();
//        }
//        return true;    
//    }
//    
//    public boolean isValidBST3(TreeNode root) {
//        if (root == null) return true; 
//        int boundL = Integer.MAX_VALUE;
//        if (root.left != null) {
//             boundL = Math.min(boundL,root.left.val);
//            if (root.left.val>= root.val || !isValidL(root.left,boundL))
//                return false; 
//        }
//        
//                int boundR = Integer.MIN_VALUE;
//               if (root.right != null) {
//             boundR = Math.max(boundR,root.right.val);
//            if (root.right.val<= root.val || !isValidR(root.right,boundR))
//                return false; 
//        }
//        return true; 
//    }
//    private boolean isValidL(TreeNode root, int bound){
//        if (root == null) return true; 
//        while(root.left !=null){
//             if (root.left.val>= bound)
//                return false; 
//        }
//        while(root.right !=null){
//             if (root.right.val>= bound)
//                return false; 
//        }
//        return true; 
//    }
//        private boolean isValidR(TreeNode root, int bound){
//        if (root == null) return true; 
//        while(root.left !=null){
//             if (root.left.val<= bound)
//                return false; 
//        }
//        while(root.right !=null){
//             if (root.right.val<= bound)
//                return false; 
//        }
//        return true; 
//    }
        
        public List<Integer> inorderTraversal(TreeNode root) {
        	List<Integer> list = new ArrayList<>();
        	if (root == null) return list;
        	Stack<TreeNode> stack = new Stack<>();
        	while (root != null|| !stack.empty()) {
        		while (root != null) {
            		stack.push(root);
            		root = root.left;
            	}
            	root = stack.pop();
            	list.add(root.val);
            	root = root.right; 
        	}    	
        	return list;       	
        }
        
        public int kthSmallest(TreeNode root, int k) {
        	if (root == null) return (Integer) null;
        	Stack<TreeNode> stack = new Stack<>();
        	while (root != null|| !stack.empty()) {
        		while (root != null) {
            		stack.push(root);
            		root = root.left;
            	}
            	root = stack.pop();
            	k = k-1; 
            	if (k==0) return root.val;
            	root = root.right; 
        	}    	
        	return root.val;
        }
// This is my own implementation done on May 25th.        
//        public boolean isValidBST(TreeNode root) {
//            int k = 0; 
//            int value = Integer.MIN_VALUE; 
//           	if (root == null || (root.left==null && root.right == null)) return true;
//            	Stack<TreeNode> stack = new Stack<>();
//            	while (root != null|| !stack.empty()) {
//            		while (root != null) {
//                		stack.push(root);
//                		root = root.left;
//                	}
//                	root = stack.pop();
//                    if ((k==0 && root.val<value) || (k>0 && root.val <= value))
//                        return false;
//                    k = k+1; 
//                	value = root.val; 
//                	root = root.right; 
//            	}    	
//            	return true; 
//    }

// 借鉴的discussion soln      
//        public boolean isValidBST(TreeNode root) {
//            if (root == null) return true;
//            Stack<TreeNode> stack = new Stack<>();
//            TreeNode pre = null;
//            while (root != null || !stack.isEmpty()) {
//                while (root != null) {
//                    stack.push(root);
//                    root = root.left;
//                }
//                root = stack.pop();
//                if(pre != null && root.val <= pre.val) return false;
//                pre = root;
//                root = root.right;
//            }
//            return true;
//    }
        // This is more like recursion + hashmap (this soln doesn't cooperates hash)
        public TreeNode buildTree(int[] preorder, int[] inorder) {
        	return helper(0,0,inorder.length-1, preorder,inorder);
        }
        public TreeNode helper(int preStart, int inStart, int inEnd, int[] preorder, int[] inorder) {
        	if (preStart > preorder.length-1 || inStart > inEnd) {
        		return null;
        	}
        	TreeNode root = new TreeNode(preorder[preStart]);
        	int inIndex = 0;
        	for (int i= inStart; i<= inEnd; i++) {
        		if(inorder[i]==root.val) {
        			inIndex= i;
        		}
        	}
        	root.left = helper(preStart+1, inStart, inIndex-1, preorder, inorder);
        	root.right = helper(preStart + inIndex - inStart +1, inIndex+1, inEnd, preorder, inorder);
        	return root; 
        }
        
        // 207 Course Schedule BFS !! 
        public boolean canFinish(int numCourses, int[][] prerequisites) {
        	int[][] matrix = new int[numCourses][numCourses];
        	int[] indegree= new int[numCourses];
        	
        	for(int i=0; i< prerequisites.length; i++) {
        		int ready = prerequisites[i][0];
        		int pre = prerequisites[i][1];
        		if (matrix[pre][ready]==0)
        			indegree[ready]++;
        		matrix[pre][ready] = 1; 
        		}
        	int count = 0; 
        	Queue<Integer> queue= new LinkedList();
        	// First offer all the classes that do not require any prerequisites 
        	for(int i =0; i<indegree.length;i++) {
        		if(indegree[i]==0) queue.offer(i);
        	}
        	while(!queue.isEmpty()) {
        		int course = queue.poll();
        		count ++;
        		for(int i=0; i<numCourses; i++) {
        			if(matrix[course][i] !=0 ) {
        				if(--indegree[i]==0)
        					queue.offer(i);
        			}
        		}
        	}
        	return count == numCourses; 
        	
        }
        // Binary Tree Non-Recursive Solution. Very good intuition and visualization 
        public void flatten(TreeNode root) {
          // If emptry return null
        	TreeNode cur= root; 
        	while (cur != null) {
        		if(cur.left !=null) {
        			TreeNode last = cur.left;
        			while (last.right != null) last = last.right;
        			last.right = cur.right;
        			cur.right = cur.left;
        			cur.left = null;
        		}
        		cur = cur.right;
        	} 
        }
        
        // Time Complexity O(nlogn) 一开始想到的就是先sort 但其实还是hash比较实用
        public List<Integer> topKFrequent(int[] nums, int k) {
            List<Integer>[] bucket = new List[nums.length+1];
            Map<Integer, Integer> frequencyMap = new HashMap<Integer, Integer>();
            for(int n :nums) {
            	frequencyMap.put(n, frequencyMap.getOrDefault(n, 0)+1);
            }
            for(int key:frequencyMap.keySet()) {
            	int frequency = frequencyMap.get(key);
            	if (bucket[frequency]==null) {
            		bucket[frequency]=new ArrayList<>();
            	}
            	bucket[frequency].add(key);
            }
            List<Integer> res = new ArrayList<>();
            
            for(int pos = bucket.length-1; pos >=0 && res.size() <k; pos--) {
            	if(bucket[pos] != null) {
            		res.addAll(bucket[pos]);
            	}
            }
            return res; 
        }
        
        public int findKthLarget(int[] nums, int k) {
			return k;
        }
        
      	
 	   public int coinChange(int[] coins, int amount) {
 		   // dp[i] <=> the smallest # of coins to reach goal i is dp[i]
 	       int[] dp = new int[amount+1];
 	       dp[0]=0;

 	       for(int i = 1; i<dp.length;i++) {
 	    	   dp[i] = Integer.MAX_VALUE;
 	    	   for(int coin: coins) {
 	    		  if(i-coin >=0 && 1+dp[i-coin]>=0 && 1+dp[i-coin]<dp[i]) {
 	 	    		  dp[i] =1+dp[i-coin];
 	 	    	   }
 	    	   }  
 	       }
 	       return dp[amount];
 	    }
 	// TODO: This quite makes sense
// 	   private int coinChangeHelper(int[] coins, int rem, int[] count) {
// 		   return 0; 
// 	   }
// 	   
// 	    public String countAndSay(int n) {
//	    	StringBuilder curr=new StringBuilder("1");
//	    	StringBuilder prev;
//	    	int count;
//	    	char say;
//	        for (int i=1;i<n;i++){
//	        	prev=curr;
//	 	        curr=new StringBuilder();       
//	 	        count=1;
//	 	        say=prev.charAt(0);
//	 	        
//	 	        for (int j=1,len=prev.length();j<len;j++){
//	 	        	if (prev.charAt(j)!=say){
//	 	        		curr.append(count).append(say);
//	 	        		count=1;
//	 	        		say=prev.charAt(j);
//	 	        	}
//	 	        	else count++;
//	 	        }
//	 	        curr.append(count).append(say);
//	        }	       	        
//	        return curr.toString();
//        
//    }
 	   
 	   public int maxProfit(int[] prices) {
 	     int n = prices.length;
 	     int[] sell = new int[n];
 	     int[] buy = new int[n];
 	     sell[0] = 0;
 	     buy[0] = -prices[0];
 	     sell[1] = Math.max(prices[1]-prices[0],0);
 	     buy[1] = Math.max(-prices[1],buy[0]);
 	     for(int i = 2; i < n; i++) {
 	    	 sell[i] = Math.max(buy[i-1]+ prices[i], sell[i-1]);
 	 	     buy[i] = Math.max(sell[i-2]-prices[i], buy[i-1]); 
 	     }
 	     return sell[n-1];
 	    }
 	   
 	    public int threeSumClosest(int[] nums, int target) {
 	        int n = nums.length;
 	    	int result = nums[0] + nums[1] + nums[n-1];
 	        Arrays.sort(nums);
 	        for (int i = 0; i< n-2;i++) {
 	        	int start = i+1, end = n-1; 
 	 	        while(start<end) {
 	 	        	int sum = nums[i] + nums[start] + nums[end]; 
 	 	        	if(sum > target) {
 	 	        		end--;
 	 	        	}else {
 	 	        		start ++ ; 
 	 	        	}
 	 	        	if(Math.abs(sum - target) < Math.abs(result -target)) {
 	 	 	        	result = sum; 
 	 	 	        }
 	 	        }   
 	        }
 	        return result; 
 	    }
 	   
 	    public int[][] reconstructQueue(int[][] people){
 	    	if (people == null || people.length ==0 || people[0].length ==0)
 	    		return new int[0][0];
 	    	
 	        Arrays.sort(people, new Comparator<int[]>(){
 	        	public int compare(int[] a, int[] b) {
 	        		if (b[0] == a[0]) return a[1]-b[1];
 	        		return b[0]-a[0];
 	        	}
 	        });
 	        
 	        int n = people.length;
 	        ArrayList<int[]> tmp = new ArrayList<>();
 	        for(int i =0; i<n; i++) {
 	        	tmp.add(people[i][1],new int[] {people[i][0],people[i][1]});
 	        }
 	        
 	        int[][] res = new int[people.length][2];
 	        int i = 0;
 	        for(int[] k: tmp) {
 	        	res[i][0] = k[0];
 	        	res[i++][1] = k[1];
 	        }
 	        return res;
 	    }
 	    
// 	    public int countSubstrings(String s) {
// 	        return 0; 
// 	    }
 	    
 	    
 	    public String addBinary(String a, String b) {
 	        int carry = 0;
// 	        ArrayList<Integer> ans = new ArrayList<Integer>();
 	        StringBuilder sb = new StringBuilder();
 	        int i = 0;
 	        while(i< a.length() && i<b.length()) {
 	        	int digit = a.charAt(i)-'0' + b.charAt(i)-'0';
 	        	sb.append(digit%10); 
 	        	carry = digit/10; 
 	        }
 	        return sb.toString();
 	    }
 // TODO： 	    
// 	    public TreeNode lowestCommonAncestor(TreeNode root, TreeNode p, TreeNode q) {
// 	        
// 	    }
 	    
 	    public int findKthLargest(int[] nums, int k) {
 	    	return 0;
 	    	
 	    }
 	    
 	    public void merge(int[] nums1, int m, int[] nums2, int n) {
 	         helper( nums1, m, nums2, n, 0,0);
 	    }
 	    private int[] helper(int[] nums1, int m, int[] nums2, int n, int idx1, int idx2){
 	        if(m==0){
 	            for (int i =0; i< n; i++){
 	                nums1[i] = nums2[i];
 	            }
 	            return nums1; 
 	        }
 	        if(n == 0){
 	            return nums1; 
 	        }
 	        if(idx1>=nums1.length){
 	            for(int i= idx2; i< nums2.length; i++){
 	                nums1[m+idx2] = nums2[idx2];
 	            }
 	            return nums1; 
 	        }
 	        if(idx2>=nums2.length){
 	            return nums1; 
 	        }
 	        if(nums1[idx1]<= nums2[idx2]){
 	            return helper(nums1, m, nums2, n, idx1+1, idx2); 
 	        }else{
 	            nums1[idx1] = nums2[idx2];
 	            return helper(nums1, m, nums2, n, idx1, idx2+1); 
 	        }
 	    }
 
// 	    
// 	    private void shuffle(int a[]) {
// 	    	final Random random  = new Random();
// 	    	for (int ind = 1; ind < a.length ; ind ++) {
// 	    		final int r= random.nextInt(ind+1);
// 	    		exch(a,ind,r);
// 	    	}
// 	    }
 	    
 	    
 	    public int strStr(String haystack, String needle) {
 	        if(needle.length() == 0 ) return 0; 
 	        if(haystack.length()== 0) return -1;
 	        int start = 0;
 	        while(haystack.length() - start>=needle.length()){
 	            int idx = haystack.indexOf(needle.charAt(0), start);
 	            System.out.println("Start is "+ start);
 	            System.out.println("String is "+ haystack.substring(start)); 
 	            System.out.println("Idx is "+ idx); 
 	            int joy = idx; 
 	            if  (idx == -1 || idx < start) return -1;
 	            boolean test = true; 
 	            for(int i = 0; i< needle.length(); i++){
 	                if(idx+i > haystack.length()-1) return -1;
 	                if( haystack.charAt(idx+i) !=needle.charAt(i) ){
 	                    start = joy+1; 
 	                    System.out.println("Update start is "+ start);
 	                    test = false; 
 	                    System.out.println("Now" + i + " bool is " + test); 
 	                    break;
 	                }
 	                     
 	            }
 	            if (test)  return idx; 
 	        }
 	        return -1;

 	    }    
 	    
 	    public boolean isValidSudoku(char[][] board) {
 	    	for(int i = 0; i<9; i++) {
 	    		HashSet<Character> rows = new HashSet<Character>();
 	    		HashSet<Character> columns = new HashSet<Character>();
 	    		HashSet<Character> cube = new HashSet<Character>();

 	    			for(int j = 0; j<9; j++) {
 	    				char cur = board[i][j]; 
 	    				if( cur != '.' && ! rows.add(cur)) return false;
 	    				if( cur != '.' && ! columns.add(cur)) return false;
 	    				int RowIndex = 3*(i/3);
 	    				int ColIndex = 3*(i%3);
 	    				if(board[RowIndex + j/3][ColIndex + j%3] != '.' && !cube.add(board[RowIndex + j/3][ColIndex + j%3]))
 	    					return false; 		
 	    			}
 	    	}
			return true;	        
 	    }
 	    
 	    public int mySqrt(int x) {
 	        int tr = x/10;
 	        while(!((tr*tr== x) || (tr*tr<x &&(tr+1)*(tr+1)>x ))){
 	            if(tr*tr==x) {
 	                return tr; 
 	            }else if(tr*tr>x){
 	                if((tr-1)*(tr-1)<=x){
 	                    return tr-1;
 	                }else {
 	                    tr = tr-2;
 	                }
 	            }else{
 	                if((tr+1)*(tr+1)>x){
 	                    return tr;
 	                }else if((tr+1)*(tr+1)==x){
 	                    return tr+1; 
 	                }else{
 	                    tr += 2;
 	                }
 	            }
 	        }
 	        return tr;             
 	            
 	    }
 	    	    
 	   private boolean isPalindrome(String str) {
 		    int left = 0;
 		    int right = str.length() - 1;
 		    while (left <= right) {
 		        if (str.charAt(left++) !=  str.charAt(right--)) return false;
 		    }
 		    return true;
 		} 
 	   
 
 	   
// 	   // Check str[l,r] is a Palindrome (both sides inclusivce)
// 	   private boolean isPalindrome(String str, int l , int r) {
//		    int left = l;
//		    int right = r;
//		    while (left <= right) {
//		        if (str.charAt(left++) !=  str.charAt(right--)) return false;
//		    }
//		    return true;
//		} 
 	    
// 	    public boolean isPalindrome(String s) {
//            if(s== "") return true; 
//            if(s.length()==0 || s.length() ==1) return true;
// 	        s = s.toUpperCase(); 
// 	        int l = 0, r = s.length()-1;
// 	        while(l<r){
//// 	        	System.out.println("Current l index is: "+ l + " r index is "+ r);
// 	while(l<=s.length()-1 && (s.charAt(l) == ' ' || s.charAt(l) -'A' >25 || s.charAt(l)<48 || (s.charAt(l) <65 &&s.charAt(l)>57 ))){l++;}
// 	             	        	if(l == s.length()) return true;
//
// 	 	        	while(r>=0 && (s.charAt(r) == ' ' || s.charAt(r) -'A' >25 || s.charAt(r)<48 || (s.charAt(r) <65 &&s.charAt(r)>57) )){r--;}
// 	            if(s.charAt(l) != s.charAt(r)) return false;
// 	            l++;
// 	            r--;
// 	        }
// 	        return true; 
// 	    }
 	    
 	    public int[] anagramMappings(int[] A, int[] B) {
 	    	int[] rst = new int[A.length];
 	        HashMap<Integer,Integer> dic = new HashMap<>();
 	        for(int i = 0; i< B.length; i++) {
 	        	dic.put(B[i], i);
 	        }
 	        for(int i =0 ; i<A.length;i++) {
 	        	rst[i] = dic.get(A[i]);
 	        }
 	        
 	        return rst;
 	    }
 	    
 	
	    
	    public void printIntMatrix(int[][] in) {
 	    	for(int[] i:in) {
 	    		System.out.println(Arrays.toString(i));
 	    		System.out.println();
 	    	}
 	    }
	    
 	    public String reverseVowels(String s) {
 	        HashMap<Integer,Character> dic = new LinkedHashMap<>();
 	        for(int i =0 ; i<s.length() ; i++) {
 	        	if(s.charAt(i) == 'a' || s.charAt(i) == 'e'||s.charAt(i) == 'i'|| s.charAt(i) == 'o'|| s.charAt(i) == 'u'|| s.charAt(i) == 'A' || s.charAt(i) == 'E'||s.charAt(i) == 'I'|| s.charAt(i) == 'O'|| s.charAt(i) == 'U') {
 	        		dic.put(i, s.charAt(i));	        	
 	        	}
 	 
 	        }
     		System.out.println("Dic is "+dic);

 	        char[] ch = new char[dic.size()]; 
 	        int idx = ch.length-1; 
			for( int i :dic.keySet()) {
 	        	ch[idx] = dic.get(i);
 	        	idx--;
 	        }
			char[] rst = new char[s.length()];
	        int idxx = 0; 
			for( int i :dic.keySet()) {
 	        	rst[i] = ch[idxx];
 	        	idxx++;
 	        }
			for(int i =0 ;i<rst.length; i++) {
				if(rst[i]=='\u0000') {
					rst[i] = s.charAt(i);
				}
			}
			return new String(rst); 	        
 	    }
 	    
 	    
 	    public int firstUniqChar(String s) {
 	        HashMap<Character,Integer> dic = new LinkedHashMap<>();
 	        for(int i =0; i<s.length(); i++) {
 	        	if(!dic.containsKey(s.charAt(i))) {
 	        		dic.put(s.charAt(i), 1);
 	        	}else {
 	        		dic.put(s.charAt(i), dic.get(s.charAt(i))+1);
 	        	}
 	        }
// 	        System.out.println("dic is "+ dic);
 	        	char x = ' ';
// 	        	System.out.println(dic.keySet());
 	        for(Character c:dic.keySet()) {
// 	        	System.out.println("C is "+c);
// 	        	System.out.println("c is "+c);
 	        	if(dic.get(c)==1) {
 	        		x = c; 
 	        		break;
 	        	}
 	        }
// 	        System.out.println("x is "+ x);
            for(int i =0; i<s.length(); i++){
                if (s.charAt(i) == x)
                    return i; 
            }
 	        return -1;
 	    }
 	    
 	    
 	    public int guess(int n) {
 	    	if(n==1000) return 0;
 	    	else if(n>1000) return -1; 
 	    	else return 1; 

 	    }
 	    public int guessNumber(int n) {
 	        int gue = n/2;
 	        while(true){
 	            if(guess(gue) == 0) return gue;
 	            else if (guess(gue) > 0){
 	                if(gue+(n-gue)/2 > gue){
 	                    gue = gue+(n-gue)/2;
 	                }else{
 	                    gue +=1; 
 	                }
 	                
 	                guess(gue);
 	            }else{
 	                if(gue/2 <gue){
 	                    gue =gue/2;
 	                }else{
 	                    gue -=1;
 	                }
 	                guess(gue);
 	            } 
 	        }
 	    }

 	    public List<List<Integer>> palindromePairs(String[] words) {
 	        List<List<Integer>> rst = new ArrayList<>();
 	        for(int i =0 ; i< words.length; i++){
 	            for(int j = i+1; j<words.length ; j++){
 	                if(words[i].length() == 0 || words[j].length() == 0) {
 	                    if(words[i]== "" && words[j]== "") {
 	                        List<Integer> put = Arrays.asList(i,j);
 	                        rst.add(put);
 	                    }else if(words[i]== "" && isPalindrome(words[j])){
// 	                    	System.out.println("I'm here 1");
 	                        List<Integer> put = Arrays.asList(i,j);
 	                        if(!rst.contains(put)) {
 	 	                        rst.add(put);
 	                        }
 	                        List<Integer> put2 = Arrays.asList(j,i);
 	                        if(!rst.contains(put2)) {
 	 	                        rst.add(put2);
 	                        }
 	                    }
 	                    if(words[j]== "" && isPalindrome(words[i])){
// 	                    	System.out.println("I'm here 2");
 	                         List<Integer> put = Arrays.asList(j,i);
 	                         rst.add(put);
 	                        if(!rst.contains(put)) {
 	 	                        rst.add(put);
 	                        }
 	                        List<Integer> put2 = Arrays.asList(i,j);
 	                        if(!rst.contains(put2)) {
 	 	                        rst.add(put2);
 	                        }
 	                     }
 	                }else{
 	                	System.out.println("I'm here ");
 	                if(words[i].charAt(0)!=words[j].charAt(words[j].length()-1) && words[i].charAt(words[i].length()-1)!=words[j].charAt(0)){
 	                    continue;
 	                }
 	                if (words[i].charAt(0)==words[j].charAt(words[j].length()-1)) {
 	                    String test= words[i].concat(words[j]);
 	                    if(isPalindrome(test)){
 	                        List<Integer> put = Arrays.asList(i,j);
 	                        rst.add(put);
 	                    }
 	                }
 	                 
 	                if (words[i].charAt(words[i].length()-1)==words[j].charAt(0)) {
 	                    String test= words[j].concat(words[i]);
 	                    if(isPalindrome(test)){
 	                        List<Integer> put = Arrays.asList(j,i);
 	                        rst.add(put);
 	                    }
 	                }
 	                    }
 	            }
 	        }
 	        return rst;  
 	    }

 	    
// 	    public String removeDuplicateLetters(String s) {
//            if(s.length() == 0 || s.length() ==1) return s; 
// 	    	char[] chars = s.toCharArray();
//
// 	    	Arrays.sort(chars);
// 	    	System.out.println("chars is " + chars);
// 	    	HashMap<Character,Integer> dic= new LinkedHashMap<>(); 
// 	    	for(char c : chars) {
// 	    		dic.put(c, s.indexOf(c)); 
// 	    	}
// 	    	System.out.println(dic);
// 	    	Set<Character> mySet = dic.keySet();
// 	    	System.out.println("set is "+ mySet);
// 	    	int i =0;
// 	    	boolean check = true;
// 	    	for(; i<s.length(); i++) {
// 	    		System.out.println("i is "+i);
// 	    		System.out.println("chars[i] is "+chars[i]);
//
// 	    		for(char j : mySet ) {
// 	    			System.out.println("j is "+ j );
// 	    			if(j!=chars[i]){
// 	    				continue;
// 	    			}else if( s.indexOf((int)j, dic.get(chars[i]))!=-1) {
// 	    				continue;
// 	    			}else{
// 	    				check = false;
// 	    				System.out.println("Check "+ check);
// 	    				break;
// 	    			}
// 	    			
// 	    		}
// 	    		if(check) break;
// 	    	}
// 	    	
// 	    	StringBuilder rst =  new StringBuilder();
// 	    	rst.append(chars[i]);
// 	    	HashSet<Character> fine = new HashSet<>();
// 	    	fine.add(chars[i]);
// 	    	HashMap<Integer, Character> hope = new HashMap<>();
// 	    	for(char x:mySet) {
// 	    		if(x!=chars[i]) {
// 	    			hope.put(s.indexOf((int)x,i),x );
// 	    		}
// 	    	}
// 	    	
// 	    	for(int j = i+1; j<chars.length; j++) {
//// 	    		if(chars[j]!= chars[i]) {
// 	    		if(!fine.contains(chars[j])) {
// 	    			rst.append(chars[j]);
// 	    			fine.add(chars[j]);
// 	    		}
// 	    	}
//			return rst.toString();
// 	   
// 	    }
 	    
 	    public String removeDuplicateLetters(String s) {
 	    	HashMap<Character, Integer> dic = new HashMap<>();
 	    	for(int i= 0; i< s.length(); i++) {
 	    		if(dic.containsKey(s.charAt(i))) {
 	    			dic.put(s.charAt(i), dic.get(s.charAt(i))+1);
 	    		}else {
 	    			dic.put(s.charAt(i), 1);
 	    		}
 	    	}
 	    	int pos = 0; 
 	    	for(int i = 0; i<s.length() ; i++) {
 	    		if(s.charAt(i)<s.charAt(pos)) {
 	    			pos = i; 
 	    		}
 	    		dic.put(s.charAt(i),dic.get(s.charAt(i))-1);
 	    		if(dic.get(s.charAt(i)) ==0) {
 	    			break;
 	    		}
 	    	}
			return s.length() == 0? "":s.charAt(pos) + removeDuplicateLetters(s.substring(pos+1).replaceAll(""+s.charAt(pos),""));
 	    	
 	    }
 	    public int minMeetingRooms(Interval[] intervals) {
 	        /** Sort by the start time. 
 	        *
 	        */
// 	    	if(intervals.length == 0) return 0; 
// 	    	Arrays.sort(intervals, (Interval a, Interval b) -> a.start - b.start);
// 	    	int count = 1;
// 	    	ArrayList<Integer> time = new ArrayList<>();
// 	    	for(int i =0; i<intervals.length; i++) {
// 	    		if(intervals[i].end>intervals[i+1].start) {
// 	    			count ++; 
// 	    			time.add(intervals[i].end);
// 	    		}else {
// 	    			
// 	    		}
// 	    	}
//			return m;
 	       if (intervals == null || intervals.length == 0)
 	          return 0;
 	         
 	      // Sort the intervals by start time
 	      Arrays.sort(intervals, (Interval a, Interval b) -> a.start - b.start);
 	      
 	     // Use a min heap to track the minimum end time of merged intervals
 	     PriorityQueue<Interval> heap = new PriorityQueue<Interval>(intervals.length,(Interval a, Interval b) -> a.end - b.end);
 	     
 	     heap.offer(intervals[0]);
 	     
 	     for(int i = 1; i<intervals.length; i++) {
 	    	 Interval interval = heap.peek();
 	    	 
 	    	 if(intervals[i].start >= interval.end) {
 	    		 interval.end = intervals[i].end;
 	    	 }else {
 	    		 heap.offer(intervals[i]);
 	    	 }
 	     }
			return heap.size();
 	    }
 	    
// 	    public int ladderLength(String beginWord, String endWord, List<String> wordList) {
// 	    	String now = beginWord;
//            String prev = beginWord; 
//	    	int min = Integer.MAX_VALUE;
// 	    	while(now != endWord) {
// 	    		System.out.println("Prev is "+ prev);
// 	    		System.out.println("Now is "+ now);
// 	    		System.out.println("min is "+ min);
//
// 	    		// cur stores the possible next-stage string from now string 
// 	    		List<String> cur = new ArrayList<>();
//	    		char[] char1 = now.toCharArray(); 
//	    		Arrays.sort(char1);
//	    		System.out.println("wordList is "+ wordList);
// 	    		for(String s: wordList) {
// 	    			System.out.println("s is " +s);
// 	    			if(now.length()!= s.length() || now == s || s == prev) continue; 
// 	    			int count = 0;
// 		    		char[] char2 = s.toCharArray(); 
// 		    		Arrays.sort(char2);
//
// 	    			for(int i = 0; i< s.length(); i++) {
// 	    				if(char1[i] != char2[i])
// 	    					count++; 
// 	    			}
// 	    			if(count == 1) cur.add(s);
// 	    			System.out.println("Cur is " + cur);
// 	 	    	}
//	    		prev = now; 
// 	    		if(cur.isEmpty()) return 0; 
// 	    		for(String str: cur) {
// 	    			int recur = ladderLength(str, endWord,wordList) ; 
// 	    			if(recur == 0) return 0; 
// 	    			if(recur +1 < min) {
// 	    				min = recur +1; 
// 	    			}
// 	    		}
// 	    		return min; 
// 	    	}
// 	    	
//			return min;
// 	        
// 	    }
 	    
 	    
 	   public int ladderLength(String start, String end, List<String> lst) {
 		  // Use queue to help BFS
 		  Queue<String> queue = new LinkedList<String>();
 		  queue.add(start);
 		  queue.add(null);
 		  
 		  // Mark visited word
 		  Set<String> visited = new HashSet<String>();
 		  visited.add(start);
 		  
 		  int level = 1;
 		  
 		  while (!queue.isEmpty()) {
 		    String str = queue.poll();
		    System.out.println("The polled str is "+ str);

 		    if (str != null) {
 		    	System.out.println("I'm here");
 		    	
 		      // Modify str's each character (so word distance is 1)
 		      for (int i = 0; i < str.length(); i++) {
 		        char[] chars = str.toCharArray();
 		        
 		        for (char c = 'a'; c <= 'z'; c++) {
 		          chars[i] = c;
 		          
 		          String word = new String(chars);
 		          
 		          // Found the end word
 		          if (word.equals(end)) return level + 1;
 		          
 		          // Put it to the queue
 		          if (lst.contains(word) && !visited.contains(word)) {
 		            queue.add(word);
 	 		    	System.out.println("The added word is "+ word);

 		            visited.add(word);
 		          }
 		        }
 		      }
 		    } else {
 		      level++;
 		      System.out.println("level is "+ level);
 		      if (!queue.isEmpty()) { 
 	 		      System.out.println("Empty again");
 		        queue.add(null);
 		      }
 		    }
 		  }
 		  
 		  return 0;
 		}
 	/**
 	 *    
 	 * @param args
 	 */
 	   
 	   
 	  public int lengthLongestPath(String input) {
   	Deque<Integer> stack = new ArrayDeque<>();
   	stack.push(0);
   	int maxLen = 0; 
   	for(String s:input.split("\n")){
   		System.out.println("s is "+ s);
   		int lev = s.lastIndexOf("\t")+1;   // number of "\t"
   		System.out.println("Lev is "+lev);
//   		System.out.println("Level is "+ lev);
   		while(lev+1 < stack.size()) stack.pop();  // find parent
   		int len = stack.peek() + s.length() - lev +1;  // remove "/t", add"/"
   		System.out.println("Before push "+stack);

   		stack.push(len);
   		System.out.println("After push "+stack);
   		// check if it is file
   		if(s.contains(".")) maxLen = Math.max(maxLen, len-1);
   	}
   	return maxLen; 
}
 	  
 	  
 	  
 	    public List<String> summaryRanges(int[] nums) {
 	        Arrays.sort(nums);
 	        List<String> rst = new ArrayList<>(); 
 	        if(nums.length ==0) return rst; 
 	        if(nums.length == 1) {
 	            rst.add(""+nums[0]);
 	            return rst;
 	        } 
 	       if(nums.length == 2) {
 	    	  if(nums[1] == nums[0]+1){
  	             rst.add("" + nums[0] + "->" + nums[1]);
  	        }else{
  	            rst.add(""+nums[nums.length -2]);
  	            rst.add(""+nums[nums.length -1]);
  	        }
	            return rst;
	        }
 	        int lastIndex = 0; // keep track of the start position of the new array 
 	        int i = 0; 
 	        while(i<nums.length-1){
 	        	if(i == nums.length -2){
 	                break; 
 	 
 	                // else if(i== lastIndex) {
 	                // rst.add(""+nums[i]);
 	             }
 	            if(nums[i+1] != nums[i]+1){
 	                rst.add(""+nums[i]);
// 	                System.out.println("I am here "+rst);
 	                lastIndex = i+1;
 	                i++;
 	                continue;
 	            }else {
 	            	while(i< nums.length -1 && nums[i+1] == nums[i]+1){
 	 	                i++; 
 	 	            }
// 	            	System.out.println("Before break");

 	            	if(i== nums.length -1) {
// 	 	            	System.out.println(" breaking");

 	            		break;
 	            	} 
 	            	else{
// 	 	            	System.out.println("After break");
 	                rst.add("" + nums[lastIndex] + "->" + nums[i]);
                     lastIndex = i+1; 
                     i++;
 	            	}
 	            }            
 	        }
 	        if(nums[nums.length -1] == nums[nums.length -2]+1){
 	             rst.add("" + nums[lastIndex] + "->" + nums[nums.length -1]);
 	        }else if(lastIndex == nums.length -2){
 	        	rst.add("" + nums[nums.length-2]);
 	            rst.add(""+nums[nums.length -1]);
 	        }
 	        else{
 	        	rst.add("" +  nums[lastIndex] + "->" + nums[nums.length-2]);
 	            rst.add(""+nums[nums.length -1]);
 	        }
 	        return rst;
 	    }
 	    
 	    
 	    // TODO: Need recovery 
 	   public String minWindow(String s, String t) {
 		   String rst = ""; 
 		   ArrayList<Character> bCheck = new ArrayList<>(); 
 		   for(int i = 0; i< t.length(); i++){
 			   if(!bCheck.contains(t.charAt(i))) {
 				  bCheck.add(t.charAt(i));
 			   }
 		   }
 		   
 		   int num = bCheck.size();
 		   ArrayList<Integer> idxDic = new ArrayList<>(); 
 		   bCheck.sort((Character a, Character b )-> a-b);
 		   System.out.println("bCheck is "+bCheck);

 		   for(int i = 0; i< s.length(); i++){
 			   if(bCheck.contains(s.charAt(i))) {
 				   idxDic.add(i);
 			   }
 		   }
 		   System.out.println("idxDic is "+idxDic);

 		  int j;
 		   for(int i : idxDic) {
 			   if(i<idxDic.size()-num) {
 				   System.out.println("num is "+ num);
 				   ArrayList<Character> really = new ArrayList<>();
 				   really.add(s.charAt(i));
 				   for( j = i+1; j<num; j++) {
 	 				   System.out.println("j is "+ j);
 	 				   System.out.println("get j is "+ idxDic.get(j));

 					   char wow = s.charAt(idxDic.get(j)); 
 					   if(!really.contains(wow)) {
 	 					   really.add(wow);
 					   }
 				   }
 	 		   really.sort((Character a, Character b )-> a-b);
 	 		   System.out.println("really is "+ really);
 			   if(really.equals(bCheck)) {
 				   System.out.println("Here");
 				   int curLength = idxDic.get(j) - i; 
 				   System.out.println("curLength "+ curLength);
 				   if(rst.length() == 0 || ( rst.length() != 0 && curLength < rst.length())) {
 					   rst = s.substring(i, idxDic.get(j));
 				   }
 			   }
 			   }
 		   }
 	        return rst; 
 	    }
 	    
 	    
// 	    public String minWindow(String s, String t) {
// 	        if(s == null || t == null)
// 	            return "";
// 	        
// 	        int[] map = new int[256];
// 	        int count = t.length();
// 	        for(char c : t.toCharArray()){
// 	            map[c]++;
// 	        }        
// 	        
// 	        int maxLength = Integer.MAX_VALUE;
// 	        int head = 0;
// 	        int start = 0, end = 0;
// 	        while(end < s.length()){
// 	            if(map[s.charAt(end++)]-- > 0){
// 	                count--;
// 	            }
// 	            
// 	            while(count == 0){
// 	                int len = end - start;
// 	                if(len < maxLength){
// 	                    maxLength = len;
// 	                    head = start;
// 	                }
// 	                
// 	                if(count == 0 && map[s.charAt(start++)]++ == 0)
// 	                    count++;
// 	                
// 	            }
// 	        }
// 	        
// 	        return maxLength == Integer.MAX_VALUE ? "" : s.substring(head, head + maxLength);
// 	    }
 	
//	        int lastEnd =0; 

// 	    public List<String> wordBreak(String s, List<String> wordDict) {
// 	    	List<String> rst= new ArrayList<>();  
// 	    	if(s.length() ==0 ) return rst;
// 	    	if(wordDict.contains(s)) rst.add(s);
// 	        HashSet<String> dic = new HashSet<>();
// 	        for(String str: wordDict){
// 	            dic.add(str);
// 	        }
//// 	        System.out.println("dic is "+ dic);
// 	        for(int i =1; i<=s.length(); i++){
//	            	System.out.println("i is "+i);
// 	            if(dic.contains(s.substring(0, i))){
// 	            	System.out.println("Substring " + s.substring(0, i));
//// 	            	lastEnd= i;
//// 	            	System.out.println("lastEnd " + lastEnd);
// 	            	
// 	            	List<String> last = wordBreak(s.substring(i),wordDict); 
//// 	            	System.out.println("last "+last );
//// 	            	System.out.println(s.substring(i));
// 	            	if(i == s.length() || last.size()!= 0) {
//// 	 	            	System.out.println("Here1");
// 	 	            	if(i== s.length()) { 
// 	 	            		last.add(s.substring(0,i) + " ");
// 	 	            		return last;
// 	 	            	}
// 	            		for(String each: last ) {
//// 	 	 	            	System.out.println("each "+each);
//
// 	 	            		each = s.substring(0,i) + " " +  each; 
// 	 	            	}
// 	 	                rst = last ; 
//
// 	            	}
// 	            	
// 	            }
// 	        }
// 	        for(String each:rst) {
// 	        	each.trim();
// 	        }
//			return rst;
// 	    }
 	    
	    HashMap<String, List<String>>  my = new HashMap<>();
	    public List<String> wordBreak(String s, List<String> wordDict) {
	    	if(my.containsKey(s)) return my.get(s);
 	    	List<String> rst= new ArrayList<>();  
 	        HashSet<String> dic = new HashSet<>();
 	        for(String str: wordDict){
 	            dic.add(str);
 	        }
 	        for(int i =1; i<=s.length(); i++){

 	        	if( i == s.length() && dic.contains(s.substring(0, i)) ) {
 	        		rst.add(s.substring(0,i));
 	        	}else if(dic.contains(s.substring(0, i))){
 	            	List<String> last = new ArrayList<String>();
 	            	last.addAll(wordBreak(s.substring(i),wordDict)); 
 	            	if( last.size()!= 0) {

 	            		for(int x =0 ; x<last.size() ; x++) {
 	            			String current = s.substring(0,i) + " " + last.get(x);
 	            			last.set(x, current);
 	 	            	}
 	 	                rst.addAll(last)  ; 
 	 	                my.put(s, rst);
 	            	} 	            	
 	            }
 	        }
 	       
			return rst;
 	    }
	    
	    
	    public int[] dailyTemperatures(int[] temperatures) {
	    	int[] rst = new int[temperatures.length];  
	    	Stack<Integer> fun = new Stack<>();
	    	if(temperatures.length== 0) {
	    		return rst;
	    	}
	    	fun.push(0);

	    	for(int i = 0; i< temperatures.length; i++) {
	    		if(fun.isEmpty() || temperatures[i]<=temperatures[fun.peek()]) {
	    		}else {
	    			while(!fun.isEmpty() && temperatures[fun.peek()]<temperatures[i]) {
	    				rst[fun.peek()] = i-fun.peek();
                        fun.pop();
	    			}
	    		}
    			fun.push(i);
	    	}
	    	rst[temperatures.length-1] = 0; 
			return rst;
	        
	    }
//	    public List<String> wordBreak(String s, List<String> wordDict) {
// 	    	List<String> rst= new ArrayList<>();  
// 	        HashSet<String> dic = new HashSet<>();
// 	        for(String str: wordDict){
// 	            dic.add(str);
// 	        }
// 	        for(int i =1; i<=s.length(); i++){
//	            	System.out.println("i "+ i);
//
// 	        	if( i == s.length() && dic.contains(s.substring(0, i)) ) {
// 	        		rst.add(s.substring(0,i));
// 	        	}else if(dic.contains(s.substring(0, i))){
// 	        		System.out.println("Hi");
// 	            	List<String> last = new ArrayList<String>();
// 	            	last.addAll(wordBreak(s.substring(i),wordDict)); 
//// 	            	System.out.println("last "+ last);
// 	            	if( last.size()!= 0) {
//
// 	            		for(int x =0 ; x<last.size() ; x++) {
//// 	            			System.out.println("new i is "+i);
//// 	            			System.out.println("current s   "+ s);
// 	            			String current = s.substring(0,i) + " " + last.get(x);
// 	            			last.set(x, current);	            			
// 	 	            	}
//// 	 	            	System.out.println("last' "+ last);
// 	 	                rst.addAll(last)  ; 
// 	            	} 	            	
// 	            }
// 	        }
// 	       
//			return rst;
// 	    }
 	   
//
// 	  public List<String> wordBreak(String s, List<String> wordDict) {
// 		  return DFS(s, wordDict, new HashMap<String, LinkedList<String>>());
// 	  }
// 	  
// 	  private List<String> DFS(String s, List<String> wordDict, HashMap<String, LinkedList<String>> map){
//		  if(map.containsKey(s))
//			  return map.get(s);
//		  
//		  LinkedList<String> res = new LinkedList<String>();
//		  if(s.length() ==0) {
//			  res.add("");
//			  return res; 
//		  }
//		  for(String word: wordDict) {
//			  if(s.startsWith(word)) {
//				  List<String> sublist = DFS(s.substring(word.length()), wordDict, map);
//				  for(String sub: sublist)
//					  res.add(word + (sub.isEmpty()? "" : " ") +sub);
//			  }
//		  }
//		  map.put(s, res);
// 		  return res;  
// 	  }
//	    
	    
	    /** 224 Basic Calculator 
	     * things to pay attention 
	     * 1. pre-process: first get rid of blank space
	     * 2. priority: () > * > +- 
	     * 3. - need to notice 
	     * 4. Stack is the intuitive structure 
	     * 5. when can we pop values from the stack 
	     * eg. 2-1*3 we cannot get 2-1 and then *3 
	     * (what about the returned value is negative/ out of bound)
	    */ 
	    public int calculate(String s) {
	    	Stack<Integer> rst = new Stack<>();
	    	int result = 0; // accumulated value 
	    	int number =0; // current value 
	    	int sign = 1; // default sign is 1; this simplifies the - 
	    	for(int i =0; i<s.length(); i++) {
	    		char c = s.charAt(i);
	    		if(Character.isDigit(c)){
	    			number = 10* number + (int)(c-'0');
	    		}else if(c=='+') {
	    			result += sign* number;
	    			number = 0; // update number 
	    			sign  = 1;
	    		}else if(c=='-') {
	    			result += sign * number; 
	    			number = 0;
	    			sign = -1; 
	    		}else if (c== '(') {
	    			rst.push(result); // that's smart 
	    			rst.push(sign); 
	    			sign = 1;
	    			result = 0;
	    		}else if(c==')') {
	    			result += sign * number; // this sign is the current sign inside the ()
	    			number = 0; // always remember to update number once it got used 
	    			result *= rst.pop(); // smart design 
	    			result += rst.pop();
	    		}};
//    				break;
//	    		case : rst.push('+');
//    				break;
//	    		case '-': rst.push('-');
//    				break;
//    			default: if(rst.peek() == '+') {
//    				rst.pop();
//    				int now = (int)s.charAt(i)+ (int)rst.pop(); 
//    				rst.push((char)now);
//    			}else if(rst.peek() == '-') {
//    				rst.pop();
//    				int now = (int)rst.pop()- (int)s.charAt(i); 
//    				rst.push((char)now);
//    			}else {
//    				rst.push(s.charAt(i));
//    			}
//    				break;
//	    		}
//	    	}
	    	if(number !=0) result += sign * number; // last update 
			return result;	        
	    }
	    
	    
//	    public int maxCoins(int[] nums) {
//	    	int max = Integer.MIN_VALUE;
//	    	HashMap<List<Integer>,Integer> dic = new HashMap<>(); 
//	    	dic.put(new ArrayList<Integer>(), 0); 
//	    	if(dic.containsKey(Arrays.asList(nums))) return dic.get(Arrays.asList(nums));
//	    	if(nums.length==0) {
//	    		return 0;
//	    	}
//	    	if(nums.length == 1) {
//	    		dic.put(Arrays.stream(nums).boxed().collect(Collectors.toList()),nums[0]);
//	    		return nums[0];
//	    	}
//	    	List<Integer> arr = Arrays.stream(nums).boxed().collect(Collectors.toList());
//	    	for(int i =0; i<nums.length;i++) {
//	    		int rst=0;
//	    		if(i==0) {
//	    			rst += nums[0]*nums[1];
//	    		}else if(i==nums.length-1) {
//	    			rst += nums[nums.length-1]*nums[nums.length-2];
//	    		}else {
//	    			rst += nums[i-1]*nums[i]*nums[i+1];
//	    		}
//	    		int save = nums[i];
//	    		arr.remove(i); 
//	    		rst += maxCoins(toArray()));
//	    		if() { 
//	    			
//	    		}
//	    		arr.add(i, save);
//	    	}
//	    	
//			return m;
//	        
//	    }
	    public int maxCoins(int[] nums) {
	    	ArrayList<Integer> conv = new ArrayList<>();
	    	for(int i:nums) {
	    		conv.add(i);
	    	}
	    	return maxCoins(conv);
	    }
	    
	    public int maxCoins(ArrayList<Integer> nums) {
	    	int max = Integer.MIN_VALUE;
	    	HashMap<List<Integer>,Integer> dic = new HashMap<>(); 
	    	dic.put(new ArrayList<Integer>(), 0); 
	    	if(dic.containsKey(nums)) return dic.get(nums);
	    	if(nums.size()==0) {
	    		return 0;
	    	}
	    	if(nums.size() == 1) {
	    		dic.put(nums,nums.get(0));
	    		return nums.get(0);
	    	}
	    	for(int i =0; i<nums.size();i++) {
	    		int rst=0;
	    		if(i==0) {
	    			rst += nums.get(0)*nums.get(1);
	    		}else if(i==nums.size()-1) {
	    			rst += nums.get(nums.size()-1)*nums.get(nums.size()-2);
	    		}else {
	    			rst += nums.get(i-1)*nums.get(i)*nums.get(i+1);
	    		}
	    		int save = nums.get(i);
	    		nums.remove(i); 
	    		rst += maxCoins(nums);
	    		if(rst>max) { 
	    			max = rst;
	    		}
	    		dic.put(nums, rst);
	    		nums.add(i, save);
	    	}	    	
			return max;	        
	    }
	   
	    
 	    public boolean canWinNim(int n) {
 	    	boolean[] dp=new boolean[n+1];
 	    	if(n<=3) return true;
 	    	dp[0] = true;
 	    	dp[1] =true;
 	    	dp[2] =true;
 	    	dp[3] =true;
 	    	for (int i = 4; i<n+1; i++) {
 	    		dp[i] = (!dp[i-1] || !dp[i-2] || !dp[i-3]);
 	    	}
 	    	return (dp[n]);   	
 	    }
 	    
 	    /**
 	     * The important thing is to consider all edge cases, 
 	     * including: negative integer, possible overflow, etc.
 	     * Use HashMap to store a remainder and its associated index while doing the division
 	     * so that whenever a same remainder comes up, we know there is a repeating 
 	     * fractional part.
 	     * 这道题真的超级考验细心 真的非常地厉害了
 	     * @param numerator
 	     * @param denominator
 	     * @return
 	     */
 	    
 	    public String fractionToDecimal(int numerator, int denominator) {
 	    	if(numerator == 0) {
 	    		return "0";
 	    	}
 	    	StringBuilder rst = new StringBuilder();
 	    	rst.append(((numerator>0) ^(denominator>0)) ? "-" : "");
 	    	long num = Math.abs((long)numerator);
 	    	long den = Math.abs((long)denominator);
 	    	
 	    	// integer part
 	    	rst.append(num/den);
 	    	num %=den;
 	    	if(num==0) {
 	    		return rst.toString();
 	    	}   	
 	    	// fraction part
 	    	rst.append('.');
 	    	HashMap<Long, Integer> map = new HashMap<>();
 	    	map.put(num, rst.length());
 	    	while(num!=0) {
 	    		num*=10;
 	    		rst.append(num/den);
 	    		num %= den;
 	    		if(map.containsKey(num)) {
 	    			int index = map.get(num);
 	    			rst.insert(index, '(');
 	    			rst.append(')');
 	    			break;
 	    		}else {
 	    			map.put(num, rst.length());
 	    		}
 	    	}
	    		return rst.toString();
 	    }
 	    
 	    /**
 	     * DFS. Once passed a pt, mark it as 0.
 	     * @param board
 	     * @return
 	     */
 	   public int countBattleships(char[][] board) {
 		   int count = 0; 
 		   int row = board.length;
 		   if(row == 0) return 0;
 		   int col = board[0].length;
 		   if(col ==0) return 0;
 		   for(int i = 0; i<row; i++) {
 			   for(int j = 0; j< col; j++) {
 				   if(board[i][j] == 'X') {
 					   int c = j;
 					   while(c+1 <col && board[i][c+1] == 'X') {
 						   board[i][c+1] = '.';
 						   c++;
 					   }
 					   int r = i; 
 					   while(r+1<row && board[r+1][j]=='X') {
 						   board[r+1][j] = '.';
 						   r++;
 					   }
 					   count ++;
 					   System.out.println("i "+i);
 					   System.out.println("j "+j);
 					   System.out.println("count "+count);
 				   }				   
 			   }
 		   }
 		   return count;
	    }
 	    
 	    
 	    public String convert(String s, int numRows) {
 	        int mod = 2*numRows -2; 
 	        StringBuilder rst = new StringBuilder();
 	        int end =numRows/2;
 	        if(numRows %2 ==0) {
 	        	end -=1 ;
 	        }
 	        int easy = 0;
 	        while(easy<s.length()) {
	        		rst.append(s.charAt(easy));
	        		easy+=mod;
 	        }
 	        for(int m = 1; m<= end; m++ ) {
 	        	int c = m; 
 	        	int pow = 0;
 	        	while(c<s.length()) {
	 	        	rst.append(s.charAt(c));
 	 	        	if(pow%2==0) {
 	 	        		c+=mod-2*m;
 	 	        	}else {
 	 	        		c+=2*m;
 	 	        	}
	 	        		pow++;
 	        	}
 	        }
 	        if(numRows%2==1) {
 	        	int extra = numRows/2;
 	        	int start = extra;
 	        	while(start<s.length()) {
	        		rst.append(s.charAt(start));
	        		start+=mod;
 	        }
 	        }
 	        return rst.toString();
 	    }
 	    
 	    
 	    public int maxProfit1(int[] prices) {
 	        int end = prices.length;
 	        int[] dp = new int[end];
 	        return helper(prices, 0, dp);
 	    }
 	    
 	   private int helper(int[] prices, int start, int[] dp){
 	        int end = prices.length;
 	        // store the maximum profit we get from day start to the end 
 	        if(end-start == 2 && prices[end-1]>prices[start]){
 	            dp[0] = prices[end-1]-prices[start];
 	            return prices[end-1]-prices[start];
 	        }else if(end-start <= 2){
 	            return 0;
 	        }else{
 	            dp[end-1] = 0;
 	            for (int i = end-2; i>=start; i--){
 	                for(int j = end-1; j>=i; j--){
 	                	int curr = 0;
 	                	if(j==end-1) {
 	                		 curr = prices[j] - prices[i];
 	                	}else {
 	                		curr  = prices[j] - prices[i] + dp[j+1];
 	                	}
 	                    if(curr > dp[i]){
 	                        dp[i] = curr;
 	                    }
 	                }
 	            }
 	        }
 	        return dp[0];
 	    }
 	   
	   public List<String> fizzBuzz(int n) {
		   List<String> rst= new ArrayList<>(); 
	        for(int i=0; i<n; i++) {
	        	if(i%3==0 && i%5 == 0) {
	        		rst.add("FizzBuzz");
	        	}else if(i%3==0) {
	        		rst.add("Fizz");
	        	}else if (i%5==0){
	        		rst.add("Buzz");
	        	}else {
	        		rst.add(Integer.toString(i));
	        	}
	        }
	        return rst; 
	    }
	   
	    public String reverseWords(String s) {
	    	String[] words = s.split("\\s+");
	    	System.out.println("words "+words);
	    	StringBuilder rst = new StringBuilder();
	    	for(int i = words.length-1; i>=0; i--) {
		    	System.out.println("i "+i);
	    		rst.append(words[i]).append(" ");
		    	System.out.println("rst "+rst);
	    	}
			return rst.toString().trim();
	        
	    }
	
	    public int removeElement(int[] nums, int val) {
	    	/** Start is the start location of to change
	    	 * i is the start pioneer location to find non-val 
	    	 * substitutes for start if necessary
	    	 * 
	    	 */
	        int start = 0;
	        int i = 0 ; 
	        int song = nums.length-1;
	        int len =nums.length ;
	        while(i<len){
	        	if(nums[i]!=val) {
	        		start ++;
	        		i++;
	        		continue;
	        	}
                song--;
	            while(i<len && nums[i]==val){                
	                i++;
	            }
                
                if(i==len){
                    if(nums[i-1] == val) return song-1;
                    else return song;
                } 
	            nums[start] = nums[i];

	            start++;
	            i++;
	        }
	        return song;
	    }
	    
	    public ListNode oddEvenList(ListNode head) {
	    	if(head == null || head.next == null) return head; 
	    	ListNode evenStart = head.next;
	    	ListNode cur = head; 
    		ListNode tmp = cur.next.next;
	    	while(tmp != null && tmp.next != null) {
	    		cur.next.next = tmp.next;
	    		cur.next = tmp;
	    		cur = tmp;
	    		tmp = tmp.next.next;
	    	}
	    	cur.next.next= null;
	    	if(tmp == null) {
	    		cur.next = evenStart; 
	    	}else {
	    		cur.next = tmp;
	    		tmp.next = evenStart;
	    	}  	
	    	
			return head;        
	    }
	    
	    
	    public void printLinkedList (ListNode head) {
	    	if(head == null) System.out.println("");
	    	while(head != null) {
	    		System.out.print(head.val + " ");
	    		head = head.next;
	    	}
	    }
	    
	    
//	    public ListNode insertionSortList(ListNode head) {
//	        ListNode rst = head;
//	      
//	        while(head!= null){
//	             ListNode put = head;
//	             head = head.next;
//	            
//	            ListNode cur= rst;
//	            while(cur.next!= null && put.val>cur.next.val){
//	                cur = cur.next; 
//	            }
//	            if(cur.next != null){
//	                put.next = cur.next;
//	                cur.next =put; 
//	            }else{
//	                if(put.val<cur.val){
//	                    put.next = cur;
//	                    cur = put;
//	                }else{
//	                    cur.next = put; 
//	                }
//	            }
//	        }
//	        
//	        return rst; 
//	    }
	    
//	    public ListNode insertionSortList(ListNode head) {
//	        if(head == null || head.next == null) return head;
//	        ListNode rst = head;
//	        head = head.next;
//	        rst.next = null; 
//	        
//	        while(head!= null){
//	            ListNode put = head;
//	            
//	            ListNode cur= rst;
//	            while(cur.next!= null && put.val>cur.next.val){
//	                cur = cur.next; 
//	            }
//	            if(cur.next != null){
//	                put.next = cur.next;
//	                cur.next =put; 
//	                if(rst.next == null){
//	                    rst = cur; 
//	                }
//	            }else{
//	                if(put.val<cur.val){
//	                    put.next = cur;
//	                    cur = put;
//	                }else{
//	                    cur.next = put; 
//	                    put.next = null; 
//	                }
//	            }
//	            head = head.next;
//	        }       
//	        return rst; 
//	    }
	    
	    public ListNode insertionSortList(ListNode head) {
	        if(head == null || head.next == null) return head;
	        ListNode rst = head;
	        
	        head = head.next;
            rst.next = null; 
            int count =0;
	System.out.println(head.val);

	        while(head != null){
	        	count++;
	            ListNode put = head;

	        	System.out.println("hi");

	            System.out.println("put "+put.val);
	            ListNode cur= rst;
	            System.out.print("cur ");
	            printLinkedList(cur);
	            System.out.print("rst ");
	            printLinkedList(rst);
	            head = head.next;

                if(put.val<cur.val){
    	            System.out.println("head1 "+head.val);

//                    head = head.next; 
    	            System.out.println("head2 "+head.val);

                    put.next = rst;
                    System.out.print("put ");
    	            printLinkedList(put);
    	            System.out.println("Period");
                    rst = put; 
    	            System.out.print("rst ");
                    printLinkedList(rst);
                    continue; 
                }
	            while( cur.next!= null && put.val>cur.next.val){
	                cur = cur.next; 
	            }
	            if(cur.next != null){
	                put.next = cur.next;
	                cur.next =put; 
//	                if(rst.next == null){
//	                    rst = cur; 
//	                }
	            }else
                {
	                // if(put.val<cur.val){
	                //     put.next = cur;
	                // if(cur == rst) rst = put;
	                // }else{
	                    cur.next = put; 
	                    put.next = null; 
	                }
                	            System.out.println("rst ");
                                printLinkedList(rst);
                                
                	            System.out.println("FANSILE ");
                	            System.out.print("head ");
                	            printLinkedList(head);
                	            
                   	            System.out.print("cur ");
                	            printLinkedList(cur);
                	            
                	            System.out.println("rst ");
                                printLinkedList(rst);

                                
	            }
	               
	        return rst; 
	    }
	    
	    
	    public int trailingZeroes(int n) {
	        // Keep track of # of #s divisible by 2 and # of #s divisible by 5, return the smaller one of the two. But the trick is the latter is always less than for former. So just need to track the later. Also notice that when some number i is 25/125... they are not simply contribute a single 5 
	        int k = 0;
	        int power5 = 1; 
	        int count5 = 0;
	        if(n>5){
	            while(power5 * 5 < n){
	                power5 =5 * power5;
	                k++;
	            }
	        }
	        System.out.println("k "+k);
	        int extra = k*(k+1)/2 - k; 
	        
	        for(int i = 1; i<=n;i++){
	            if(i%5==0) count5++;
	        }
	        return count5+extra;
	    }
	
	    
        public int largestPalindrome(int n) {
        // if input is 1 then max is 9 
        if(n == 1){
            return 9;
        }
        
        // if n = 3 then upperBound = 999 and lowerBound = 99
        int upperBound = (int) Math.pow(10, n) - 1, lowerBound = upperBound / 10;
        long maxNumber = (long) upperBound * (long) upperBound;
        
        // represents the first half of the maximum assumed palindrom.
        // e.g. if n = 3 then maxNumber = 999 x 999 = 998001 so firstHalf = 998
        int firstHalf = (int)(maxNumber / (long) Math.pow(10, n));
        
        boolean palindromFound = false;
        long palindrom = 0;
        
        while (!palindromFound) {
            // creates maximum assumed palindrom
            // e.g. if n = 3 first time the maximum assumed palindrom will be 998 899
            palindrom = createPalindrom(firstHalf);
            
            // here i and palindrom/i forms the two factor of assumed palindrom
            for (long i = upperBound; upperBound > lowerBound; i--) {
                // if n= 3 none of the factor of palindrom  can be more than 999 or less than square root of assumed palindrom 
                if (palindrom / i > maxNumber || i * i < palindrom) {
                    break;
                }
                
                // if two factors found, where both of them are n-digits,
                if (palindrom % i == 0) {
                    palindromFound = true;
                    break;
                }
            }

            firstHalf--;
        }

        return (int) (palindrom % 1337);
    }
        
        
        public boolean isAnagram(String s, String t) {
            char[] char1= s.toCharArray();
            char[] char2 = t.toCharArray();
            Arrays.sort(char1);
            Arrays.sort(char2);

            for(int i = 0; i< s.length(); i++){
                if(char1[i] !=char2[i]) return false; 
            }
               return true;
         }

    private long createPalindrom(long num) {
        String str = num + new StringBuilder().append(num).reverse().toString();
        return Long.parseLong(str);
    }
 	   
    
    public TreeNode sortedArrayToBST(int[] nums) {
    	return helpersortedArrayToBST(nums, 0, nums.length -1);
        
    }
    
    private TreeNode helpersortedArrayToBST(int[] nums, int l, int h) {
    	if(nums.length == 0) return null; 
        if(nums.length == 1) return(new TreeNode(nums[0]));
        if(l>h) return null;
        if(l == h) {
        	return(new TreeNode(nums[l]));
        }
        int mid = (l+h)/2;
        TreeNode rst= new TreeNode(nums[mid]);
        rst.left = helpersortedArrayToBST(nums,l,mid-1);
        rst.right = helpersortedArrayToBST(nums,mid+1,h);
        return rst; 
    }
    
   public String reverseWords2(String s) {
   		String[] words = s.split("\\s+");
   	 for(int i =0; i<words.length; i++) {
   			words[i]= revSingleWord(words[i]);
   		}
 	   StringBuilder rst = new StringBuilder();
 	   for(int i =0; i<words.length; i++) {
 		   rst.append(words[i]).append(" ");
 	   }
		return rst.toString();

    }

   private String revSingleWord(String s) {
	   StringBuilder rst = new StringBuilder();
	   for(int i = s.length()-1; i>=0; i--) {
		   rst.append(s.charAt(i));
	   }
	   return rst.toString().trim();
   }
   
   
   public String largestNumber(int[] nums) {
       StringBuilder rst = new StringBuilder(); 
       List<Integer> arr = Arrays.stream(nums).boxed().collect(Collectors.toList());
       while(!arr.isEmpty()){
           int msb = -1; 
           int idx = -1;
           for(int i =0; i<arr.size(); i++){
               if(getMSB(arr.get(i))>msb){
                   msb= getMSB(arr.get(i)); 
                   idx = i;
//               }else if (getMSB(arr.get(i))==msb){
//                   int see = Wrest(getMSB(arr.get(idx),getMSB(arr.get(i))); 
//                   if(see == 0){
//                       // test1
//                       arr.remove(idx);
//                       int t1= 
//                       // test2
//                   }else if(see >0){
//                       idx = i; 
//                   }
               }
           }
           rst.append(arr.get(idx));
           arr.remove(idx);
       }
       return rst.toString();
   }
   
   private int getMSB(int n){
       while(n/10 !=0){
           n = n/10;
       }
       return n; 
   }
   
   // Know that the msb of n1 & n2 are the same, return the one with a higher msb as we continue to look
   // -1 means the first one, 1 means the second one, 0 means can't tell 
   private int Wrest(int n1, int n2){
       
       int msb1;
       int msb2;

       while (n1 != 0 && n2 != 0 ){
           msb1 = getMSB(n1);
           msb2 = getMSB(n2);
           if(msb1>msb2) return -1; 
           else if(msb1<msb2) return 1;
           else {
               n1 = n1/10;
               n2 = n2/10; 
           }
       }
       return 0; 
   }
   /** Given a binary tree, determine if it is height-balanced.
    * Where height-balanced is defined as the depth of the two 
    * subtrees of every node never differ by more than 1.
    * Feel like implement BFS 
    * @param root
    * @return
    */
   public boolean isBalanced(TreeNode root) {
	   if(root == null) return true; 
	   int depth = 0; 
	   Queue<TreeNode> que = new LinkedList<>();
	   que.add(root);
	   int level = 1;
	   while(!que.isEmpty()) {
		   int size = que.size();
		   for(int i = 0; i<size; i++) {
			   if(que.peek().left == null && que.peek().right == null) {
				   if(depth ==0 )  {depth = level;}
				   else if(level>depth+1 || level<depth-1) {return false; }
			   }else {
				   if(que.peek().left != null) {
					   que.add(que.peek().left);
				   }
				   if(que.peek().right != null) {
					   que.add(que.peek().right);
				   }
			   }
			   que.poll();
		   }
		   level++;
	   }
	return true;
       
   }
   
   public int evalRPN(String[] tokens) {
       Stack<Integer> rst = new Stack<>();
       for(int i =0;i<tokens.length ; i++){
           if(tokens[i] == "+") {
        	   int a = rst.pop();
        	   int b = rst.pop();
        	   rst.push(a+b);
           }else if(tokens[i] == "-") {
        	   int a = rst.pop();
        	   int b = rst.pop();
        	   rst.push(b-a);
           }else if(tokens[i] == "*") {
        	   int a = rst.pop();
        	   int b = rst.pop();
        	   rst.push(a*b);
           }else if(tokens[i] == "/") {
        	   int a = rst.pop();
        	   int b = rst.pop();
        	   rst.push(b/a);
           }else {
        	   int num = Integer.parseInt(tokens[i]); 
        	   rst.push(num);
           }
       }
	return rst.peek();
   }
   
   int time = -1;

   public int minMutation(String start, String end, String[] bank) {
	 if(bank == null || bank.length ==0 || start.length() != end.length()) return -1; 
	 helper(start, end, bank, new boolean[bank.length], 0);
	 return time;    	   	     
   }
   public void helper(String start, String end, String[] bank, boolean[] isVisit, int t) {
	     if(start.equals(end)) {
	    	 if(time != -1)
	    		 time = Math.min(time, t);
	    	 else time = t;
	    	 return; 
	     }
	     
	     for(int i =0; i< bank.length; i++) {
	    	 if(isVisit[i]) continue;
	    	 if(singleDiff(start,bank[i])) {
	    		 isVisit[i] = true;
	    		 helper(bank[i],end,bank,isVisit, t+1);
	    		 isVisit[i] = false; 
	    	 }
	     }
   }
   
   private boolean singleDiff (String a, String b){
       int count = 0; 
       for(int i=0;i<a.length(); i++) {
    	   if(a.charAt(i)!= b.charAt(i))
    		   count ++;
       }
       return count ==1; 
   }
   
   
   /** The problem is overflow, should try to avoid doing +/-
    * Simply compare the numbers and a more suitable data structure is needed
    * @param nums
    * @param k
    * @param t
    * @return
    */
   public boolean containsNearbyAlmostDuplicate(int[] nums, int k, int t) {
       int len= nums.length;
       if(len <2) return false;
       if(k==0) return false; 
       if(k<0) return false; 
       if(t == Integer.MAX_VALUE) return true; 
       // i is in [0,len-1) for distinct 
       for (int i =0; i< len-1; i++){
           int bd = Math.min(i+k,len-1);
           if (check(i,bd,t,nums)) return true; 
       }
       return false;
   }
   
   private boolean check(int l, int r, int t, int[] A){
       int dif = Integer.MAX_VALUE;
       for(int j = l ; j< r; j++){
           for(int k = j+1; k<=r; k++){
               int cur = Math.abs(A[j]-A[k]); 
               if(cur<=t) return true; 
               if(cur<dif)
                   dif = cur;
           }
       }
       return false; 
   }
   
   // public ListNode reverseList(ListNode head) {
   //     if(head == null || head.next == null){
   //         return head;
   //     }else{
   //         ListNode newHead = reverseList(head.next);
   //         ListNode tail = newHead;
   //         while(tail.next != null){
   //             tail = tail.next;
   //         }
   //         tail.next = head;
   //         head.next = null; 
   //                 return newHead; 
   //     }     
   // }
   public ListNode reverseList(ListNode head) {
       if(head == null || head.next == null){
           return head;
       }else{
           ListNode nT = head.next;
           ListNode tmp = nT.next;
           
           // Get rid of the first
           nT.next = head;
           head.next = null;
           while(tmp!= null){
                ListNode nH = tmp;
                tmp = tmp.next;
                nH.next = nT;
                nT = nH; 
           }
           return nT; 
       }
   }
   

//   public boolean containsNearbyAlmostDuplicate(int[] nums, int p, int t) {
//       if(nums.length <2) return false;
////       if(p == 0) return true; 
//          if(t<0) return false; 
//      int[] max;
//      if(t>nums.length){
//          max = new int[nums.length]; 
//      }else{
//          max= new int[nums.length- t+1]; 
//      }
//      for(int i =0; i<max.length; i++){
//          int m = Integer.MIN_VALUE;
//          for(int j = i; j< Math.min(nums.length, i+t); j++){
//              for(int k =j+1; k< Math.min(nums.length, i+t); k++){
//                  int now =Math.abs(nums[k]- nums[j]); 
//                  if(now>m){
//                      m = now; 
//                  }
//              }
//          }
//          max[i] = m; 
//                     if(max[max.length-1] == Integer.MIN_VALUE) max[max.length-1] = Integer.MAX_VALUE;
//      }
//      for(int i = 0; i< max.length; i++){
//          if(max[i]<=p) return true; 
//      }
//      return false; 
//  }
   
   public int minimumTotal(List<List<Integer>> A) {
       int size= A.size(); 
       if(size == 0) return 0; 
       List<List<Integer>> DP = new ArrayList<List<Integer>>(); 
       DP.add(A.get(size-1));
       for(int i = size-2; i>=0; i--){
           ArrayList<Integer> line = new ArrayList<>(); 
           int hori = A.get(i).size(); 
           for(int j =0; j< hori; j++){
               line.add(A.get(i).get(j)+Math.min(DP.get(0).get(j),DP.get(0).get(j+1) ));
           }
           DP.add(0,line);
       }
       return DP.get(0).get(0);
   }
   
   public TreeNode constructMaximumBinaryTree(int[] nums) {
	   if(nums.length == 0) return null;
	   if(nums.length == 1) return new TreeNode(nums[0]); 
	 return HconstructMaximumBinaryTree(nums,0,nums.length-1);       
   }
   public TreeNode HconstructMaximumBinaryTree(int[] nums, int start, int end) {
	   if(start == end) return new TreeNode(nums[start]); 
	   if(start>end) return null;
	   int[] max = getMax(nums, start, end);
	   TreeNode root = new TreeNode(max[1]);
	   root.left = HconstructMaximumBinaryTree(nums,start,max[0]-1); 
	   root.right = HconstructMaximumBinaryTree(nums,max[0]+1,end);
	   return root;
	   }
   
   // [start, end] return [i, max] where i is the index of the max value in the array 
   private int[] getMax(int[] arr, int start, int end) {
	   int max = Integer.MIN_VALUE; 
	   int i = start;
	   int idx = -1; 
	   for(; i<=end; i++) {
		   if(arr[i]>max) {
			   max= arr[i];
			   idx = i; 
		   }
	   }
	   return new int[]{idx, max}; 
   }
   
   public int canCompleteCircuit(int[] gas, int[] cost) {
	   ArrayList<Integer> startCandi = new ArrayList<>();
	   for(int i =0; i<gas.length ; i++) {
		   if(gas[i]<= cost[i])
			   startCandi.add(i);
	   }
	   for(int i: startCandi) {
		   int cur= i;
		   int box = 0; 
		   while(box + gas[i]> cost[i] && cur != i) {
			   cur = (cur+1)% gas.length;
			   box = box + gas[i] -cost[i];
		   }
		   if(cur == i) return i; 
	   }
	   return -1;
   }
   
   
   public boolean isIsomorphic(String s, String t) {
       if(s.length()!=t.length()) return false;
       HashMap<Character, Character> dic = new HashMap<>();
       HashSet<Character> forbidden = new HashSet<>();
       for(int i=0; i< s.length(); i++){
           if(dic.containsKey(s.charAt(i))){
               if(t.charAt(i)!=dic.get(s.charAt(i))) return false;
           }else{
        	   if(forbidden.contains(t.charAt(i))) return false;
               dic.put(s.charAt(i),t.charAt(i));
               forbidden.add(t.charAt(i));
//               System.out.println("dic "+dic);
           }
       }
       return true;
   }
   
   
   public List<Integer> lexicalOrder(int n) {
	   List<Integer> rst = new ArrayList<>();
	   int msb = 0; 
       int cur = n;
       while(cur%10 != 0) {
    	   cur/=10;
    	   msb++;
       }
	   while(cur/10 !=0 ) {
		   cur /=10;
		   msb++;
	   }
       int i = 0 ;
       while(i<= msb) {
    	   System.out.println("i "+ i);
    	   int upper = Math.min(2*powerBase10(i)-1, n);
    	   System.out.println("upper "+ upper);

    	   for(int j =powerBase10(i); j<= upper; j++ ) {
    		   rst.add(j);
    		   System.out.println("j "+j);
    	   }
    	   i++;
       }
       for(int x = 2; x<=9; x++) {
    	   i =0; 
    	      while(i< msb) {
    	    	   int upper = Math.min((x+1)*powerBase10(i)-1, n);
    	    	   for(int j =x*powerBase10(i); j<= upper; j++ ) {
    	    		   rst.add(j);
    	    	   }
    	    	   i++;
    	       }       
    	}
    	   return rst; 
   }
   // return 10^i 
   private int powerBase10 (int i) {
	   if(i== 0) return 1;
	   return 10 * powerBase10(i-1); 
   }
   
   
   public int reverseBits(int n) {
//       int j = Integer.highestOneBit(n);
//       int move = 32-powerOf2(j);
	   int move =	Integer.numberOfLeadingZeros(n);
       int rst = 0; 
       while(n>0){
           if((n & 1) ==1){
               rst <<=1;
               rst |=1;
           }else{
               rst <<=1;
           }
           n>>=1; 
       }
       for(int i=0; i<=move; i++){
          rst<<=1; 
       }
       return rst; 
   }
   
//   private int powerOf2(int n) {
//	   if(n==1) return 0;
//	   return 1+powerOf2(n/2);
//   }
//   
   
   public int maxProfit(int k, int[] prices) {
	return k;

   }

   Comparator<Integer> conve = (Integer a, Integer b) -> b-a;
   public String largestNumber1(int[] nums) {
       String[] dic = new String[nums.length];
//       Arrays.sort((Integer[])nums,conve);
       for(int i =0; i< dic.length; i++) {
    	   dic[i] = Integer.toString(nums[i]);
       }
       Arrays.sort(dic, (String a, String b) -> b.compareTo(a));
       StringBuilder rst = new StringBuilder();
       for(String s: dic) {
    	   rst.append(s);
       }
	return rst.toString();
   }
   
   public int candy(int[] r) {
       int rst = 0;
       int len = r.length; 
       if(len <=1) return rst+len; 
       int[] dic = new int[len];
       for(int i=0; i<len ; i++){
           dic[i] = 1;
       }

       int lastMax = 0; 
       if(r[0]>r[1]) dic[0] ++;
        else lastMax++;
     
       System.out.println();
       for(int i =1; i<=len-2; i++){
    	   System.out.println("i "+i);
           if(r[i]<=r[i-1] && r[i]<=r[i+1]) {
        	   lastMax=i+1;
               continue;
           }else if(r[i]<=r[i-1] && r[i]>r[i+1]){
//                               for(int j = lastMax; j<=i;j++){
//                                                   dic[j]++;
//                               }
               for(int j = i; j>= lastMax;j--){
                   dic[j] = Math.max(dic[j], dic[j+1]+1);
               }
           }else if (r[i]>r[i-1] && r[i]<=r[i+1]){
             dic[i] = dic[i-1]+1;
            lastMax =i+1;
         }else{
            dic[i]= Math.max(dic[i-1]+1, dic[i+1]+1);
         } 
       }
       if(r[len-1]>r[len-2]){
           for(int j = lastMax; j<=len-1;j++){
                                                   dic[j]= dic[j-1]+1;
                               }
       }
       System.out.println("End");
       System.out.println("lastMax "+lastMax + " dic ");
       System.out.println();
       for(int i=0; i<len ; i++){
           rst+= dic[i];
       }
       return rst; 
   }
   
//   private List<List<Integer>> rst = new ArrayList<>();
//   public List<List<Integer>> pathSum(TreeNode root, int sum) {
//	   ArrayList<Integer> cur = new ArrayList<>();
//	   helperPSum(root,sum,cur);
//	   return rst;   
//   }
//   public void helperPSum(TreeNode root, int sum, ArrayList<Integer> cur) {
//       if(root == null || root.val>sum) return;
//       if(root.left == null && root.right ==null){
//           if(root.val == sum) {
//        	   cur.add(root.val);
//        	   rst.add(cur);
//           }
//           return; 
//       }
//       cur.add(root.val);
//       helperPSum(root.left, sum-root.val, cur);
//       helperPSum(root.right, sum-root.val, cur);
//   }

   
   private List<List<Integer>> rst = new ArrayList<>();
   public List<List<Integer>> pathSum(TreeNode root, int sum) {
	   ArrayList<Integer> cur = new ArrayList<>();
	   helperPSum(root,sum,cur);
	   return rst;   
   }
   public void helperPSum(TreeNode root, int sum, ArrayList<Integer> cur) {
       if(root == null) {
           //  if(!cur.isEmpty())
           // cur.remove(cur.size()-1);
           return;
       }
               	   cur.add(root.val);

       if(root.left == null && root.right ==null &&root.val == sum){
           
        	   rst.add(new ArrayList<Integer>(cur));
               cur.remove(cur.size()-1);
                return; 
       }else{
       helperPSum(root.left, sum-root.val, cur);
       helperPSum(root.right, sum-root.val, cur);
       }
      cur.remove(cur.size()-1);
   }
//   public int candy(int[] r) {
//       int rst = 0;
//       int len = r.length; 
//       if(len <=1) return rst+len; 
//       int[] dic = new int[len];
//       for(int i=0; i<len ; i++){
//           dic[i] = 1;
//       }
//       int lastMax = 0; 
//       if(r[0]>r[1]) {dic[0] ++;
//                     }else{lastMax++;} 
//       for(int i =1; i<=len-2; i++){
//           if(r[i]<=r[i-1] && r[i]<=r[i+1]) {
//               lastMax+=2;
//               continue;
//           }else if(r[i]<=r[i-1] && r[i]>r[i+1]){
//                               for(int j = lastMax; j<=i;j++){
//                                                   dic[j]++;
//                               }
//           }else if (r[i]>r[i-1] && r[i]<=r[i+1]){
//               dic[i]++;
//              lastMax +=2;
//           }else{
//              dic[i]++;
//              lastMax ++;
//           } 
//           //r[i]>r[i-1] && r[i]>r[i+1]
//           // if(r[i]>r[r-1])  
//       }
//       if(r[len-1]>r[len-2]){
//           for(int j = lastMax; j<=len-1;j++){
//                                                   dic[j]++;
//                               }
//       }
//       for(int i=0; i<len ; i++){
//           rst+= dic[i];
//       }
//       return rst; 
//   }
   
   
//   public List<Integer> rightSideView(TreeNode root) {
//       List<Integer> rst = new ArrayList<>();
//       if(root == null) return rst; 
//       rst.add(root.val);
//       if(root.right!=null) {
//    	   rst.add(root.right.val);
//       }else {
//    	   rst.add(root.left.val);
//       }
//       while()
//   }
   
   
   public int sumNumbers(TreeNode root) {
       if(root == null) return 0;
       if(root.val == 0) return 0;
       if(root.left == null && root.right == null) return root.val;
       int left= sumNumbers(root.left);
       int right = sumNumbers(root.right);
       int l1 = 0;
       int l2 = 0;
       if(left !=0){
           String par1 = Integer.toString(root.val) + Integer.toString(left);
           l1 = Integer.parseInt(par1);
       }
       if(right !=0){
           String par2 = Integer.toString(root.val) + Integer.toString(right);
           l2 = Integer.parseInt(par2);
       }
       

       return l1+l2;
   }
   
   public void connect(TreeLinkNode root) {
       if(root == null) return;
       Queue<TreeLinkNode> que = new LinkedList<>();
//       que.add(root);
       root.next = null; 
       if(root.left != null) que.add(root.left);
       if(root.right!=null) que.add(root.right);
       int len = que.size();
       while(!que.isEmpty()) {
//    	    cur = null;
    	   TreeLinkNode[] lst = new TreeLinkNode[len];
    	   for(int i =0; i<len; i++) {
    		   TreeLinkNode cur = que.poll();
    		   lst[i] = cur; 
//               if(que.peek() == null) {
//                   cur.next= null;
//                   return; 
//               }else{
//            	   cur.next = que.peek();
//               }  
               if(cur.left != null)  que.add(cur.left);
        	   if(cur.right != null) que.add(cur.right);
    	   }
    	   for(int i = 0; i<len-1; i++) {
    		   lst[i].next = lst[i+1];
    	   }
    	   lst[len-1].next = null;
//    	   TreeLinkNode cur = que.poll();
//    	   cur.next = null;
    	   len = que.size();
       }
       
   }
   
   
   public TreeNode sortedListToBST(ListNode head) {
       if(head == null) return null;
       if(head.next == null) return new TreeNode(head.val);
       ListNode pre = head; 
       ListNode l = head;
       ListNode r = head; 
       while(r != null && r.next !=null){
           if(l!=head){
               pre=pre.next; 
           }
           l = l.next;
           r=r.next.next;
       }
       System.out.println("l "+l.val);
       System.out.println("pre "+pre.val);

       TreeNode root= new TreeNode(l.val);
       pre.next = null;
       root.left = sortedListToBST(head);
       root.right = sortedListToBST(l.next);
       return root;
   }
   
   public int findCircleNum(int[][] M) {
       int n = M.length; 
       int groupNum = 0; 
       int count =0;
       HashMap<Integer, Integer> dic = new HashMap<>();
       for(int i =0; i< n; i++){
           for(int j = i; j<n;j++){
        	   if(M[i][j] ==1) {
        		   if(!dic.containsKey(i) && !dic.containsKey(j)){
                       dic.put(i,groupNum);
                       dic.put(j,groupNum);
                       groupNum++;
                       count ++;
                   }else if(!dic.containsKey(i) && dic.containsKey(j)){
                       int group = dic.get(j);
                       dic.put(i,group);
                   }else if(dic.containsKey(i) && !dic.containsKey(j)){
                       int group = dic.get(i);
                       dic.put(j,group);
                   }else{
                       if(dic.get(i)!=dic.get(j)){
                           int obse = dic.get(j);
                           int mergeTo = dic.get(i);
                           Set<Integer> keys = dic.keySet();
                           for(int k: keys){
                               if(dic.get(k)==obse){
                                   dic.replace(k,mergeTo);
                               }
                           }
                           count --;
                       }

                   }
        	   }
               
           }
       }
       return count; 
   }
   
   
//   public int numDistinct(String s, String t) {
//       //dp[i][j] represents # of distinct representations in s[i,end] of t[j,end] (both sides inclusive)
//       if(s.length() ==0) return 1;
//       if(t.length()==0) return 0; 
//       int[] dp = new int[t.length()];
//       
//       // fill in the order from the end to the front
//       for(int i = s.length()-1; i>=0; i--){
//           for(int j = t.length()-1; j>=0;j--){
//               if(s.charAt(i) == t.charAt(j)){
//                   if(i == s.length()-1){
//                       dp[j]=1;
//                   }else{
//                           if(j==t.length()-1){
//                               dp[j] +=1;
//                           }else{
//                               dp[j]= dp[j]+dp[j+1];
//                           }
//                       }
//               }
//        	   System.out.println("new");
//               printIntArray(dp);
//
//               for(int p =0; p<t.length(); p++) {
//                   printIntArray(dp);
//                   System.out.println();
//               }             
//           }
//       }
//       return dp[0];
//   }
   
   public int numDistinct(String s, String t) {
       //dp[i][j] represents # of distinct representations in s[i,end] of t[j,end] (both sides inclusive)
       if(s.length() ==0) return 0;
       if(t.length()==0) return 0; 
       int[][] dp = new int[t.length()][s.length()];
       
       // fill in the order from the end to the front
       int count = 0 ;
       for(int j = t.length()-1; j>=0;j--){
           if(s.charAt(s.length()-1) == t.charAt(j) && count ==0){
        	   count ++;
        	   dp[j][s.length()-1]=1;
           }
       }
       for(int i = s.length()-2; i>=0; i--){
           for(int j = t.length()-1; j>=0;j--){
               if(s.charAt(i) == t.charAt(j)){
                   if(i == s.length()-1){
                       dp[j][i] =1;
                   }else{
                           if(j==t.length()-1){
                               dp[j][i] =dp[j][i+1]+1;
                           }else{
                               dp[j][i] = dp[j][i+1] + dp[j+1][i+1];
                           }
                       }
               }else {
            	   if(i != s.length()-1)
            	   dp[j][i] = dp[j][i+1];
               }
           }
       }
       return dp[0][0];
   }
   
   
   
   // return sum of l1, l2, and the carry 
   public ListNode helperAT(ListNode l1, ListNode l2, int carry) {
       if(l1 == null && l2 == null) return new ListNode(carry);
       if(l1 ==null && l2 != null) {
           if(carry ==0) return l2;
           return helperAT(l2, new ListNode(carry), 0);
       }
       if(l1 != null && l2 == null) {
           if(carry ==0) return l1;
           return helperAT(l1, new ListNode(carry), 0);
       }
       ListNode last1 = l1;
       ListNode last2 = l2; 
       // Base case 
       if(last1.next == null && last2.next == null){
           int rst = l1.val + l2.val +carry;
           if(rst<10){
                return new ListNode(rst);
           }else{
               ListNode one = new ListNode(1); 
               ListNode two = new ListNode(rst%10);
               one.next = two;
               return one; 
           }
       }
       if(last1.next == null) {
    	   l1 = null; 
       }else {
    	   ListNode end1= l1;

    	   while(last1.next != null){
             if(last1 !=l1){
                 end1= end1.next;
             }
               last1 = last1.next;
           }
           end1.next = null;
       }
       
       if(last2.next == null) {
    	   l2 = null;
       }else {
       
       ListNode end2 = l2;

        while(last2.next != null){
         if(last2 !=l2){
             end2= end2.next;
         }
           last2 = last2.next;
       }
        
           end2.next = null;
       }
       
       ListNode one = null;
       int count = last1.val + last2.val + carry;
       if( count>=10){
           System.out.println(">=10");
           System.out.println("l1");
           printLinkedList(l1);
           System.out.println("l2");
           printLinkedList(l2);
            one = helperAT(l1,l2,1);
            System.out.println("one");
            printLinkedList(one);
       }else{
           System.out.println("<10");
           System.out.println("l1");
            printLinkedList(l1);
            System.out.println("l2");
            printLinkedList(l2);
            one = helperAT(l1,l2,0);
            System.out.println("one");
            printLinkedList(one);
       }
       ListNode two = new ListNode(count %10);
       System.out.println("two");
       printLinkedList(two);
       ListNode oneEnd = findLast(one);
       oneEnd.next = two;
       return one; 
       
   }
   
   public ListNode findLast(ListNode head) {
	   if(head == null) return null;
	   ListNode cur = head;
	   while(cur.next != null){
		   cur =cur.next;
	   }
	   return cur; 
   }
   public ListNode addTwoNumbers1(ListNode l1, ListNode l2) {
       return helperAT(l1,l1,0);
   }
   
   
   public ListNode addTwoNumbers(ListNode l1, ListNode l2) {
	   if(l1 == null && l2 == null) return null;
	   if(l1 == null && l2 != null) return l2;
	   if(l1 != null && l2 == null) return l1;

//       Deque<Integer> que1= new LinkedList<>();
//       Deque<Integer> que2= new LinkedList<>();

	   Stack<Integer> st1 = new Stack<>();
	   Stack<Integer> st2 = new Stack<>();
	   while(l1 != null) {
		   st1.add(l1.val);
		   l1 = l1.next;
	   }
	   int n1 = 0;
	   while(!st1.isEmpty()) {
		   n1 = n1*10 + st1.pop();
	   }
	   while(l2 != null) {
		   st2.add(l2.val);
		   l2 = l2.next;
	   }
	   int n2 = 0;
	   while(!st2.isEmpty()) {
		   n2 = n2*10 + st2.pop();
	   }
	   int rst = n1 + n2;
	   ListNode head = new ListNode(rst%10);
	   ListNode cur = head;

	   while(rst/10 != 0) {
		   rst /=10;
		   cur.next = new ListNode(rst%10);
		   cur = cur.next; 
	   }
	   return head; 
		  
   }
   
   public List<TreeNode> generateTrees(int n) {
       // List<TreeNode> dp[][] = new ArrayList<TreeNode>[n][n];
       // for(int i =1; i<=n; i++){
       //     dp[i][i] = new TreeNode(i) ; 
       // }
       List<List<List<TreeNode>>> dp = new ArrayList<List<List<TreeNode>>>(n);
       // Fill in the table 
       // i is left bound, j is right bound, both sides are inclusive
       for(int i =0; i<n; i++){
    	   List<TreeNode> same =  new ArrayList<>();
    	   same.add(new TreeNode(i));
           System.out.println("Here 0");
           dp.get(i).set(i, same);
           System.out.println("Here 1");
           for(int j =i+1; j<n; j++){
//               if(i>j) {
//            	   List<TreeNode> greater =  new ArrayList<>();
//            	   greater.add(null);
//            	   dp.get(i).set(j,greater);
//               } 
               System.out.println("Here 2");

               if(dp.get(i).get(j) != null) continue;
               else{
                   List<TreeNode> list = new ArrayList<>();
                   for(int k = i; k<= j; k++){
                       List<TreeNode> left = dp.get(i).get(k-1);
                       if(left == null) left.add(null);
                       List<TreeNode> right = dp.get(k+1).get(j);
                       if(right == null) right.add(null);
                       for(TreeNode lnode: left){
                           for(TreeNode rnode : right){
                               TreeNode root= new TreeNode(k);
                               root.left = lnode;
                               root.right = rnode;
                               list.add(root); 
                           }
                       }   
                   }
                   System.out.println("Here 3");

                   dp.get(i).set(j, list);
                   System.out.println("Here 4");

               }
           }
       }
       return dp.get(n-1).get(n-1);
   }
   
   public String longestPalindrome(String s) {
       int rstLen = 0; 
       int len = s.length();
       String rst = "";
       for(int i = 0; i<len; i++ ){
           for(int j = i+1; j<=len; j++ ){
               if(isPalindrome(s.substring(i,j))&& j-i>rstLen){
                       rst = s.substring(i,j);
                       rstLen = j-i;
               }
           }
       }
                  return rst; 
   }
   
   public int integerBreak(int n) {
       int[] dp = new int[1+n];
       dp[0] = 0;
       dp[1] = 0;
       dp[2] = 1;
       if(n == 2) return 1; 
       for(int i = 3; i<n+1 ; i++){
           int put = 0 ; 
           for(int j = i/2; j>1; j--){
               if(Math.max(j,dp[j]) * Math.max(i-j, dp[i-j])>put)
                   put = Math.max(j,dp[j]) * Math.max(i-j, dp[i-j]);
           }
           dp[i] = put; 
       }
       return dp[n];
   }
   
   public int lengthOfLongestSubstringKDistinct(String s, int k) {
       if(s==null || k==0) return 0;
      
       int len = s.length();
       if(len<=k) return len; 
       int rst = 0; 
       for(int i=0; i<len-rst; i++){
           HashMap<Character, Integer> dic = new HashMap<>();
           int total = k; 
           int end= i;
           while(total!= -1 && end!= len){
               if(!dic.containsKey(s.charAt(end))){
                   dic.put(s.charAt(end),end);
                   total--;
               }
               end++;
//               System.out.println("i "+i);
//               System.out.println("total "+total);
//               System.out.println("end "+end);
           }
           if(total==-1){
        	   if(end-1-i>rst) rst = end-1 -i;            
        	   }else if(end - i > rst) rst =end -i;  
       }
       
       return rst; 
   }
   
   
   public int nthUglyNumber(int n) {
       if(n<=6) return n;
       int[] arr = new int[n];
       for(int i =0; i<6; i++){
           arr[i] = i+1;
       }
       int two = 0;
       int three = 0;
       int five = 0; 
       for(int i =6; i<n; i++){
           while(arr[two]*2<=arr[i-1]) two++;
           System.out.println("two "+ two);
           while(arr[three]*3<=arr[i-1]) three++;
           System.out.println("three "+ three);
           while(arr[five]*5<=arr[i-1]) five++;
           System.out.println("five "+ five);
           int put = Math.min(2*arr[two],3*arr[three]);
           arr[i] = Math.min(put,5*arr[five]);
       }
       return arr[n-1];
   }
   
//   public int numDecodings(String s) {
//       if(s==null) return 0;
//       int len = s.length();
//       int[] dic = new int[len+2];
//       dic[len+1] = 1; 
//       dic[len] = 1;
//       for(int i = len-1; i>=1; i--){
//    	   if((s.charAt(i-1) == '1'&& s.charAt(i) - '0'>0) || (s.charAt(i-1)=='2' && s.charAt(i) - '0'<=6 && s.charAt(i) - '0'>0)){
//               dic[i] = dic[i+1] + dic[i+2];
//    	   }else{
//               dic[i] = dic[i+1];
//    	   }
////    	   
////               if(s.charAt(i-1)=='1' ||s.charAt(i-1)=='2'  ){
////                   dic[i] = dic[i+1];
////               }else return 0; 
////           }else 
////           }
//       }
//       printIntArray(dic);
//   	System.out.println();
//
//       return dic[1];
//   }
   
   public int numDecodings(String s) {
       if(s==null) return 0;
       if(s.charAt(0) == '0') return 0;
       int len = s.length();
       if(len ==1) return 1;
       int[] dic = new int[len+1];
       dic[len] = 1; 
       if(s.charAt(len-1) == '0') dic[len-1] = 0; 
           else dic[len-1] = 1; 
       for(int i = len-2; i>=0; i--){
           if(s.charAt(i) == '0'){
               dic[i] = 0;
           }else if( (s.charAt(i)=='1' && s.charAt(i+1) - '0'>0 )||(s.charAt(i)=='2' && s.charAt(i+1) - '0'<=6 && s.charAt(i+1) - '0'>0)){
                   dic[i] = dic[i+1] + dic[i+2];
           }else {
                           dic[i] = dic[i+1];
           }
       }
       return dic[0];
   }
   
   public boolean checkPossibility(int[] A) {
       if(A==null) return false; 
	   int[] nums = A;
       int count = 0;
       
       int len = nums.length;
       if(len <=2) return true; 
       int i = 0; 
       if(nums[0]>nums[1]) {
    	   count ++; 
       }
       i++;
       for(; i<len-1 ;i++){
           if(nums[i]<nums[i-1] || nums[i]>nums[i+1]) {
        	   if(nums[i]<nums[i-1] && nums[i]>nums[i+1]) return false; 
        	   else if(nums[i]<nums[i-1]) nums[i] = nums[i-1];
        	   else {
        		   if(nums[i+1]>=nums[i-1]) nums[i] = nums[i+1];
        		   else nums[i+1] = nums[i];
        	   }
        	   count ++;
        	   i++;  
           }
           System.out.println("i "+ i);

           System.out.println("Count "+ count);

       }
       if(nums[len-1]<nums[len-2]) count ++; 
       System.out.println("Count "+ count);
       return count <=1; 
   }
   // 错在bijection 还是要用hashmap
//   public boolean wordPattern(String pattern, String str) {
//       if(pattern == null || str == null) return false;
//       int len1 = pattern.length();
//       if(len1 ==0){
//           if(str.length() == 0) return true;
//           return false; 
//       }
//       String[] dic = new String[256]; 
//       String[] lst = str.split("\\s+");
//       printStringArray(lst);
//       int len2 = lst.length;
//       if(len1 != len2) return false; 
//       for(int i = 0 ; i<len1; i++){
//           if(dic[pattern.charAt(i)] == null){
//        	   System.out.println("null");
//               dic[pattern.charAt(i)] = lst[i];
//           }else{
//        	   System.out.println("here");
//               printStringArray(dic);
//               System.out.println();
//               if(!lst[i].equals(dic[pattern.charAt(i)])) return false; 
//           }
//       }
//       return true; 
//   }
 public boolean wordPattern(String pattern, String str) {
  if(pattern == null || str == null) return false;
  int len1 = pattern.length();
  if(len1 ==0){
  if(str.length() == 0) return true;
  return false; 
  }
  String[] lst = str.split("\\s+");
  int len2 = lst.length;
  if(len1 != len2) return false;
  HashMap<Character, String> dic = new HashMap<>();
  HashSet<String> dupli = new HashSet<>();
  for(int i = 0 ; i<len1; i++){
  if(!dic.containsKey(pattern.charAt(i))){
	  if(dupli.contains(lst[i])) return false; 
	  dic.put(pattern.charAt(i), lst[i]);
	  dupli.add(lst[i]);
  }else{
  if(!lst[i].equals(dic.get(pattern.charAt(i)))) return false; 
  }
  }
  return true;	   
    }
 
 public int lengthOfLastWord(String s) {
     if(s==null) return 0; 
     String[] splited = s.split("\\s+");
//     System.out.println("splited.length " + splited.length);
     if(splited.length == 0 || (splited.length == 1 && splited[0].length() == 0)) return 0;
     return splited[splited.length-1].length();
 }
 
 public ListNode partition(ListNode head, int x) {
     if(head == null || head.next == null) return head; 
     ListNode preCur = head; 
     ListNode cur = head;
     ListNode pre= head;
     ListNode pio = head.next;
     System.out.println("Why");
     printLinkedList(head);
     if(cur.val>=x){
         System.out.println("Why");
         printLinkedList(head);
         while(cur != null && cur.val >=x){
             if(cur != head)  preCur = preCur.next; 
             cur =cur.next;
             System.out.println("Why2");
             printLinkedList(head);
         }
         
         System.out.println("Why3");
         printLinkedList(head);
         if(cur == null) return head;
         System.out.println("Bad");
         preCur.next = cur.next;
         cur.next = head; 
         head = cur;
         pre = preCur;
         pio = preCur.next;
     }
     
     while(pio!=null){
         if(pio.val>=x){
             pre = pio;
             pio = pio.next; 
         }else{
        	 if(pre == cur) {
        		 pre = pio;
                 pio = pio.next; 
                 cur = cur.next;
        	 }else { 
             pre.next = pio.next;
             pio.next = cur.next;
             cur.next = pio;
             pio = pre.next;
             cur = cur.next;
        	 }
         }
     }
     return head; 
 }
 
 public List<List<Integer>> combine(int n, int k) {
	// Trivial non-sense
     if(k<0) return new ArrayList<>();
     if(n<1) return new ArrayList<>();
	 List<List<Integer>> rstt = new ArrayList<>();
	 combine(rstt, new ArrayList<Integer>(), 1,n,k);
	 return rstt; 
 }
 
 public static void combine(List<List<Integer>> rstt, List<Integer> cur, int start, int n, int k) {
     if(k==0) {
    	 rstt.add(new ArrayList<Integer>(cur));
    	 return;
     }
     for(int i = start; i<=n; i++) {
    	 cur.add(n);
    	 combine(rstt,cur,i+1,n,k-1);
         cur.remove(cur.size()-1);
     }   
 }
 

// List<List<Integer>> appen = combine(n-1,k-1);
// if(appen.size() != 0){
//     for(List<Integer> each: appen){
//         cur.addAll(each);
//         rstt.add(cur);
//     }
// }
// List<List<Integer>> appen2 = combine(n-1,k);
// if(appen2.size() != 0){
//     for(List<Integer> each: appen2){
//         cur.addAll(each);
//         rstt.add(cur);
//     }
// }
// return rstt; 
 
 public int maxDistToClosest(int[] s) {
     int rst = 0 ; 
     // l is the leftmost 1 position
     int l = 0;
     // r is the rightmost 1 position
     int r = 0;
     // Find the start position 
     while(l<s.length && s[l]!=1) l++;
     rst = l;
     while(l<s.length){
         // Invariant, r>l
    	 r = l+1;
         while(r<s.length && s[r]!=1) r++;
         // there is no longer candidate interval; we are safe to return  
         if(r>=s.length){
             rst = Math.max(rst,s.length-1-l);
             return rst; 
         }
         // now [l,r] is the interval that we could look at 
         int cur =  (r+l)/2-l;
         rst = Math.max(rst,cur);
         l = r; 
     }
     return rst; 
 }
 /** 这个思路非常简单的解法的问题就在于 如果max是多个重复 很难keep track of all the indexes
  * 还是乖乖用BFS是正道
  * @param matrix
  * @return
  */
// public List<int[]> pacificAtlantic(int[][] matrix) {
//	 if(matrix == null) 	return null;
//	 int m = matrix.length; 
//	 if(m== 0 )	return null;
//	 int n = matrix[0].length;
//	 if(n==0) 	return null;
//	 List<int[]> rst = new ArrayList<>();
//	 HashMap<Integer, Integer> dic = new HashMap<>();
//	 for(int i =0; i<m; i++) {
//		 int max = matrix[i][0]; 
//		 int idx = 0; 
//		 for (int j =0; j<n; j++) {
//			 if(matrix[i][j]>max) {
//				 max = matrix[i][j];
//				 idx = j; 
//			 }
//		 }
//		 dic.put(i, idx);
//	 }
//	 System.out.println(dic);
//	 
//	 for(int j =0; j<n; j++) {
//		 int max = matrix[0][j]; 
//		 int idx = 0; 
//		 for (int i =0; i<m; i++) {
//			 if(matrix[i][j]>max) {
//				 max = matrix[i][j];
//				 idx = i; 
//			 }
//		 }
////		 if(dic.containsKey(idx)) {
//			 if(dic.get(idx) != j) {
//				 int[] put = new int[2];
//				 put[0] = idx;
//				 put[1] = j;
//				 rst.add(put);
//			 }
////		 }		 
//	 }
//     
//	 Set<Integer> key = dic.keySet();
//	 Iterator<Integer> iter = key.iterator();
//	 while(iter.hasNext()) {
//		 int[] put = new int[2];
//		 int hey = iter.next();
//		 put[0] = hey;
//		 put[1] = dic.get(hey);
//		 rst.add(put);
//	 }
//	 return rst;
// }
 
 public List<int[]> pacificAtlantic(int[][] matrix) {
     if(matrix == null) 	return null;
     List<int[]> rst = new ArrayList<>();
     int m = matrix.length;
     if(m== 0)	return rst;
     int n = matrix[0].length; 
     if(n== 0)	return rst;
     int[][] visited1 = new int[m][n];
     int[][] visited2 = new int[m][n];
      for(int i=0;i<n;i++) {
 	 visited1 = dfs(matrix, visited1, 0,i);
      visited2= dfs(matrix, visited2, m-1,i);
  }
  for(int j = 0; j<m; j++) {
 	 visited1 = dfs(matrix, visited1, j,0);
      visited2= dfs(matrix, visited2, j,n-1);
  }
     for(int i =0; i<m; i++) {
		    for (int j =0; j<n; j++) {
			    if(visited1[i][j]==1 &&visited2[i][j]==1) {
				    int[] put = new int[]{i,j};
                 rst.add(put);
			    }
		    }
	    }
     return rst; 
 }
 
     private int[][] dfs(int[][] matrix, int[][] visited, int x, int y){
     int[][] dirs = {{-1,0},{1,0},{0,-1},{0,1}};
     int m = matrix.length;
     int n = matrix[0].length; 
     
     Queue<int[]> que1 = new LinkedList<>();
     que1.add(new int[]{x,y});
     visited[x][y]=1;
     while(!que1.isEmpty()){
         int[] cell = que1.poll();
         for(int[] d: dirs){
             int r= cell[0] + d[0];
             int c = cell[1]+d[1];
             if(r<0 || r>=m || c<0 || c>=n || matrix[r][c]<matrix[cell[0]][cell[1]] || visited[r][c] ==1) continue;
             que1.add(new int[] {r,c});
             visited[r][c] = 1; 
         }
     }
     return visited;
     }
   // 错在没有记录visited 很危险啊啊啊啊啊！！！！ 这个还没有infinite loop  
//     public int longestIncreasingPath(int[][] matrix) {
//         if(matrix == null) 	return 0;
//         int m = matrix.length;
//         if(m== 0)	return 0;
//         int n = matrix[0].length; 
//         if(n== 0)	return 0;
//         int rst = 0; 
//         int idx1 = 0; 
//         int idx2 = 0; 
//            for(int i =0; i<m; i++) {
// 		    for (int j =0; j<n; j++) {
// 			    if(dfs(matrix, i,j)>rst) {
// 				    rst= dfs(matrix, i,j);
// 				    idx1= i;
// 				    idx2= j;
// 			    }
// 		    }
// 	    }
// 		    System.out.println("i "+idx1);
// 		    System.out.println("j" +idx2);
//
//         return rst; 
//     }
//     
//        
//     
//         private int dfs(int[][] matrix, int x, int y){
//         int[][] dirs = {{-1,0},{1,0},{0,-1},{0,1}};
//         int m = matrix.length;
//         int n = matrix[0].length; 
//    	 int[][] visited = new int[m][n];
//         int len = 0; 
//         Queue<int[]> que1 = new LinkedList<>();
//         que1.add(new int[]{x,y});
//         len++; 
//         visited[x][y] = 1; 
//         while(!que1.isEmpty()){
//             int[] cell = que1.poll();
//             for(int[] d: dirs){
//                 int r= cell[0] + d[0];
//                 int c = cell[1]+d[1];
//                 if(r<0 || r>=m || c<0 || c>=n || matrix[r][c]<=matrix[cell[0]][cell[1]]) {
//                     continue;
//                 }
//                 que1.add(new int[] {r,c});
//                 len++;
//                 visited[r][c] = 1; 
//             }
//         }
//         System.out.println("visited");
//         printIntMatrix(visited);
//         return len;
//         }
 
     public int longestIncreasingPath(int[][] matrix) {
         if(matrix == null) 	return 0;
         int m = matrix.length;
         if(m== 0)	return 0;
         int n = matrix[0].length; 
         if(n== 0)	return 0;
         int rst = 0; 
	     int[][] cache = new int[m][n];
            for(int i =0; i<m; i++) {
 		    for (int j =0; j<n; j++) {
                 int len = dfs(matrix, i,j,m,n,cache);
 			     rst = Math.max(rst,len);
 		    }
 	    }
         return rst; 
     }
     
        
     
         private int dfs(int[][] matrix, int x, int y, int m, int n, int[][] cache){
        	 if(cache[x][y] != 0) return cache[x][y];
         int[][] dirs = {{-1,0},{1,0},{0,-1},{0,1}};
         int len = 1; 
             for(int[] d: dirs){
                 int r= x+ d[0];
                 int c = y+d[1];
                 if(r<0 || r>=m || c<0 || c>=n || matrix[r][c]<=matrix[x][y] ) {
                     continue;
                 }
                 len = Math.max(len, 1+dfs(matrix, r, c, m,n,cache));
             }
             cache[x][y] = len; 
             return len;
         }
         
         public int coinChange1(int[] coins, int amount) {
             if(amount <0) return -1; 
             if(amount ==0) return 0; 
             int rst = 0; 
             Queue<Integer> bfs = new LinkedList<>();
             HashSet<Integer> visited  = new HashSet<>();
             bfs.add(amount); 
             visited.add(amount);
             while(!bfs.isEmpty()){
                 rst++;
                 int size = bfs.size();
                 for(int j=0 ; j<size; j++){
                     int cur = bfs.poll(); 
                 for(int i :coins){
                     if(cur-i==0) return rst; 
                     else if (cur-i<0 || visited.contains(cur-i)){
                         continue;
                     }else{
                        bfs.add(cur-i); 
                        visited.add(cur-i);
                     }
                 }
                 }
                 
             }
             return -1; 
         }
         
         private Stack<TreeNode> stack = new Stack<>();
         private List<Integer> rst1 = new ArrayList<>();

         public List<Integer> inorderTraversal1(TreeNode root) {
             if(root == null) return rst1; 
             helper(root);
             return rst1; 
         }
         // Add all the left side values to the rst, 
         private void helper(TreeNode root){
             if(root!= null && root.left == null && root.right ==null) {
                 rst1.add(root.val);
                 return; 
             }
                 TreeNode cur = root;
             while(cur.left!= null || !stack.empty()){
                 stack.push(cur);
                 while(cur.left!= null){
                     cur = cur.left; 
                     stack.push(cur);
                 }
                  while(!stack.isEmpty()){
                     TreeNode ri = stack.pop();
                     rst1.add(ri.val);
                     if(ri.right!=null) {
                         cur = ri.right;
                         break;
                     } 
                 }
              }       
         }
         public int mySqrt1(int x) {
        	    int l =1; 
        	    int r = x ;
        	        while(l<=r){
        	            int mid = l + (r-l)/2;
        	            if(mid*mid == x) return mid;
        	            else if (mid*mid > x) r = mid-1;
        	            else{
        	                if((mid+1)*(mid+1)>x) return mid;
        	                else l = mid+1; 
        	            }
        	        }
        	        return r; 
        	    }
         
         
         // private int rst = 0; 
         public int sumNumbers1(TreeNode root) {
             if(root == null) return 0;
             return helper(root,0);
         }
         // Pre: root is not null 
         // Return the sum value from the currrent root to the leaf 
         private int helper(TreeNode root, int add){
             if(root.left== null && root.right == null) {
                 return add*10 + root.val;
             }
             add = add*10 + root.val; 
             int left = 0; 
             if(root.left != null){
                 left=helper(root.left, add);
             }
             int right = 0; 
             if(root.right!= null){
                  right=helper(root.right, add);
             }
             return left+right; 
         }
         
         public int numSquares(int n) {
             Queue<Integer> que = new LinkedList<>();
             HashSet<Integer> dic = new HashSet<>(); 
             ArrayList<Integer> check = new ArrayList<>(); 
             que.add(n);
             int i = 1; 
             while(i*i<=n){
                 if(i*i == n) return 1; 
                 check.add(i);
                 que.add(n-i*i);
                 dic.add(n-i*i);
                 i++;
             }
             int count = 1; 
             while(!que.isEmpty()) {
            	 int size = que.size();
            	 count++; 
            	 for(int j = 0;j<size; j++) {
            		 int cur = que.poll();
            		 int k = 1;
            		 while(k*k<=cur) {
            			 if(cur-k*k ==0) return count; 
            			 else if(cur/k<k) continue; 
            			 else {
            				 if(!dic.contains(cur-k*k)) {
            					 dic.add(cur-k*k);
            					 que.add(cur-k*k);
            				 }
            			 }
            			 k++;
            		 }
            	 }
             }
			return -1;
         }
         
         public int removeDuplicates(int[] n) {
             if(n == null) return 0;
             int len = n.length; 
             if(len == 0) return 0; 
             int cur = 1;
             int pio = 1; 
             int local = 0; 
             while(pio<len){
                if(cur == 0 || n[pio]!=n[cur-1]) {
                	n[cur] = n[pio];
                	cur++;
                	pio++;
                	local = 0; 
                }else {
                	local++;
                	if(local!=2){
                		n[cur] = n[pio];
                		cur++;
                		pio++;
                	}else {
                		while(pio<len && n[pio] == n[cur-1]) pio++;
                		if(pio==len) return cur; 
                		n[cur] = n[pio];
                		cur++;
                		pio++;
                		local = 0; 
                	}
                }
             }
             for(int i =0; i<cur; i++) {
            	 System.out.print(n[i] + " ");
             }
             return cur; 
         }
         
         public int[] nextGreaterElement(int[] nums1, int[] nums2) {
             int len1 = nums1.length; 
             int len2 = nums2.length;
             int[] rst = new int[len1];
             for(int i =0; i<len1; i++){
                 int idx = binarySearch(nums2,nums1[i]);
                 for(int j = idx+1; j <len2; j++){
                     if(nums2[j]>nums1[i]){
                         rst[i] = nums2[j]; 
                         break;
                     }
                 }
                 if(rst[i] == 0) rst[i] = -1;
             }
             return rst; 
         }
         
 /** Pre: promise that i is in nums2
      * Pos: return the index of i in nums 
      */ 
        private int binarySearch(int[] arr, int value) {
        	  int k=0;
              for(int i=0;i<arr.length;i++){

                  if(arr[i]==value){
                      k=i;
                      break;
                  }
              }
          return k;
		}
        // 772. Basic Calculator III       
//        public int calculate(String s) {
//            if(s== null) return 0; 
//            int len = s.length();
//            if(len == 0) return 0; 
//            Stack<Character> stack = new Stack<>(); 
//            for(int i = 0; i< len; i++){
//                char c= s.charAt(i);
//                if(c=='(' || c=='+'|| c== '-' || c=='*' || c=='/') stack.push(c);
//                else if (c==' ') continue;
//                else if(c==')'){
//                    
//                }else if(Character.isDigit(c)){
//                    stack.push(c); 
//
//                    // char next = stack.peek(); 
//                	
//                }
//            }
//        }
       
//        // Pre: the given stack does not contains any ')' but at least one '('
//        // Pos: return the value contained inside the stack till we see the first '('
//        private int eval(Stack<Character> stack){
//            while(!stack.isEmpty() && stack.peek() != '('){
//                int i2 = getNum(stack); 
//                if(stack.peek() == '(') return i2;
//                else {
//                    char operator = stack.pop();
//                    int i1 = getNum(stack);
//                    if(operator == '+') {
//                    	stack.push(i1+i2);
//                    }else if(operator == '-') {
//                    	stack.push(i1-i2);
//                    }else if(operator == '+') {
//                    	stack.push(i1+i2);
//                    }else if(operator == '+') {
//                    	stack.push(i1+i2);
//                    }
//                    	
//                }
//                       
//                    }
//        }
//        
//        private int getNum(Stack<Character> stack) {
//        	int i2 = Character.getNumericValue(stack.pop()); 
//            if(!Character.isDigit(stack.peek())) return i2;
//            else {
//            	int time = 10; 
//                while(!stack.isEmpty() && Character.isDigit(stack.peek())) {
//                	i2 = Character.getNumericValue(stack.pop())*time + i2;
//                	time *=10 ; 
//                }
//            }
//			return i2;
//        }

//        public int[] nextGreaterElements(int[] A) {
//            Deque<Integer> deque = new LinkedList<>();
//            int len = A.length; 
//            if(len == 1) return new int[]{-1}; 
//            int[] rst = new int[len]; 
//            for(int i = len-2; i>=0; i--){
//                if(A[i]<deque.peekFirst()){
//                    rst[i] = deque.peekFirst(); 
//                    deque.addLast(A[i]);
//                }else{
//                    
//                }
//            }
//        }
        
        public String countOfAtoms(String f) {
            int len = f.length(); 
            if(len ==1) return f;
            HashMap<String,Integer> dic = new HashMap<>();
            int i =0;
            Stack<String> stack = new Stack<>(); 
            while( i<len){
                char c = f.charAt(i);
                // start of a new element 
                if(Character.isUpperCase(c)){
                    i =  helperCOA(f,i,dic);
                }else if(c == '('){
                    int j = helperCOA(f,i+1,dic,stack);
                    while(f.charAt(j)!=')'){
                        j = helperCOA(f,j,dic,stack); 
                    }
                    if(Character.isDigit(f.charAt(j))){
                    int num = Character.getNumericValue(f.charAt(j)); 
                    j++; 
                    while(j<len && Character.isDigit(f.charAt(j))){
                        num = num*10 + Character.getNumericValue(f.charAt(j));
                        j++; 
                    }
                        while(!stack.isEmpty()){
                            String put = stack.pop();
                            dic.put(put,dic.get(put)*num);
                        }
                       i = j; 
                }
            }
            }
            Set<String> keys = dic.keySet();
            Object[] key = keys.toArray();
            Arrays.sort(key);
            StringBuilder ret = new StringBuilder();
            for(int p =0; p<key.length; p++) {
            	ret.append(key[p]);
            	if(dic.get(key[p])!=1) ret.append(dic.get(key[p]));
            }
            return ret.toString();
        }
        // Input i: the char at index i is UpperCase letter, return the next index that we should look at ()
        // till we get to next upperCase letter/ '('
        private int helperCOA(String f, int i, HashMap<String, Integer> dic){
            int j = i+1; 
            int num = 0; 
                    while(j<f.length() && Character.isLowerCase(f.charAt(j))){
                        j++; 
                    }
                    String ele = f.substring(i,j);
             if(j == f.length()) {
            	 dic.put(ele, 1);
            	 return f.length();
             }   
            if(!Character.isDigit(f.charAt(j))){
                num =1;
            }else{
                 num = Character.getNumericValue(f.charAt(j)); 
                    j++; 
                    while(Character.isDigit(f.charAt(j))){
                        num = num*10 + f.charAt(j++)-'0';
                    }
                }
            dic.put(ele, dic.getOrDefault(ele, 0)+num);
//            if(dic.containsKey(ele)){
//                dic.put(ele, dic.get(ele)+num);
//            }else{
//                dic.put(ele,num);
//            }
            return j;
        }
        // Input i: the char at index i is '(', return the next index that we should look at (after its corresponding ')')
         private int helperCOA(String f, int i, HashMap<String, Integer> dic, Stack<String> stack){
            int j = i;
            //Sanity check 
            System.out.println("j "+j);
            if(!Character.isUpperCase(f.charAt(j))) {
            	System.out.println("Sanity Check is BAD");
            	return Integer.MIN_VALUE;
            }else {
            int num = 0; 
            j++;
                    while(j<f.length() && Character.isLowerCase(f.charAt(j))){
                        j++; 
                    }
                    String ele = f.substring(i,j);
                    stack.push(ele); 
            if(!Character.isDigit(f.charAt(j))){
                num =1;
            }else{
                 num = Character.getNumericValue(f.charAt(j)); 
                    j++; 
                    while(j<f.length() && Character.isDigit(f.charAt(j))){
                        num = num*10 + Character.getNumericValue(f.charAt(j));
                        j++; 
                    }
                } 
            if(dic.containsKey(ele)){
                dic.put(ele, dic.get(ele)+num);
            }else{
                dic.put(ele,num);
            }
            return j;
            }
        }
         
         public int minCostII(int[][] costs) {
             if(costs == null) return 0;
             int h = costs.length;
             if(h==0) return 0; 
             int c = costs[0].length;
             if(c==0) return 0; 
             int[][] dp = new int[h][c];
             int min = Integer.MAX_VALUE; 
             int minIdx = -1; 
             int secondMin = Integer.MAX_VALUE; 
             int sMIdx = -1; 
             for(int i =0; i<c; i++){
                 dp[0][i] = costs[0][i]; 
                 if(dp[0][i]<min){
                     secondMin = min;
                     sMIdx = minIdx; 
                     min = dp[0][i];
                     minIdx = i;
                 }else if(dp[0][i]>=min &&dp[0][i]<secondMin){
                     secondMin = dp[0][i];
                     sMIdx = i; 
                 }
             }
             for(int j=1; j<h; j++){
                 int lmin= Integer.MAX_VALUE; 
                 int lminIdx = -1;
                 int lsecondMin= Integer.MAX_VALUE;
                 int lsMIdx =-1; 
                 for (int i = 0; i< c; i++){
                    
                     if(i!=minIdx){
                         dp[j][i] = costs[j][i] + min;
                     }else{
                         dp[j][i] = costs[j][i] + secondMin;
                     }
                      if(dp[j][i]<lmin){
                     lsecondMin = lmin;
                     lsMIdx = lminIdx; 
                     lmin = dp[j][i];
                     lminIdx = i;
                 }else if(dp[j][i]>=lmin &&dp[j][i]<lsecondMin){
                     lsecondMin = dp[j][i];
                     lsMIdx = i; 
                 }
                     
                 }
                 min = lmin;
                 minIdx=lminIdx;
                 secondMin=lsecondMin;
                 sMIdx=lsMIdx;
             }
          return min; 
             
         }
         
         public int findNthDigit(int n) {
             if(n<=9) return n; 
             StringBuilder sb = new StringBuilder();
             // i is the next num that we should append to sb
             int i =1; 
             for(; i<10; i++)  sb.append(i); 
             while(sb.length()<n){
                 String s = Integer.toString(i); 
                 if(sb.length()+s.length()>=n){
                     return Character.getNumericValue(s.charAt(n-sb.length()-1));
                 }else{
                     sb.append(s);
                     i++;
                 }
             }
             return -1; 
         }
        
         
         // BFS Here we use int[]{x,y} to represent states in position x and with y liters of gas in the tank  ; for instance, {10,50} represents we are currently at position 10 with 50 liters of gas in the tank  
//         public int minRefuelStops(int target, int startFuel, int[][] stations) {
//             // trival case where we should return 
//             int len = stations.length; 
//             // no help 
//             if(len==0){
//                 if(startFuel>=target) return 0;
//                 return -1; 
//             }
//             // We first do us a favor by sorting the stations base on their position from low -> high 
//             Arrays.sort(stations, (int[] a, int[] b) -> a[0]-b[0]);
//             int s = 0; 
//             // We promise that no duplicates will be visited by only visiting stations that are farther 
//             Queue<int[]> bfs = new LinkedList<>();
////             HashSet<Integer> mark = new HashSet<>(); 
//             bfs.add(new int[] {0,startFuel});
////             while(stations[s][0]<=0) {
////            	 mark.add(s);
////            	 s++; 
////             }
//             int stops = 0; 
//             while(!bfs.isEmpty()){
//            	 int size = bfs.size();
//            	 for(int i =0; i<size; i++) {
//	                 int[] state = bfs.poll();
//	                 int reachable = state[0] + state[1]; 
//	                 System.out.println("reachable "+reachable);
//	                 if(reachable>=target) return stops;  
//	                 int search =0; 
//	                 while(search<len) {
//	                	 if( stations[search][0]>state[0] &&stations[search][0]<=reachable)
//	                	 bfs.add(new int[] {stations[search][0],state[1]-(stations[search][0]-state[0])+stations[search][1]});
//	                	 if(stations[search][0]>reachable) break; 
//	                	 search++; 
//	                 }
//            	 }
//            	 stops++;
//             } 
//             return -1; 
//         }
         
         // BFS Here we use int[]{x,y} to represent states in position x and with y liters of gas in the tank  ; for instance, {10,50} represents we are currently at position 10 with 50 liters of gas in the tank  
         public int minRefuelStops(int target, int startFuel, int[][] stations) {
             // trival case where we should return 
             int len = stations.length; 
             // no help 
             if(len==0){
                 if(startFuel>=target) return 0;
                 return -1; 
             }
             // We first do us a favor by sorting the stations base on their position from low -> high 
             Arrays.sort(stations, (int[] a, int[] b) -> a[0]-b[0]);
//             int s = 0; 
             // We promise that no duplicates will be visited by only visiting stations that are farther 
             Queue<int[]> bfs = new LinkedList<>();
             HashSet<Integer> mark = new HashSet<>(); 
             bfs.add(new int[] {0,startFuel});
             int stops = 0; 
             while(!bfs.isEmpty()){
            	 int size = bfs.size();
            	 for(int i =0; i<size; i++) {
	                 int[] state = bfs.poll();
	                 int reachable = state[0] + state[1]; 
	                 System.out.println("reachable "+reachable);
	                 if(reachable>=target) return stops;  
	                 int search =0; 
	                 while(search<len) {
	                	 if(stations[search][0]<=state[0]) search++; 
	                	 else if (stations[search][0]>state[0] && stations[search][0]<=reachable) {
		                	 bfs.add(new int[] {stations[search][0],state[1]-(stations[search][0]-state[0])+stations[search][1]});
//		                	 mark.add(search); 
		                	 search ++; 
		                 }else break;
	                 }
            	 }
            	 
                 stops++;
             } 
             return -1; 
         }
         
//       // BFS Here we use int[]{x,y} to represent states in position x and with y liters of gas in the tank  ; for instance, {10,50} represents we are currently at position 10 with 50 liters of gas in the tank  
//       public int minRefuelStops(int target, int startFuel, int[][] stations) {
//           // trival case where we should return 
//           int len = stations.length; 
//           // no help 
//           if(len==0){
//               if(startFuel>=target) return 0;
//               return -1; 
//           }
//           // We first do us a favor by sorting the stations base on their position from low -> high 
//           Arrays.sort(stations, (int[] a, int[] b) -> a[0]-b[0]);
//           int s = 0; 
//           // We promise that no duplicates will be visited by only visiting stations that are farther 
//           Queue<int[]> bfs = new LinkedList<>();
//           HashSet<Integer> mark = new HashSet<>(); 
//           bfs.add(new int[] {0,startFuel});
//           int stops = 0; 
//           while(!bfs.isEmpty()){
//          	 int size = bfs.size();
//          	 for(int i =0; i<size; i++) {
//	                 int[] state = bfs.poll();
//	                 int reachable = state[0] + state[1]; 
//	                 System.out.println("reachable "+reachable);
//	                 if(reachable>=target) return stops;  
//	                 int search =s; 
//	                 while(search<len && !mark.contains(search) &&stations[search][0]<=reachable) {
////	                	 mark.add(search);
//	                	 bfs.add(new int[] {stations[search][0],state[1]-(stations[search][0]-state[0])+stations[search][1]});
//	                	 search++; 
//	                 }
//          	 }
//          	 stops++;
//           } 
//           return -1; 
//       }
         
//         public String simplifyPath(String path) {
//        	 Deque<String> rst = new LinkedList<>(); 
//        	 String arr[] = path.split("/"); 
//             int len = arr.length; 
//             int i = 0; 
//             while(i<len){
//            	 if(arr[i] == ".") i++;
//            	 else if (arr[i] == "..") {
//            		 if(!rst.isEmpty()) rst.pop(); 
//            	 }else if(arr[i] == " ") i++;
//            	 else {
//            		 rst.push(arr[i]);
//            	 }
//             }
//             StringBuilder ret = new StringBuilder();
//             ret.append("/");
//             while(! rst.isEmpty()) {
//            	 ret.append(c)
//             }
//         }
         public String simplifyPath(String path) {
        	 Deque<String> rst = new LinkedList<>(); 
        	 String arr[] = path.split("/"); 
             int len = arr.length; 
             int i = 0; 
             while(i<len){
            	 if(arr[i].equals(".") || arr[i].equals("") || arr[i].equals(" ")) {
            		 i++;
            	 }
            	 else if (arr[i].equals("..") ) {
            		 if(!rst.isEmpty()) rst.removeFirst(); 
            		 i++;
            	 }else {
            		 rst.addFirst(arr[i]);
            		 i++;
            	 }
             }
             StringBuilder ret = new StringBuilder();
             if(rst.isEmpty()) {
                 ret.append("/");
                 return ret.toString(); 
             }
             while(! rst.isEmpty()) {
            	 ret.append("/");
            	 ret.append(rst.removeLast());
             }
             return ret.toString();
         }
        
         public int dominantIndex(int[] nums) {
             int len = nums.length; 
             if(len == 1) return 0; 
             Deque<Integer> deque = new LinkedList<>();
             deque.addFirst(0); 
             for(int i =1; i<len ; i++){
                 if(deque.isEmpty()) deque.addFirst(i);
                 else{
                     int l = deque.peekFirst();
                     int r = deque.peekLast();
                     if(nums[i]>= 2*nums[r]){
                         deque = new LinkedList<>();
                         deque.addFirst(i);
                     }else if (2*nums[i]<=nums[r] || nums[i] == nums[r] || nums[i] == nums[l]) continue;
                     else if(nums[i]>nums[r]){
                         deque.addLast(i); 
                         while(2*nums[deque.peekFirst()]<=nums[i]) deque.removeFirst(); 
                     }else {
                         deque.addFirst(i);
                     }
//                      else if(if(nums[i]<nums[l])){
                         
//                      }
                     // 其实最后一个不用else if 直接用else就可以了 因为它就是最后一种情况了
                 }
             }
             if(deque.size() >1 || deque.isEmpty()) return -1;
             else return deque.peekFirst(); 
         }
         
//         public String minWindow1(String s, String t) {
//             if(s == null || (s.length() <t.length())) return ""; 
//             int[] ddic = new int[256];
//             for(int i=0; i< t.length(); i++){
//                 ddic[t.charAt(i)] +=1; 
//             }
//             int l = 0;
//             int r= 0; 
//             int count = t.length(); 
//             String rst = ""; 
//             // [l,r)
//             while(l<s.length()){
//                 int[] dic = ddic.clone(); 
//                 while(l<s.length() && dic[s.charAt(l)]==0) l++;
//                 if(l==s.length()) {
//                	 System.out.println("Here1");
//                	 return rst; 
//                 }
//                 dic[s.charAt(l)]--;
//                 count --; 
//                 System.out.println("l: "+l);
//                 r = l+1;
//                 while(count >0 && r<s.length()){
//                     if(dic[s.charAt(r)]>0) {
//                         dic[s.charAt(r++)]--;
//                         count --;
//                     }else r++; 
//                 }
//                 System.out.println("r: "+r);
//
////                 if(r==s.length()){
////                	 System.out.println("Here2");
////                     if(count ==0) {
////                         if(rst.equals("") || r-l<rst.length())  return s.substring(l,r); 
////                     }
////                     else return rst; 
////                 }
//                 if(count ==0 && (rst.equals("") || r-l<rst.length())) {
//                	 System.out.println("why");
//                	 rst= s.substring(l,r); 
//                 }
//                 count = t.length(); 
//                 l = l+1;
//                 System.out.println("l2: "+l);
//             }
//             
//             return rst;
//         }
         
         public String minWindow1(String s, String t) {
             if(s == null || (s.length() <t.length())) return ""; 
             int[] ddic = new int[256];
             for(int i=0; i< t.length(); i++){
                 ddic[t.charAt(i)] +=1; 
             }
             int l = 0;
             int r= 0; 
             int count = t.length(); 
             String rst = ""; 
             // [l,r)
             while(r<s.length()){
                 // if(count == 0){
                     // if(rst == "" || r-l<rst.length()) rst = s.substring(l,r);   
                 // }
                 int[] dic = ddic.clone(); 
                 while(l<s.length() && dic[s.charAt(l)]==0) l++;
                 if(l==s.length()) return rst; 
                 dic[s.charAt(l)]--;
                 count --; 
                 r = l+1;
                 while(count >0 && r<s.length()){
                     if(dic[s.charAt(r)]>0) {
                         dic[s.charAt(r++)]--;
                         count --;
                     }else dic[s.charAt(r++)]--;
                 }
                 if(r==s.length()){
                     while(count ==0) {
                         if(rst.equals("") || r-l<rst.length())  rst= s.substring(l,r);
                         System.out.println("l before : "+l);
                         System.out.println("0:" +dic[s.charAt(l)]);
                         if(dic[s.charAt(l++)]++==0) {
                        	 System.out.println("safe");
                        	 count ++; 
                         }
                         System.out.println("0:" +dic[s.charAt(l)]);
                         System.out.println("l after: "+l);
                     }
                     break;
                 }
                 if(rst.equals("") || r-l<rst.length()) rst= s.substring(l,r); 
                 count = t.length(); 
                 l = l+1;
             }
             return rst;
         }
         
         public int compareVersion(String version1, String version2) {
             String[] arr1 = version1.split("\\.");
             String[] arr2 = version2.split("\\.");
             int i = 0; 
             while(i<arr1.length && i<arr2.length){
//                 if(arr1[i].compareTo(arr2[i])>0) return 1;
//                 else if (arr1[i].compareTo(arr2[i])<0) return -1;
            	 if(Integer.parseInt(arr1[i])>Integer.parseInt(arr2[i])) return 1;
                 else if (Integer.parseInt(arr1[i])<Integer.parseInt(arr2[i])) return -1;
                 else i++;
             }
             if(i== arr1.length && i== arr2.length) return 0;
             if(i== arr1.length) {
            	 for(int k=i; k<arr2.length; k++) {
            		 if(Integer.parseInt(arr2[k])!=0) return -1;
            	 }
            	 return 0; 
             }  
             if(i== arr2.length) {
            	 for(int k=i; k<arr1.length; k++) {
            		 if(Integer.parseInt(arr1[k])!=0) return 1;
            	 }
            	 return 0; 
             }
             return 0; 
         }
         
         public int kthGrammar(int n, int k) {
             return helperKG(n,k-1); 
         }
         
         private int helperKG(int n, int k){
             if(n==1) return 0;
             if(n==2) return k;
             int N =1; 
             for(int i =0; i<n-1; i++){
                  N*=2; 
             }
             if(k<N/2){
                 return helperKG(n-1,k);
             }else{
                 return (1- helperKG(n-1,k-N/2)); 
             }
         }
         
         public void nextPermutation(int[] nums) {
             if(nums== null || nums.length <= 1) return;
             int len = nums.length; 
             int r= len-1;
             while(r>0 ){
                 if(nums[r]>nums[r-1]){
                     for(int i = len-1; i>=r; i--){
                         if(nums[i]>nums[r-1]){
                             int tep = nums[r-1];
                             nums[r-1] = nums[i];
                             nums[i] = tep; 
                             return; 
                         }
                     }
                     
                  }else r--; 
             }
             if(r==0) {
                 Arrays.sort(nums);
                 return; 
             }
         }
         
         public List<String> summaryRanges1(int[] nums) {
             List<String> rst = new ArrayList<>(); 
             Stack<Integer> time = new Stack<>();
             int len = nums.length;
             if(len ==0) return rst; 
//             if(len ==1) {
//                 rst.add("" + nums[0] );
//                 return rst; 
//             }
             time.add(nums[0]);
             for(int i =1; i<len; i++){
                 if(nums[i]==time.peek()+1){
                     if(time.size()>1) time.pop();
                     time.add(nums[i]);
                 }else if(nums[i]>time.peek()+1){
                     int end = time.pop();
                     if(time.isEmpty()){
                         rst.add(""+end); 
                     }else{
                         int start = time.pop();
                         rst.add(""+ start +"->" + end); 
                     }
                     time.add(nums[i]);
                 }
                 // last one is actually else 
             }
             if(!time.isEmpty()){
            	 System.out.println("hi");
                 int end = time.pop();
                 if(time.isEmpty()){
                         rst.add(""+end); 
                     }else{
                     int start = time.pop();
                     rst.add(""+ start +"->" + end); 
                 }
             }
             return rst; 
         }
         
   	  public String getPermutation(int n, int k) {
	        if(n==1) return "1"; 
	        List<String> dic = new ArrayList<>();
	        StringBuilder sb = new StringBuilder();
	        backTrack(dic,n,sb);
	        dic.sort((String a, String b)-> a.compareTo(b));
	        return dic.get(k-1); 
	    }
	    
	    private void backTrack(List<String> dic,int n, StringBuilder sb){
	        if(sb.length() == n) dic.add(sb.toString());
	        else{
	            for(int i =1; i<=n; i++){
	                if(sb.toString().indexOf((char)'0'+i)==-1){
	                    sb.append((char)'0'+i);
	                    backTrack(dic,n,sb); 
	                    sb.deleteCharAt(sb.length()-1);
	                }
	            }
	        }
	    }
         
		public ArrayList<Integer> maxset(ArrayList<Integer> a) {
		    if(a == null || a.size()<1) return a; 
		    ArrayList<Integer> rst = new ArrayList<>();
		    int l = 0;
		    long max = 0; 
		    int r  = 0; 
		    long sum = 0; 
		    while(r<a.size()){
		        if(a.get(r)>=0) r++; 
		        else{
		            sum = 0; 
		            for(int i =l; i<r; i++){
		            sum += a.get(i); 
		            }
		        if(sum>=max) {
		            if(sum>max || (sum==max && r-l>rst.size())){
		                rst = new ArrayList<>();
		                for(int i =l; i<r; i++){
		                rst.add( a.get(i)); 
		                }
		             }
		            max = sum;
		        }
		        l = r+1; 
		        r++;
		        }
		    }
		    sum = 0; 
		            for(int i =l; i<r; i++){
		            sum += a.get(i); 
		            }
		        if(sum>=max) {
		        	System.out.println("why");
		            if(sum>max || (sum==max && r-l>rst.size())){
			        	System.out.println("why??");
		                rst = new ArrayList<>();
		                for(int i =l; i<r; i++){
		                rst.add(a.get(i)); 
		                }
		             }
		            max = sum;
		        }
		    return rst; 
		}
		public long reverse(long a) {
		    long ans = 0; 
		    for(int i =0; i<32; i++){
		        if((a & 1) == 1){
		            ans+=1; 
		            
		        }
		        ans <<=1;
		            a>>=1;
		    }
		    return ans;
		    
		}
		public class mySolution {
			public ArrayList<String> prettyJSON(String a) {
			    ArrayList<String> rst = new ArrayList<>(); 
			    if(a== null || a.length()==0) return rst;
			    int idx =0;
			    int ident = 0; 
			    while(idx<a.length()){
			        char c = a.charAt(idx);
			        if(c=='{' || c=='['){
			            StringBuilder id= new StringBuilder();
			            for(int i=0; i<ident; i++){
			                 id.append("\t");
			            }
			            id.append(c);
			            rst.add(id.toString()); 
			        }else if(c=='}' || c==']'){
			            StringBuilder id= new StringBuilder();
			            for(int i=0; i<ident-1; i++){
			                 id.append("\t");
			            }
			            id.append(c);
			            rst.add(id.toString()); 
			        }else{
			            StringBuilder id = new StringBuilder();
			            for(int i=0; i<ident; i++){
			                 id.append("\t");
			            }
			            while(idx<a.length()){
			                char cur = a.charAt(idx);
			                if(cur!='{' && cur!='}' && cur!='[' &&cur!=']'){
			                    id.append(cur); 
			                    idx++; 
			                }else break; 
			            }
			            rst.add(id.toString()); 
			        }
			    }
			    return rst; 
			}
		}
		
		public class Solution {
			public ArrayList<String> prettyJSON(String A) {
			    ArrayList<String> res = new ArrayList<>();
			    StringBuilder str = new StringBuilder();
			    int n = A.length();
			    int tabs = 0;
			    
			    for (int i = 0; i < n; ) {
			        
			        i = skipSpace(A, i);
			        
			        if (i >= n)
			            break;
			        
			        str = new StringBuilder();
			        char c = A.charAt(i);
			        
			        if (delimiter(c)) {
			            
			            if (isOpenBracket(c)) {
		    	            for (int j = 0; j < tabs; j++)
		    	                str.append("\t");	                
			                tabs++;
			            } else if (isClosedBracket(c)) {
			                tabs--;
		    	            for (int j = 0; j < tabs; j++)
			                    str.append("\t");
			            }
			            
			            str.append(c);
			            i++;
			            
		    	        if (i < n && canAdd(A.charAt(i))) {
		    	            str.append(A.charAt(i));
		    	            i++;
			            }
			            
			            res.add(str.toString());
			            
			            continue;
			        }
			        
			        while (i < n && !delimiter(A.charAt(i))) {
			            str.append(A.charAt(i));
			            i++;
			        }
			        
			        if (i < n && canAdd(A.charAt(i))) {
			            str.append(A.charAt(i));
			            i++;
			        }
			        
			        StringBuilder strB = new StringBuilder();
			        
			        for (int j = 0; j < tabs; j++)
			            strB.append("\t");
			        
			        strB.append(str);
			        res.add(strB.toString());
			    }
			    
			    return res;
			    
			}
			
			public boolean canAdd(char c) {
			    if (c == ',' || c == ':')
			        return true;
			        
			    return false;
			}
			
			public boolean delimiter(char c) {
			    if (c == ',' || isOpenBracket(c) || isClosedBracket(c))
			        return true;
			    return false;
			}
			
			public boolean isOpenBracket(char c) {
			    if (c == '[' || c == '{')
			        return true;
			    return false;
			}
			
			public boolean isClosedBracket(char c) {
			    if (c == ']' || c == '}')
			        return true;
			    return false;
			}
			
			public int skipSpace(String A, int i) {
			    int n = A.length();
		        while (i < n && A.charAt(i) == ' ')
		            i++;
		        return i;
			}
			
			
		}
		
	
	    
	
	    
		   // Complete the reverse function below.
	    static String reverse1(String expr) {
	        if(expr == null || expr.length() <=1) return expr;
	        int len = expr.length(); 
	        Stack<String> stack = new Stack<>();
	        int i = 0; 
	        while(i<len){
	            char c = expr.charAt(i); 
	            if(c=='/' || c=='*' || c=='+') {
	            	stack.push(Character.toString(c)); 
	            	i++; 
	            }else if (c=='-' && (i==0 || Character.isDigit(expr.charAt(i-1)) || expr.charAt(i-1)=='.')){
//	            	if(i-1>=0 && ) {
	            		stack.push(Character.toString(c));
	            		i++; 
//	            	}
	            }else {
	                StringBuilder cur = new StringBuilder(); 
	                if(c=='-') {
	                	cur.append(c); 
	                	i++;
	                }
	                while(i<len && (Character.isDigit(expr.charAt(i)) || expr.charAt(i)=='.')){
	                    cur.append(expr.charAt(i)); 
	                    i++; 
	                }
	                stack.push(cur.toString());
	            }
	        } 
	        StringBuilder sb = new StringBuilder();
	        while(!stack.isEmpty()) {
	        	sb.append(stack.pop());
	        }
	        return sb.toString();

	    }
	    
//      // Complete the reverse function below.
//      static String reverse(String expr) {
//          if(expr == null || expr.length() <=1) return expr;
//          int len = expr.length(); 
//          Stack<String> stack = new Stack<>();
//          int i = 0; 
//          while(i<len){
//              char c = expr.charAt(i); 
//              if(c=='/' || c=='*' || c=='+') {
//                  stack.push(Character.toString(c)); 
//                  i++; 
//              }else if (c=='-'){
//                  if(i-1>=0 && Character.isDigit(expr.charAt(i-1))) {
//                      stack.push(Character.toString(c));
//                      i++; 
//                  }
//              }else{
//                  StringBuilder cur = new StringBuilder(); 
//                  if(c=='-') {
//                      cur.append(c); 
//                      i++;
//                  }
//                  while(i<len && (Character.isDigit(expr.charAt(i)) || expr.charAt(i)=='.')){
//                      cur.append(expr.charAt(i)); 
//                      i++; 
//                  }
//                  stack.push(cur.toString());
//              }
//          } 
//          StringBuilder sb = new StringBuilder();
//          while(!stack.isEmpty()) {
//              sb.append(stack.pop());
//          }
//          return sb.toString();


        // Complete the reverse function below.
      static String reverse(String expr) {
          if(expr == null || expr.length() ==0) return expr;
          int len = expr.length(); 
          Stack<String> stack = new Stack<>();
          int i = 0; 
          while(i<len){
              char c = expr.charAt(i); 
              if(c=='/' || c=='*' || c=='+') {
                  stack.push(Character.toString(c)); 
                  i++; 
              }else if (c=='-' && i!=0 && (Character.isDigit(expr.charAt(i-1)) || expr.charAt(i-1)=='.')){
//                  if(i-1>=0 && ) {
                      stack.push(Character.toString(c));
                      i++; 
//                  }
              }else {
                  StringBuilder cur = new StringBuilder(); 
                  if(c=='-') {
                      cur.append(c); 
                      i++;
                  }
                  while(i<len && (Character.isDigit(expr.charAt(i)) || expr.charAt(i)=='.')){
                      cur.append(expr.charAt(i)); 
                      i++; 
                  }
                  stack.push(cur.toString());
              }
          } 
          StringBuilder sb = new StringBuilder();
          while(!stack.isEmpty()) {
              sb.append(stack.pop());
          }
          return sb.toString();

      }
      
      // Complete the latestStudent function below.
	    static String latestStudent(List<List<String>> attenData) {
	    	// Date to avg lateness 
	    	HashMap<String, Integer> dic = new HashMap<>();
	    	HashMap<String, Integer> count = new HashMap<>();
	    	for(int i=0; i<attenData.size(); i++) {
	    		String key =attenData.get(i).get(0); 
	    		if(dic.containsKey(key)) {
	    			dic.put(key, dic.get(key)+Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2)))); 
	    		}else{
	    			dic.put(key, Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))));
	    		}
	    		if(count.containsKey(key)) {
	    			count.put(key, count.get(key)+1);
	    		}else count.put(key,1); 
	    	}
	    	for(String key: dic.keySet()) {
	    		dic.put(key, dic.get(key)/count.get(key));
	    	}
	    	
	    	// students' name to their relative lateness
	    	HashMap<String, Integer> students = new HashMap<>();
	    	for(int i=0; i<attenData.size(); i++) {
	    		String key = attenData.get(i).get(1); 
	    		if(students.containsKey(key)) {
	    			students.put(key, dic.get(key)+Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))-dic.get(attenData.get(i).get(0)))); 
	    		}else{
	    			students.put(key, Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))-dic.get(attenData.get(i).get(0))));
	    		}
	    	}
	    	
	    	int max = 0; 
	    	String rst = ""; 
	    	for(String key: students.keySet()) {
	    		if(students.get(key)>max) {
	    			max = students.get(key);
	    			rst = key; 
	    		}else if(students.get(key)==max) {
	    			if(key.compareTo(rst)<0) rst = key; 
	    		}
	    	}
	    	return rst; 
	    }
	    
//	      // Complete the latestStudent function below.
//        static String latestStudent(List<List<String>> attenData) {
//            // Date to avg lateness 
//            if(attenData==null || attenData.size()==0) return "";
//            HashMap<String, Integer> dic = new HashMap<>();
//            HashMap<String, Integer> count = new HashMap<>();
//            for(int i=0; i<attenData.size(); i++) {
//                String key =attenData.get(i).get(0); 
//                if(dic.containsKey(key)) {
//                    dic.put(key, dic.get(key)+Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2)))); 
//                }else{
//                    dic.put(key, Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))));
//                }
//                if(count.containsKey(key)) {
//                    count.put(key, count.get(key)+1);
//                }else count.put(key,1); 
//            }
//            for(String key: dic.keySet()) {
//                dic.put(key, dic.get(key)/count.get(key));
//            }
//            
//            // students' name to their relative lateness
//            HashMap<String, Integer> students = new HashMap<>();
//            for(int i=0; i<attenData.size(); i++) {
//                String key = attenData.get(i).get(1); 
//                if(students.containsKey(key)) {
//                    students.put(key, students.get(key)+Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))-dic.get(attenData.get(i).get(0)))); 
//                }else{
//                    students.put(key, Math.max(0, Integer.parseInt(attenData.get(i).get(3))-Integer.parseInt(attenData.get(i).get(2))-dic.get(attenData.get(i).get(0))));
//                }
//            }
//            
//            int max = 0; 
//            String rst = ""; 
//            for(String key: students.keySet()) {
//                if(students.get(key)>max) {
//                    max = students.get(key);
//                    rst = key; 
//                }else if(students.get(key)==max) {
//                    if(key.compareTo(rst)<0) rst = key; 
//                }
//            }
//            return rst; 
//        }
	    
	    public int myAtoi(String str) {
	        long rst = 0;
	        if(str==null) return (int)rst; 
	        int i =0 ; 
	        while(i<str.length()&&str.charAt(i)==' ') i++;
	        int sign = 1; 
	        boolean valid = false; 
	        while(i<str.length()){
	            char c = str.charAt(i); 
	            if(!valid && c == '-'){
	                sign = -1; 
	                valid = true;
	            }else if (!valid && c == '+'){
	                sign = 1; 
	                valid = true;
	            }else if(Character.isDigit(c)) {
	                rst = rst*10 + Character.getNumericValue(c); 
	                System.out.println("rst is "+rst);
	                int intR = (int)rst; 
	    	        if(intR != rst) {
	    	        	System.out.println("here");
	    	            if(sign >0) return Integer.MAX_VALUE;
	    	            else return Integer.MIN_VALUE;
	    	        }
	                valid = true;
	            }else{
	                if(valid ==false) return 0; 
	                else break; 
	            }
	            i++;    
	        }
	        return sign* (int)rst; 
	    }
	    
	    public int divide(int dividend, int divisor) {
			if(dividend==Integer.MIN_VALUE && divisor==-1) return Integer.MAX_VALUE;
	        if(dividend > 0 && divisor > 0) return divideHelper(-dividend, -divisor);
	        else if(dividend > 0) return -divideHelper(-dividend,divisor);
	        else if(divisor > 0) return -divideHelper(dividend,-divisor);
	        else return divideHelper(dividend, divisor);
	    }
	    
	    private int divideHelper(int dividend, int divisor){
	        // base case
	        if(divisor < dividend) return 0;
	        // get highest digit of divisor
	        int cur = 0, res = 0;
	        while((divisor << cur) >= dividend && divisor << cur < 0 && cur < 31) cur++;
	        res = dividend - (divisor << cur-1);
	        if(res > divisor) return 1 << cur-1;
	        return (1 << cur-1)+divide(res, divisor);
	    }
	    
	    public int firstMissingPositive(int[] nums) {
	        if(nums.length == 0) return 1; 
	        int max = nums[0];
	        ArrayList<Integer> check = new ArrayList<>();
	        for(int i =1; i<max ; i++){
	            check.add(i); 
	        }
	        for(int i = 1; i<nums.length; i++){
	            if(nums[i]<=0) continue;
	            else if(nums[i]<=max){
	                if(!check.isEmpty())
	                    check.remove(nums[i]); 
	            }else{
	                for(int j = max+1; j<nums[i]; j++){
	                    check.add(j);
	                }
	                max  = nums[i]; 
	            }
	        }
	        if(check.isEmpty()) return max+1;
	        else return check.get(0); 
	    }
	    
	    //这段代码的目标错了
//	    List<List<Integer>> rst = new ArrayList<List<Integer>>();
//	    public List<List<Integer>> permuteUnique(int[] nums) {
//	        Arrays.sort(nums); 
//	        backtrack(rst, new ArrayList<Integer>(), nums, 0);
//	        return rst; 
//	    }
//	    
//	    private void backtrack(List<List<Integer>> rst, List<Integer> cur, int[] nums, int start){
//	        if(start>nums.length) return; 
//	        if(cur.size() == nums.length){
//	            rst.add(new ArrayList<>(cur));
//	            int dup = start; 
//	            while(start< nums.length && nums[dup] == nums[start]) start++; 
//	            return; 
//	        }
//	        int i =start; 
//	        while(i<nums.length){
//	            cur.add(nums[i]);
//	            backtrack(rst, cur, nums, start+1);
//	            cur.remove(cur.size()-1); 
//	            i++; 
//	        }
//	    }
	    
	    public int minPathSum(int[][] grid) {
	        if(grid.length == 0 || grid[0].length == 0) return 0; 
	        int[][] dp = new int[grid.length][grid[0].length]; 
	        dp[grid.length-1][grid[0].length-1] = grid[grid.length-1][grid[0].length-1];
	        for(int i=grid.length-2; i>=0; i--){
	            dp[i][grid[0].length-1] = dp[i+1][grid[0].length-1] + grid[i][grid[0].length-1];
	        }
	        for(int i = grid[0].length-2; i>=0; i--){
	            dp[grid.length-1][i] = dp[grid.length-1][i+1] + grid[grid.length-1][i]; 
	        }
	        
	        for(int i = grid.length-2; i>=0; i--){
	            for(int j = grid[0].length-2; j>=0; j--){
	                dp[i][j] = Math.min(dp[i+1][j],dp[i][j+1]) + grid[i][j];
	            }
	        }
	    	printIntMatrix(dp);
	        return grid[0][0]; 
	    }
	    
	    public int numIslands(char[][] grid) {
	        int rst= 0;
	        for(int i = 0; i<grid.length; i++){
	            for(int j = 0; j<grid[0].length; j++){
	                if(grid[i][j]=='1'){
	                    rst++;
	                    grid[i][j] = '0'; 
	                    grid = dfs(grid, i, j); 
	                }
	            }
	        }
	        return rst; 
	        
	    }
	    
	    int[][] dir= {{1,1},{1,-1},{-1,-1},{-1,1}}; 
	    private char[][] dfs(char[][] grid, int x, int y){
	        
	        for(int[] d:dir){
	            if(x+d[0] < 0 || x+d[0]>=grid.length || y+d[1]<0||y+d[1] >= grid[0].length || grid[x+d[0]][y+d[1]] == 0) continue;
	            else{
	                grid[x+d[0]][y+d[1]] = '0';
	                grid = dfs(grid,x+d[0], y+d[1]); 
	            }
	        }
	        return grid; 
	    }
	    private void printMatrix(int[][] grid) {
	    	for(int i =0; i<grid.length; i++) {
	    		for(int j =0; j<grid[0].length; j++) {
	    			System.out.print(grid[i][j] + " ");
	    		}
	    		System.out.println();
	    	}
	    }
	    
	    public int superpalindromesInRange(String L, String R) {
	        int rst = 0; 
	        Long l = sqrt(Long.parseLong(L)); 
	        Long r = sqrt(Long.parseLong(R)); 
	        if(l != Long.parseLong(L)/l) l++;

	        for(Long i = l ; i<=r; i++){
	            if(isPN(i) && isPN(i*i)) rst++; 
	            
	        }
	        return rst; 
	    }
	    
	    public Long sqrt(Long x) {
	    if (x == 0l)
	        return 0l;
	    Long left = 1l, right = Long.MAX_VALUE;
	    while (true) {
	        Long mid = left + (right - left)/2l;
	        if (mid > x/mid) {
	            right = mid - 1;
	        } else {
	            if (mid + 1l > x/(mid + 1l))
	                return mid;
	            left = mid + 1l;
	        }
	        }
	    }

	    
	    public boolean isPN(Long x) {
//	    if (x<0l || (x!=0l && x%10l==0l)) return false;
	    if  (x!=0l && x%10l==0l) return false;

	    Long rev = 0l;
	    while (x>rev){
	    	rev = rev*10l + x%10l;
	    	x = x/10l;
	    }
	    return (x.equals(rev) || x.equals(rev/10l));
	    }
	    
//	    public int superpalindromesInRange(String L, String R) {
//	        int l = (int) Math.sqrt(Long.parseLong(L));
//	        int r = (int) Math.sqrt(Long.parseLong(R));
//	      
//	        List<Integer> allP = new ArrayList<>();
//	        help(allP,l,r);
//	        int cnt = 0;
//	        for(int p : allP){
//	            long k = (long)p*p;
//	            if(isP(k)) {
//	            	System.out.println("p "+p);
//	            	System.out.println("k "+k);
//	            	cnt++;
//	            	
//	            } 
//	        }
//	        return cnt;
//	    }
	   
	    
//	    int createP(int input, int b, int isOdd) { 
//	        long n = input; 
//	        long palin = input; 
//	        if (isOdd == 1) 
//	            n /= b; 
//	        while (n > 0) { 
//	            palin = palin * b + (n % b); 
//	            n /= b; 
//	        } 
//	        return (int)palin; 
//	    } 
//
//	    void help(List<Integer> allP,int k,int n) { 
//	        int number; 
//	        for (int j = 0; j < 2; j++) { 
//	            int i = 1; 
//	            while ((number = createP(i, 10, j % 2)) <= n) { 
//	                if(number >= k)
//	                    allP.add(number); 
//	                i++; 
//	            } 
//	        }
//	    } 
	    
	    
	        public int[] exclusiveTime(int n, List<String> logs) {
	            int[] rst = new int[n];
	            Stack<Integer> stack = new Stack<>();
	            if(logs == null || logs.size() == 0) return rst; 
	            String s0 = logs.get(0);  
	            String[] str0 = s0.split(":");            
	            stack.push(Integer.parseInt(str0[2])); 
                boolean start = false; 
	            for(int i =1; i<logs.size(); i++){
	                String prev = logs.get(i-1);
	                String[] pstr = prev.split(":");            
	                int pidx = Integer.parseInt(pstr[0]); 
	                String s = logs.get(i);
	                String[] str = s.split(":");  
	                int idx = Integer.parseInt(str[0]); 
	               
	                if(str[1].equals("start")){
	                	start = false; 
	                    // if(!stack.isEmpty())
	                	if(pstr[1].equals("start"))
	                		rst[pidx] += Integer.parseInt(str[2]) - stack.peek(); 
	                	else 
	                		rst[pidx] += Integer.parseInt(str[2]) - stack.peek()-1; 
    		
	                    stack.push(Integer.parseInt(str[2])); 
	                }else{
	                    if(!start){
	                        rst[idx] += Integer.parseInt(str[2]) - stack.peek()+1; 
	        	            System.out.println(i + " "+Arrays.toString(rst)); 
	                        start = true; 
	                        stack.push(Integer.parseInt(str[2]));
	                    }else{
	                    	System.out.println("why "+ Integer.parseInt(str[2]) + "& "+stack.peek() );
	                        rst[idx] += Integer.parseInt(str[2]) - stack.peek();
	        	            System.out.println(i + " "+Arrays.toString(rst)); 

	                        stack.push(Integer.parseInt(str[2]));
	                    }
	                }
	            }
	            System.out.println(Arrays.toString(rst)); 
	            return rst;    
	        }
	        
	        public int compress(char[] chars) {
	            // l is the idx of the next writing position 
	            int l = 0;
	            // k is the idx of the 1st char that we are currently compressing 
	            int k = 0; 
	            // r is the pioneer idx
	            int r = 0;
	            while(r<chars.length){
	                while(r < chars.length && chars[r] == chars[k]) r++;
	                if(r-k==1){
	                    chars[l] = chars[k];
	                    l+=1; 
	                    k= r; 
	                }else{
	                    chars[l] = chars[k];
	                    l+=1;
	                    int num = r-k;
	                    System.out.println("num "+num);
	                        String str= Integer.toString(num);
		                    System.out.println("str "+str);
	                        char[] arr = str.toCharArray(); 
	                        System.out.println("arr "+Arrays.toString(arr));
	                        for(int i = 0; i<arr.length; i++){
	                            chars[i+l] = arr[i]; 
	                            
	                        }
	                        l+=arr.length; 
	                    k = r; 
	                }
	            }
	            // if(l >= chars.length-1) return l; 
	            // if(r-k==1) {
	            //     chars[l++] = chars[k]; 
	            // }else{
	            //         chars[l++] = chars[k]; 
	            //         int num = r-k;
	            //         while(num>0){
	            //             chars[l++] = (char)('0'+ num%10); 
	            //             num /=10; 
	            //         }
	            // }
            	System.out.println("sol");
	            for(int o=0; o<l; o++) {
	            	System.out.print(chars[o] + " ");
	            }
	            return l; 
	        }
	        
	   private int[] please(int[] house, int[] store) {
	        // write your code in Java SE 8
	        int[] rst = new int[house.length]; 
		        Arrays.sort(house);
		        
		        Arrays.sort(store);
		        for(int i =0; i< house.length; i++) {
		        	int idx = Arrays.binarySearch(store, house[i]); 
		        	if(idx<0) {
		        		idx = -(idx+1);
		        		if(idx == 0) rst[i] = store[0];
		        		else if(idx == store.length-1) rst[i] = store[store.length-1]; 
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

	    
	    public int smallestRangeII(int[] A, int K) {
	        if(A.length<=1) return 0;
	        Arrays.sort(A);
	        int min = A[0];
	        int max = A[A.length-1];
	        if(max-min == 2*K) return 0;
	        int dif = Math.min(Math.abs(max-min-2*K), max-min);
	        for(int i =1; i<A.length-1; i++){
	            if(A[i]-2*K<min && A[i]+2*K>max){
	                int why1 = Math.min(Math.abs(max-(A[i]-K) - 2*K), max-A[i]-K);
	                System.out.println("why1 "+why1);
	                int why2 = Math.min(Math.abs(A[i]+K - min - 2*K), A[i]+K - min );
	                System.out.println("why2"+why2);
	                int w = Math.min(why1, why2);
	                System.out.println("w "+w);
	                dif = Math.max(dif, w); 
	                System.out.println("dif "+dif);

	            }
	        }
	        return dif; 
	    }
	    
//	    public int partitionDisjoint(int[] A) {
//	        int l = 0, r = A.length -1;
//	        int max = Integer.MIN_VALUE, min = Integer.MAX_VALUE; 
//	            
//	        while(max<=min && l <=r){
//	            System.out.println("l "+l);
//	            System.out.println("r "+r);
//	            if(A[l]>max) max = A[l];
//	            if(A[r]<min) min = A[r];
//	            l++;
//	            r--; 
//	      
//
//	            // if(max>min){
//	            //     r--; 
//	            // }else{
//	            //     l++; 
//	            // }
//	        }
//	        if(r+2>=A.length) return 1; 
//	        return r+2; 
//	        
//	    }
	    public int partitionDisjoint(int[] A) {
	        int l = 0, r = A.length -1;
	        int[] maxArr = new int[A.length];
	        int[] minArr = new int[A.length]; 
	        int max = A[0], min = A[A.length -1];
	        maxArr[0] = max;
	        minArr[A.length-1] = min;
	        for(int i =1; i<A.length; i++){
	            if(A[i]>max){
	                max = A[i];
	            }
	                   maxArr[i] = max;
	        }
	        for(int i = A.length-2; i>=0; i--){
	            if(A[i] < min) min = A[i];
	            minArr[i] = min; 
	        }
		        for(int i=0; i<A.length-1; i++){
	                if(maxArr[i]<=minArr[i+1]) return i+1; 
	            }
//	 	        while(max<=min && l <=r){
//	 	            if(A[l]>max) max = A[l];
//	 	            if(A[r]<min) min = A[r];
//	 	            l++;
//	 	            r--; 
		      

//	 	            // if(max>min){
//	 	            //     r--; 
//	 	            // }else{
//	 	            //     l++; 
//	 	            // }
//	 	        }
//	 	        if(r+2>=A.length) return 1; 
//	 	        return r+2;
	        return 1; 
	    }
	      
	    private void printSet(Set<Character> A) {
	    	for(char c:A) {
	    		System.out.print(c +" ");
	    	}
	    }
	    
	    public List<String> wordSubsets(String[] A, String[] B) {
	        HashMap<Character, Integer> Map = new HashMap<>(); 
	        for(String s:B){
	            HashMap<Character, Integer> map = new HashMap<>(); 
	            for(int i= 0; i< s.length(); i++){
	                if(map.containsKey(s.charAt(i))){
	                    map.put(s.charAt(i),map.get(s.charAt(i))+1);
	                }else{
	                    map.put(s.charAt(i),1);
	                }
	            }
	            for(char c: map.keySet()){
	                if(Map.containsKey(c)){
	                    Map.put(c,Math.max(Map.get(c),map.get(c))); 
	                }else{
	                    Map.put(c, map.get(c)); 
	                }
	            }               
	        }
	        
	        System.out.println("Map");
	        printSet(Map.keySet());
	        
	        List<String> rst  = new ArrayList<>();
	        for(String s: A){
	            HashMap<Character, Integer> Amap = new HashMap<>(); 
	            for(int i= 0; i< s.length(); i++){
	                if(Amap.containsKey(s.charAt(i))){
	                    Amap.put(s.charAt(i),Amap.get(s.charAt(i))+1);
	                }else{
	                    Amap.put(s.charAt(i),1);
	                }
	            }
	            System.out.println("Amap");
		        printSet(Amap.keySet());
//	            if(Amap.keySet().size() != Map.keySet().size()) continue; 
	            boolean add = true; 
	            for(char c: Map.keySet()){
	                if(!Amap.containsKey(c) || Amap.get(c)<Map.get(c)) {
	                    add = false;
	                    break;
	                }
	            }
	            if(add) rst.add(s);
	        }
	      
	        return rst; 
	    }
	    
	    
	    public int maxSubarraySumCircular(int[] A) {
	        int total =0, maxSum = -30000, curMax = 0, minSum = 30000, curMin =0;
	        for(int a: A){
	            curMax = Math.max(curMax+a, a);
	            maxSum = Math.max(maxSum, curMax); 
	            curMin = Math.min(curMin+a,a);
	            minSum = Math.min(minSum , curMin);
	            total +=a; 
	        }
//	        return Math.max(maxSum, total - minSum); 
	        return maxSum > 0 ? Math.max(maxSum, total - minSum ): maxSum; 
	    }
	       
//	    public double mincostToHireWorkers(int[] quality, int[] wage, int K) {
//	        // First need to calculate quality/wage array and get the first K ones
//	        double[] rate = new double[wage.length];
//	        HashMap<Double, Integer> map = new HashMap<>();
//	        System.out.println("quality " +Arrays.toString(quality) );
//	        System.out.println("wage " +Arrays.toString(wage) );
//
//	        for(int i =0 ; i< wage.length; i++){
//	            rate[i] = (double)quality[i]/(double)wage[i]; 
//	            System.out.println("this item "+rate[i]);
//	            map.put(rate[i],i); 
//	        }
//	        System.out.println("beforeSort " +Arrays.toString(rate) );
//	        
//	        Arrays.sort(rate); 
//	        System.out.println("afterSort"+Arrays.toString(rate) );
//
//	        double rst = 0; 
//	        double see = rate[wage.length-K]; 
//	        System.out.println("see "+ see);
//	        for(int i =wage.length-1; i>=wage.length-K; i--){
//	            // rst += wage[map.get(rate[i])] * rate[i]/see;
//	            rst += quality[map.get(rate[i])]/see;
//	        }
//	        return rst; 
//	    }
	    
	    
	    public double mincostToHireWorkers(int[] q, int[] w, int K) {
	        double[][] workers = new double[q.length][2];
	        for (int i = 0; i < q.length; ++i)
	            workers[i] = new double[]{(double)(w[i]) / q[i], (double)q[i]};
	        Arrays.sort(workers, (a, b) -> Double.compare(a[0], b[0]));
//	        for(int i =0; i<=1; i++) {
//	        	System.out.println("afterSort"+Arrays.toString(workers[i]) );
//	        }
	        
	        double res = Double.MAX_VALUE, qsum = 0;
	        PriorityQueue<Double> pq = new PriorityQueue<>();
	        for (double[] worker: workers) {
	            qsum += worker[1];
	            pq.add(-worker[1]);
	            if (pq.size() > K) qsum += pq.poll();
	            System.out.println("worker0 "+ worker[0]);
	            if (pq.size() == K) res = Math.min(res, qsum * worker[0]);
	        }
	        return res;
	    }
	    
	    
	    public boolean isPossible(int[] nums) {
	        Map<Integer, Integer> freq = new HashMap<>(), appendfreq = new HashMap<>();
	        for (int i : nums) freq.put(i, freq.getOrDefault(i,0) + 1);
	        for (int i : nums) {
	            if (freq.get(i) == 0) continue;
	            else if (appendfreq.getOrDefault(i,0) > 0) {
	                appendfreq.put(i, appendfreq.get(i) - 1);
	                appendfreq.put(i+1, appendfreq.getOrDefault(i+1,0) + 1);
	            }   
	            else if (freq.getOrDefault(i+1,0) > 0 && freq.getOrDefault(i+2,0) > 0) {
	                freq.put(i+1, freq.get(i+1) - 1);
	                freq.put(i+2, freq.get(i+2) - 1);
	                appendfreq.put(i+3, appendfreq.getOrDefault(i+3,0) + 1);
	            }
	            else return false;
	            freq.put(i, freq.get(i) - 1);
	        }
	        return true;
	    }
	    
	    
	    
	    public int reverseBits2(int n) {
	        int rst = 0; 
	        int i = 0; 
	        int cur = 0; 
	        while(n>0){
	            cur = n & 1; 
	            n >>>= 1; 
	            rst <<= 1;
	            rst += cur; 
	            i++;
	        }
	        for(int j = i; j<32; j++){
	            rst <<= 1;
	        }
	        return rst; 
	    }
	    
	    
	    public int minAddToMakeValid(String s) {
	        if(s.length() == 0) return 0;
	        if(s.equals("(") || s.equals(")")) return 1;
	        int len = s.length();
	        int[][] dp = new int[len][len];
	        for(int i= 0; i<len; i++){
	            dp[i][i] = 1; 
	        }
	        // for(int i =0; i<len; i++){
	        //     for(int j = i+1; j<len; j++){
	        //         int second = Integer.MAX_VALUE; 
	        //         if(s.charAt(i) == '(' && s.charAt(j) == ')'){
	        //             second = dp[i+1][j-1]; 
	        //         }
	        //         int first = Integer.MAX_VALUE;
	        //         for(int k = i; k<j; k++){
	        //             if(dp[i][k] + dp[k+1][j] < first) 
	        //                 first= dp[i][k] + dp[k+1][j]; 
	        //         }
	        //         dp[i][j] = Math.min(first, second); 
	        //     }
	        // }
//	        for(int j =1; j<len; j++){
//	            for(int i = 0; i<len-j; i++){
////	                int second = Integer.MAX_VALUE; 
//	                if(s.charAt(i) == '(' && s.charAt(i+j) == ')'){
//	                    dp[i][i+j] = dp[i+1][i+j-1]; 
//	                    continue; 
//	                }
//	                int first = Integer.MAX_VALUE;
//	                for(int k = i; k<i+j; k++){
//	                    if(dp[i][k] + dp[k+1][i+j] < first) 
//	                        first= dp[i][k] + dp[k+1][i+j]; 
//	                }
//	                dp[i][i+j] = first; 
//	            }
//	        }
            for(int j =1; j<len; j++){
	            for(int i = 0; i<len-j; i++){
	                int second = Integer.MAX_VALUE; 
	                if(s.charAt(i) == '(' && s.charAt(i+j) == ')'){
	                    second = dp[i+1][i+j-1]; 
	                }
	                
	                int first = Math.min(1+ dp[i+1][i+j], 1+dp[i][i+j-1]); 
	                if(i+2 < len && s.charAt(i) == '(' && s.charAt(i+1) == ')') first = Math.min(first, dp[i+2][i+j]);
	                if(i+j-2>=0 && s.charAt(i+j) == ')' && s.charAt(i+j-1) == '(') first = Math.min(first, dp[i][i+j-2]);
	                	
//	                		Integer.MAX_VALUE;
//	                for(int k = i; k<i+j; k++){
//	                    if(dp[i][k] + dp[k+1][i+j] < first) 
//	                        first= dp[i][k] + dp[k+1][i+j]; 
//	                }
	                dp[i][i+j] = Math.min(first, second); 
	            }
	        }
//	        printMatrix(dp);
	        return dp[0][len-1]; 
	    }
	    
	    
	    public int minAddToMakeValid1(String s) {
	    	Deque<Character> stack = new LinkedList<>();
	    	for(int i=0; i<s.length(); i++) {
	    		char c =s.charAt(i); 
	    		if( c=='(')
	    			stack.push(c);
	    		else {
	    			if(stack.peek() == '(')
	    				stack.pop();
	    			else if(stack.peekLast() == '(')
	    				stack.pollLast();
	    		}
	    	}
	    	return stack.size();
	    }
	    
	    
	    public int wordsTyping(String[] sentence, int rows, int cols) {
	        int rst = 0;
	        int idx = 0; 
	        int len = sentence.length; 
	        for(int i = 0; i< rows; i++){
	            // for(int j = 0; j< cols; j++){
	                int left= cols;
	                while(sentence[idx].length() <= left){
	                    left -= (sentence[idx].length()+1);
//	                    idx ++;
	                    idx = (idx+1)%len; 
	                    rst++; 
	                }
//	    	        System.out.println("rst "+ rst);
	                if(idx == len) return rst/len; 
	            // }
	        }
	        System.out.println("rst "+ rst);
	        return rst/len; 
	    }
	    
	    
	    public int[] threeEqualParts(int[] A) {
	        // 先数1的数量 如果%3 ！=0 return false 尤其注意最后一个1的位置
	        int len = A.length;
	        int one = 0;
	        // last 1 
	        int pos =40000; 
	        for(int i = len-1; i>=0; i--){
	            if(A[i] == 1){
	                if(pos > 30000) pos = i; 
	                one++; 
	            } 
	        }
	        if(one == 0) return new int[]{0,len-1}; 
	        if(one %3 != 0) return new int[]{-1,-1}; 
	        int note = len-1-pos;

	        int div = one/3; 
	        StringBuilder sb = new StringBuilder();
	        int firstOne = -1; 
	        boolean start =false; 
	        int count = 0; 
	        int i = 0;
	        while(i<len && count < div){
	            if(start){
	                sb.append(A[i]); 
	                if(A[i]==1){
	                    count ++; 
	                }
	            }
	            if(A[i] == 1 && firstOne == -1){
	                firstOne = i; 
	                start = true;
	                sb.append(A[i]);
	                count ++; 
	            } 
	            i++; 
	        }
	        int k = 0; 

	        while(i<len && k<note){
	            if(A[i] ==1) return new int[]{-1,-1};
	            sb.append(A[i]);
	            i++; 
	            k++; 
	        }
//	        System.out.println("i "+ i);

	        if(k<note) return new int[]{-1,-1}; 
	        String s1 = sb.toString(); 
//	        System.out.println("s1 "+ s1);
	        
	        int ah = 0; 
	        int j = i; 
	        while(j<len && A[j] == 0){
	            j++;
	            ah++; 
	        } 
	        String s2 = helper(A, note, i, div); 
//	        System.out.println("s2 "+ s2);
	        
	        if(s2.equals("")) return new int[]{-1,-1};
	        if(!s1.equals(s2)) return new int[]{-1,-1};
//	        System.out.println(i+ah+s1.length());
	        String s3 = helper(A, note, i+ah+s1.length(), div);
//	        System.out.println("s3 "+ s3);
	        if(s3.equals("")) return new int[]{-1,-1};
	        if(!s3.equals(s2)) return new int[]{-1,-1};

	        return new int[]{i-1, i+ah+s1.length()}; 
	    }
	    
	    private String helper(int[] A, int note, int i, int div){
	        int len = A.length;
	        StringBuilder sb = new StringBuilder();
	        int firstOne = -1; 
	        boolean start =false; 
	        int count = 0; 

	        while(i<len && count < div){
	            if(start){
	                sb.append(A[i]); 
	                if(A[i]==1){
	                    count ++; 
	                }
	            }
	            if(A[i] == 1 && firstOne == -1){
	                firstOne = i; 
	                start = true;
	                sb.append(A[i]);
	                count ++; 
	            } 
	            i++; 
	        }
	        int k = 0; 
	        while(i<len && k<note){
	            if(A[i] ==1) return "";
	            sb.append(A[i]);
	            i++; 
	            k++;
	        }
	        if(k<note) return ""; 
	        return sb.toString();
	    }
	    public int numSubarraysWithSum(int[] A, int s) {
	        int len = A.length;
	        int[] sum = new int[len];
	        // sum[i] = sum of A[0] ... A[i]
	        for(int i=0; i<len; i++){
	            if(i==0){
	                sum[i] = A[i]; 
	            }else{
	                sum[i] = sum[i-1] + A[i]; 
	            }
	        }
	        System.out.println(Arrays.toString(sum));
	        int rst = 0; 
	        // [l,r] 0
	        int r = 0; 
	        for(int l = 0; l < len; l++){
	        	r= Math.max(r, l);
	            if(l==0){
	                while(r<len && sum[r] < s) r++;
	                if(r==len) return rst; 
	                int cur = r; 
	                while(cur<len && sum[cur] ==s) cur++;
	                rst +=cur- r;
	                
	            }else{
	                if(r==len){
	                        if(sum[len-1] - sum[l-1] == s) {
		                    	if(A[l] == 0) {
		                    		rst +=len-l; 
		                    	}else {
		                    		rst+=1; 
		                    	}
		                    } 
	                }else{
	                    while(r<len && sum[r]-sum[l-1] < s) r++;
	                    if(r==len) return rst; 
	                    int cur = r; 
	                    while(cur<len && sum[cur]-sum[l-1] ==s) cur++;
	                    rst +=cur -r; 
	                }
	            }
	        } 
	        return rst; 
	    }
	    
	    
	    public int minFallingPathSum(int[][] A) {
	        int n = A.length; 
	        int rst = Integer.MAX_VALUE; 
	        int[][] dp = new int[n][n];
	        
	        for(int i = 0; i< n; i++){
	            dp[n-1][i] = A[n-1][i]; 
	        }
	        for(int i = n-2; i>=0; i--){
	            for(int j =0; j<n; j++){
	                int cur = Integer.MAX_VALUE; 
	                for(int k = Math.max(0,j-1); k<= Math.min(n-1, j+1); k++){
	                    cur = Math.min(cur, dp[i+1][k]); 
	                }
	                dp[i][j] = A[i][j] + cur; 
	            }
	        }
	        printMatrix(dp);
	        for(int i = 0; i< n; i++){
	            rst = Math.min(rst, dp[0][i]); 
	        }
	        return rst; 
	    }
	    
	    
	    public List<String> traceDisease(List<String> info){
	    	List<String> rst = new ArrayList<>();
	    	int ppl = info.size();
//	    	Set<String> set = new HashSet<>(); 
	    	int idx = 0; 
	    	Map<String, Integer> map = new HashMap<>(); 
 	    	String names = ""; 
	    	for(String s: info) {
	    		String[] arr = s.split("\\s");
	    		names += arr[0] + " ";
	    		for(int i =1; i< arr.length; i++) {
	    			if(!map.containsKey(arr[i])) {
	    				map.put(arr[i], idx);
	    				idx++; 
	    			}
	    		}
	    	}
	    	names = names.trim();
	    	int locs = map.size(); 
	    	boolean[][] table1 = new boolean[locs][365]; 
	    	String[][] table2 = new String[ppl][365];
	    	
	    	
	    	int col = 0; 
	    	for(int i= 0; i< ppl; i++) {
	    		String[] arr = (info.get(i).split("\\s")); 
	    		table2[i][0] = arr[1]; 
	    		if(!arr[1].equals("HEALTHY")) {
	    			table1[map.get(arr[2])][0] = true; 
	    		}
	    	}

	    	while(col<364 && !stop(table2, col)) {
	    		col++; 
	    		
		    	for(int i= 0; i< ppl; i++) {
		    		int todayLoc = -1; 
		    		int lastLoc =-1; 
		    		String[] arr = (info.get(i).split("\\s")); 
	    			int len = arr.length; 
	    			if(col>=len-2) {
//	    				todayLoc = (col-1-2)%(len-2)+2;
	    				todayLoc = col%(len-2) +2; 
	    			}else {
	    				todayLoc = 2+col; 
	    			}
	    			
	    			if(col-1 >= len-2) {
//	    				lastLoc = (col-2-2)%(len-2)+2; 
	    				lastLoc = (col-1)%(len-2) +2; 
	    			}else {
	    				lastLoc = 2+col-1; 
	    			}
//	    			System.out.println("i "+i + " lastLoc" + arr[lastLoc] + " todayLoc "+ arr[todayLoc]);
	    			
	    			
		    		String prev = table2[i][col-1]; 
		    		if(prev.equals("SICK"))
	    				table2[i][col] = "RECOVERING";
	    			else if (prev.equals("RECOVERING")) {
	    				table2[i][col] = "HEALTHY";
	    			}else {
	    				if(table1[map.get(arr[lastLoc])][col-1]) {
	    					table2[i][col] = "SICK";
	    				}else
	    					table2[i][col] = "HEALTHY";
	    			}
		    		
		    		if(!table2[i][col].equals("HEALTHY")) {
		    			table1[map.get(arr[todayLoc])][col] = true; 
		    		}
		    	}
	    	}
	    	
	    	if(col == 365) {
	    		rst = outPut(table2, 364);
	    		rst.add(0, names); 
	    		rst.add("365");
	    		System.out.println("result "+365);
	    		return rst; 
	    	}else {
	    		rst= outPut(table2, col); 
	    		rst.add(0, names); 
	    		rst.add(col+1+"");
	    		int o = col+1; 
	    		System.out.println(""+col+1 );
	    		return rst; 
	    	}
	    }
	    
	    // end is inclusive 
 	    private List<String> outPut(String[][] table2, int end){
 	    	int n = table2.length; 
 	    	List<String> rst = new ArrayList<>();
 	    	for(int i = 0; i<=end; i++) {
 	    		String s = ""; 
 	    		for(int j = 0; j<n; j++) {
 	    			s += table2[j][i] + " ";
 	    		}
 	    		s.trim();
 	    		rst.add(s);
 	    	}
 	    	System.out.println(rst.size()-end);
 	    	return rst; 
 	    }
 	    
	    private boolean stop(String[][] table2, int col) {
	    	int m = table2.length;
	    	for(int i=0; i<m; i++) {
	    		if(!table2[i][col].equals("HEALTHY")) return false;
	    	}
	    	return true; 
	    }
	    
	    
	    public void printMatrix2(int[][] dp) {
	    	for (int i = 0; i < dp.length; i++) {
		        for (int j = 0; j < dp[i].length; j++) {
		            System.out.print(dp[i][j] + " ");
		        }
		        System.out.println();
		    }
	    }
	    
	    public void printMatrix3(boolean[][] matrix) {
	    	for (int i = 0; i < matrix.length; i++) {
		        for (int j = 0; j < matrix[i].length; j++) {
		            System.out.print(matrix[i][j] + " ");
		        }
		        System.out.println();
		    }
	    }    
	    
	    
	    
//	    public int knightDialer(int n) {
//	        if(n==1) return 10; 
//	        n++;
//	        int[] zero = new int[n]; 
//	        int[] one = new int[n];
//	        int[] two = new int[n];
//	        int[] four = new int[n]; 
//	        zero[1] = 1; 
//	        one[1] = 1;
//	        two[1] = 1;
//	        four[1] = 1; 
//	        for(int i = 2; i<=n-1; i++){
//	            zero[i] = 2*four[i-1];
//	            one[i] = four[i-1] + two[i-1];
//	            two[i] = 2*one[i-1];
//	            four[i] = 2 * one[i-1] + zero[i-1];
//	        }
//	        return 4*one[n-1] + 2*four[n-1] + 2*two[n-1] + zero[n-1];
//	    }
	    
	    public int knightDialer(int n) {
	        int mod = 1000000007; 
	        if(n==1) return 10; 
	        n++;
	        long[] zero = new long[n]; 
	        long[] one = new long[n];
	        long[] two = new long[n];
	        long[] four = new long[n]; 
	        zero[1] = 1; 
	        one[1] = 1;
	        two[1] = 1;
	        four[1] = 1; 
	        for(int i = 2; i<=n-1; i++){
	            zero[i] = (2*four[i-1])%mod;
	            one[i] = (four[i-1] + two[i-1])%mod;
	            two[i] = (2*one[i-1])%mod;
	            four[i] = ( (2 * one[i-1])%mod + (zero[i-1]%mod))%mod;
	        }
	        long l = ((4*one[n-1])%mod + (2*four[n-1])%mod + (2*two[n-1])%mod + zero[n-1])%mod;
	        return Math.toIntExact(l);
	    }
	    
	    private void printStringSet(Set<String> set) {
	    	for(String s: set) {
	    		System.out.print(s + " ");
	    	}
//	    	return;
	    }
	    
	    int[][] dirs = {{-1,0},{1,0},{0,-1},{0,1}};
	    public int shortestBridge(int[][] A) {
	        Set<String> i1 = new HashSet<String>();
	        Set<String> i2 = new HashSet<String>();
	        int m = A.length; 
	        int n = A[0].length; 
	        boolean stop = false; 
	        for(int i =0; i< m; i++){
	            if(stop) break; 
	            for(int j = 0; j< n; j++){
	                if(A[i][j] == 1){
	                    if(i1.size() ==0){
	                        i1 = isLand(i,j, A, i1); 
	                    }else if(i1.size() > 0 && i1.contains(i+"->"+j)){
	                        continue;
	                    }else{
	                        i2 = isLand(i,j, A, i2);
	                        stop = true; 
	                        break; 
	                    }
	                }
	            }
	        }
	        
	        
//		    System.out.println("s1 "+ i1.size());
//		    printStringSet(i1); 
//		    
//		    System.out.println("s2 "+ i2.size());
//		    printStringSet(i2); 
		    
	        int rst = Integer.MAX_VALUE;
	        
	        
	        // now do a BFS
	        // if(i1.size() < i2.size()){
	        
//	        int why = bfs(0,0,i1,i2,A); 
//	        System.out.println(why);
	            for(String cur: i1){
	                int[] arr = helper(cur);
	                int ii1 = arr[0], ii2 = arr[1];
	                int dis = bfs(ii1, ii2, i1, i2, A);
	                if(dis == 1)  System.out.println("ii1 "+ ii1 + "ii2 "+ ii2);
	                if(dis < rst) rst = dis; 
	            }
	        return rst; 
	    }
	    
	    // find the nearest distance from (x,y) in i1 to anypoint in i2 
	    private int bfs(int x, int y, Set<String> i1, Set<String> i2, int[][] A){
	        int m = A.length; 
	        int n = A[0].length; 
	        int rst = 0; 
	        // boolean mid = true; 
	        Set<String> visited = new HashSet<>(); 
	        Queue<String> q = new LinkedList<>(); 
	        q.add(x+"->"+y);
	        while(!q.isEmpty()){
	            int size = q.size();
	            for(int i =0; i<size; i++){
	                String cur = q.poll(); 

		            
	                int[] arr = helper(cur);
	                int ii1 = arr[0], ii2 = arr[1];
	                
	                
	                visited.add(ii1+"->"+ii2); 
	                
	                for(int[] d: dirs){
	                    int xx = ii1+ d[0], yy = ii2+d[1];
	                    if(xx<0 || yy< 0 || xx>m-1 || yy> n-1 || i1.contains(xx+"->"+yy) || visited.contains(xx+"->"+yy)) continue;
	                    if(i2.contains(xx+"->"+yy)) {
	                    	printStringSet(visited);
	                    	System.out.println("really weird"); 
	                    	return rst; 
	                    }
	                    q.add(xx+"->"+yy); 
	                }
	            }
	            rst++; 
	        }
	        
	        return Integer.MAX_VALUE; 
	    }
	    
//	    private boolean mid(Set<String> isL, int x, int y, int[][] A) {
//	    	int m = A.length, n = A[0].length; 
//	    	for(int[] d: dirs){
//	            int xx = x+ d[0], yy = y+d[1];
//	            if(xx<0 || yy< 0 || xx>m-1 || yy> n-1 || isL.contains(xx+"->"+yy)) continue;
//	            return false; 
//	        }
//	    	return true; 
//	    }
	    
	    private int[] helper(String cur){
	        int i1 = 0, i2= 0; 
	        int i =0; 
	        while(cur.charAt(i) != '-'){
	            char c = cur.charAt(i);
	            i1 = i1 * 10 + (c-'0'); 
	            i++; 
	        }
	        i+=2;
	        while(i<cur.length()){
	            char c = cur.charAt(i);
	            i2 = i2 * 10 + (c-'0'); 
	            i++; 
	        }
	        return new int[]{i1,i2}; 
	    }

	    private Set<String> isLand(int x, int y, int[][] A, Set<String> isL){
	        int m = A.length; 
	        int n = A[0].length; 
	        isL.add(x+"->"+y);
	        for(int[] d: dirs){
	            int xx = x+ d[0], yy = y+d[1];
	            if(xx<0 || yy< 0 || xx>m-1 || yy> n-1 || A[xx][yy] ==0 || isL.contains(xx+"->"+yy)) continue;
	            isL = isLand(xx,yy,A,isL); 
	        }
	        return isL; 
	    }
	    
	    
	    public int shortestDistance(int[][] maze, int[] start, int[] destination) {
	        Set<String> visited = new HashSet<>(); 
	        return dfs(maze, start[0], start[1], destination, visited, 0); 
	    }
	    
	    private int dfs(int[][] maze, int x, int y, int[] dest, Set<String> visited, int steps){
	        int rst = Integer.MAX_VALUE;
	        int m = maze.length, n = maze[0].length; 
	        if(x<0 || x>=m || y<0 || y>=n || visited.contains(x+"->"+y)) return Integer.MAX_VALUE; 
	        if(x ==dest[0] && y==dest[1]) return steps; 
	        visited.add(x+"->"+y);
	        for(int[] d: dirs){
	            int step =steps; 
	            int xx = x+ d[0], yy = y+d[1];  
	           
	            if(xx<0 || xx>=m || yy<0 || yy>=n || maze[xx][yy] == 1|| visited.contains(xx+"->"+yy)) continue;
	            step++; 
	            while(!(xx<0 || xx>=m || yy<0 || yy>=n || maze[xx][yy] == 1)){
	                xx += d[0];
	                yy +=d[1]; 
	                step++; 
	            }
	            xx -= d[0];
	            yy -=d[1]; 
	            step--; 
//	            if(visited.contains(xx+"->"+yy)) continue;
	            int now = dfs(maze, xx, yy, dest, visited, step);
	            if(now == 9) {
	            	System.out.println("xx "+xx + "yy "+yy +"step "+ step );
	            }
	            if( now < rst) {
	                rst = now; 
	            }
	        }
	        return rst; 
	    }
	    
	    private int[] helper1(String cur){
	        int i1 = 0, i2= 0; 
	        int i =0; 
	        while(cur.charAt(i) != '-'){
	            char c = cur.charAt(i);
	            i1 = i1 * 10 + (c-'0'); 
	            i++; 
	        }
	        i+=2;
	        while(i<cur.length()){
	            char c = cur.charAt(i);
	            i2 = i2 * 10 + (c-'0'); 
	            i++; 
	        }
	        return new int[]{i1,i2}; 
	    }
//        public String shortestSuperstring(String[] A) {
//	    	StringBuilder sb = new StringBuilder(); 
//	    	int min = Integer.MAX_VALUE;
//	    	String rst = sb.toString(); 
//	        for(int i=0; i<A.length; i++){
//	            HashSet<Integer> set = new HashSet<>();
//	            for(int k=0 ; k<A.length; k++) set.add(k);
//	            set.remove(i); 
//	            sb = new StringBuilder(); 
//	            sb.append(A[i]);
//	            int cur = i; 
//	            for(int j = 1; j<A.length; j++) {
//	            	int[] idx = findNext(set,A,A[cur]);
//	            	int next; 
//	            	System.out.println("cur "+ A[cur]);
//	            	if(idx[0] != -1) System.out.println(A[idx[1]]);
//	            	if(idx[0] == -1) {
//	            		next = set.iterator().next(); 
//	            		sb.append(A[next]);
//	            	}else {
//	            		next = idx[0];
//	            		sb.append(A[next].substring(idx[1]));
//	            	}
//            		cur = next; 
//            		set.remove(next);
//	            }
////	            System.out.println(j);
//	            System.out.println(sb.toString());
//	            if(sb.length() < min) {
//	            	min = sb.length();
//	            	rst = sb.toString(); 
//	            }
//	        }
//	        return rst; 
//	    }
//	    
//	    // return [max index, max length]
//	    private int[] findNext(HashSet<Integer> set, String[] A, String cur){
//	        int max = 0, maxIdx = -1;
//	        for(int i:set){
//	            int num = overlap(cur, A[i]); 
//	            if(num > max){
//	                max = num;
//	                maxIdx = i; 
//	            }
//	        }
//	        if(max ==0){
//	            return new int[]{-1,-1}; 
//	        }else{
//	            return new int[]{maxIdx, max}; 
//	        }
//	    }
//    
//    // find the maximum length of the tail of s1 that overlaps with the head of s2. 
//    private int overlap(String s1, String s2){
//        int l1 = s1.length();
//        int l2 = s2.length();
//        // i is the length of the overlapping substring
//        for(int i = Math.min(l1,l2)-1; i>=0; i--){
//            if(s1.substring(l1-i).equals(s2.substring(0,i))) return i; 
//        }
//        return 0; 
//    }
    
	    public String shortestSuperstring(String[] A) {
	    	HashSet<Integer> set = new HashSet<>();
            for(int k=0 ; k<A.length; k++) set.add(k);
            return helperShortest(A, set);
	    }
	    
	    public String helperShortest(String[] A, HashSet<Integer> copy) {
	    	StringBuilder sb = new StringBuilder();
	    	
	    	int min = Integer.MAX_VALUE;
	    	String rst = sb.toString(); 
	    	Iterator<Integer> iter = copy.iterator(); 
	    	int[] start = new int[A.length];
	    	int ix = 0; 
	    	while(iter.hasNext()) {
	    		start[ix] = iter.next();
	    		ix++; 
	    	}
	    	
	    	int len = copy.size(); 
	        for(int s=0; s<len; s++){
	        	String check = ""; 
	        	HashSet<Integer> set = (HashSet<Integer>)copy.clone();
	        	int cur = start[s]; 
	            set.remove(cur); 
	            sb = new StringBuilder(); 
	            sb.append(A[cur]);
	            int size = set.size();
		        for(int j = 0; j<size; j++) {

	            	int[] idx = findNext(set,A,A[cur]);
	            	int next; 
	            	
	            	if(idx[0] == -1) {
	            		check = sb.toString() + helperShortest(A, set);
	            		System.out.println("break " + sb.toString() + ";" + helperShortest(A, set));
	            		break; 
//	            		next = set.iterator().next(); 
//	            		sb.append(A[next]);
	            	}else {
	            		next = idx[0];
	            		sb.append(A[next].substring(idx[1]));
	            	}
	            	System.out.println("cur "+ A[cur]);
	            	if(idx[0] != -1) System.out.println(A[next]);
            		cur = next; 
            		set.remove(next);
	            }
		        if(!check.equals("")) check = sb.toString(); 
	            System.out.println(sb.toString());
	            if(check.length() < min) {
	            	min = check.length();
	            	rst = check;
	            }
	        }
	        return rst; 
	    }
	    
	    // return [max index, max length]
	    private int[] findNext(HashSet<Integer> set, String[] A, String cur){
	        int max = 0, maxIdx = -1;
	        for(int i:set){
	            int num = overlap(cur, A[i]); 
	            if(num > max){
	                max = num;
	                maxIdx = i; 
	            }
	        }
	        if(max ==0){
	            return new int[]{-1,-1}; 
	        }else{
	            return new int[]{maxIdx, max}; 
	        }
	    }
	    // find the maximum length of the tail of s1 that overlaps with the head of s2. 
	    private int overlap(String s1, String s2){
	        int l1 = s1.length();
	        int l2 = s2.length();
	        // i is the length of the overlapping substring
	        for(int i = Math.min(l1,l2)-1; i>=0; i--){
	            if(s1.substring(l1-i).equals(s2.substring(0,i))) return i; 
	        }
	        return 0; 
	    }
	    
	    public String largestTimeFromDigits(int[] A) {
	        int max = 24*60-1; 
	        String rst =""; 
	        List<List<Integer>> lst = permuteUnique(A);
	        int validMax = -1;
	        for(List<Integer> l:lst){
	            if(l.get(2)>=6 || l.get(0)>3 || (l.get(0) == 2 && l.get(1)>4) ) continue; 
	            int time = (l.get(0)*10 + l.get(1))*60 + l.get(2)*10+l.get(3);
	            if(time<=max && time > validMax){
	                validMax = time;
	                rst = ""+l.get(0)+l.get(1)+":"+l.get(2)+l.get(3);
	            }
	        }
	        if(validMax ==-1) return "";
	        else return rst; 
	    }
	    
	    
	    public List<List<Integer>> permuteUnique(int[] nums) {
	        List<List<Integer>> res = new ArrayList<List<Integer>>();
	        if(nums==null || nums.length==0) return res;
	        boolean[] used = new boolean[nums.length];
	        List<Integer> list = new ArrayList<Integer>();
	        Arrays.sort(nums);
	        dfs(nums, used, list, res);
	        return res;
	    }

	    public void dfs(int[] nums, boolean[] used, List<Integer> list, List<List<Integer>> res){
	        if(list.size()==nums.length){
	            res.add(new ArrayList<Integer>(list));
	            return;
	        }
	        for(int i=0;i<nums.length;i++){
	            if(used[i]) continue;
	            if(i>0 &&nums[i-1]==nums[i] && !used[i-1]) continue;
	            used[i]=true;
	            list.add(nums[i]);
	            dfs(nums,used,list,res);
	            used[i]=false;
	            list.remove(list.size()-1);
	        }
	    }
	    
	    public int largestComponentSize(int[] A) {
	        int n = A.length; 
	        int max = 0; 
	        Arrays.sort(A);
	        UnionFind uf = new UnionFind(n);
	        for (int i = 0; i < n - 1; i++) {
	            for (int j = i + 1; j < n; j++) {
	                if (!primeTo(A[i], A[j])) uf.union(i, j);
	            }
	        }
	        
	        uf.clean();
	        
	        HashMap<Integer, Integer> map = new HashMap<>(); 
	        int[] parent = uf.parent();
	        int[] realParent = new int[parent.length];
	        for(int i=0; i<parent.length; i++) {
	        	realParent[i] = A[parent[i]];
	        }
	        System.out.println(Arrays.toString(A));
	        System.out.println(Arrays.toString(realParent));
	        for(int i = 0; i< parent.length; i++){
	            int val = map.getOrDefault(parent[i],0);
	            map.put(parent[i], val+1); 
	        }
	        for(int i: map.keySet()){
	            if(map.get(i) > max)
	                max = map.get(i); 
	        }
	        return max;
	    }
	    
	    // If relatively prime return true 
	    private boolean primeTo(int a, int b) {
	        int t;
	        while(b != 0){
	            t = a;
	            a = b;
	            b = t%b;
	        }
	        return a==1;
	    }
	    
	    class UnionFind {
	        private int[] parent, rank;
	        // index to its group size 
	   
	        public UnionFind(int n) {
	            parent = new int[n];
	            rank = new int[n];
	            for (int i = 0; i < n; i++) {
	                parent[i] = i;
	            }
	        }
	        
	        public int find(int p) {
	        	while (p != parent[p]) {
	                parent[p] = parent[parent[p]]; // path compression by halving
	                p = parent[p];
	            }
	            return p;
	        }
	        
	        public void union(int p, int q) {
	            int rootP = find(p);
	            int rootQ = find(q);
	            if (rootP == rootQ) return;
	            if (rank[rootQ] > rank[rootP]) {
	                parent[rootP] = rootQ;
	            }
	            else {
	                parent[rootQ] = rootP;
	                if (rank[rootP] == rank[rootQ]) {
	                    rank[rootP]++;
	                }
	            }
	        }
	        
	        public void clean() {
	        	for(int i= 0; i<parent.length; i++) {
	        		find(i);
	        	}
	        }
	        
	        public int[] parent() {
	            return parent;
	        }
	    }
	    
	    
	    public boolean isCousins(TreeNode root, int x, int y) {
	        if(root==null) return false; 
	        if(root.left != null && root.right != null && (root.left.val == x && root.right.val== y  || root.left.val == y && root.right.val== x)) return false; 
	        
	        if (isCousins(root.left, x, y) || isCousins(root.right,x,y)) 
	        {
	        	System.out.println("here1"+ isCousins(root.left, x, y)  +"or "+ isCousins(root.right, x, y) );
	        	return true; 
	        }
	        if((depth(root.left, x) == depth(root.right, y) && depth(root.left, x) != -1) || (depth(root.left, y) == depth(root.right, x)&& depth(root.left, y) != -1) ) 
	        {
	        	System.out.println("here2");
	        	return true; 
	        }
	        return false; 
	    }
	    
	    private int depth(TreeNode root, int x){
	        if(root == null) return -1; 
            if(root.val ==x) return 1; 
            int l = depth(root.left, x);
            int r = depth(root.right, x); 
            if(l == r && l ==-1) return -1; 
	        return 1+Math.max(l,r);
	    }
	    
	    public int orangesRotting(int[][] matrix) {
	        int rst = 0; 
	        int m = matrix.length;
	        int n = matrix[0].length;
	        for(int i = 0; i<m; i++){
	            for(int j =0; j<n; j++){
	                if(matrix[i][j] == 1){
	                    int cur = BFS(matrix, new int[m][n],i, j);
	                    if(cur > rst){
//	                    	System.out.println("Why "+ i + " & " + j);
	                        rst = cur; 
	                    }
	                }
	            }
	        }
	        if(rst == Integer.MAX_VALUE) return -1;
	        return rst; 
	    }
	    
	    // given a 1 at position (x,y), find the Manhattan distance to its closest 2. 
	    private int BFS(int[][] matrix, int[][] visited, int x, int y){
	        int m = matrix.length;
	        int n = matrix[0].length; 
//	        if(matrix[x][y] == 2) return 0; 
	        int[][] dirs = {{-1,0},{1,0},{0,-1},{0,1}};
	        int dis = 0; 
 	         Queue<int[]> que1 = new LinkedList<>();
	         que1.add(new int[]{x,y});
	         visited[x][y]=1;
	         while(!que1.isEmpty()){
	        	 dis++; 
	        	 int size = que1.size();
	        	 for(int i = 0; i<size; i++) {
	        		 int[] cell = que1.poll();
//		             System.out.println(dis);
		             for(int[] d: dirs){
		                 int r= cell[0] + d[0];
		                 int c = cell[1]+d[1];
		                 if(r<0 || r>=m || c<0 || c>=n || matrix[r][c] == 0 || visited[r][c]==1) continue;
//		                 System.out.println("r & c" + r + "," +c);
		                 if(matrix[r][c] == 2) return dis; 
		                 que1.add(new int[] {r,c});
		                 visited[r][c] = 1;
		             }
	        	 }
	             
	         }
	         return Integer.MAX_VALUE;
	    }  
	    public int numSquarefulPerms(int[] A) {
	        List<List<Integer>> perm = permuteUnique (A); 
	        Map<Integer, Boolean> hasmap = new HashMap<>(); 
	        int result = 0; 
	        for(List<Integer> l : perm){
	            if(check(l,hasmap)) result++;
	        }
	        return result; 
	    }
	    
	    private boolean check(List<Integer> l, Map<Integer, Boolean> map){
	        for(int i =0; i<l.size()-1; i++){
	            int whatIf = l.get(i) + l.get(i+1);
	            boolean cur = false; 
	            if(map.containsKey(whatIf))  cur = map.get(whatIf);
	            map.put(whatIf, sqrt(whatIf));
	            if(cur) continue;
	            else return false; 
	        }
	        return true; 
	    }
	    
	    
	    private boolean sqrt(int x) {
	        if (x == 0)
	            return true;
	        int left = 1, right = Integer.MAX_VALUE;
	        while (true) {
	            int mid = left + (right - left)/2;
	            if(mid * mid == x) return true; 
	            if (mid > x/mid) {
	                right = mid - 1;
	            } else {
	                if (mid + 1 > x/(mid + 1))
	                    return false;
	                left = mid + 1;
	            }
	        }
	    }

	    private void backtrack(List<List<Integer>> list, List<Integer> tempList, int [] nums, boolean [] used){
	        if(tempList.size() == nums.length){
	            list.add(new ArrayList<>(tempList));
	        } else{
	            for(int i = 0; i < nums.length; i++){
	                if(used[i] || i > 0 && nums[i] == nums[i-1] && !used[i - 1]) continue;
	                used[i] = true; 
	                tempList.add(nums[i]);
	                backtrack(list, tempList, nums, used);
	                used[i] = false; 
	                tempList.remove(tempList.size() - 1);
	            }
	        }
	    }
	    
	  
	    public int findJudge(int N, int[][] trust) {
	    	if (trust == null || trust.length ==0 ){
                if(N==1) return 1; 
                return -1; 
            } 
	    	if (trust == null || trust.length ==0 ) return -1; 
	    	int len = trust.length; 
	    	HashMap<Integer,Integer> table = new HashMap(); 
	    	ArrayList<Integer> arr = new ArrayList(); 
	    	for (int i =0; i< len; i++) {
	    		if (table.containsKey(trust[i][0])) {
	    			table.put(trust[i][0], -200);
	    		}else {
	    			table.put(trust[i][0], -200);
	    		}
	    		
                int cnt = -1; 
                if (table.containsKey(trust[i][1])){
                    cnt = table.get(trust[i][1])+1;
                    table.put(trust[i][1], cnt);
                }else{
                    table.put(trust[i][1],1);
                    cnt = 1; 
                }
                arr.add(trust[i][1]);
	    	}
	    	for(int i =0; i< arr.size(); i++) {
	    		if(table.containsKey(arr.get(i)) && table.get(arr.get(i)) == N-1) return arr.get(i); 
	    	}
	    	return -1; 
	    }
	    class Grid{
	    	private int m, n;
	    	public Grid(int m, int n) {
	    		this.m = m;
	    		this.n = n; 
	    	}
	    	
	    	@Override 
	    	public boolean equals(Object o) {
	    		if(o == this) {
	    			return true; 
	    		}
	    		if(!(o instanceof Grid)) {
	    			return false; 
	    		}
	    		Grid c = (Grid) o; 
		    	
		    	return (c.m == this.m && c.n == this.n); 
	    	}    	
	    }
	    
	    public int[] gridIllumination(int N, int[][] lamps, int[][] queries) {
	    	int len = queries.length; 
	    	int[] rst = new int[len]; 
	    	HashMap<Grid, ArrayList<Grid>> map = new HashMap<>();
	    	for(int i =0; i<lamps.length; i++) {
	    		Grid cur = new Grid(lamps[i][0],lamps[i][1]);
	    		for(int m =0 ; m< N; m++) {
	    			ArrayList<Grid> curList;
	    			if(map.containsKey(new Grid(lamps[i][0],m))) {
	    				curList = map.get(new Grid(lamps[i][0],m));
	    	    		curList.add(cur);
	    			}else {
	    				curList = new ArrayList<Grid>();
	    			}
	    			map.put(new Grid(lamps[i][0],m),curList); 
	    		}
	    		for(int m =0 ; m< N; m++) {
	    			ArrayList<Grid> curList;
	    			if(map.containsKey(new Grid(m,lamps[i][1]))) {
	    				curList = map.get(new Grid(m,lamps[i][1]));
	    	    		curList.add(cur);
	    			}else {
	    				curList = new ArrayList<Grid>();
	    			}
	    			map.put(new Grid(m, lamps[i][1]),curList); 
	    		}
	    	    int[][] dirs = {{-1,-1},{1,1},{1,-1},{-1,1}};
	    	    
	    	    for(int j =0 ; j< 4; j++) {
	    	    	int mm = lamps[i][0], nn = lamps[i][1];
	    	    	mm += dirs[j][0]; 
	    	    	nn += dirs[j][1];
	    	    	while(mm>=0 && mm <N && nn >=0 && nn<N) {
	    	    		ArrayList<Grid> curList = map.get(new Grid(mm, nn));
	    	    		curList.add(cur);
		    			map.put(new Grid(mm, nn),curList); 
		    		} 
		    	}
	    	}
	    	for(int i =0 ; i<queries.length; i++) {
	    		if(map.containsKey(new Grid(queries[i][0],queries[i][0]))){
	    			ArrayList<Grid> curList= map.get(new Grid(queries[i][0],queries[i][0]));
	    			for(int k = 0; k<curList.size(); k++) {
	    				Grid grid = curList.get(k); 
	    				int m = grid.m; 
	    				int n = grid.n; 
	    				if(lamps[m][n] ==1) rst[i] = 1; 
	    				lamps[m][n] = 0 ;
	    			} 
	    		}
	    	}
	        return rst; 
	    }
	    
	    public List<String> commonChars(String[] A) {
	        List<String> rst = new ArrayList<>();
	        HashMap<Character, Integer> map = new HashMap<>(); 
	        for(int i = 0; i<A[0].length(); i++) {
	        	char c= A[0].charAt(i);
	        	if(map.containsKey(c)) {
	        		map.put(c, map.get(c)+1);
                }else {
                	map.put(c,1);
                }
	        }
	        for(int i = 1; i<A.length; i++) {
	            int[] arr = new int[26]; 
	            for(int j=0; j<A[i].length(); j++)
	                arr[A[i].charAt(j)-'a'] ++; 
	            for(int k = 0; k<26; k++) {
	                char c = (char)('a'+k); 
	                if(map.containsKey(c)) {
	                    if(map.get(c)>arr[k])
	                        map.put(c, arr[k]); 
	                }
	            }
	        }
	        for(int k = 0; k<26; k++) {
	            char c = (char)('a'+k); 
	            if(map.containsKey(c)) {
	                for(int i = 0; i<map.get(c); i++)
	                    rst.add(""+c);
	            }
	        }
	        return rst;
	    }
	    
	    public boolean isValid(String S) {
	    	int len = S.length(); 
	    	if(len == 0) return true; 
	        if(len<3 || S.charAt(0) != 'a' || S.charAt(len-1)!= 'c') return false; 
	        boolean rst = (S.charAt(0) == 'a' && S.charAt(1) == 'b' && S.charAt(2) == 'c' && isValid(S.substring(3)) ||
	        		(S.charAt(0) == 'a' && S.charAt(1) == 'b' && S.charAt(len-1) == 'c' && isValid(S.substring(2,len-1))) || 
	        		(S.charAt(0) == 'a' && S.charAt(len-2) == 'b' && S.charAt(len-1) == 'c' && isValid(S.substring(1,len-2))) || 
	        		S.charAt(len-3) == 'a' && S.charAt(len-2) == 'b' && S.charAt(len-1) == 'c' && isValid(S.substring(0,len-3)));
	        if (rst) 
	        	return true; 
	        else {
	        	for(int i = 1; i<len-1; i++) {
	        		if(S.charAt(i) =='b' && isValid(S.substring(1,i)) &&  isValid(S.substring(i+1, len-1))) return true; 
	        	}
	        	return false; 
	        }
	    }
	    
	    public int longestOnes(int[] nums, int K) {
	        int max = 0, zero = 0, k = K; // flip at most k zero
	        for (int l = 0, h = 0; h < nums.length; h++) {
	            if (nums[h] == 0)                                           
	                zero++;
	            while (zero > k)
	                if (nums[l++] == 0)
	                    zero--;                                     
	            max = Math.max(max, h - l + 1);
	        }                                                               
	        return max; 
	    }
	    
	    public int mergeStones(int[] stones, int k) {
	    	int len = stones.length; 
	    	if(len == 1) return 0; 
	    	if(len<k) return -1; 
	    	int sum = Integer.MAX_VALUE; 
	    	ArrayList<Integer> idx = new ArrayList<>();  
	    	for(int i =0; i<= len- k; i++) {
	    		int curSum = 0; 
	    		for(int j = 0; j<k; j++) {
	    			curSum += stones[i+j]; 
	    		}
	    		if(curSum < sum) {
	    			sum = curSum; 
	    			idx = new ArrayList<>();
	    			idx.add(i); 
	    		}else if(curSum == sum) {
	    			idx.add(i);
	    		}
	    	}
	    	if(len-k+1 ==1) return sum; 
	    	int rst = Integer.MAX_VALUE; 
	    	int[] arr = new int[len-k+1]; 
	    	for(int p: idx) {
	    		for(int i = 0; i<p; i++) {
		    		arr[i] = stones[i]; 
		    	}
		    	arr[p] = sum; 
		    	for(int j = p+1; j<arr.length; j++) {
		    		arr[j] = stones[j+k-1]; 
		    	}
		    	int next = mergeStones(arr, k); 
		    	if(next == -1) return -1;
		    	rst = Math.min(rst, next + sum);
	    	}
	    	return rst; 
	    }
	    
	    
	    public String minWindow2(String s, String t) {
	        if(s == null || t== null || t.length() == 0 || (s.length() <t.length())) return ""; 
	        int[] dic = new int[256];
	        for(int i=0; i< t.length(); i++){
	            dic[t.charAt(i)] +=1; 
	        }
	        int l = 0;
	        int r= 0; 
	        int count = t.length(); 
	        String rst = ""; 
	        // [l,r)
	        while(r<s.length()){
	            // if(dic[s.charAt(r++)]-->0) count --; 
	            if(dic[s.charAt(r)]>0) count --; 
	            dic[s.charAt(r)]--;
	            r++; 
	            while(count ==0) {
	                if(rst.equals("") || r-l<rst.length())  rst= s.substring(l,r); 
	                // if(dic[s.charAt(l++)]++ ==0) count ++; 
	                if(dic[s.charAt(l)]==0) count ++; 
	                dic[s.charAt(l)]++;
	                l++; 
	            }
	        }
	        return rst;
	    }
	    
	    public int largestSumAfterKNegations(int[] A, int K) {
	    	// greedy 
	        int len = A.length;
	        PriorityQueue<Integer> pq = new PriorityQueue<Integer>(); 
	        int rst = 0; 
	        int sum = 0; 
	        for(int i =0; i<len; i++) {
	        	pq.add(A[i]); 
	        	sum += A[i]; 
	        }
	        for(int i =0; i<K; i++) {
	        	int cur = pq.poll(); 
	        	sum = sum -cur - cur; 
	        	pq.add(-cur); 
	        }
	        return sum; 
	       
	    }
	    public int helper(int N){
	        if(N <= 0) return 0; 
	    	if(N ==1) return 1;
	    	if(N == 2) return 2; 
	    	if(N == 3) return 6;
	    	if(N ==4) return 7; 
	        return (N * (N-1))/(N-2) + (N-3); 
	    }
	    HashMap<Integer, Integer> map = new HashMap<>(); 
	    public int clumsy(int N) {
	            if(N<=0) return 0; 
	        	map.put(1, 1); 
		    	map.put(2,2);
		    	map.put(3,6);
		    	map.put(4, 7);
		    	if(map.containsKey(N)) return map.get(N); 
		    	int rst = 0;
		    	if(N <=7) 
		    		rst = N * (N-1) / (N-2) + (N-3) - helper(N-4); 
		    	else
		    		rst = (N * (N-1) / (N-2) + (N-3) - 2*helper(N-4) + 2*(N-7)+ clumsy(N-4)); 
		    	map.put(N , rst); 
		    	return rst;
	    }
//	    public int minDominoRotations(int[] A, int[] B) {
//	        int[] dic = new int[6];
//	        int len = A.length; 
//	        for(int i = 0; i< len; i++){
//	            dic[A[i]-1] ++; 
//	            dic[B[i]-1] ++;
//	        }
//	        int max = 0; 
//	        int num = -1; 
//	        for(int i =0; i<6; i++){
//	            if(dic[i] > max){
//	                max = dic[i];
//	                num = i+1; 
//	            } 
//	        }
//	        if(max < len) return -1; 
//	        int rst = len; 
//	        int a = 0; 
//	        for(int i= 0; i< len; i++){
//	            if(A[i] == num){
//	                a ++; 
//	                rst --; 
//	                if(B[i] == num)
//	                    rst --; 
//	            }
//	        }
//	        if(rst < len -a) return -1; 
//	        return Math.min(a, len-a); 
//	    }
	    
	    public int minDominoRotations(int[] A, int[] B) {
	        int[] dic = new int[6];
	        int len = A.length; 
	        for(int i = 0; i< len; i++){
	            dic[A[i]-1] ++; 
	            dic[B[i]-1] ++;
	        }
	        int max = 0; 
	        int num = -1; 
	        for(int i =0; i<6; i++){
	            if(dic[i] > max){
	                max = dic[i];
	                num = i+1; 
	            } 
	        }
	        if(max < len) return -1; 
	        int rst = max; 
	        int a = 0; 
	        for(int i= 0; i< len; i++){
	            if(A[i] == num){
	                a ++; 
	                rst --; 
	                if(B[i] == num)
	                    rst --; 
	            }
	        }
	        if(rst < len -a) return -1; 
	        int d1 = Math.min(len-a, len-rst);
	        int d2 = 0;
	        for(int i =0; i< len; i++){
	            if(B[i] != num) d2++; 
	        }
	        return Math.min(d1,d2);
//	        System.out.println("len "+ len);
//	        System.out.println("max "+ max);
//	        System.out.println("a "+ a);
//	        System.out.println("rst "+ rst);
	    }
	    
	    
	    public TreeNode bstFromPreorder(int[] preorder) {
	        return helper2(preorder, 0, preorder.length-1);
	    }
	    
	    
	    private TreeNode helper2(int[] preorder, int s, int t) {
	    	if(s>t) return null; 
	    	if(s == t) return new TreeNode(preorder[s]); 
	    	int root = preorder[s]; 
	    	int i = s+1;
	    	while( i<=t) {
	    		if(preorder[i] > root) break; 
	    		i++; 
	    	}
	    	TreeNode rt = new TreeNode(root);
	    	rt.left = helper2(preorder, s+1, i-1);
	    	rt.right = helper2(preorder, i, t);
	    	return rt;
	    }
	    
	    
	    public int bitwiseComplement(int n) {
//	    	int cur = 0; 
//	    	int pow = 0; 
//	        while(n>0) {
//	        	if(pow == 0) {
//	        		cur = 1-cur;
//	        	}else {
//	        		cur = (int) (1-cur + (n%2)* Math.pow(10, pow));
//	        	}
//	        	pow++; 
//	        	n = n/2; 
//	        }
//	        return cur; 
	    	if(n==0) return 1; 
	        int nbits = (int)(Math.log(n) / Math.log(2) + 1); 
	        return (int) (Math.pow(2, nbits)-1 - n); 
	    }
//	    
//	    public int numPairsDivisibleBy60(int[] time) {
//	    	int count = 0;
//	    	for(int i =0; i< time.length; i++) {
//	    		for(int j = i+1; j< time.length; j++) {
//	    			if((time[i] + time[j])%60 ==0) count ++; 
//	    		}
//	    	}
//	        return count; 
//	    }
	    
	    public int numPairsDivisibleBy60(int[] time) {
	    	int rst = 0; 
	    	int[] count = new int[60];
	    	for(int i =0; i< time.length; i++) {
	    		count[time[i]%60] ++; 
	    	}
	    	for(int i =0; i< time.length; i++) {
	    		if(time[i] == 30) rst += count[30] -1; 
	    		rst += count[60-time[i]%60];
	    		count[time[i]%60] --; 
	    	}
	        return rst; 
	    }
	    
	    
//	    public int shipWithinDays(int[] weights, int D) {
//	    	// dp[i][j] represents cost of moving weights[i,end] in j-1 days 
//	    	int len = weights.length; 
//	    	if(D == 1) {
//	    		int sum = 0; 
//	    		for(int i =0; i< len; i++)
//	    			sum += weights[i]; 
//	    		return sum;
//	    	} 
//	        int[][] dp = new int[len][D];
//	        for(int j = len-1; j>=0 ; j--) {
//	        	if(j == len-1) dp[j][1] = weights[j]; 
//	        	else dp[j][1] = dp[j+1][1] + weights[j]; 
//	        }
//	        for(int i = 2; i< D; i++) {
//	        	for(int j = len-1; j>=0 ; j--) {
//	        		int fill = Integer.MAX_VALUE;
//	        		int sum = 0; 
//	        		for(int k = j; k<len-1; k++) {
//	        			sum += weights[k];
//	        			if(Math.max(sum, dp[k+1][i-1])<fill) fill = Math.max(sum, dp[k+1][i-1]); 
//	        		}
//	        		dp[j][i] = fill; 
//	        	}
//	        }
//	        int rst = Integer.MAX_VALUE;
//	        int sum = 0; 
//	        for(int i=0; i<len-1; i++) {
//	        	sum += weights[i]; 
//	        	if(Math.max(sum, dp[i+1][D-1])<rst) rst = Math.max(sum, dp[i+1][D-1]);
//	        }
//	        return rst; 
//	    }
	    
	    
	    public int shipWithinDays(int[] weights, int D) {
	    	int len = weights.length; 
	    	int[] sum = new int[len];
	    	sum[0] = weights[0]; 
	    	for(int i =1; i<len; i++) {
	    		sum[i] = sum[i-1] + weights[i];  
	    	}
	    	System.out.println(Arrays.toString(sum));
	    	if(D == 1) return sum[len-1];
	    	int avg; 
	    	if(sum[len-1]%D ==0) avg = sum[len-1]/D; 
	    	else avg = sum[len-1]/D+1;
//	    	System.out.println("avg "+ avg);
	    	int idx = Arrays.binarySearch(sum, avg);
    		if(idx < 0 ) {
    			if(-idx == sum.length +1 ) return sum[len-1];
    			else idx = -idx -1; 
    		}
	    	while(checkDays(sum, avg) > D) {
//		    	System.out.println("idx "+ idx);
	    		avg = sum[idx]; 
	    		idx++;
	    	}
	    	return avg; 
	    }
	    
	    private int checkDays(int[] sum, int goal) {
	    	int idx = -1;
	    	int goalDup = goal; 
	    	int prev = -1; 
	    	while(idx < sum.length-1) {
	    		int next = Arrays.binarySearch(sum, goal);
	    		System.out.println("next "+next);
	    		if(next >= 0 ) {
	    			idx = next; 
	    		}else{
	    			if(-next == sum.length +1 ) return goal/goalDup; 
	    			idx = -next -2; 
	    		}
	    		if(prev == idx) return Integer.MAX_VALUE;
	    		prev = idx;
	    		System.out.println("idx "+idx);
	    		goal += goalDup; 
	    	}
	    	return goal/goalDup; 
	    }
//	    private String fN; 
//	    private String lN; 
//	    public practice(String firstName, String lastName) {
//	    	this.fN = firstName; 
//	    	this.lN = lastName; 
//	    }
	    
//		public static void main(String[] args) {
//			int a = 5; 
//			practice hi = new practice();
//			int[] weights = {1,2,3,4,5,6,7,8,9,10}; 
//			String[] arr = {"lov"};
//			int[] sum = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55};
//			String wu = "1"; 
////			System.out.println(hi.checkDays(sum, 6));
////			System.out.println(hi.shipWithinDays(weights, 10));
////			System.out.println(arr);
////			String[] as = new String[10];
////			Object[] ao = new Object[10];
////			
////			ao = as; 
////			ao[0] = 2110;
////			String s = as[0]; 
////			ArrayList<String> ls= new ArrayList<String>();
////			ArrayList<Object> lo= new ArrayList<Object>();
////			lo= ls; //Suppose this is legal
////			lo.add(2110); //Type-checks: Integer subtype Object
////			String s = ls.get(0); //Type-checks: ls is a List<String>
//		}
	    
//	    public static void main(String[] pars) { 
//	        try { 
//	           System.out.println("try-block 0"); 
//	           throw new ArithmeticException("fake exception 1"); 
//	        } catch (ArithmeticException ae) { 
//	           System.out.println("catch-block 0"); 
//	           System.out.println(ae); 
//	           throw new ArithmeticException("fake exception 2"); 
//	        } 
//	    }
	    
	    public boolean canThreePartsEqualSum(int[] A) {
	        int sum = 0; 
	        for(int i =0; i<A.length; i++) {
	        	sum += A[i]; 
	        }
	        System.out.println("sum "+ sum);
	        if (sum%3!=0) return false; 
	        int eachPart = sum/3;
	        int i = 0; 
	        // sum of [0,i)
	        int curSum = 0; 
	        while(i<A.length && curSum!=eachPart) {
	        	curSum += A[i]; 
	        	i++;
	        }
	        if(curSum!=eachPart) return false; 
	        curSum = 0; 
	        while(i<A.length && curSum!=eachPart) {
	        	curSum += A[i]; 
	        	i++;
	        }
	        if(curSum!=eachPart) return false; 
	        curSum = 0; 
	        while(i<A.length && curSum!=eachPart) {
	        	curSum += A[i]; 
	        	i++;
	        }
	        if(curSum!=eachPart) return false; 

	        return true;     
	    }
	    
//	    public int smallestRepunitDivByK(int K) {
//	        if(K%2==0) return -1; 
//	        int cur = 1; 
//	        int len = 1; 
//	        while(cur<1111111111 && cur%K !=0) {
////		        System.out.println(""+ cur);
//	        	cur = cur*10 + 1; 
//	        	len++; 
//	        }
////	        System.out.println(""+ cur);
//	        if(cur%K!=0) {
//	        	if(1111111111%K!=0)
//	        		return -1; 
//	        	else return 10;
//	        }else return len; 
//	    }
//	    
//	    public int maxScoreSightseeingPair(int[] A) {
//	    	// Check[i] stores the min value of A[j] +j where j ranges from [i...end]
//	    	int[] check = new int[A.length];
//	    	check[A.length-1] = A[A.length-1] + A.length -1; 
//	    	for(int i = A.length-2; i>=0; i--) {
//	    		if(A[i]+i < check[i+1]) {
//	    			check[i] = A[i]+i;
//	    		}else check[i] = check[i+1];
//	    	}
////	    	
////	        int left = 0; 
////	        int right =-1; 
//	        int max = Integer.MIN_VALUE; 
//	        for(int i = 0; i<A.length; i++) {
//	        	if(A[i] + i - check[i+1]>max) {
//	        		max = A[i] + i - check[i+1]; 
//	        	}
//	        }
//	        // right is the index where index 0 should be mapped to 
//	        int newLeft= left+1; 
//	        while(A[newLeft] - A[left] - (newLeft - left)<=0) newLeft++; 
//	        while(left < A.length -1){
//	        	if(right == A.length -1) {
//	        		max = Math.max(max, A[A.length -1] + A[left] + left - (A.length-1));
//	        		continue; 
//	        	}
//	        	while(right<A.length) {
//	        		int cur = A[left] + A[right] + left -right; 
//	        		if(cur>= max) {
//	        			max = cur; 
//	        		}
//	        	}
//	            
////	        }
//	        return max; 
//	    }
	    
	    public int maxScoreSightseeingPair(int[] A) {
	    	// Check[i] stores the min value of A[j] +j where j ranges from [i...end]
	    	int[] check = new int[A.length];
	    	check[A.length-1] = -A[A.length-1] + A.length -1; 
	    	for(int i = A.length-2; i>=0; i--) {
	    		if(-A[i]+i < check[i+1]) {
	    			check[i] = -A[i]+i;
	    		}else check[i] = check[i+1];
	    	}
	        int max = Integer.MIN_VALUE; 
	        for(int i = 0; i<A.length-1; i++) {
	        	if(A[i] + i - check[i+1]>max) {
	        		max = A[i] + i - check[i+1]; 
	        	}
	        }
        return max; 
    }
	    
    public int smallestRepunitDivByK(int K) {
    	int rem = 1;
    	int count = 1;
    	int bound = K-1;
    	int repeat = 0;
    	
    	if(K==1) return 1;
		while(count>0){
			count++;
			int remainderN = rem*10 + 1;
			remainderN = remainderN % K;
			if(remainderN == 0){
				break;
			}else if(remainderN == rem || remainderN == 1 || repeat == 2){
				count = -1;
				break;
			}else if (remainderN == bound){
				repeat++;
			}
    		rem = remainderN;
    	}
    	return count;
    }
	    
	    public boolean queryString(String S, int N) {
	    	for (int i=1; i<=N; i++){
	    		String stuff = Integer.toBinaryString(i);
	    		if(!S.contains(stuff)) {
	    			return false; 
	    		}
	    	}
	    	return true;
	    }
	    
	    // backTracking是肯定的
	    public List<List<Integer>> findSubsequences(int[] nums) {
	        // 1. sort the array 
	        Arrays.sort(nums); 
	        List<List<Integer>> rst = new ArrayList<>();
	        // 2. iterate 
	        List<Integer> cur = new ArrayList<>();
	        return backTrack(nums, 0, rst, cur); 
	    }
	    
        /** This is the backtrack part: pick or not pick 
	     * the condition to add it to the rest is to check the length >= 2 
	    */
	    private List<List<Integer>> backTrack(int[] nums, int start, List<List<Integer>> rst, List<Integer> cur){
	    	if(cur.size() >= 2) rst.add(new ArrayList<>(cur)); 
//	    	if(start >= nums.length) return rst; 
	    	for(int i = start; i< nums.length; i++) {
	    		cur.add(nums[i]); 
		    	backTrack(nums, i+1, rst, cur); 
		    	cur.remove(cur.size()-1); 
	    	}
	    	return rst; 
	    }
	    
	    public List<String> addOperators(String num, int target){
	    	List<String> rst = new ArrayList<String>();
	    	if(num == null || num.length() == 0) return rst; 
	    	helperAddOp(rst, "", num, target, 0,0,0);
	    	return rst;  
	    }
	    
	    public void helperAddOp(List<String> rst, String path, String num, int target, int pos, long eval, long multed) {
	    	if(pos == num.length()) {
	    		if(target == eval)
	    			rst.add(path);
	    		return; 
	    	}
	    	for(int i = pos; i< num.length(); i++) {
	    		// 注意这里是pos == ‘0’ 跟i无关 就是只要i越过pos就不行
	    		if(i != pos && num.charAt(pos) == '0') break; 
	    		long cur = Long.parseLong(num.substring(pos,  i+1));
	    		if(pos == 0) {
	    			helperAddOp(rst, path + "+" + cur, num, target, i+1, eval + cur, cur); 
	    			helperAddOp(rst, path + "-" + cur, num, target, i+1, eval - cur, -cur); 
	    			// 这个处理很smart
	    			helperAddOp(rst, path + "*" + cur, num, target, i+1, eval - multed + multed * cur,  multed * cur); 
	    		}
	    	}
	    }
	    
	    // 扫雷 In general, I like DFS better. 
	    public char[][] updateBoard(char[][] board, int[] click){
	    	int m = board.length, n = board[0].length;
	    	int row = click[0], col = click[1];
	    	if(board[row][col]=='M') { // Mine 
	    		board[row][col] = 'X'; 
	    	}else { // Empty // Get the number of mines first
	    		int count = 0; 
	    		for(int i = -1; i<2; i++) {
	    			for(int j = -1; j<2; j++) {
	    				if(i ==0 && j==0) continue; 
	    				int r = row+i , c = col + j;
	    				if (c<0 || r >= m || c < 0 || c>=n) continue; 
	    				if(board[r][c] =='M' || board[r][c] == 'X') count ++; 
	    			}
	    		}
	    		if(count > 0) { // if it is not a 'B', stop further DFS 
	    			board[row][col] = (char)(count +'0'); 
	    		}else { // Contnue DFS to adjacent cells 
	    			board[row][col] = 'B';
	    			for(int i = -1; i<2; i++){
	    				for(int j = -1; j<2; j++) {
	    					if (i == 0 && j == 0) continue;
	    					int r = row+i , c = col + j;
		    				if (c<0 || r >= m || c < 0 || c>=n) continue; 
		    				if(board[r][c] =='E') updateBoard(board, new int[]{r,c}); 
	    				}
	    			}
	    		}
	    	}
	    	return board;
	    }
	    
	    // 801 Find Eventual Safe State (蔓延)
	    public List<Integer> eventualSafeNodes(int[][] graph){
	    	List<Integer> res = new ArrayList<>(); 
	    	if(graph == null || graph.length == 0) return res; 
	    	
	    	int nodeCount = graph.length; 
	    	int[] color = new int[nodeCount];
	    	
	    	for(int i =0; i< nodeCount; i++) {
	    		if(dfs(graph, i, color)) res.add(i); 
	    	}
	    	return res; 
	    }
	    
	    private boolean dfs(int[][] graph, int start, int[] color) {
	    	if(color[start] != 0) return color[start] ==1; 
	    	color[start] =2;
	    	
	    	for(int newNode: graph[start]) {
	    		if(!dfs(graph, newNode, color)) return false; 
	    	}
	    	color[start] =1; 
	    	return true; 
	    }
	    // 785 Is Graph Bipartite
	    public boolean isBipartite(int[][] graph) {
	    	int n = graph.length;
	    	int[] colors = new int[n];
	    	Arrays.fill(colors,-1);
	    	
	    	for(int i =0; i< n; i++) {
	    		if(colors[i] == -1 && !validColor(graph, colors, 0,1))
	    			return false; 
	    	}
	    	return true; 
	    }
	    
	    public boolean validColor(int[][] graph, int[] colors, int color, int node) {
	    	if(colors[node] != -1)
	    		return colors[node] == color; 
	    	colors[node] = color; 
	    	for(int next: graph[node]) {
	    		if(!validColor(graph, colors, 1-color, next))
	    			return false; 
	    	}
	    	return true; 
	    }
	    
	    // 788 Rotated Digits 
//	    Not sure what this questions is exactly about, but ...
	    public int rotatedDigits(int N) {
	    	int[] dp = new int[N+1];
	    	int count = 0; 
	    	for(int i =0; i<= N; i++) {
	    		if(i<10) {
	    			if(i==0 || i==1 || i==8) dp[i] =1; 
	    			else if (i == 2 || i == 5 || i == 6 || i ==9) {
	    				dp[i] = 2; 
	    				count ++;
	    			}
	    		}else {
	    			int a = dp[i/10], b = dp[i%10];
	    			if(a == 1 && b==1) dp[i] = 1; 
	    			else if (a >=1 && b >=1) {
	    				dp[i] = 2; 
	    				count ++;
	    			}
	    		}
	    	}
	    	return count; 
	    }
	    
	    // 140 Word Break II Hmmmm Interesting Easier than word ladder II 
	    public List<String> wordBreak(String s, Set<String> wordDict){
	    	return DFS(s, wordDict, new HashMap<String, LinkedList<String>>());
	    }
	    // DFS function returns an array including all substrings derived from s. 
	    private List<String> DFS(String s, Set<String> wordDict, HashMap<String, LinkedList<String>> map){
	    	if(map.containsKey(s))
	    		return map.get(s);
	    	LinkedList<String> res = new LinkedList<>(); 
	    	if(s.length() ==0) {
	    		res.add("");
	    		return res; 
	    	}
	    	for(String word: wordDict) {
	    		if(s.startsWith(word)) {
	    			List<String> sublist = DFS(s.substring(word.length()), wordDict, map);
	    			for(String sub:sublist)
	    				res.add(word + (sub.isEmpty()? "":" ") + sub);
	    		}
	    	}
	    	map.put(s, res);
	    	return res; 
	    }
	    
	    // 839 Similar String Groups 
	    public int numSimilarGroups(String[] A) {
	    	if(A.length <2) return A.length;
	    	int res = 0; 
	    	for(int i = 0; i< A.length; i++) {
	    		if(A[i] == null) continue; 
	    		String str = A[i];
	    		A[i] = null; 
	    		res ++;
	    		dfs2(A,str);
	    	}
	    	return res; 
	    }
	    // srsly bad naming 
	    private void dfs2(String[] arr, String str) {
	    	for(int i =0; i<arr.length; i++) {
	    		if(arr[i] == null) continue; 
	    		if(helper2(str, arr[i])) { // both string str and arr[i] belong in same group 
	    			String s = arr[i];
	    			arr[i] = null; 
	    			dfs2(arr,s);
	    		}
	    	}
	    }
	    
	    private boolean helper2(String s, String t) {
	    	int res = 0, i = 0; 
	    	while(res <=2 && i< s.length()) {
	    		if(s.charAt(i) != t.charAt(i)) res ++ ;
	    		i++; 
	    	}
	    	return res == 2 || res ==0;
	    }
	    
	    // 97. Interleaving String Wait Wait so what is this question again? 
	    public boolean isInterLeave(String s1, String s2, String s3) {
	    	if(s3.length() != s1.length() + s2.length())
	    		return false; 
	    	boolean[][] table = new boolean[s1.length() +1][s2.length() +1];
	    	for(int i =0; i <s1.length() +1; i++) {
	    		for(int j = 0; j<s2.length(); j++) {
	    			if(i ==0 && j == 0)
		    			table[i][j] = true; 
		    		else if (i ==0)
		    			table[i][j] = table[i][j-1] && s2.charAt(j-1) == s3.charAt(i+j-1);
		    		else if (j == 0)
		    			table[i][j] = table[i-1][j] && s1.charAt(i-1) == s3.charAt(i+j-1);
		    		else
		    			table[i][j] = (table[i-1][j] && s1.charAt(i-1) == s3.charAt(i+j-1)) || (table[i][j-1] && s2.charAt(j-1) == s3.charAt(i+j-1)); 
	    		}
	    	}	
	    	return table[s1.length()][s2.length()];
	    }
	    
	    // hmmm 这道题跟那个1、2、3涂色finally safe state很像
	    class Graph{
	    	private final int V;
	    	private final List<List<Integer>> adj;
	    	
	    	public Graph(int V) {
	    		this.V = V; 
	    		adj = new ArrayList<>(V);
	    		for(int i =0; i<V; i++)
	    			adj.add(new LinkedList<>());
	    	}
	    }
	    
//	    private boolean isCyclicUtil(int i, boolean[] visited, boolean[] recStack) {
//	    	// Mark the current node as visited and part of recursion stack 
//	    	if(recStack[i])
//	    		return true;
//	    	if(visited[i])
//	    		return false; 
//	    	visited[i] = true; 
//	    	recStack[i] = true; 
//	    	// 这里很有意思的地方是backtrack给自己的孩子传的东西是带着自己的特质的 
//	    	// 但是回到sibling层的时候传给sibling的时候是要还原父母传给自己的状态的
//	    	List<Integer> children = adj.get(i); 
//	    	for(Integer c: children) {
//	    		if(isCyclicUtil(c, visited, recStack))
//	    			return true; 
//	    	}
//	    	recStack[i] = false; 
//	    	return false; 
//	    }
//	    
//	    private void addedge(int source, int dest) {
//	    	adj.get(source).add(dest); 
//	    }
//	    
//	    // returns true if the graph contains a cycle, else false 
//	    private boolean isCyclic() {
//	    	// Mark all the vertices as not visited and not part of recursion stack
//	    	boolean[] visited = new boolean[V]; 
//	    	boolean[] recStack = new boolean[V];
//	    	
//	    	// Call the recursion helper function to detect cycle in different DFS trees 
//	    	// detect cycle in different DFS trees 
//	    	for(int i =0; i<V; i++) {
//	    		if(isCyclicUtil(i, visited, recStack))
//	    			return true; 
//	    	}
//	    	return false; 
//	    }
	    
	    // Union Find Hmmm where are these two functions?
//	    int isCycle(Graph graph) {
//	    	// allocate memory for creating V subsets
//	    	int parent[] = new int[graph.V]; 
//	    	
//	    	// Initialize all subsets as signel element sets 
//	    	for(int i =0; i<graph.V; i++)
//	    		parent[i] = -1; 
//	    	// Integrate through all edges of graph, find subset of both
//	    	// vertices of every edge, if both subsets are same, then there 
//	    	// is cycle in graph 
//	    	for(int i =0; i<graph.hashCode(); i++) {
//	    		int x = graph.find(parent, graph.edge[i].sec);
//	    		int y = graph.find(parent, graph.edge[i].dest);
//	    		 
//	    		if(x==y) return 1; 
//	    		graph.Union(parent, x, y); 
//	    	}
//	    	return 0; 
//	    }
	    public List<Integer> numIslands2(int m, int n, int[][] positions){
	    	List<Integer> rst = new ArrayList<>(); 
	    	if(m<=0 || n<=0) return rst; 
	    	
	    	int count = 0; // number of islands 
	    	int[] roots = new int[m*n]; // one island = one tree; 
	    	Arrays.fill(roots, -1);
	    	
	    	for(int[] p: positions) {
	    		int root = n * p[0] + p[1]; // assume new point is isoloated island 
	    		roots[root] = root; // add new island 
	    		count ++;  
	    		for(int[] dir: dirs) {
	    			int x = p[0] + dir[0];
	    			int y = p[1] + dir[1];
	    			int nb = n * x + y;
	    			if(x<0 || x>=m || y<0 || roots[nb] ==-1) continue; 
	    			int rootNb = findIsland(roots, nb); 
	    			if(root != rootNb) { // if neighbor is in another island 
	    				roots[root] = rootNb; // union tow islands 
	    				root = rootNb; // current tree root = joined tree root 
	    				count --; 
	    			}
	    		}
	    		rst.add(count);
	    	}
	    	return rst; 
	    }
	    // 684 redundant connection   
	    public int[] findRedundantConnection(int[][] edges) {
	    	int[] parent = new int[2001];
	    	for(int i = 0; i< parent.length; i++) parent[i] = i; 
	    	for(int[] edge: edges) {
	    		int f = edge[0], t = edge[1];
	    		if(find(parent, f) == find(parent, t)) return edge; 
	    		else parent[find(parent, f)] = find(parent, t);
	    	}
	    	return new int[2];
	    }
	    // why this is not instead while loop? 
	    private int find(int[] parent, int f) {
	    	if(f != parent[f])
	    		parent[f] = find(parent, parent[f]);
	    	return parent[f];
	    }
	    
	    private int findIsland(int[] roots, int id) {
	    	while(id != roots[id]){
	    		roots[id] = roots[roots[id]];// only one line added
	    		id = roots[id]; 
	    	}
	    	return id;
	    }
	    
	   
	    /**
	     * I believe that must be better methods than Map<Integer, Map<Integer, Integer>>
	     * @param pars
	     */
	    public int kthSmallest(int[][] matrix, int k) {
	    	return 0; 
	    }
	    

	    class Node {
	        public int val;
	        public Node next;

	        public Node() {}

	        public Node(int _val,Node _next) {
	            val = _val;
	            next = _next;
	        }
	    }
	    
	    /**
	     * 708. Insert into a Cyclic Sorted List
	     * For this question, the only problem is that the last link
	     * can only add sth inside of it iff we have traversed all the other 
	     * links and they all cannot work 
	     * @param pars
	     */
	    public Node insert(Node head, int insertVal) {
			return head;
	        
	    }
	    
	    
	    public int[][] allCellsDistOrder(int R, int C, int r0, int c0) {
	        int[][] rst = new int[R*C][2];
	        int cnt = 1; 
	        int maxDis = Math.max(Math.max(Math.max(r0+c0, R-1-r0+C-1-c0), r0 + C-1-c0), R-1-r0+ c0);
	        rst[0] = new int[] {r0, c0};
	        for(int dis = 1; dis<=maxDis; dis++){
	            for(int i =0; i<= dis; i++){
	                int j = dis-i;
	                if(r0+i<R && r0+i >=0 && c0+j<C && c0+j>=0){
	                    rst[cnt] = new int[]{r0+i, c0+j};
	                    cnt++; 
	                }
	                if(i != 0 && r0-i<R && r0-i >=0 && c0+j<C && c0+j>=0){
	                    rst[cnt] = new int[]{r0-i, c0+j};
	                    cnt++; 
	                }
	                    
	                if(j!= 0 && r0+i<R && r0+i >=0 && c0-j<C && c0-j>=0){
	                    rst[cnt] = new int[]{r0+i, c0-j}; 
	                    cnt++;
	                }
	                    
	                if(i!= 0 && j!=0 && r0-i<R && r0-i >=0 && c0-j<C && c0-j>=0){
	                    rst[cnt] = new int[]{r0-i, c0-j};
	                    cnt++; 
	                }
	                     
	            }
	        }
	        return rst; 
	    }
	    
	    String st = "";
	    HashSet<String> set = new HashSet<>(); 
	    
	    public int[] numMovesStones(int a, int b, int c) {
	        int min = Math.min(c, Math.min(a,b));
	        int max = Math.max(c, Math.max(a,b));
	        int mid = a + b + c - min - max; 
	        if(mid + 1 == max && mid -1 == min) return new int[]{0,0}; 
	        int ret = max - mid + 1 + mid - min+1; 
	        if(min+1 == mid || mid+1 == max) return new int[]{1, ret}; 
	        return new int[]{2, ret}; 
	    }
	    
		    HashSet<String> visited = new HashSet<>();
		    public int[][] colorBorder(int[][] grid, int r0, int c0, int color) {
		        // DFS or BFS 
		        int[][] newGrid = new int[grid.length][grid[0].length];
		        for(int i = 0; i<grid.length; i++) {
		        	for(int j = 0; j< grid[0].length; j++) {
		        		newGrid[i][j] = grid[i][j]; 
		        	}
		        }
			     
			     int curColor = grid[r0][c0]; 
			     return colorHelper(grid, r0, c0, color,curColor, newGrid); 
		    }
		    
		    private int[][] colorHelper(int[][] grid, int r0, int c0, int color, int curCol,int[][] newGrid){
		    	int m = grid.length;
		    	int n = grid[0].length;
	            visited.add(r0+"," + c0);
	            if(r0 == 0 || r0 == m-1 || c0== 0 || c0== n-1 || ((r0 + 1 >=0 && r0+1 < m && grid[r0+1][c0] != curCol) || 
		            		(r0 - 1 >=0 && r0-1 < m && grid[r0-1][c0] != curCol) || 
		            		(c0 - 1 >=0 && c0-1 < n && grid[r0][c0-1] != curCol) ||  
		            		((c0 + 1 >=0 && c0+1 < n && grid[r0][c0+1] != curCol)))) {
	            		newGrid[r0][c0] = color; 
		            }
		        for(int[] dir: dirs) {
		        	int[] cur = new int[]{r0 + dir[0], c0 + dir[1]}; 
		        	if(visited.contains(cur[0]+ ","+cur[1]) || cur[0] < 0 || cur[0] > m-1 || cur[1] < 0 || cur[1] > n-1 || grid[cur[0]][cur[1]]!=curCol)
		        		continue; 
		            // if(cur[0] == 0 || cur[0] == m-1 || cur[1] == 0 || cur[1] == n-1 || ((cur[0] + 1 >=0 && cur[0]+1 < m && grid[cur[0]+1][cur[1]] != curCol) || 
		            // 		(cur[0] - 1 >=0 && cur[0]-1 < m && grid[cur[0]-1][cur[1]] != curCol) || 
		            // 		(cur[1] - 1 >=0 && cur[1]-1 < n && grid[cur[0]][cur[1]-1] != curCol) ||  
		            // 		((cur[1] + 1 >=0 && cur[1]+1 < n && grid[cur[0]][cur[1]+1] != curCol)))) {
		            // 	grid[cur[0]][cur[1]] = color; 
		            // }
		        	newGrid = colorHelper(grid, cur[0], cur[1], color, curCol,newGrid); 
		        }
		        return newGrid; 
		   }
//	    HashSet<String> visited = new HashSet<>();
//	    public int[][] colorBorder(int[][] grid, int r0, int c0, int color) {
//	        // DFS or BFS 
//	    	  
//		     visited.add(r0+"," + c0);
//		     int curColor = grid[r0][c0]; 
//		     return colorHelper(grid, r0, c0, color,curColor); 
//	    }
//	    
//	    private int[][] colorHelper(int[][] grid, int r0, int c0, int color, int curCol){
//	    	int m = grid.length;
//	    	int n = grid[0].length; 
//	        for(int[] dir: dirs) {
//	        	int[] cur = new int[]{r0 + dir[0], c0 + dir[1]}; 
//	        	if(visited.contains(cur[0]+ ","+cur[1]) || cur[0] < 0 || cur[0] > m-1 || cur[1] < 0 || cur[1] > n-1 || grid[cur[0]][cur[1]]!=curCol)
//	        		continue; 
//	            if(cur[0] == 0 || cur[0] == m-1 || cur[1] == 0 || cur[1] == n-1 || ((cur[0] + 1 >=0 && cur[0]+1 < m && grid[cur[0]+1][cur[1]] != curCol) || 
//	            		(cur[0] - 1 >=0 && cur[0]-1 < m && grid[cur[0]-1][cur[1]] != curCol) || 
//	            		(cur[1] - 1 >=0 && cur[1]-1 < n && grid[cur[0]][cur[1]-1] != curCol) ||  
//	            		((cur[1] + 1 >=0 && cur[1]+1 < n && grid[cur[0]][cur[1]+1] != curCol)))) {
//	            	grid[cur[0]][cur[1]] = color; 
//	            }
//	            grid = colorHelper(grid, cur[0], cur[1], color, curCol); 
//	        }
//	        return grid; 
//	   }
	    HashMap<String, Integer> mapHmax = new HashMap<>(); 
	    public int HmaxUncrossed(int[] A, int[] B, int stA, int stB) {
	    	if (mapHmax.containsKey(stA+ ","+stB)) return mapHmax.get(stA+ ","+stB); 
	    	if(stA == A.length || stB == B.length) return 0; 
	    	int rst= 0; 
	        int lenB = B.length; 
	        for(int i =stB; i< lenB; i++) {
	        	for(int j =stA; j<A.length; j++) {
	        		if(A[j] == B[i]) {
	        			rst = Math.max(rst, Math.max(1+ HmaxUncrossed(A, B, j+1, i+1), HmaxUncrossed(A, B, j, i+1))); 
//	        			System.out.println("1 is " + HmaxUncrossed(A, B, stA+1, stB+1) + "stA is" + stA + "stB is" + stB);
//	        			System.out.println("2 is "+HmaxUncrossed(A, B, stA, stB+1) + "stA is" + stA + "stB is" + stB);
	        		}
	        	}
	        }
	        mapHmax.put(stA+ ","+stB, rst); 
	        return rst; 
	    }
	    
	    public int maxUncrossedLines(int[] A, int[] B) {
	    	return HmaxUncrossed(A,B,0,0); 
	    }
	    
	    public boolean isEscapePossible(int[][] blocked, int[] source, int[] target) {
	        return HisEscapePossible(blocked, source, target); 
	    }
	    
	    private boolean HisEscapePossible(int[][] blocked, int[] source, int[] target) {
	    	if(visited.contains(source[0] + ","  + source[1])) return false; 
			for(int[] dir: dirs) {
				  int x = source[0] + dir[0], y = source[1]+dir[1]; 
				  if(!visited.contains(x+","+y)) {
					  
				  }
			}
			return false; 
	    }
	   
	    public TreeNode bstToGst(TreeNode root) {
	    	if(root == null) return null; 
	    	TreeNode right = bstToGst(root.right);  
	    	TreeNode leftM = leftMost(right); 
	    	int diff = leftM.val;
	    	
	    	TreeNode newRoot = new TreeNode(root.val + diff);
	    	newRoot.right = right; 
	    	TreeNode left = bstToGst(root.left); 
	    	
	    	newRoot.left = helperBstToGst(left, newRoot.val); 
	        return newRoot; 
	    }
	    
	    private TreeNode leftMost(TreeNode root) {
	    	if(root == null) return null; 
	    	TreeNode search = root; 
	    	while(search.left!= null) {
	    		search = search.left; 
	    	}
	    	if(search.right == null) return search;
	    	return leftMost(search.right); 
	    }
	    private TreeNode helperBstToGst(TreeNode root, int add) {
	    	if(root == null) return root; 
	    	TreeNode right = helperBstToGst(root.right,  add); 
	    	TreeNode left = helperBstToGst(root.left,  add); 
	    	TreeNode newRoot = new TreeNode(root.val + add);
	    	newRoot.left = left; 
	    	newRoot.right = right; 
	    	return newRoot; 
	    }
	    
	    public int lastStoneWeight(int[] stones) {
	    	Arrays.sort(stones);
	    	ArrayList<Integer> arr = new ArrayList<>(); 
	    	for(int i = stones.length-1; i>=0; i--) {
	    		arr.add(stones[i]); 
	    	}
	    	System.out.println(Arrays.toString(stones));
	    	while(arr.size() > 1) {
	    		int comp = arr.get(0) - arr.get(1);
	    		System.out.println("0: "+ arr.get(0) + " 1: "+ arr.get(1));
	    		if(comp == 0) {
	    			arr.remove(0);
	    			arr.remove(0);
	    		}
	    		else if(comp <0) {
	    			System.out.println("< branch");
	    			arr.add(-comp);
	    			arr.remove(0);
	    			arr.remove(0); 
	    			Collections.sort(arr, Collections.reverseOrder());
		    	}else {
		    		arr.add(comp);
		    		arr.remove(0);
		    		arr.remove(0);
	    			Collections.sort(arr, Collections.reverseOrder());
		    	}
	    		System.out.println(arr.size());
	    	} 
	    	if(arr.size() == 0) return 0; 
	    	return arr.get(0); 
	    }
	    
	    public String removeDuplicates(String S) {
	    	if(S.length() == 1) return S; 
	    	for(int i=0; i< S.length()-1; i++) {
	        	if(S.charAt(i) == S.charAt(i+1)) {
	        		String st1 = removeDuplicates(S.substring(0, i));
	        		String st2 = removeDuplicates(S.substring(i+1));
	        		while(st1.length() > 0 && st2.length() > 0 && 
	        				st1.charAt(st1.length()-1) == st2.charAt(0)) {
	        			st1 = st1.substring(0, st1.length()-1); 
	        			st2 = st2.substring(1); 
	        		}
	        		return st1 + st2; 
	        	}
	        }
	    	return S; 
	    }
//	        char[] arr = S.toCharArray();
//	        
//	        boolean con = false; 
//	        for(int i=0; i< arr.length-1; i++) {
//	        	if(arr[i] == arr[i+1]) {
//	        		
//	        	}
//	        }
	    
//	    HashMap<String, Integer> lSCmap = new HashMap<>(); 
//	    public int longestStrChain(String[] words) {
//	    	HashMap<String, ArrayList<String>> map = new HashMap<>(); 
//	    	Arrays.sort(words);  
////	    	for(int i = 0; i<words.length; i++) {
////	    		words[i]; 
////	    	}
//	        for(int i = 0; i<words.length; i++) {
//	        	map.put(words[i], new ArrayList<>()); 
//	        	for(int j = i+1; j<words.length; j++) {
//	        		if(words[j].length() == words[i].length() +1) {
//	        			for(int k = 0; k<words[j].length(); k++) {
//	        				if(words[i].equals(words[j].substring(0, k) + words[j].substring(k+1))) {
//	        					map.get(words[i]).add(words[j]); 
//	        				}
//	        			}
//	        		}
//	        	}
//	        }
//	        int max = 1;
//	        for(String s: map.keySet()) {
//	        	max = lSChelper(s, map, max); 
//	        }
//	        return max; 
//	    }
//	    
//	    
//	    private int lSChelper(String s, HashMap<String, ArrayList<String>> map, int max) {
//	    	int len  = 0; 
//	    	int cur = 1; 
//	    	for(String str: map.get(s)) {
//	    		len = lSChelper(str, map, 1);
//	    		if(1+len > cur) cur = 1+len; 
//	    		if(1+len > max) max = len;
//	    	}
//	    	lSCmap.put(s, cur); 
//	    	return max; 
//	    }
	    
	       HashMap<String, Integer> lSCmap = new HashMap<>(); 
		    public int longestStrChain(String[] words) {
		    	HashMap<String, ArrayList<String>> map = new HashMap<>(); 
		    	Arrays.sort(words, Comparator.comparing(String::length).thenComparing(String::compareTo));
//		    	Arrays.sort((String s1, String s2) -> );  
//		    	for(int i = 0; i<words.length; i++) {
//		    		words[i]; 
//		    	}
		        for(int i = 0; i<words.length; i++) {
		        	map.put(words[i], new ArrayList<>()); 
		        	for(int j = i+1; j<words.length; j++) {
		        		if(words[j].length() == words[i].length() +1) {
		        			for(int k = 0; k<words[j].length(); k++) {
		        				if(words[i].equals(words[j].substring(0, k) + words[j].substring(k+1))) {
		        					map.get(words[i]).add(words[j]); 
		        				}
		        			}
		        		}
		        	}
		        }
		        int max = 1;
		        for(String s: map.keySet()) {
		        	int cur = lSChelper(s, map);
	                if(max < cur) max = cur; 
		        }
		        return max; 
		    }
		    
		    
		    private int lSChelper(String s, HashMap<String, ArrayList<String>> map) {
	            if(lSCmap.containsKey(s)) return lSCmap.get(s); 
		    	int len  = 0; 
		    	int cur = 1; 
		    	for(String str: map.get(s)) {
		    		len = lSChelper(str, map);
		    		if(1+len > cur) cur = 1+len; 
		    		// if(1+len > max) max = len;
		    	}
		    	lSCmap.put(s, cur); 
		    	return cur; 
		    }
	    
		    public boolean isRobotBounded(String instructions) {
		    	int[] loc = new int[]{0,0}; 
		    	HashSet<String> set = new HashSet<>();
		    	set.add("0,0"); 
		        int[] dir = new int[]{0,1}; 
		        for(int i =0; i< instructions.length(); i++) {
		        	char cur = instructions.charAt(i);
		        	if(cur == 'G') {
		        		loc[0] += dir[0];
		        		loc[1] += dir[1]; 
		        	}else if(cur == 'L') {
		        		if(loc[0]==0) {
		        			loc[0] = -loc[1]; 
		        			loc[1] = 0; 
		        		}else {
		        			loc[1] = loc[0]; 
		        			loc[0] = 0;
		        		}
		        	}else{
		        		if(loc[0]==0) {
		        			loc[0] = loc[1]; 
		        			loc[1] = 0; 
		        		}else {
		        			loc[1] = -loc[0]; 
		        			loc[0] = 0;
		        		}
		        	}
		        	String newKey = loc[0] + "," + loc[1]; 
		        	if(set.contains(newKey)) return true; 
		        	set.add(newKey); 
		        }
		        return false; 
		    }
		    
		    
		    class TimeMap {
		        private Map<String, TreeMap<Integer, String>> map;
		        /** Initialize your data structure here. */
		        public TimeMap() {
		            map = new HashMap<>(); 
		        }
		        
		        public void set(String key, String value, int timestamp) {
		            if(map.containsKey(key)){
		                map.put(key, new TreeMap<>()); 
		            }else{
		                map.get(key).put(timestamp, value); 
		            }
		        }
		        
		        public String get(String key, int timestamp) {
		            TreeMap<Integer, String> treeMap = map.get(key);
		            if(treeMap == null)
		            	return "";
		            Integer floor = treeMap.floorKey(timestamp);
		            if(floor==null)
		            	return ""; 
		            return treeMap.get(floor); 
		        }
		    }

		    /**
		     * Your TimeMap object will be instantiated and called as such:
		     * TimeMap obj = new TimeMap();
		     * obj.set(key,value,timestamp);
		     * String param_2 = obj.get(key,timestamp);
		     */
		 private static int sanity(int[] array) {
			 System.out.println(Arrays.toString(array));
			 array[0] = array[0]+1; 
			  System.out.println(Arrays.toString(array));
			  return 1; 
		 }
		 class Position{
			 int num; 
			 String date; 
			 
			 public Position(int num,  String date) {
				 this.num = num;
				 this.date= date; 
			 }
			 
	
		 }
//		 
//		 private List<Integer> palantir (List<List<Position>> input){
//			 return null; 
//		 }
		int num = 0; 
	    public int numTilePossibilities(String tiles) {
	    	char[] arr = tiles.toCharArray();
	    	Arrays.sort(arr);
	        bbackTrack(arr, 0); 
	        return num; 
	    }
	    
	    private void bbackTrack(char[] arr, int index) {
	    	if(index == arr.length) return; 
	    	num++; 
	    	for(int i = index; i< arr.length; i++) {
	    		if(i > index && arr[i] == arr[i-1]) continue;
	    		
	    		bbackTrack(arr, index+1); 
	    	}
	    }
	    
	    
	    public String smallestSubsequence(String text) {
	    	int len = text.length();
	        int[] cnt = new int[text.length()];
	        cnt[len-1] = 1; 
	        HashSet<Character> set = new HashSet<>();
	        set.add(text.charAt(len-1));
	        for(int i = len-2; i>=0 ; i--) {
	        	if(set.contains(text.charAt(i)))
	        		cnt[i] = cnt[i+1]; 
	        	else {
	        		set.add(text.charAt(i));
	        		cnt[i] = cnt[i+1]+1; 
	        	}	
	        }
//	        int sum = set.size(); 
	        int[][] table = new int[26][len];
	        for(int i = len-1; i>=0; i--) {
	        	if(i == len-1)
	        		table[text.charAt(i)-'a'][i] = 1; 
	        	else {
	        		for(char j = 'a'; j<='z'; j++) {
	        			if(j == text.charAt(i))
	    	        		table[text.charAt(i)-'a'][i] = table[text.charAt(i)-'a'][i+1]+1; 
	        			else
	        				table[j-'a'][i] = table[j-'a'][i+1]; 
	        		}
	        	}	
	        }
	        HashSet<Integer> used = new HashSet<>(); 
	        StringBuilder sb= new StringBuilder(); 
	        int cur = set.size();
//	        int i =0; 
	        printMatrix(table); 
	        for(int i = 0; i< len; i++) {
	        	if(cur == 0) return sb.toString(); 
//	        	while(i < len && (i == 0 || cnt[i] == cnt[i-1])) {
//		        	i++; 
//		        }
//	        	if(i!=0) i--;
	        	System.out.println("i "+i); 
	        	System.out.println("cur " +cur);
	        	System.out.println("cnti " +cnt[i]); 
//	        	System.out.println("len " +len);
//	        	System.out.println(Arrays.toString(cnt));
	        	if(cnt[i] == cur && (i == len-1 || (i+1 < len && cnt[i+1] !=cur)) ) {
	        		System.out.println("true");
	        		for(int j =0; j<26; j++) {
	        			if(!used.contains(j) && table[j][i] >0 && len - i - cnt[i] >= table[j][i]) {
	        				used.add(j); 
	        				System.out.println("j " +j);
	        				sb.append((char)('a'+j)); 
	        				System.out.println("sb "+sb.toString());
	        				cur --;
	        				break; 
	        			}
	        		}
	        	}
	        }
	        return "hi" + sb.toString(); 
	    }
	    public TreeNode sufficientSubset(TreeNode root, int limit) {
	    	return sufficientSubset( root, (long)limit); 
	    }
	    public TreeNode sufficientSubset(TreeNode root, long limit) {
	        if(root == null) return root;
	        if(root.left == null && root.right == null){
	            if(root.val < limit) return null;
	            return root; 
	        }
	        TreeNode left = sufficientSubset(root.left, (long)(limit- root.val));
	        TreeNode right = sufficientSubset(root.right, (long)(limit- root.val)); 
	        root.left = left;
	        root.right = right; 
	        if(left == null && right == null && root.val < limit) return null;
	        return root; 
	    }
	    
	    
//	    public List<String> invalidTransactions(String[] transactions) {
//	        List<String> rst = new ArrayList<>(); 
//	        Map<String, List<String>> map = new HashMap<>();
//	        String previousName = ""; 
//	        String previousLocation = ""; 
//	        String previousString = ""; 
//	        int previousTime = 0;
//	        for(int i = 0; i<transactions.length; i++) {
//	        	String st = transactions[i]; 
//	        	int index1 = st.indexOf(',');
//	        	int index2 = st.indexOf(',', index1+1);
//	        	int index3 = st.indexOf(',', index2+1);
//	        	String name = st.substring(0,index1);
//	        	String rest = st.substring(index1+1);
//	        	if(map.containsKey(rest)) {
//	        		List<String> curList = map.get(name);
//	        		curList.add(rest); 
//	        		map.put(name, curList);
//	        	}else {
//	        		List<String> curList = new ArrayList<>();
//	        		curList.add(rest);
//	        		map.put(name, curList); 
//	        	}
//	        }
//	        
//	        Set<String> keySet = map.keySet(); 
//	        Iterator<String> iter = keySet.iterator();
//            while(iter.hasNext()){
//                String st = iter.next();
//                List<String> rest = map.get(st); 
//                rest.sort((String s1, String s2) -> Integer.parseInt(s1.substring(0,  s1.indexOf(','))).compareTo(Interger.parseInt(s2.substring(0,  s2.indexOf(',')))));
//                for(int i = 0; i<rest.size(); i++) {
//    	        	String st = rest.get(i); 
//	                int index1 = st.indexOf(',');
//		        	int index2 = st.indexOf(',', index1+1);
//		        	//int index3 = st.indexOf(',', index2+1);
//		        	String time = st.substring(0,index1);
//		        	int time = Integer.parseInt(st.substring(index1+1, index2));
//	        	int money = Integer.parseInt(st.substring(index2+1, index3));
//	        	String location = st.substring(index3+1); 
//                if(name.equals(previousName) && location.equals(previousLocation) && time - previousTime >= 60) ) {
//	        		rst.add(st);
//	        		rst.add(previousString); 
//	        	}else if(money > 1000) {
//	        		rst.add(st); 
//	        	}
//	        	previousString = st; 
//	        	previousName = name; 
//	        	previousLocation = location; 
//	        	previousTime = time; 
//            }
//
//	         
//	        Arrays.sort(transactions, (String s1, String s2) -> s1.substring(0, s1.indexOf(',')).compareTo(s2.substring(0, s2.indexOf(','))) != 0 ? s1.substring(0, s1.indexOf(',')).compareTo(s2.substring(0, s2.indexOf(','))) : (Integer.parseInt(s1.substring(s1.indexOf(',')+1, s1.indexOf(',',s1.indexOf(',')+1)).compareTo(Integer.parseInt(s2.substring(s2.indexOf(',')+1, s2.indexOf(',',s2.indexOf(',')+1) ))))));
//	        System.out.println(Arrays.toString(transactions));
//	        for(int i = 0; i<transactions.length; i++) {
//	        	String st = transactions[i]; 
//	        	int index1 = st.indexOf(',');
//	        	int index2 = st.indexOf(',', index1+1);
//	        	int index3 = st.indexOf(',', index2+1);
//	        	String name = st.substring(0,index1);
//	        	int time = Integer.parseInt(st.substring(index1+1, index2));
//	        	int money = Integer.parseInt(st.substring(index2+1, index3));
//	        	String location = st.substring(index3+1); 
//	        	
//	        	if((i>0 && name.equals(previousName) && location.equals(previousLocation) && time - previousTime >= 60) ) {
//	        		rst.add(st);
//	        		rst.add(previousString); 
//	        	}else if(money > 1000) {
//	        		rst.add(st); 
//	        	}
//	        	previousString = st; 
//	        	previousName = name; 
//	        	previousLocation = location; 
//	        	previousTime = time; 
//	        }
//	        return rst; 
//	    }
	    
	    
	    public ListNode removeZeroSumSublists(ListNode head) {
	        if(head == null) return null; 
	        if(head.val == 0) return null; 
	        if(head.val != 0 && head.next == null) return head; 
	        Map<Integer, Integer> map = new HashMap<>();  // value to index 
	        ListNode start = head; 
	        int sum = head.val; 
	        map.put(sum,0); 
	        start = head.next; 
	        int min = -1; 
	        int max = 10000;
	        int len = 0; 
	        int index = 1; 
	        while( start.next != null ){
	            sum += start.val; 
	            if(sum == 0 ){
	                // head = start.next; 
	                // map = new HashMap<>(); 
	                min = 0; 
	                max = index; 
	            }else if(map.containsKey(sum)){
	                int prev = map.get(sum);
	                if(index-prev > len){
	                    min = prev;
	                    max = index;
	                    len = index-prev;
	                }
	                
	            }
	            index++;
	        }
	        
	        ListNode newS = head; 
	        ListNode newE = head; 
	        for(int i = 0; i< min; i++){
	            newS = newS.next; 
	        }
	        for(int i = 0; i< max; i++){
	            newE = newE.next; 
	        }
	        newS.next = newE.next; 
	        return head; 
	    }
	    
	    public int maximumSum1(int[] arr) {
	        int n = arr.length; 
	        int[][] dp = new int[n][n];
	        int[] start = new int[n]; 
	        for(int i =0; i<n; i++){
	            int startMax = i; 
	            for(int j = i; j<n; j++){
	                if(i ==j){
	                    dp[i][j] = arr[i];
	                }else{
	                    dp[i][j] = arr[j] + dp[i][j-1]; 
	                }
	                if(dp[i][j] > dp[i][startMax])
	                    startMax = j; 
	            }
	            start[i] = startMax; 
	        }
	        int[] end = new int[n];
	        for(int j = 0; j<n ; j++){
	            int endMin = j; 
	            for(int i = 0; i<= j; i++){
	                if(dp[i][j] > dp[endMin][j])
	                    endMin = i; 
	            }
	            end[j] = endMin; 
	        }
	        int result = Integer.MIN_VALUE; 
	        for(int i =0; i< n; i++){
	            if(dp[i][start[i]] > result) result = dp[i][start[i]]; 
	            if(i == n-1 && dp[end[i]][i] > result) result = dp[end[i]][i];
	            if(i == 1) {
	            	System.out.println("dp1");
	            	System.out.println(dp[i+1][start[i+1]]);
	            	System.out.println("dp2");
	            	System.out.println(dp[end[i-1]][i-1]);
	            }
	            if( i>=1 && i+ 1 < n && dp[i+1][start[i+1]] + dp[end[i-1]][i-1] > result ) result =dp[i+1][start[i+1]] + dp[end[i-1]][i-1]; 
	        }
	        printMatrix2(dp); 
	        System.out.println("start");
	        System.out.println(Arrays.toString(start));
	        System.out.println("end");
	        System.out.println(Arrays.toString(end));
	        return result; 
	    }
	    
		 
	    public int maxNumberOfBalloons(String text) {
	    	if(text == null) return 0; 
//	        int[] arr = new int[5];
	        int numB = 0, numA = 0, numL = 0, numO = 0, numN = 0; 
	        for(int i =0; i<text.length(); i++) {
	        	char c = text.charAt(i); 
	        	if(c == 'b') {
	        		numB ++;
	        	}else if(c == 'a') {
	        		numA ++; 
	        	}else if(c=='l') {
	        		numL ++; 
	        	}else if(c=='o') {
	        		numO++;
	        	}else if(c == 'n') {
	        		numN++;
	        	}
	        }
	        int rst = Math.min(numA, numB); 
	        rst = Math.min(rst,  numN);
	        rst = Math.min(rst, numL/2);
	        return Math.min(rst, numO/2);
	    }
	    
	    public String reverseParentheses(String s) {
	    	StringBuilder sb = new StringBuilder(); 
	    	if(s == null) return null;
	    	int i = 0; 
//	    	if(s.charAt(0) == '(') return reversePHelper(s.substring(1), 1); 
	    	while(i < s.length() && s.charAt(i) != '(') {
	    		sb.append(s.charAt(i));
	    		i++;
	    	}
	    	System.out.println(sb);
	    	if(i == s.length()) return sb.toString();
	    	String inner = reversePHelper(s.substring(i+1), 1); 
	    	System.out.println("inner");
	    	System.out.println(inner);
	    	sb.append(inner);
	    	int newStart = i+inner.length() + 3;
	    	System.out.println(newStart);;
	    	if( newStart == s.length()) return sb.toString();
//	    	while(newStart<s.length() && s.charAt(newStart)== ')') {
//	    		newStart ++; 
//	    	}
	    	StringBuilder sb2 = new StringBuilder(); 
	    	sb2.append(s.substring(newStart, s.length()-2)); 
	    	if(s.charAt(0) != '(') {
	    		sb2.append(s.charAt(s.length()-1));
	    		sb.append(sb2);
	    	}else {
	    		sb2 = sb2.reverse();
	    		sb.insert(0,  sb2);
	    	}
	    	
	        return sb.toString(); 
	    }
	    // dir = 0 => ; dir = 1 <=
	    // Prereq: s starts with {
	    // return the processed String of {}
	    private String reversePHelper(String s, int dir) {
	    	int i = 0; 
	    	StringBuilder sb = new StringBuilder(); 
	    	while(s.charAt(i) != ')') {
	    		if(s.charAt(i) == '(') {
	    			String inner = reversePHelper(s.substring(i+1), 1); 
	    			sb.append(inner); 
	    			i += inner.length() + 3; 
	    		}else {
	    			sb.append(s.charAt(i));
		    		i++;
	    		}
	    	}
	    	System.out.println("dir");
	    	System.out.println(dir);
//	    	if(dir == 0) return sb.toString();
//	    	else return sb.reverse().toString();
	    	return sb.reverse().toString();
	    }
	    
	    
//	    public int kConcatenationMaxSum(int[] arr, int k) {
//	        int n  = arr.length; 
//	        int[] dp = new int[n];
//	        dp[n-1] = arr[n-1];
//	        
//	        if(k == 1) {
//	        	int rst = Math.max(0, arr[n-1]);
//	        	for(int i = n-2; i>=0; i--) {
//	        		dp[i] = Math.max(arr[i]+dp[i+1], arr[i]); 
//	        		rst = Math.max(rst, dp[i]);
//	        	}
//	        	return rst;
//	        }else {
//	        	int[] newArr = new int[2*n];
//	        	for(int i =0; i<n; i++) {
//	        		newArr[i] = arr[i];
//	        	}
//	        	for(int i =0; i<n; i++) {
//	        		newArr[n+i] = arr[i];
//	        	}
//	        }
//	    }
	    
//	    public static List<String> funWithAnagrams(List<String> s) {
//	        // Write your code here
//	    	
//	    }
//	    public int[] sf1(int n) {
//	    	LinkedList<Integer> res = new LinkedList<>(); 
//	    	int mask = 1; 
//	    	int count = 0; 
//	    	for(int i = 0; i < 32; i++){ //找到1的位置存下来part
//	    		if((n & mask)!=0){
//	    			count++;
//	    			res.addFirst(i);
//	    		}
//	    		mask <<= 1;
//	    	}
//	    	res.addFirst(count);
////	    	return res.toArray(new LinkedList<Integer>()); 
//	    }
	    
	    public static List<String> funWithAnagrams(List<String> s) {
	        // Write your code here
	    	HashSet<String> set = new HashSet<>();
	    	List<String> rst = new ArrayList<>(); 
	    	for(int i =0 ;i < s.size(); i++) {
	    		String cur = s.get(i);
	    		String processed = helperFun(cur);
	    		if(!set.contains(processed)) {
	    			rst.add(cur);
	    			set.add(processed);
	    		}
	    	}
	    	Collections.sort(rst);
//	    	rst.sort((String a, String b) => a.compareTo(b));
	    	return rst; 
	    }
	    
	    private static String helperFun(String s) {
	    	char[] arr = s.toCharArray(); 
	    	Arrays.sort(arr);
	    	return String.valueOf(arr);
	    }
	    
	    
	    public static int maxStep(int n, int k) {
	        // Write your code here
	    	int sum =0; 
	    	int step = 1; 
	    	while(sum < k ) {
	    		sum += step; 
	    		step ++; 
	    	}
	    	int rst = (1+n) * n /2;  
	    	if(sum == k) {
	    		return rst-1; 
	    	}else {
	    		return rst; 
	    	}

	    }
	    
	    public int flowers(int[] A, int Full) {
	    	int step = 0; 
	    	int capacity = Full; 
	    	int c = -1; 
	    	int n = A.length;
	    	while(c < n) {
	    		while(c + 1 < n && A[c+1] <= capacity) {
	    			capacity = capacity - A[c+1];
	    			c++; 
	    			step++;
	    		}
	    		if( c == n-1 ) return step; 
	    		step += 2 * (c+1); 
	    		capacity = Full; 
	    	}
	    	return step; 
	    }
	    
	    public int domino(int[] A, int[] B) {
	    	int n = A.length; 
	    	HashMap<Integer, Integer> map = new HashMap<>();
	    	int target = -1; 
	    	for(int i = 0; i< n; i++) {
	    		if(A[i] == B[i]) {
	    			map.put(A[i], map.getOrDefault(A[i], 0) + 1); 
	    			if(map.get(A[i]) >= n) target = A[i]; 
	    		}else {
	    			map.put(A[i], map.getOrDefault(A[i], 0) + 1); 
	    			if(map.get(A[i]) >= n) target = A[i]; 
	    			map.put(B[i], map.getOrDefault(B[i], 0) + 1); 
	    			if(map.get(B[i]) >= n) target = B[i]; 
	    		}
	    	}
	    	if(target == -1) return -1; 
	    	int getInA = 0; 
	    	int getInB = 0; 
	    	for(int i = 0; i< n; i++) {
	    		if(A[i] != target) getInA ++;
	    		if(B[i] != target) getInB ++; 
	    	}
	    	return Math.min(getInA, getInB);
	    }
	    
	    public String[] printFrame(int n){
	    	String[] rst = new String[n];
	    	String firstAndLast = "";
	    	for(int i = 0; i< n; i++)
	    		firstAndLast += "*";
	    	rst[0] = firstAndLast; 
	    	rst[n-1] = firstAndLast; 
	    	String other = "*";
	    	for(int i = 1; i< n-1; i++) 
	    		other+=" ";
	    	other += "*"; 
	    	for(int i = 1; i< n-1; i++) {
	    		rst[i] = other; 
	    	}
	    	return rst; 
	    }
	    
	    public String concatSum(String str1, String str2) {
	    	int l1 = str1.length()-1;
	    	int l2 = str2.length()-1; 
	    	String rst = "";
	    	while(l1 >=0 && l2 >=0) {
	    		int n1 = Character.getNumericValue(str1.charAt(l1));
	    		int n2 = Character.getNumericValue(str2.charAt(l2));
	    		int curSum = n1 + n2; 
	    		String cur = ""+curSum; 
	    		rst = cur + rst;
	    		l1 --; 
	    		l2 --; 
	    	}
	    	if(l1 >=0) {
	    		rst = str1.substring(0, l1+1) + rst; 
	    	}
	    	if(l2 >=0) {
	    		rst = str2.substring(0, l2+1) + rst; 
	    	}
	    	return rst;
	    }
	    
	    public int findSumOfMaxSumNum(int[][] m, int n) {
	    	int max = 0; 
	    	int row = m.length;
	    	int col = m[0].length; 
	    	for(int i = 0; i<= row - n; i++) {
	    		for(int j = 0; j<= col - n; j++) {
	    			int curSum = 0; 
	    			for(int k1 = 0; k1 < n; k1++) {
	    				for(int k2 = 0; k2 < n; k2++) {
	    					curSum += m[i+k1][j+k2]; 
	    				}
	    			}
	    			if(curSum > max) max = curSum; 
	    		}
	    	}
	    	return max; 
	    }
	    
	    private static int ribbon(int[] arr, int k) {
	        int hi = 0;
	        for (int i : arr) hi = Math.max(hi, i);
	        int lo = 1;
	        int res = 0;
	        while (lo <= hi) {
	          int mid = (lo + hi)/2;
	          int curr = 0;
	          for (int i : arr) curr += i/mid;
	          if (curr >= k) {
	            res = Math.max(res, mid);
	            lo = mid + 1;
	          } else hi = mid - 1;

	        }
	        return res;
	      }
	        
	   

	    
	    public int maxLen(int[] a, int[] b) {
	    	HashSet<Integer> setA = new HashSet<>(); 
	    	HashSet<Integer> setB = new HashSet<>(); 
	    	for(int i = 0; i<a.length; i++)
	    		setA.add(a[i]);
	    	for(int i = 0; i<b.length; i++)
	    		setB.add(b[i]);
	    	int max = Integer.MIN_VALUE;
	    	int minDiff = Integer.MAX_VALUE;
	    	int len1 = a.length; 
	    	int len2 = b.length; 
	    	for(int i = 0; i+1< len1; i++) {
	    		minDiff = Math.min(minDiff, a[i+1]-a[i]);
	    	}
	    	ArrayList<Integer> possible = factors(minDiff);
//	    	System.out.println(possible);
	    	for(int i = 0; i< possible.size(); i++) {
	    		int diff = possible.get(i);
	    		ArrayList<Integer> curList = new ArrayList<>(); 
	    		curList.add(a[0]);
	    		int curLast = a[0];
	    		int curFirst = a[0];
	    		int count = 1; 
	    		while(setB.contains(curFirst - diff)){
	    			curList.add(0,curFirst - diff);
	    			curFirst -= diff; 
	    		}
	    		while(setA.contains(curLast + diff) || setB.contains(curLast + diff)){
	    			if(setA.contains(curLast + diff))
	    				count ++; 
	    			curList.add(curLast + diff);
	    			curLast += diff; 
	    		}
	    		
	    		if(count == len1 && curList.size() > max) max = curList.size(); 
	    		
	    	}
	    	if(max != Integer.MIN_VALUE) return max; 
	    	return -1;
	    }
	    
	    public ArrayList<Integer> factors(int n){ 
	    	ArrayList<Integer> rst = new ArrayList<>(); 
	        
	        for (int i = 1; i <= Math.sqrt(n); i ++) { 
	            if (n % i == 0) { 
	                rst.add(i);
	                rst.add(n/i);
	            } 
	        } 
	  
	        return rst; 
	    } 
//	    
//	    public boolean stringTransfer(String a, String b){
//	    	
//	    }
	    public int candy(int numPiles, int[] candyPiles, int numHours) {
			   
			   int len = candyPiles.length; 
			   int high = candyPiles[0];
			   for(int i = 0; i<len; i++) {
				   if(candyPiles[i] > high)
					   high = candyPiles[i]; 
			   }
			   int low = 1;  
			   while(low < high) {
				   int mid = low + (high-low)/2; 
				   int total = 0; 
				   for(int i = 0; i<len; i++) {
					   
					   total += candyPiles[i]/mid;
					   if(candyPiles[i]%mid != 0) 
						   total+=1; 
					   
				   }
				   if(total <= numHours)
					   high = mid; 
				   else 
					   low = mid+1; 
			   }
			   return high;
		   }
	    
	    public static List<String> costsOfNodes(List<String> lines) {
	        // Write your code here
	            
	            List<String> rst = new ArrayList<>(); 
	            int n = lines.size(); 
	            HashMap<String, ArrayList<String>> map = new HashMap<>(); 
	            for(int i = 0; i<n; i++){
	                String line = lines.get(i); 
	                String[] arr = line.split(",");
	                for(int j = 1; j<arr.length; j++){
	                    if(map.containsKey(arr[j])){
	                        ArrayList<String> lst = map.get(arr[j]);
	                        lst.add(arr[0]);
	                        map.put(arr[j], lst);
	                    }else{
	                        ArrayList<String> lst = new ArrayList<>();
	                        lst.add(arr[0]);
	                        map.put(arr[j], lst);
	                    }
	                }   
	            }

	            Map<String, Integer> res = bfs(map);
	            Iterator<String> ite = res.keySet().iterator();
	            while(ite.hasNext()) {
	            	String key = ite.next();
	            	rst.add(key + "," + res.get(key));
	            }
	            return rst; 
	             
	        }
	    
	    public static Map<String, Integer> bfs(HashMap<String, ArrayList<String>> map2){
	    	Iterator<String> iter = map2.keySet().iterator();
	    	Map<String, Integer> result = new HashMap<>(); 
	    	while(iter.hasNext()) {
	    		String curr = iter.next();
	    		HashSet<String> all = new HashSet<>();
	    		travel(map2,curr,all);
	    		result.put(curr, all.size());
	    	}
	    	return result;
	    }
	    
	    public static void travel(HashMap<String, ArrayList<String>> map2, String key, Set<String> seen) {
	    	List<String> now = map2.get(key);
	    	if(now.isEmpty()) return;
	    	for(String i: now) {
	    		seen.add(i);
	    		travel(map2,i, seen);
	    	}
	    	return;
	    }
	    
	    public void testTest(int[] arr) {
	    	arr[0] = 1; 
	    }
	    
//	    public int lengthLongestPath
	    public class Order{
	    	int price; 
	    	int quantity; 
	    	String side; 
	    	
	    	public Order(int price, int quantity, String side) {
	    		this.price = price;
	    		this.quantity = quantity; 
	    		this.side = side; 
	    	}
	    }
	    
	    public int processTransactions(List<Order> stream) {
	    	int numShares = 0; 
	    	HashMap<Integer, Integer> buyMap = new HashMap<>();
	    	HashMap<Integer, Integer> sellMap = new HashMap<>();
	    	PriorityQueue<Integer> buyQ = new PriorityQueue<>(Collections.reverseOrder()); // higher better 
	    	PriorityQueue<Integer> sellQ = new PriorityQueue<>(); // lower better 
	    	for(int i =0; i<stream.size(); i++) {
	    		Order order = stream.get(i); 
	    		String side = order.side; 
	    		int quantity = order.quantity;
	    		int price = order.price; 
	    		boolean canContinue = true; 		
	    		if(side == "Buy") {
    				while(quantity > 0 && canContinue) {
	    				if(!sellQ.isEmpty()) {
	    					int curBestPrice = sellQ.peek(); 
		    				if(curBestPrice > price) { // do nothing 
		    					canContinue = false; 
		    				}else{
		    					if(sellMap.get(curBestPrice) > quantity) {
		    						numShares += quantity; 
		    						sellMap.put(curBestPrice, sellMap.get(curBestPrice) - quantity); 
		    					}else {
		    						numShares += sellMap.get(curBestPrice);
		    						quantity -= sellMap.get(curBestPrice);
		    						sellQ.remove(curBestPrice); 
		    						sellMap.remove(curBestPrice); 
		    					}
		    				}
	    				}else canContinue = false; 
    				}
    				if(quantity > 0) {
    					buyMap.put(price, buyMap.getOrDefault(price, 0) + quantity); 
    					buyQ.add(price); 
    				}
    			}else {
    				while(quantity > 0 && canContinue) {
	    				if(!buyQ.isEmpty()) {
	    					int curBestPrice = buyQ.peek(); 
		    				if(curBestPrice < price) { // do nothing 
		    					canContinue = false; 
		    				}else{
		    					if(buyMap.get(curBestPrice) > quantity) {
		    						numShares += quantity; 
		    						buyMap.put(curBestPrice, buyMap.get(curBestPrice) - quantity); 
		    					}else {
		    						numShares += buyMap.get(curBestPrice);
		    						quantity -= buyMap.get(curBestPrice);
		    						buyQ.remove(curBestPrice); 
		    						buyMap.remove(curBestPrice); 
		    					}
		    				}
	    				}else 
	    					canContinue = false; 
    				}
    				if(quantity > 0) {
    					sellMap.put(price, sellMap.getOrDefault(price, 0) + quantity); 
    					sellQ.add(price); 
    				}
    			}
	    	}
	    	return numShares; 
	    }
	    
	   
	    
	    public static void main(String[] args){ 
	    	practice hi = new practice();
	    	practice.Order o1 = new practice().new Order(160, 5, "Buy");
	    	practice.Order o2 = new practice().new Order(158, 10, "Sell");
	    	practice.Order o3 = new practice().new Order(165, 8, "Buy");
	    	practice.Order o4 = new practice().new Order(156, 3, "Buy");
	    	List<Order> stream = new ArrayList<>(); 
	    	stream.add(o1); 
	    	stream.add(o2); 
	    	stream.add(o3); 
	    	stream.add(o4); 
	    	int num = hi.processTransactions(stream);
	    	System.out.print(num);
	    }
}
