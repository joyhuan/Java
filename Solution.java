//import java.util.HashMap;
//import java.util.Set;
import java.util.*;

public class Solution{
    public int numberOfPairs(int[] a, int k ){
      if(a==null) return 0;
      // int to its index 
      HashMap<Integer, Integer> dic= new HashMap<>();
      int len = a.length;
      for(int i=0; i<len; i++){
        if(!dic.containsKey(a[i])) dic.put(a[i],i);
      }
      int count = 0; 
      // can choose itself 
      Set<Integer> keys =dic.keySet();  
      for(int i:keys){
        if(i<=k/2 && dic.containsKey(k-i)){
          count ++; 
//          dic.remove(k-a[i]); 
        } 
      }
      return count; 
    }
    
      public static void main(String[] args){
        Solution hi = new Solution(); 
        int[] arr = {1, 2, 3, 4, 5, 5, 5, 5, 6, 7, 8, 9, 1, 9};
        System.out.println(hi.numberOfPairs(arr, 10));
    }
  }