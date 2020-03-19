/** Contains static methods for dealing with int arrays
    --sorting, searching, etc.-- together with the methods 
  that they use. Included, in this order, are:<br>

  linearSearch: linear search an array<br>

  binarySearch: binary search a sorted array<br>

  min: to find the minimum of b[h..k]<br>

  selectionSort: to sort an array, which calls min<br>

  selectionSort1: to sort an array, which finds min in place<br>

  insertValue: to insert a value into its sorted position in b[h..k-1]<br>

  insertionSort: to sort an array<br>

  partition0: as done in class<br>

  partition: used in quicksort, done a bit more efficiently<br>

  medianOf3: swap median of b[h], b[(h+k)/2, b[k] into b[k]<br>

  copy: return a copy of b[h..k]<br>

  merge: merge sorted b[h..e] and b[e+1..k] into b[h..k]<br>

  mergeSort: to sort b[h..k] recursively<br>

  quickSort0: the basic quicksort algorithm<br>

  quickSort: quicksort0 changed to fix inefficiencies<br>

  Dutch National Flag<br>

  toString: create a String that contains elements of an array
  
  @Author David Gries

 */
public class IntArrays {

    /** = index of x in b ---or the length of b if x is not in b. */
    public static int linearSearch(int[] b, int x) {
        // invariant: x is not in b[0..i-1]
        int i;
        for (i= 0; i != b.length && x != b[i]; i= i+1) {}
        return i;
    }

    /** Assume virtual elements b[-1] = -infinity and b.[b.length] = +infinity.<br>
       Return a value i that satisfies R: b[i] <= x < b[i+1] */
    public static int binarySearch(int[] b, int x) {
        int i= -1; int j= b.length;
        // {P:b[i] <= x < b[j] and -1 <= i < j <= b.length}, i.e.

        //          0--------i------------j-------- b.length
        // post: b |       |<=x|    ?   |>x|      |
        //          -------------------------------
        while (j != i+1) {
            int e= (i+j)/2;
            // {-1 <= i < e < j <= b.length}
            if (b[e] <= x) i= e;
            else j= e;
        }
        return i;
    }

    /** Return the position of the minimum value of b[h..k].<br>
        Precondition: h < k. */
    public static int min(int[] b, int h, int k) {
        int p= h; // will contain index of minimum
        int i= h;
        // {inv: b[p] is the minimum of b[h..i]}, i.e.
        //
        //    h-------p------------i---------k
        // b | b[p] is min of this  |    ?    |
        //    --------------------------------

        while (i!= k) {
            i= i+1;
            if (b[i] < b[p])
            {p= i;}
        }
        return p;
    }

    /** Sort b --put its elements in ascending order. */
    public static void selectionSort(int[] b) {
        int j= 0;
        // {inv P: b[0..j-1] is sorted and b[0..j-1] <= b[j..]}, i.e.

        //          0---------------j--------------- b.length
        // inv : b |  sorted, <=   |    >=          |
        //          --------------------------------

        while (j != b.length) {
            int p= min(b, j, b.length-1);
            // {b[p] is minimum of b[j..b.length-1]}
            // Swap b[j] and b[p]
            int t= b[j]; b[j]= b[p]; b[p]= t;

            j= j+1;
        }
    }


    /** Sort b --put its elements in ascending order. */
    public static void selectionSort1(int[] b) {
        int j= 0;
        // {inv P: b[0..j-1] is sorted and b[0..j-1] <= b[j..]}

        //          0---------------j--------------- b.length
        // inv : b |  sorted, <=   |    >=          |
        //          --------------------------------

        while (j != b.length) {
            // Put into p the index of smallest element in b[j..]
            int p= j; 
            for (int i= j+1; i != b.length; i++) {
                if (b[i] < b[p])  p= i;
            }

            // Swap b[j] and b[p]
            int t= b[j]; b[j]= b[p]; b[p]= t;
            j= j+1;
        }
    }


    /** Precondition: b[h..k-1] is sorted, and b[k] contains a value.<br>
        Sort b[h..k] */
    public static void insertValue(int[] b, int h, int k) {
        int v= b[k];
        int i= k;
        /* inv P: (1) Placing v in b[i] makes b[h..k] a
         permutation of its initial value
         (2) b[h..k] with b[i] removed is initial b[h..k-1]
         (3) v < b[i+1..k]
         */
        while ((i != h) && v < b[i-1]) {
            b[i]= b[i-1];
            i= i-1;
        }
        b[i]= v;
    }

    /** Sort b[h..k] --put its elements in ascending order. */
    public static void insertionSort(int[] b, int h, int k) {
        // inv: h <= j <= k+1  and  b[h..j-1] is sorted

        //          h---------------j-------------k
        // inv : b |  sorted,      |     ?         |
        //          -------------------------------

        for (int j= h; j <= k; j++) {
            // Sort b[h..j], given that b[h..j-1] is sorted
            insertValue(b,h,j);
        } 
    }

    /** b[h..k] has at least three elements.<br>
      Let x be the value initially in b[h].<br>
      Permute b[h..k] and return integer j satisfying R:<br><br>

         b[h..j-1] <= b[j] = x <= b[j+1..k]
     */
    public static int partition0(int[] b, int h, int k) {
        // 
        //          h---------------------------k
        // pre:  b |x|     ?                     |  for some x
        //          -----------------------------
        // 
        //          h---------------------------k
        // post: b |   <= x      |x|     >= x    |
        //          -----------------------------
        int j= h;
        int i= k;
        // {inv P: b[h+1..i-1] <= b[h] = x <= b[j+1..k]}, i.e.
        //  
        //          h---------j-------i------------k
        // post: b |   <= x  |x|  ?    |    >= x    |
        //          --------------------------------
        while (j < i) {
            if (b[j+1] <= b[j]) {
                int temp= b[j]; b[j]= b[j+1]; b[j+1]= temp;
                j= j+1;
            }
            else {
                int temp= b[j+1]; b[j+1]= b[i]; b[i]= temp;
                i= i-1;
            }
        }

        return j;
    }

    /** b[h..k] has at least three elements.<br>
        Let x be the value initially in b[h].<br>
        Permute b[h..k] and return integer j satisfying R:<br><br>

      b[h..j-1] <= b[j] = x <= b[j+1..k]
     */
    public static int partition(int[] b, int h, int k) {
        // {Q: Let x be the value initially in b[h]}
        int j;
        // Truthify R1: b[h+1..j] <= b[h] = x <= b[j+1..k];
        int i= h+1; j= k;
        // {inv P: b[h+1..i-1] <= b[h] = x <= b[j+1..k]}, i.e.
        //  
        //
        //    h---------i------j----------k
        // b |x| <= x  |    ?   |  >= x    |
        //    -----------------------------
        while (i <= j) {
            if (b[i] < b[h]) i= i+1;
            else if (b[j] > b[h]) j= j-1;
            else {// {b[j] < x < b[i]}
                int t1= b[i]; b[i]= b[j]; b[j]= t1;
                i= i+1; j= j-1;
            }
        }
        int t= b[h]; b[h]= b[j]; b[j]= t;
        // {R}
        return j;
    }

    /** Permute b[h], b[(h+k)/2], and b[k] to put their median in b[h]. */
    public static void medianOf3(int[] b, int h, int k) {
        int e= (h+k)/2;  // index of middle element of array
        int m;           // index of median
        // Return if b[h] is median; else store index of median in m
        if (b[h] <= b[e]) {
            if (b[h] >= b[k]) return;
            // {b[h] is smallest of the three}
            if (b[e] <= b[k]) m= e;
            else m= k;
        }
        else {
            if (b[h] <= b[k]) return;
            // {b[h] is largest of the three}
            if (b[e] <= b[k]) m= k;
            else m= e;
        }
        int t= b[h]; b[h]= b[m]; b[m]= t;
    }

    /** Return a copy of array segment b[h..k]. */
    public static int[] copy(int[] b, int h, int k) {
        int[] c= new int[k+1-h];
        // inv: b[h..i-1] has been copied to c[0..i-h-1]
        for (int i= h; i != k+1; i++) 
        {c[i-h]= b[i];}
        return c;
    }


    /** Segments b[h..e] and b[e+1..k] are already sorted.<br>
        Permute their values so that b[h..k] is sorted.
     */
    public static void merge (int b[], int h, int e, int k) {
        int[] c= copy(b,h,e);
        // {c is a copy of original b[h..e]}
        int i= h; int j= e+1; int m= 0;
        /* inv: b[h..i-1] contains its final, sorted values
         b[j..k] remains to be transferred
         c[m..e-h] remains to be transferred
         */
        for (i= h; i != k+1; i++) {
            if (j <= k && (m > e-h || b[j] <= c[m])) {
                b[i]= b[j]; j= j+1;
            } else {
                b[i]= c[m]; m= m+1;
            }
        }
    }


    /** Sort b[h..k]. */
    public static void mergeSort(int[] b, int h, int k) {
        if (h >= k) return;
        
        int e= (h+k)/2;
        mergeSort(b, h, e);   // Sort b[h..e]
        mergeSort(b, e+1, k); // Sort b[e+1:k]
        merge(b, h, e, k);    // Merge the 2 segments
    }

    /** Sort b[h..k] */
    public static void quickSort0(int[] b, int h, int k) {
        if (k+1-h < 2) return;

        int j= partition(b,h,k);
        // b[h..j-1] <= b[j] <= b[j+1..k]
        quickSort(b,h,j-1);
        quickSort(b,j+1,k);
    }

    /** Sort b[h..k]. */
    public static void quickSort(int[] b, int h, int k) {
        int h1= h;
        int k1= k;
        // invariant: b[h..k] is sorted if b[h1..k1] is 
        while(true) {         
            if (k1-h1 < 10) {
                insertionSort(b,h1,k1);
                return;
            }

            medianOf3(b, h1, k1);
            // {b[h1] is between b[(h1+k1)/2] and b[k1]}

            int j= partition(b,h1,k1);
            // b[h1..j-1] <= b[j] <= b[j+1..k1], i.e.
            //
            //    h1--------j--------k1
            // b |   <= x  |x|  >= x   |   for some x
            //    ---------------------
            if (j-h1 <= k1-j) {
                quickSort(b,h1,j-1);  //sort smaller segment recursively
                // b[j+1..k1] remains to be sorted
                h1= j+1;
            }
            else {
                quickSort(b,j+1,k1); //sort smaller segment recursively
                // b[h1..j-1] remains to be sorted
                k1= j-1;
            }
        } 
    }

    /** Swap the values of b so that the negative ones are
      first, then the 0's, and finally the positive ones.
      In the original problem, the negative values are red
      balls, 0's are white balls, positive values are blue balls*/
    public static void DutchNationalFlag(int[] b) {
        //        0------------------------------------ b.length
        // pre: b|              ?                      |
        //        -------------------------------------
        // 
        //        0------------------------------------ b.length
        // post:b|   <0    |    = 0     |    >0        |
        //        -------------------------------------

        int h= 0;
        int k= 0;
        int t= b.length-1;

        //        0-------h-------k------t------------- b.length
        // inv :b|  <0   |  = 0  |   ?    |    >0      |
        //        -------------------------------------

        while (k <= t) {
            if (b[k] < 0) {
                int temp= b[k]; b[k]= b[h]; b[h]= temp;
                h= h+1;  k= k+1;
            }
            else if (b[k] == 0) {
                k= k+1;
            }
            else {
                int temp= b[k]; b[k]= b[t]; b[t]= temp;
                t= t-1; 
            }
        }

    }

    /** Return the values of b, separated by ", " and with
      the whole list delimited by "[" and "]". */
    public static String toString(int[] b) {
        String res= "[";
        // inv: res contains b[0..k-1], with "[" at
        //      beginning and values separated by ", "
        for (int k= 0; k != b.length; k= k+1) {
            if (k != 0)
                res= res + ", ";
            res= res + b[k];
        }
        return res + "]";
    }

}
