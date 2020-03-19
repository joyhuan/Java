/** Contains class contains static methods for sorting and search using a
 Comparable[]. If they do not. If the array elements are not consistent, a
 ClassCastException it thrown. For example, having
 one element be a String and another an Integer won't work.<br><br>

 Included, in this order, are:<br>
 
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

  mergesort: to sort b[h..k] recursively<br>

  quicksort0: the basic quicksort algorithm<br>

  quicksort: quicksort0 changed to fix inefficiencies<br>

  toString: create a String that contains elements of an array<br>
  
  @Author David Gries

 */
public class ComparableArrays {
    
    /** = index of x in b ---or the length of b if x is not in b. */
    public static int linearSearch(Comparable[] b, Comparable x) {
        // invariant: x is not in b[0..i-1]
        int i;
        for (i= 0; i != b.length && x.compareTo(b[i]) != 0; i= i+1) {}
        return i;
    }

    /** Assume virtual elements b[-1] = -infinity and b.[b.length] = +infinity.
     <br>Return a value i that satisfies R: b[i] <= x < b[i+1] */
    public static int binarySearch(Comparable[] b, Comparable x) {
        int i= -1; int j= b.length;
        // inv: b[0i] <= x < b[j] and -1 <= i < j <= b.length, i.e.

        //          0--------i----------j---------- b.length
        //       b |       |<=x|    ?  |>x|       |
        //          -------------------------------
        while (j != i+1) {
            int e= (i+j)/2;
            // -1 <= i < e < j <= b.length
            if (x.compareTo(b[e]) >= 0) i= e;
            else j= e;
        }
        return i;
    }

    
    /** = the position of the minimum value of b[h..k].<br>
        Precondition: h < k. */
    public static int min(Comparable[] b, int h, int k) {
        int p= h; 
        int i= h;
        // inv: b[p] is the minimum of b[h..i], i.e.
        //
        //    h-------p------------i---------k
        // b | b[p] is min of this  |    ?    |
        //    --------------------------------

        while (i!= k) {
            i= i+1;
            if (b[p].compareTo(b[i]) > 0) {
                p= i;
            }
        }
        return p;
    }
    
    /** Sort b --put its elements in ascending order. */
    public static void selectionSort(Comparable[] b) {
        int j= 0;
        // inv P:  b[0..j-1] is sorted and b[0..j-1] <= b[j..], i.e.

        //          0---------------j--------------- b.length
        // inv : b |  sorted, <=   |    >=          |
        //          --------------------------------

        while (j != b.length) {
            int p= min(b, j, b.length-1);
            // b[p] is minimum of b[j..b.length-1]
            // Swap b[j] and b[p]
            Comparable t= b[j]; b[j]= b[p]; b[p]= t;
            
            j= j+1;
        }
    }

    /** Sort b --put its elements in ascending order. */
    public static void selectionSort1(Comparable[] b) {
        int j= 0;
        // inv P: b[0..j-1] is sorted and b[0..j-1] <= b[j..], i.e.

        //          0---------------j--------------- b.length
        // inv : b |  sorted, <=   |    >=          |
        //          --------------------------------

        while (j != b.length) {
            // Put into p the index of smallest element in b[j..]
            int p= j;
            for (int i= j+1; i != b.length; i++) {
                if (b[i].compareTo(b[p]) < 0)  p= i;
            }
            
            // Swap b[j] and b[p]
            Comparable t= b[j]; b[j]= b[p]; b[p]= t;
            j= j+1;
        }
    }

    /** Precondition: b[h..k-1] is sorted, and b[k] contains a value.<br>
       Sort b[h..k]. */
    public static void insertValue(Comparable[] b, int h, int k) {
        Comparable v= b[k];
        int i= k;
        /* inv P: (1) Placing v in b[i] makes b[h..k] a
                      permutation of its initial value
         (2) b[h..k] with b[i] removed is initial b[h..k-1]
         (3) v < b[i+1..k]
         */
        while ((i != h) && v.compareTo(b[i-1])< 0) {
            b[i]= b[i-1];
            i= i-1;
        }
        b[i]= v;
    }


    /** Sort b[h..k] --put its elements in ascending order. */
    public static void insertionSort(Comparable[] b, int h, int k) {
        // inv: h <= j <= k+1  and  b[h..j-1] is sorted, i.e.

        //          h---------------j--------------k
        // inv : b |  sorted,      |     ?          |
        //          --------------------------------

        for (int j= h; j <= k; j= j+1) {
            // Sort b[h..j], given that b[h..j-1] is sorted
            insertValue(b,h,j);
        } 
    } 

    /** b[h..k] has at least three elements.<br>
      Let x be the value initially in b[h].<br>
      Permute b[h..k] and return integer j satisfying R:<br><br>

         b[h..j-1] <= b[j] = x <= b[j+1..k]
     */
    public static int partition0(Comparable[] b, int h, int k) {
        // 
        //          h---------------------------k
        // pre:  b |x|     ?                     |  for some x
        //          -----------------------------
        // 
        //          h-------------j-------------k
        // post: b |   <= x      |x|     >= x    |
        //          -----------------------------
        int j= h;
        Comparable x= b[h];
        int i= k;
        // inv P: b[h..j-1] <= x <= b[i+1..k], i.e.
        //  
        //          h---------j-------i------------k
        // post: b |   <= x  |x|  ?    |    >= x    |
        //          --------------------------------
        while (j < i) {
            if ( x.compareTo(b[j+1]) >= 0) {
                b[j]= b[j+1];
                j= j+1;
            }
            else {
                Comparable temp= b[j+1]; b[j+1]= b[i]; b[i]= temp;
                i= i-1;
            }
        }
        b[j]= x;
        return j;
    }

    /** Permute b[h..k] and return integer j satisfying R:<br><br>

          b[h..j-1] <= b[j] = x <= b[j+1..k]<br><br>

     where x stands for the value initially in b[h].<br>
     Precondition: b[h..k] has at least three elements. */
    public static int partition(Comparable[] b, int h, int k) {
        int j;
        // Truthify R1: b[h+1..j] <= b[h] = x <= b[j+1..k];
        int i= h+1; j= k;

        // inv P:  b[h+1..i-1] <= x <= b[j+1..k], i.e.
        //
        //    h---------i------j----------k
        // b |x| <= x  |    ?   |  >= x    |
        //    -----------------------------
        while (i <= j) {
            if (b[i].compareTo(b[h]) <= 0) {
                i= i+1;
            }
            else if (b[j].compareTo(b[h]) >= 0) {
                j= j-1;
            }
            else {// b[j] < x < b[i]
                Comparable t1= b[i]; b[i]= b[j]; b[j]= t1;
                i= i+1; j= j-1;
            }
        }
        Comparable temp= b[h]; b[h]= b[j]; b[j]= temp;
        // R
        return j;
    }

    /** Permute b[h], b[(h+k)/2], and b[k] to put their median in b[h]. */
    public static void medianOf3(Comparable[] b, int h, int k) {
        int e= (h+k)/2;  // index of middle element of array
        int m;           // index of median
        // Return if b[h] is median; else store index of median in m
        if (b[h].compareTo(b[e]) <= 0) {
            if (b[h].compareTo(b[k]) >= 0) return;
            // b[h] is smallest of the three
            if (b[k].compareTo(b[e]) <= 0)  m= k;
            else  m= e;
        }
        else {
            if (b[h].compareTo(b[k]) <= 0) return;
            // b[h] is largest of the three
            if (b[k].compareTo(b[e]) <= 0) m= e;
            else m= k;
        }
        Comparable t= b[h]; b[h]= b[m]; b[m]= t;
    }

    /** = a copy of array segment b[h..k]. */
    public static Comparable[] copy(Comparable[] b, int h, int k) {
        Comparable[] c= new Comparable[k+1-h];
        // inv: b[h..i-1] has been copied to c[0..i-h-1]
        for (int i= h; i != k+1; i= i+1) {
            c[i-h]= b[i];
        }
        return c;
    }


    /** Sort b[h..k]  --put its elements in ascending order.<br>
        Precondition: Segments b[h..e] and b[e+1..k] are already sorted.
     */
    public static void merge (Comparable b[], int h, int e, int k) {
        Comparable[] c= copy(b,h,e);
        // c is a copy of original b[h..e]
        int i= h; int j= e+1; int m= 0;
        /* inv: b[h..i-1] contains its final, sorted values
         b[j..k] remains to be transferred
         c[m..e-h] remains to be transferred
         */
        for (i= h; i != k+1; i= i+1) {
            if (j <= k && (m > e-h || b[j].compareTo(c[m]) <= 0)) {
                b[i]= b[j]; j= j+1;
            }
            else {
                b[i]= c[m]; m= m+1;
            }
        }
    }


    /** Sort b[h..k]  --put its elements in ascending order. */
    public static void mergeSort(Comparable[] b, int h, int k) {
        if (h >= k) return;

        int e= (h+k)/2;
        mergeSort(b, h, e);   // Sort b[h..e]
        mergeSort(b, e+1, k); // Sort b[e+1..k]
        merge(b, h, e, k);    // Merge the 2 segments
    }

    /** Sort b[h..k]  --put its elements in ascending order. */
    public static void quickSort0(Comparable[] b, int h, int k) {
        if (k+1-h < 2) return;

        int j= partition(b,h,k);
        // b[h..j-1] <= b[j] <= b[j+1..k]
        quickSort0(b,h,j-1);
        quickSort0(b,j+1,k);
    }

    /** Sort b[h..k]  --put its elements in ascending order. */
    public static void quickSort(Comparable[] b, int h, int k) {
        tailRecursionLoop: while (true) {
            int h1= h;
            int k1= k;
         // invariant: b[h..k] is sorted if b[h1..k1] is 
            if (k1-h1 < 10) {
                insertionSort(b,h1,k1);
                return;
            }
            medianOf3(b,h1,k1);
            // b[h1] is between b[(h1+k1)/2] and b[k1]
            
            int j= partition(b,h1,k1);
            // b[h..j-1] <= b[j] <= b[j+1..k], i.e.
            //
            //    h1---------j-------k1
            // b |   <= x  |x|  >= x   |   for some x
            //    ---------------------
            if (j-h1 <= k1-j) {
                quickSort(b,h1,j-1);
                // b[j+1..k1] remains to be sorted
                h1= j+1;
            }
            else {
                quickSort(b,j+1,k1);
                // b[h1..j-1] remains to be sorted
                k1= j-1;
            }
        }
    }

    /** Return the values of b, separated by ", " and with
     the whole list delimited by "[" and "]".
     */
    public static String toString(Comparable[] b) {
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
